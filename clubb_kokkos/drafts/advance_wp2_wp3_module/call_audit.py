"""Compare Native/draft source expressions only; no compiled model or numerical run."""
from pathlib import Path
import runpy,re,ast,json,hashlib
ns=runpy.run_path(str(Path(__file__).resolve().parent/'source_audit.py'))
Cpp=ns['Cpp'];fexpr=ns['fexpr'];statements=ns['statements'];defs=ns['defs'];source=ns['source'];HERE=ns['HERE']

def canonical(n):
 if not isinstance(n,tuple):return n
 if n[0]=='str':return ('str',n[1].lower())
 if n[0]=='v' and re.fullmatch(r'\w+_stats_\d+',n[1]):
  field=n[1].split('_stats_')[0]
  return ('bin','/',('v','one'),('v','rcond')) if field=='rcond' else ('bin','/',('v',field),('v','dt'))
 if n[0]=='call':
  name=n[1];args=tuple(canonical(x) for x in n[2])
  if name=='real1d':return args[0]
  if name=='subview':
   args=tuple(('v','all') if a==('v','all') else a for a in args)
   return ('call','section',args)
  return ('call',name,args)
 return tuple(canonical(x) if isinstance(x,tuple) else x for x in n)
reviews=[]
for d in defs:
 name=d['name'];raw=re.search(r'(?ms)^  subroutine '+name+r'\(.*?^  end subroutine '+name+r'\b',source).group()
 if name=='advance_wp2_wp3':
  start=raw.index('    if ( wp2_sponge_damp_settings')
  end=raw.index('    end if ! wp3_sponge_damp_settings%l_sponge_damping',start)+len('    end if ! wp3_sponge_damp_settings%l_sponge_damping')
  raw=raw[:start]+raw[end:]
 expected=[]
 for line in statements(raw):
  if not re.match(r'call\s+',line,re.I):continue
  if re.match(r'call\s+band_solve',line,re.I):continue
  line=re.sub(r'^call\s+','',line,flags=re.I)
  line=re.sub(r'clubb_params\(\s*:\s*,\s*(\w+)\s*\)',r'section(clubb_params,all,\1)',line,flags=re.I)
  line=re.sub(r'l_count_sample\s*=\s*','',line,flags=re.I)
  expected.append(canonical(fexpr(line)))
 names={x[1] for x in expected};actual=[]
 for m in re.finditer(r'\b(\w+)\s*\(',d['body']):
  if m.group(1).lower() not in names:continue
  pos=m.end();depth=1;end=pos
  while depth:depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  actual.append(canonical(Cpp(d['body'][m.start():end],set()).parse()))
 if expected!=actual:
  for i,(a,b) in enumerate(zip(expected,actual)):
   if a!=b:print(name,'call difference',i,'Native',str(a)[:1800],'Cpp',str(b)[:1800]);break
  print('counts',len(expected),len(actual));raise SystemExit(1)
 reviews.append({'name':name,'ordered_calllists':len(expected),'mapped_arguments':sum(len(c[2]) for c in expected)})
record={'source_review':'All Native CALLs and full argument expression trees in order, excluding4 exact driver-scope sponge budget calls and2 manually reviewed solver language-adapter calls','routines':reviews,'total_calllists':sum(x['ordered_calllists'] for x in reviews),'total_arguments':sum(x['mapped_arguments'] for x in reviews),'target_sha256':hashlib.sha256((HERE/'advance_wp2_wp3_module.cpp').read_bytes()).hexdigest(),'runtime_evidence':'None: source-only review.'}
(HERE/'call_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print('All',record['total_calllists'],'ordered calllists,',record['total_arguments'],'mapped arguments checked')
