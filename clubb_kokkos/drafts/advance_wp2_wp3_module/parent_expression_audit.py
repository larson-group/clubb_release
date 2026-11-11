"""Compare Native/draft source expressions only; no compiled model or numerical run."""
from pathlib import Path
import runpy,re,json
n=runpy.run_path(str(Path(__file__).resolve().parent/'source_audit.py'));Cpp=n['Cpp'];fexpr=n['fexpr'];reviews=[]
for d in n['defs']:
 if d['name'] not in ['advance_wp2_wp3','wp23_solve','wp23_rhs']:continue
 raw=re.search(r'(?ms)^  subroutine '+d['name']+r'\(.*?^  end subroutine '+d['name']+r'\b',n['source']).group()
 if d['name']=='advance_wp2_wp3':
  start=raw.index('    if ( wp2_sponge_damp_settings');end=raw.index('    end if ! wp3_sponge_damp_settings%l_sponge_damping',start)+len('    end if ! wp3_sponge_damp_settings%l_sponge_damping');raw=raw[:start]+raw[end:]
 decls={}
 for line in n['statements'](raw):decls.update(n['declaration'](line))
 arrays={k for k,a in decls.items() if a.dims};arrays.update(re.findall(r'\b(gr_\w+)\s*\(',d['body']));arrays.add('pdf_implicit_coefs_terms_coef_wp4_implicit')
 expected=[]
 for line in n['statements'](raw):
  m=re.match(r'^(\w+\([^()]*\))\s*=(?!=)\s*(.*)$',line)
  if not m or ':' in m.group(1) or m.group(1).startswith('wp2_min_array'):continue
  lhs,rhs=m.groups();expected.append((fexpr(lhs),fexpr(rhs)))
 names={e[0][1] for e in expected};actual=[]
 for m in re.finditer(r'(?m)^\s*(\w+\([^;{}]*?\))\s*=(?!=)\s*([^;]*);',d['body']):
  if m.group(1).split('(')[0].lower() not in names:continue
  actual.append((Cpp(m.group(1),arrays).parse(),Cpp(m.group(2),arrays).parse()))
 if expected!=actual:
  for i,(a,b) in enumerate(zip(expected,actual)):
   if a!=b:print(d['name'],'difference',i,'Native',str(a)[:2200],'Cpp',str(b)[:2200]);break
  print('length',len(expected),len(actual));raise SystemExit(1)
 reviews.append({'name':d['name'],'ordered_explicit_indexed_assignments':len(expected)})
record={'canonical_sha256':n['hashlib'].sha256(n['NATIVE'].read_bytes()).hexdigest(),'target_sha256':n['hashlib'].sha256(n['TARGET'].read_bytes()).hexdigest(),'reviews':reviews,'note':'Explicit indexed scientific/budget assignment expression trees; whole-array language adapters and wp2_min MAX/warning/ceiling retained separate manual review.'}
(n['HERE']/'parent_expression_audit.json').write_text(json.dumps(record,indent=2)+'\n');print(reviews)
