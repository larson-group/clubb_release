from pathlib import Path
import ast,re,json,sys,hashlib
R=Path(__file__).resolve().parents[3];sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'));from fortran_procedures import statements,declaration
from generate_derived_types import split_top
P=R/'clubb_kokkos/drafts/advance_helper_module/source_audit.py';s=P.read_text();s=s[s.index('class Cpp:'):s.index('reviews=[]')];s=s.replace("    assert b in [('n',2.),('n',3.)]","    if b not in [('n',2.),('n',3.)]:return ('call','pow',(a,b))")
line=next(l for l in s.splitlines() if "if name in {'max','min'}:return" in l)
indent=line[:len(line)-len(line.lstrip())]
s=s.replace(line,indent+"if name in {'max','min'}:\n"+indent+" out=args[0]\n"+indent+" for a in args[1:]:out=('if',('bin','>' if name=='max' else '<',out,a),out,a)\n"+indent+" return out")
exec(s)
def strip_comments(t):return re.sub(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*.*?\*/',lambda m:m.group() if m.group()[0] in '\"\'' else '',t,flags=re.S)
def balanced(t,start,op='(',cl=')'):
 depth=1;i=start;quote=''
 while depth:
  c=t[i]
  if quote:
   if c==quote and t[i-1]!='\\':quote=''
  elif c in '\"\'':quote=c
  else:depth+=(c==op)-(c==cl)
  i+=1
 return i

def review(module):
 D=R/'clubb_kokkos/drafts'/module;O=json.loads((D/'SOURCE_OUTLINE.json').read_text());S=(R/('src/CLUBB_core/'+module+'.F90')).read_text();T=(D/(module+('.hpp' if module=='Nc_Ncn_eqns' else '.cpp'))).read_text();C=strip_comments(T);names=O['native_inventory_and_order'];defs=[]
 for m in re.finditer(r'\b(?:void|real_t)\s+(\w+)\s*\(',C):
  if m.group(1) not in names:continue
  e=balanced(C,m.end());p=e
  while C[p].isspace():p+=1
  if C[p]!='{':continue
  end=balanced(C,p+1,'{','}');defs.append({'name':m.group(1),'args':[x.strip().split()[-1] for x in split_top(C[m.end():e-1])],'body':C[p+1:end-1]})
 defs=defs[:6] if module=='Nc_Ncn_eqns' else defs[-4:];assert [d['name'] for d in defs]==names
 arrays=set()
 for p in O['routines']:
  for a in p['arguments']:
   if a['dims']:arrays.add(a['name'].lower())
 for m in re.finditer(r'\b(?:Real[1-4]D|Bool1D)\s+(\w+)\s*\(',C):arrays.add(m.group(1).lower())
 arrays|={a+'_host' for a in arrays}
 def canon(n):
  if not isinstance(n,tuple) or not n:return n
  if n[0]=='v':
   name=n[1].replace('_host','')
   if name in {'j','k','free_axis'} or name.startswith('slice_'):return ('v','cell')
   if name=='cloud_frac_maxval':return ('call','maxval',(('call','cloud_frac',(('v','cell'),('v','cell'))),))
   if name=='ice_supersat_frac_maxval':return ('call','maxval',(('call','ice_supersat_frac',(('v','cell'),('v','cell'))),))
   return ('v',name)
  if n[0]=='call':
   name=n[1].replace('nc_ncn_eqns_detail::','').replace('_host','');args=tuple(canon(x) for x in n[2])
   if name=='any' and args[0]==('v','l_frozen_hm'):return ('v','frozen_hm_present')
   if name in {'any','all'} and 'hydromet' in str(args):return ('v','any_hydromet' if name=='any' else 'all_hydromet_below')
   if name=='real1d' and len(args)==1:return args[0] # rank1stats overload adapter
   return ('call',name,args)
  return tuple(canon(x) if isinstance(x,tuple) else x for x in n)
 def f(t):
  t=t.replace('err_info%err_code','err_code')
  t=re.sub(r'(?<!\w):','free_axis',t)
  return canon(fexpr(t))
 def c(t):
  t=t.replace('Kokkos::ALL','free_axis').replace('Kokkos::subview','subview').replace('err_info.err_code','err_code')
  node=Cpp(t,arrays).parse()
  def sub(n):
   if isinstance(n,tuple) and n and n[0]=='call' and n[1]=='subview':
    v=n[2][0][1];aa=[]
    for axis,a in enumerate(n[2][1:]):
     if v=='clubb_params' and axis==1:aa.append(a)
     elif a==('v','free_axis'):aa.append(a)
     else:assert a[:2]==('bin','-') and a[3]==('n',1.0),(t,a);aa.append(a[2])
    return ('call',v,tuple(aa))
   if isinstance(n,tuple):return tuple(sub(x) if isinstance(x,tuple) else x for x in n)
   return n
  return canon(sub(node))
 def align(n,other):
  if n[0]=='v' and n[1] in arrays and other[0]=='call' and other[1]==n[1] and all(a==('v','cell') for a in other[2]):return other
  if n[0]=='call' and other[0]=='call' and n[1]==other[1] and len(n[2])==len(other[2]):return ('call',n[1],tuple(align(a,b) for a,b in zip(n[2],other[2])))
  if n[0]=='bin' and other[0]=='bin' and n[1]==other[1]:return ('bin',n[1],align(n[2],other[2]),align(n[3],other[3]))
  if n[0]=='unary' and other[0]=='unary' and n[1]==other[1]:return ('unary',n[1],align(n[2],other[2]))
  return n
 def equal(a,b):return align(a,b)==align(b,a)
 reports=[];manual=[]
 for d,p in zip(defs,O['routines']):
  assert d['args']==[a['name'] for a in p['arguments']],p['name']
  lo,hi=p['source_lines'];raw='\n'.join(S.splitlines()[lo-1:hi]);seq=statements(raw);expected=[]
  for line in seq:
   m=re.match(r'^(\w+(?:%\w+)?(?:\([^=]*\))?)\s*=(?!=)\s*(.*)$',line)
   if not m:continue
   lhs,rhs=m.groups()
   if lhs=='err_info%err_code':manual.append([p['name'],line,'Native all-column fatal boundary']);continue
   expected.append((f(lhs),f(rhs),line))
  actual=[]
  for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',d['body']):
   try:actual.append((m.start(),c(m.group(1)),c(m.group(2))))
   except Exception:pass
  actual.sort();last=-1;missing=[]
  for lhs,rhs,line in expected:
   hit=None
   for index,(pos,a,b) in enumerate(actual):
    if index<=last:continue
    if equal(lhs,a) and equal(rhs,b):hit=index;break
   if hit is None:missing.append({'source':line,'lhs':lhs,'rhs':rhs})
   else:last=hit
  predicates=[]
  for line in seq:
   m=re.match(r'(?:if|elseif)\s*\((.*)\)\s*then$',line,re.I)
   if m and not re.search(r'stats%|err_info%err_code',m.group(1),re.I):predicates.append(f(m.group(1)))
  cp=[]
  for m in re.finditer(r'\bif\s*\(',d['body']):
   e=balanced(d['body'],m.end())
   try:cp.append(c(d['body'][m.end():e-1]))
   except Exception:pass
  last=-1;pm=[]
  for e in predicates:
   hit=next((i for i,a in enumerate(cp) if i>last and equal(e,a)),None)
   if hit is None:pm.append(e)
   else:last=hit
  reports.append({'name':p['name'],'arguments':len(p['arguments']),'ordered_assignment_trees':len(expected)-len(missing),'assignment_mismatches':missing,'ordered_predicates':len(predicates)-len(pm),'predicate_mismatches':pm})


 return locals()

def complete(record):
 module=record['module'];D=record['D'];O=record['O'];S=record['S'];T=record['T'];C=record['C'];defs=record['defs'];c=record['c'];f=record['f'];equal=record['equal'];reports=record['reports'];manual=record['manual'];names=record['names'];call_reviews=[]
 for d,p in zip(defs,O['routines']):
  lo,hi=p['source_lines'];raw='\n'.join(S.splitlines()[lo-1:hi]);seq=statements(raw);expected=[];actual=[]
  for l in seq:
   m=re.match(r'call\s+(\w+)\s*\((.*)\)\s*$',l,re.I)
   if m:expected.append((m.group(1),split_top(m.group(2))))
  callnames={n.lower() for n,a in expected}
  for m in re.finditer(r'\b(?:f2k::)?(\w+)\s*\(',d['body']):
   if m.group(1).lower() not in callnames:continue
   end=balanced(d['body'],m.end());actual.append((m.group(1),split_top(d['body'][m.end():end-1])))
  assert [n.lower() for n,a in expected]==[n.lower() for n,a in actual],(module,p['name'],'call order')
  count=0
  for (en,ea),(an,aa) in zip(expected,actual):
   assert len(ea)==len(aa),(en,ea,aa)
   for e,a in zip(ea,aa):
    if e.startswith('"'):assert ast.literal_eval(e)==ast.literal_eval(a)
    else:assert equal(f(e),c(a)),(module,p['name'],en,e,a,f(e),c(a))
    count+=1
  call_reviews.append({'name':p['name'],'ordered_calls':len(expected),'actuals':count})
  for m in re.finditer(r'\b(\w+)\s*\(',d['body']):
   if m.group(1).lower() in {n.lower() for n in names}:assert m.group(1) in names,('C++ case',m.group(1))
  assert not reports[names.index(p['name'])]['assignment_mismatches']
  assert not reports[names.index(p['name'])]['predicate_mismatches']
 kernels=[]
 for m in re.finditer(r'KOKKOS_LAMBDA\([^)]*\)\s*\{',C):
  end=balanced(C,m.end(),'{','}');body=C[m.end():end-1]
  assert not re.search(r'\b(?:gr|err_info|stats)\.|f2k::|std::cerr|std::string|to_host_space|to_execution_space|parallel_for|parallel_reduce',body),body
  assert not re.search(r'\w+_host\(',body),body
  for idx in re.findall(r'for\s*\((\w+)\s*=',body):assert re.search(r'\bint\s+'+idx+r'\b',body),(module,idx)
  kernels.append(body)
 if module=='Nc_Ncn_eqns':
  assert len(kernels)==4
  assert T.count('using Policy = Kokkos::RangePolicy<typename Field0::execution_space>;')==4
  assert T.count('typename Field0::memory_space>')>=8
  assert 'Nc_in_cloud_to_Ncnm(int' not in T
  assert not re.search('f2k::|optimization_barrier|_safe|std::erfc|fma',C)
  details=['Six device scalar functions preserve all21 assignments and10 predicates; four public array generics call only their exact scalar overloads.','Explicit output View return types retain first-input memory/execution space for homogeneous rank1/rank2 conformable Views, including host/strided fields. All scalar required arguments remain unchanged.','Both private elemental integral formulas use Native ERFC/SQRT branches and tolerances, with no JAX eager where/safe denominator substitution.','Native module flowchart/scientific derivations and original call/moment ordering retained. No new public physics helper or mutable global state.']
 else:
  # Complete diagnostics: all literal fragments retain case/content. The only
  # omitted source fatal-assignment duplicates become one common final fill.
  helper=defs[-1]['body'];assert re.findall(r'failure_code\(0\) = (\d+);',helper)==[str(n) for n in range(13)]
  assert helper.count('Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);')==1
  assert helper.index('switch (failure_code_host(0))')<helper.index('Kokkos::deep_copy(err_info.err_code')
  case6=re.search(r'case 6:\s*\{(.*?)\n  case 7:',helper,re.S).group(1);assert 'err_info.err_header_global' not in case6
  assert re.search(r'if \(\(\(precip_frac\(k - 1\) - .*?> eps\)\)',helper,re.S)
  assert T.count('Independent columns retain Native physical top-to-bottom recurrence.')==4
  for n in ['precip_fraction','component_precip_frac_weighted']:
   d=next(x for x in defs if x['name']==n)
   assert d['body'].count('gr_k_ub_zt - gr_grid_dir_indx')==2
   assert d['body'].count('k += -(gr_grid_dir_indx)')==2
  assert 'f2k::Real1D(precip_frac_tol)' in T
  assert 'all_hydromet_below &&' in T and '>= hydromet_tol(species - 1)' in T and '< hydromet_tol(species - 1)' in T
  assert not re.search(r'fp==|optimization_barrier|_safe|\bfma\s*\(',C)
  details=['Both selector branches, warm/frozen options and both physical grid directions preserved. Four top-to-bottom recurrences remain serial within independent columns; physical bounds and step keep Native1based metadata.','MAXVAL profile scans use increasing storage level and negative HUGE identity; only immutable frozen-species flags use integer reductions. No reassociated floating reduction. Compiled MAXVAL tie/nonfinite behavior remains untested.','ANY(hydromet>=tol) and ALL(hydromet<tol) are separately computed, preserving empty-species and NaN distinctions; no complement/finite-value shortcut. Species boost order and all corrective formula associations remain Native.','Full12check helper runs one ordered profile work item. First failure metadata and four values cross to host; exact source diagnostics are printed, then every error column is marked fatal. Parent still calls/checks each column in Native order, including preexisting all-column fatal status.','Case6 preserves Native absence of the global header. The final consistency check remains a one-sided residual>eps; no abs or extra finite check. Diagnostic spacing is not a numerical claim.','Native sampled rank1 stats update uses the existing explicit const-element overload; no new stats/UDT owner or device bridge. Generic staging restores four caller outputs. No JAX fp==tol identity, optimization barrier or safe denominator imported.']
 files=[D/(module+ext) for ext in (['.hpp'] if module=='Nc_Ncn_eqns' else ['.cpp','.hpp'])]
 result={'canonical_source':'src/CLUBB_core/'+module+'.F90','canonical_sha256':hashlib.sha256((R/('src/CLUBB_core/'+module+'.F90')).read_bytes()).hexdigest(),'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},'routines':reports,'ordered_CALLs':call_reviews,'kernels_reviewed':len(kernels),'manual_boundaries':manual,'review_kind':'Owner complete source review, with separate cleanup; independent full-file review pending.','manual_review':details,'remaining_lower_contracts':'None: only existing constants, Kokkos scalar math, grid/types/array owner, Native stats/debug service.','runtime_status':'No builds/tests, CPU/GPU comparisons or BFB claim. All scalar branches, thresholds, selector1/2, grid directions, stats/errors and array surfaces remain runtime-unverified.'}
 (D/'SOURCE_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
 return result
if __name__=='__main__':
 result=complete(review('precipitation_fraction'))
 print('precipitation_fraction source audit passed:',len(result['routines']),'Native routines;',result['kernels_reviewed'],'device bodies')
