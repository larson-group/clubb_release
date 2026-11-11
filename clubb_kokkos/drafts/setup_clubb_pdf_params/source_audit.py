"""Owner source audit of extended PDF contracts/formulas; no model execution."""
from pathlib import Path
import ast,hashlib,json,re,sys
H=Path(__file__).resolve().parent;R=H.parents[2];sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,declaration
from generate_derived_types import split_top
O=json.loads((H/'SOURCE_OUTLINE.json').read_text());source=(R/'src/CLUBB_core/setup_clubb_pdf_params.F90').read_text();text=(H/'setup_clubb_pdf_params.cpp').read_text()
# Protect strings while dropping comments: URLs and string slashes are retained.
def strip_comments(t):
 return re.sub(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*.*?\*/',lambda m:m.group() if m.group()[0] in '\"\'' else '',t,flags=re.S)
code=strip_comments(text)
q=(R/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();exec(q[q.index('class Cpp:'):q.index('reviews=[]')])
def balanced(t,start,opening='(',closing=')'):
 p=start;depth=1;quote=''
 while depth:
  c=t[p]
  if quote:
   if c==quote and t[p-1]!='\\':quote=''
  elif c in '\"\'':quote=c
  else:depth+=(c==opening)-(c==closing)
  p+=1
 return p

defs=[]
for m in re.finditer(r'\b(?:void|real_t)\s+(\w+)\s*\(',code):
 if m.group(1) not in O['native_inventory_and_order']:continue
 end=balanced(code,m.end());tail=code[end:].lstrip()
 if not tail.startswith('{'):continue
 start=code.index('{',end)+1;stop=balanced(code,start,'{','}')
 defs.append({'name':m.group(1),'args':[a.strip().split()[-1] for a in split_top(code[m.end():end-1])],'body':code[start:stop-1],'at':m.start()})
# Thin template adapters precede Native definitions; select concrete definitions.
core=defs[-18:];assert [d['name'] for d in core]==O['native_inventory_and_order']
cpp_arrays=set();physics=set();reviews=[]
for p in O['routines']:
 for a in p['arguments']:
  if a['dims'] and a['kind']=='real':cpp_arrays.add(a['name'].lower())
for m in re.finditer(r'\bReal([1-4])D\s+(\w+)\s*\(',code):cpp_arrays.add(m.group(2).lower())
for m in re.finditer(r'const auto\s+(\w+)\s*=\s*(?:pdf_params|precip_fracs|hydromet_pdf_params)\.',code):cpp_arrays.add(m.group(1).lower())
for x in ['hm_metadata_hydromet_tol','hm_metadata_hmp2_ip_on_hmm2_ip','err_info_err_code']:cpp_arrays.add(x)
cpp_arrays|={n+'_host' for n in cpp_arrays}
# Metadata and prescribed correlation matrices use physical PDF indices rather
# than column/level axes. All other field indices remain in ASTs, with only
# Native array-section free axes abstracted to the same independent cell.
def canon(n):
 if not isinstance(n,tuple) or not n:return n
 if n==('unary','-',('v','default_real_positive_huge')):return ('v','default_real_negative_huge')
 if n[0]=='v':
  name=n[1].replace('_host','')
  if name in {'i','k'} or name.startswith('slice_') or name=='free_axis':return ('v','cell')
  if name=='default_real_negative_huge':return ('v','default_real_negative_huge')
  return ('v',name)
 if n[0]=='call':
  name=n[1].replace('_host','');args=tuple(canon(x) for x in n[2])
  if name=='huge':return ('v','huge')
  if name=='hydromet2pdf_idx':return ('v','jpdf_idx')
  if name=='nc_in_cloud_to_ncnm' and len(args)==12:args=args[2:]
  if name=='zt2zm_api' or name=='zm2zt_api':return ('call',name,args)
  return ('call',name,args)
 return tuple(canon(x) if isinstance(x,tuple) else x for x in n)
def f(s):
 s=re.sub(r'hydromet_pdf_params\(:,:\)%(\w+)\(([^)]*)\)',r'hydromet_pdf_params_\1(free_axis,free_axis,\2)',s)
 s=re.sub(r'hydromet_pdf_params\(:,:\)%(\w+)',r'hydromet_pdf_params_\1(free_axis,free_axis)',s)
 s=s.replace('real( -huge( 0.0 ), kind = core_rknd )','default_real_negative_huge')
 s=re.sub(r'(?<!\w):', 'free_axis',s)
 return canon(fexpr(s))
def c(s):
 s=s.replace('std::numeric_limits<real_t>::max()','huge()').replace('real_t(std::numeric_limits<float>::max())','default_real_positive_huge')
 s=re.sub(r'Kokkos::subview\s*\(', 'subview(',s)
 s=s.replace('Kokkos::ALL','free_axis')
 n=Cpp(s,cpp_arrays).parse()
 def sub(t):
  if isinstance(t,tuple) and t and t[0]=='call' and t[1]=='subview':
   name=t[2][0][1];args=[]
   for arg in t[2][1:]:
    if arg==('v','free_axis'):args.append(arg)
    else:
     assert arg[:2]==('bin','-') and arg[3]==('n',1.0),(s,arg)
     args.append(arg[2])
   return ('call',name,tuple(sub(a) for a in args))
  if isinstance(t,tuple):return tuple(sub(a) if isinstance(a,tuple) else a for a in t)
  return t
 n=sub(n)
 if n==('unary','-',('v','default_real_positive_huge')):n=('v','default_real_negative_huge')
 return canon(n)
# The source AST reader collapses whole rank1/2 array sections to a name.
# Restore those array-valued leaves to explicit free axes in each assignment.
def align_whole(n,other):
 if n[0]=='v' and n[1] in cpp_arrays and other[0]=='call' and other[1]==n[1]:
  if all(x==('v','cell') for x in other[2]):return other
 if n[0]=='call' and other[0]=='call' and n[1]==other[1] and len(n[2])==len(other[2]):return ('call',n[1],tuple(align_whole(a,b) for a,b in zip(n[2],other[2])))
 if n[0]=='bin' and other[0]=='bin' and n[1]==other[1]:return ('bin',n[1],align_whole(n[2],other[2]),align_whole(n[3],other[3]))
 if n[0]=='unary' and other[0]=='unary' and n[1]==other[1]:return ('unary',n[1],align_whole(n[2],other[2]))
 return n
manual=[]
for d,p in zip(core,O['routines']):
 assert d['args']==[a['name'] for a in p['arguments']],p['name']
 lo,hi=p['source_lines'];raw='\n'.join(source.splitlines()[lo-1:hi]);expected=[]
 for l in statements(raw):
  m=re.match(r'^(\w+(?:\([^=]*\))?(?:%\w+(?:\([^=]*\))?)?)\s*=(?!=)\s*(.*)$',l)
  if not m:continue
  lhs,rhs=m.groups()
  if '"' in rhs or "'" in rhs or 'hydromet_list' in rhs or 'parameter' in lhs:continue
  if lhs.lower()=='err_info%err_code(i)':manual.append({'routine':p['name'],'assignment':l,'reason':'first-fatal device boundary manually reviewed'});continue
  if 'compute_rtp2_from_chi' in rhs and ':' in rhs:manual.append({'routine':p['name'],'assignment':l,'reason':'post-DO j scalarization manually reviewed'});continue
  try:expected.append((f(lhs),f(rhs),l))
  except Exception as e:raise AssertionError((p['name'],l,str(e)))
 actual=[]
 for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',d['body']):
  try:actual.append((m.start(),c(m.group(1)),c(m.group(2)),m.group()))
  except Exception:pass
 for m in re.finditer(r'Kokkos::deep_copy\s*\(',d['body']):
  end=balanced(d['body'],m.end());args=split_top(d['body'][m.end():end-1])
  try:actual.append((m.start(),c(args[0]),c(args[1]),'deep_copy'))
  except Exception:pass
 actual=[(a,b,t) for pos,a,b,t in sorted(actual)]
 used=set();last=-1;matched=0;missing=[]
 for lhs,rhs,line in expected:
  found=None
  for i,(al,ar,t) in enumerate(actual):
   if i in used or i<=last:continue
   if align_whole(lhs,al)==align_whole(al,lhs) and align_whole(rhs,ar)==align_whole(ar,rhs):found=i;break
  if found is None:missing.append({'source':line,'lhs':lhs,'rhs':rhs})
  else:used.add(found);last=found;matched+=1
 reviews.append({'name':p['name'],'arguments':len(p['arguments']),'assignment_trees_matched':matched,'assignment_trees_unmatched':missing})
print(json.dumps(reviews,indent=2))

# Calls: compare each complete field-specific actual argument in source order.
call_reviews=[]
for d,p in zip(core,O['routines']):
 lo,hi=p['source_lines'];raw='\n'.join(source.splitlines()[lo-1:hi]);expected=[];actual=[]
 for line in statements(raw):
  if re.match(r'call\s+',line,re.I):
   m=re.match(r'call\s+(\w+)\s*\((.*)\)\s*$',line,re.I);assert m,line
   expected.append((m.group(1),split_top(m.group(2))))
 names={n.lower() for n,args in expected}
 for m in re.finditer(r'\b(?:f2k::)?(\w+)\s*\(',d['body']):
  if m.group(1).lower() not in names:continue
  end=balanced(d['body'],m.end());actual.append((m.group(1),split_top(d['body'][m.end():end-1])))
 assert [n.lower() for n,a in expected]==[n.lower() for n,a in actual],(p['name'],'call order',expected,actual)
 count=0
 for (en,ea),(an,aa) in zip(expected,actual):
  assert len(ea)==len(aa),(p['name'],en,'arity',ea,aa)
  for index,(e,a) in enumerate(zip(ea,aa)):
   if en.lower()=='pack_hydromet_pdf_params' and index==len(ea)-1:
    assert e.replace(' ','')=='hydromet_pdf_params(:,:)' and a=='*hydromet_pdf_params';count+=1;continue
   fe=f(e);ca=c(a)
   # Case preserved independently for string literals; Fortran AST lexer folds.
   if fe[0]=='str':assert ast.literal_eval(e)==ast.literal_eval(a);count+=1;continue
   assert align_whole(fe,ca)==align_whole(ca,fe),(p['name'],en,index,e,a,fe,ca)
   count+=1
 call_reviews.append({'name':p['name'],'ordered_calls':len(expected),'actual_arguments':count})
# Branch source ASTs are retained for cell predicates. Diagnostic ANY, present,
# stats flags and C++ language boundaries are explicitly reviewed below.
branch_reviews=[]
for d,p in zip(core,O['routines']):
 lo,hi=p['source_lines'];raw='\n'.join(source.splitlines()[lo-1:hi]);expected=[]
 for l in statements(raw):
  m=re.match(r'(?:if|elseif)\s*\((.*)\)\s*then\s*$',l,re.I) or re.match(r'where\s*\((.*)\)\s*$',l,re.I)
  if not m:continue
  expr=m.group(1)
  if re.search(r'any\s*\(|present\s*\(|stats%|err_info%err_code|var_on_stats_list',expr,re.I):continue
  if p['name']=='setup_pdf_parameters_api' and 'hydromet(i,k,j)' in expr:continue # first-fatal integer flag boundary above
  expected.append(f(expr))
 actual=[]
 for m in re.finditer(r'\bif\s*\(',d['body']):
  end=balanced(d['body'],m.end())
  try:actual.append(c(d['body'][m.end():end-1]))
  except Exception:continue
 remaining=list(actual)
 for e in expected:
  hit=None
  for i,a in enumerate(remaining):
   if align_whole(e,a)==align_whole(a,e):hit=i;break
  assert hit is not None,(p['name'],'predicate',e)
  remaining=remaining[hit+1:]
 branch_reviews.append({'name':p['name'],'ordered_scientific_predicates':len(expected)})
# Device captures and ownership: every extended lambda must contain only
# component aliases and Views/scalars, never host callbacks or owning UDTs.
kernels=[]
for m in re.finditer(r'KOKKOS_LAMBDA\([^)]*\)\s*\{',code):
 end=balanced(code,m.end(),'{','}');body=code[m.end():end-1]
 assert not re.search(r'\b(?:gr|pdf_params|hm_metadata|err_info|precip_fracs|hydromet_pdf_params)\.',body),body
 assert not re.search(r'f2k::|to_host_space|to_execution_space|std::cerr|std::string|parallel_for|parallel_reduce|hydromet2pdf_idx|pdf2hydromet_idx',body),body
 assert not re.search(r'\w+_host\(',body),body
 # Loop indices must be declared inside the lambda or deliberately captured
 # from host; any mutable for-index is necessarily local and cannot be const.
 for idx in re.findall(r'\bfor\s*\(\s*(\w+)\s*=',body):assert re.search(r'\bint\s+'+idx+r'\b',body),(idx,body)
 kernels.append(body)
assert not re.search(r'\bfma\s*\(|\bthread_local\b',code)
# True Native public inventory; adapters repeat names without adding options.
header=(H/'setup_clubb_pdf_params.hpp').read_text();declarations=re.findall(r'(?m)^void\s+(\w+)\s*\(',header)
assert declarations==[n for n in O['native_inventory_and_order'] if n in O['public_inventory']]
for d in core:
 calls=[m.group(1) for m in re.finditer(r'\b(\w+)\(',d['body']) if m.group(1).lower() in {n.lower() for n in O['native_inventory_and_order']}]
 for n in calls:assert n in O['native_inventory_and_order'],('C++ case',n)
assert 'pdf_params_rt_2(rt_2_column - 1, k - 1)' in text
assert 'const int rt_2_column = hydromet_dim + 1;' in text
assert re.search(r'compute_mean_stdev\(.*?sigma_w_1,\s*sigma_w_1,',core[0]['body'],re.S)
for r in reviews:assert not r['assignment_trees_unmatched'],r
record={'canonical_source':'src/CLUBB_core/setup_clubb_pdf_params.F90','canonical_sha256':hashlib.sha256((R/'src/CLUBB_core/setup_clubb_pdf_params.F90').read_bytes()).hexdigest(),'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [H/'setup_clubb_pdf_params.cpp',H/'setup_clubb_pdf_params.hpp',H/'f2k_dependencies.hpp']},'routines':reviews,'ordered_calls':call_reviews,'ordered_predicates':branch_reviews,'kernels_reviewed':len(kernels),'manual_boundaries':manual,'review_kind':'Owner source audit plus separate cleanup; no independent full-file reviewer yet.','manual_review':['Full4,385line Native and2,303line actualJAX file read before recorded18routine outline. Native is authority at every noted difference.','Parent negative-check reduction is integer-only; host prints first j/k/i violation and marks only that column before immediate return. Precipitation debug0 return, debug2 symmetric assertions and final fatal context retain Native order.','All field-specific stats calls preserve names, argument order, guards and dynamic first-two-character TRIM naming. Normal mean stats use negative default-real HUGE widened to real_t; actual normal means retain negative real_t HUGE.','Optional hydro collection uses existing owner; no duplicate UDT/global mutable storage. Native pack updates active fields only. Lower/diagonal-only true correlation writes and undefined Native upper/optional-output conjunctions are documented rather than silently corrected.','Source scientific scalar/index expression trees preserve grouping and Native literal kinds; powers squared use direct multiplication, no FMA workaround or JAX safe denominators/gradient floors.','All six public and twelve private definitions remain in source order; array-section overloads are thin shared-owner adapters. Scalar elemental helpers are device callable. Lower matrix diagnostics remain host APIs, not callbacks inside kernels.','Tunable indices use existing zero-based parameter owner; all grid/hydrometeor/PDF metadata keep Native1based values. Column-first fields and strided matrix section pack/restore preserve existing layouts.','Post-DO j diagnostic rt_2 actual and duplicated sigma_w_1 source call explicitly retained; exact diagnostics-only Native/JAX differences recorded. No Native scientific fix.'],'pending_integrations':['shared Real4D and rank4 array_types owner extension','two declared transient Native array bridges (precip_fraction22args; Nc_in_cloud_to_Ncnm10inputs plus2shape args)','stable matrix/diagnose/corrmetadata/indexmapping/clip owners as detailed in INTEGRATION'],'runtime_status':'No builds/tests/CPU/GPU comparisons run. No BFB claim; all compiled/numerical behavior and option/error/stats paths remain pending reached-model integration.'}
(H/'SOURCE_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'routines':len(core),'assignments':sum(x['assignment_trees_matched'] for x in reviews),'calls':sum(x['ordered_calls'] for x in call_reviews),'actuals':sum(x['actual_arguments'] for x in call_reviews),'predicates':sum(x['ordered_scientific_predicates'] for x in branch_reviews),'kernels':len(kernels)}))
# Keep the stats registration guards' case and order as well as update calls.
for d,p in zip(core,O['routines']):
 lo,hi=p['source_lines'];raw='\n'.join(source.splitlines()[lo-1:hi]);expected=[];actual=[]
 for t,destination in [(raw,expected),(d['body'],actual)]:
  for m in re.finditer(r'(?:f2k::)?var_on_stats_list\s*\(',t,re.I):
   end=balanced(t,m.end());args=split_top(t[m.end():end-1]);destination.append([re.sub(r'\s+','',a) for a in args])
 assert expected==actual,(p['name'],'case-sensitive stats guards',expected,actual)
# Universal adapters may introduce no new required/optional physics argument.
# Every source Out array avoids reading unspecified caller staging values.
for p in O['routines']:
 arrays=[a for a in p['arguments'] if a['kind']=='real' and a['dims']]
 if not arrays:continue
 t=header if p['name'] in O['public_inventory'] else text
 m=re.search(r'template\s*<[^;{}]*>\s*inline\s+void\s+'+p['name']+r'\s*\([^{}]*\)\s*\{',t,re.S);assert m,p['name']
 end=balanced(t,m.end(),'{','}');body=t[m.end():end-1]
 for a in arrays:
  assert 'to_execution_space('+a['name']+(', false' if a['intent']=='out' else '')+')' in body,(p['name'],a['name'],a['intent'])
  if a['intent'] in {'out','inout'}:assert 'from_execution_space('+a['name']+', '+a['name']+'_exec);' in body

record['case_sensitive_stats_registration_guards']='All source guards, names and order checked.'
record['adapter_source_contracts']='All16 shared-owner section adapters keep original arguments, restore every Out/InOut array and avoid reading caller Out values during staging.'
(H/'SOURCE_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n')
