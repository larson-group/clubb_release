"""Review source inventory/contracts/calls and device boundaries; run no model."""
from pathlib import Path
import re,json,hashlib,sys,ast
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,find_procedure
from generate_derived_types import split_top
source=(root/'src/CLUBB_core/advance_xm_wpxp_module.F90').read_text()
target=(folder/'advance_xm_wpxp_module.cpp').read_text()
header=(folder/'advance_xm_wpxp_module.hpp').read_text()
# Remove comments before token/signature work, preserving string contents.
def clean(text):return re.sub(r'//[^\n]*','',text)
cpp=clean(target)
native_names=re.findall(r'^  subroutine (\w+)',source,re.M)
nb={n:re.search(r'^  subroutine '+n+r'\(.*?^  end subroutine '+n+r'\b',source,re.M|re.S).group() for n in native_names}
def body(text,match):
 depth,pos=1,match.end()
 while depth:depth+=(text[pos]=='{')-(text[pos]=='}');pos+=1
 return text[match.end():pos-1]
cb={};locations=[];argument_counts={}
for n in native_names:
 matches=list(re.finditer(r'\bvoid\s+'+n+r'\s*\(([^()]*)\)\s*\{',cpp))
 assert len(matches)==1,(n,len(matches))
 m=matches[0];locations.append(m.start());cb[n]=body(cpp,m)
 args=[a.split('=')[0].split()[-1].lstrip('*&') for a in split_top(m[1])]
 signature=statements(nb[n])[0];expected=split_top(signature[signature.index('(')+1:signature.rindex(')')])
 assert args==expected,(n,args,expected)
 argument_counts[n]=len(args)
 if n!='damp_coefficient':
  proc=find_procedure(n)
  for a in proc.arguments:
   declaration=next(d for d in split_top(m[1]) if re.search(r'\b'+a.name+r'\b',d))
   assert ('*' in declaration)==a.optional,(n,a.name,declaration)
assert locations==sorted(locations)
public=re.search(r'void\s+advance_xm_wpxp\s*\(([^()]*)\)',clean(header))
assert [a.split()[-1] for a in split_top(public[1])]==split_top(statements(nb['advance_xm_wpxp'])[0][statements(nb['advance_xm_wpxp'])[0].index('(')+1:-1])

all_calls=set(native_names)|{'stats_update','stats_begin_budget','stats_update_budget','stats_finalize_budget','Skx_func','calc_wp3_on_wp2','calc_turb_adv_range','diffusion_zm_lhs','diffusion_zt_lhs','term_ma_zm_lhs','term_ma_zt_lhs','xpyp_term_ta_pdf_rhs','xpyp_term_ta_pdf_lhs','xpyp_term_ta_pdf_lhs_godunov','band_solve','monotonic_turbulent_flux_limit','pos_definite_adj','fill_holes_vertical_api','clip_covar'}

def cpp_calls(text):
 result=[]
 for m in re.finditer(r'(?<!\w)(?:f2k::)?('+'|'.join(sorted(all_calls,key=len,reverse=True))+r')\s*\(',text):
  depth,pos=1,m.end()
  while depth:depth+=(text[pos]=='(')-(text[pos]==')');pos+=1
  result.append((m[1],split_top(text[m.end():pos-1])))
 return result

def native_calls(text):
 result=[];sponge=0
 for s in statements(text):
  if re.match(r'if\s*\(',s,re.I) and re.search(r'\bthen\s*$',s,re.I):
   if 'sponge_damp_settings%' in s:sponge=1;continue
   if sponge:sponge+=1
  if re.match(r'end\s*if',s,re.I) and sponge:sponge-=1;continue
  if sponge:continue
  if m:=re.match(r'call\s+(\w+)\s*\((.*)\)',s,re.I):result.append((m[1],split_top(m[2])))
 return result

# Canonicalize source arguments and scalar index conversions independently.
# Integer indices are linear expressions in native loop variables/metadata.
def balanced_at(s,start):
 depth,pos=1,start
 while depth:depth+=(s[pos]=='(')-(s[pos]==')');pos+=1
 return s[start:pos-1],pos

def polynomial(s):
 s=s.replace('gr.','gr_')
 t=ast.parse(s,mode='eval').body
 def add(a,b,scale=1):
  out=dict(a)
  for k,v in b.items():out[k]=out.get(k,0)+scale*v
  return {k:v for k,v in out.items() if v}
 def norm(n):
  if isinstance(n,ast.Constant):return {'1':n.value}
  if isinstance(n,ast.Name):return {n.id:1}
  if isinstance(n,ast.UnaryOp):return {k:(-v if isinstance(n.op,ast.USub) else v) for k,v in norm(n.operand).items()}
  if isinstance(n,ast.BinOp):
   a,b=norm(n.left),norm(n.right)
   if isinstance(n.op,ast.Add):return add(a,b)
   if isinstance(n.op,ast.Sub):return add(a,b,-1)
   if isinstance(n.op,ast.Mult):
    if set(a)<={'1'}:return {k:a.get('1',0)*v for k,v in b.items()}
    if set(b)<={'1'}:return {k:b.get('1',0)*v for k,v in a.items()}
  raise ValueError(s)
 return repr(sorted(norm(t).items()))


def argument(a,is_cpp=False):
 a=a.strip().replace('%','.')
 if not is_cpp:
  while m:=re.search(r'\*\*\s*2\b',a):
   end=m.start();start=end-1
   while a[start].isspace():start-=1
   if a[start]==')':
    depth=1;start-=1
    while depth:depth+=(a[start]==')')-(a[start]=='(');start-=1
   while start>=0 and (a[start].isalnum() or a[start] in '_.'):start-=1
   start+=1;base=a[start:end]
   a=a[:start]+'('+base+'*'+base+')'+a[m.end():]
 a=re.sub(r'_core_rknd\b','',a)
 a=re.sub(r'\breal_t\(([^()]*)\)',r'(\1)',a)
 a=a.replace('Kokkos::','')
 a=a.replace('to_host_space','HOST')
 if is_cpp:
  a=a.lstrip('&')
  # Host scalar reads and View staging have no numerical transformation.
  for function in ['HOST','to_execution_space']:
   while m:=re.search(r'\b'+function+r'\(([^()]*)\)',a):a=a[:m.start()]+m[1]+a[m.end():]
  if a=='lhs_sequence':a='subview(lhs_ta_wpsclrp,ALL,ALL,ALL,(sclr-1))'
  while m:=re.search(r'\bsubview\(',a):
   raw,end=balanced_at(a,m.end());items=split_top(raw);base=items[0];indices=[]
   for value in items[1:]:
    value=value.strip()
    if value=='ALL':indices.append(':')
    elif value.startswith('std::make_pair('):
     raw=split_top(value[value.index('(')+1:-1]);indices.append(polynomial(raw[0]+'+1')+':'+polynomial(raw[1]))
    elif value.strip('() ').startswith(('iC','ic_','ialtitude','iwpxp')):indices.append(value.strip('() '))
    else:indices.append(polynomial('('+value+')+1'))
   value=base+'['+','.join(indices)+']';a=a[:m.start()]+value+a[end:]
  # The only scalar-array arguments are the column fcor diagnostic and
  # passive-scalar thresholds. Keep native integer index meaning.
  for base in ['sclr_tol','fcor']:
   pattern=r'\b'+base+r'\('
   while m:=re.search(pattern,a):
    raw,end=balanced_at(a,m.end());a=a[:m.start()]+base+'['+polynomial('('+raw+')+1')+']'+a[end:]
 else:
  # Native explicit sections retain their argument rank and lower bounds.
  for m in reversed(list(re.finditer(r'\b([\w.]+)\(([^()]*)\)',a))):
   base=m[1];items=split_top(m[2])
   if base not in {'sclr_tol','fcor'} and not any(':' in i for i in items):continue
   indices=[]
   for value in items:
    if value.strip()==':':indices.append(':')
    elif ':' in value:
     lo,hi=value.split(':');indices.append(polynomial(lo)+':'+polynomial(hi))
    elif value.startswith(('iC','ic_','ialtitude','iwpxp')):indices.append(value)
    else:indices.append(polynomial(value))
   a=a[:m.start()]+base+'['+','.join(indices)+']'+a[m.end():]
 # Source powers in call thresholds are all integer squares. Repeated-operand
 # multiplication is the ordinary explicit target lowering.
 return re.sub(r'[\s()]','',a).lower()

call_counts={};argument_checks=0;stats_counts={};known_differences=[]
for n in native_names:
 actual=cpp_calls(cb[n]);expected=native_calls(nb[n])
 assert [x[0] for x in actual]==[x[0] for x in expected],(n,[x[0] for x in actual],[x[0] for x in expected])
 call_counts[n]=len(actual)
 stats_counts[n]=sum(name.startswith('stats_') for name,_ in actual)
 # Resolve explicit staging of array-valued stats expressions.
 stats_values={}
 for m in re.finditer(r'\b(stats_values_\d+)\([^;=]*\)\s*=\s*([^;]+);',cb[n]):
  value=m[2]
  value=re.sub(r'\b(\w+)\(i(?:,\s*k)?\)',r'\1',value)
  stats_values[m[1]]=value
 for (name,aa),(same,ee) in zip(actual,expected):
  if name=='band_solve':
   known_differences.append('xm_wpxp_solve: explicit Native-bands/C++-bands and optional-pointer API boundary')
   assert len(aa)==14
   continue
  if name in native_names:
   params=split_top(statements(nb[name])[0][statements(nb[name])[0].index('(')+1:-1]);mapped={};pos=0
   for val in ee:
    if m:=re.match(r'^(\w+)\s*=(?!=)\s*(.*)',val):mapped[m[1]]=m[2]
    else:mapped[params[pos]]=val;pos+=1
   ee=[mapped.get(param,'nullptr') for param in params]
  assert len(aa)==len(ee),(n,name,len(aa),len(ee))
  for index,(a,e) in enumerate(zip(aa,ee)):
   if name.startswith('stats_') and index==1:
    a=re.sub(r'^(?:to_host_space|f2k::Real1D)\((stats_values_\d+)\)$',r'\1',a)
    a=stats_values.get(a,a)
   assert argument(a,True)==argument(e), (n,name,index,a,e,argument(a,True),argument(e))
   argument_checks+=1

diagnostic_writes={}
for name in native_names:
 if name=='xm_wpxp_solve':continue # Its two source writes share one explicit host stream block.
 expected=[]
 for statement in statements(nb[name]):
  if match:=re.search(r'\bwrite\(fstderr,\s*\*\)\s*(.*)',statement,re.I):
   payload=re.sub(r"new_line\('c'\)",'',match[1])
   expected.append(re.findall(r'"[^"]*"',payload))
 actual=[]
 for match in re.finditer(r'message\s*<<\s*(.*?);',cb[name],re.S):
  if match[1].startswith('std::setprecision'):continue
  actual.append([value for value in re.findall(r'"(?:\\.|[^"\\])*"',match[1])
                 if value not in ['" "','"\\n"']])
 assert actual==expected,(name,len(actual),len(expected))
 diagnostic_writes[name]=len(expected)

# No owning UDT is copied into a kernel: only its component Views/scalars.
kernels=0
for m in re.finditer(r'KOKKOS_LAMBDA\([^{}]*?\)\s*\{',cpp):
 b=body(cpp,m);kernels+=1
 assert not re.search(r'\b(?:gr|err_info|pdf_implicit_coefs_terms|nu_vert_res_dep)\.',b)
 assert 'f2k::' not in b
 assert 'to_host_space(' not in b
 assert 'std::cerr' not in b
 assert 'Kokkos::parallel_for(' not in b

for filename in ['advance_xm_wpxp_module.cpp','advance_xm_wpxp_module.hpp','f2k_dependencies.hpp']:
 text=(folder/filename).read_text();c=clean(text);c=re.sub(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'','',c)
 for left,right in [('(',')'),('[',']'),('{','}')]:
  depth=0
  for char in c:depth+=(char==left)-(char==right);assert depth>=0,(filename,left)
  assert depth==0,(filename,left,depth)
 assert all(line.rstrip()==line and '\t' not in line for line in text.splitlines())
 assert all(len(line)<=100 for line in text.splitlines())
assert 'Kokkos::deep_copy(xm_tndcy_wpxp_cl, zero)' in target
assert 'l_clipping_needed(column) = needed' in target
assert 'Kokkos::MemoryTraits<Kokkos::Unmanaged>' in target
assert 'Native BFB claim' in target
assert 'fma(' not in cpp
assert 'host(j, column)' in target
assert 'j < host.extent(0)' in target
assert 'f2k::Real1D(stats_values_' in target
record={
 'status':'Full supported source draft; not compiled or runtime validated',
 'source':'src/CLUBB_core/advance_xm_wpxp_module.F90; all 6472 lines read',
 'routine_order':native_names,'argument_counts':argument_counts,
 'all_signature_names_order_optionals':True,'source_order_call_counts':call_counts,
 'all_non_solver_call_arguments_checked':argument_checks,
 'stats_calls_by_routine':stats_counts,'device_kernel_bodies_checked':kernels,
 'native_diagnostic_write_labels_in_order':diagnostic_writes,
 'character_diagnostics':'Existing owner shape(character,column); helper iterates extent0',
 'no_udt_capture_host_access_or_nested_launch_in_kernel':True,
 'supported_scope_exclusion':'Three sponge boundaries and eight associated stats calls; require upstream rejection matching JAX',
 'source_definedness_choice':'Clipping correction initializes per-column masks false and nonclipping-column tacl stats zero, matching JAX',
 'native_contract_edges':['Explicit scalar TA source sequence-associates rank3 LHS into rank2 RHS output and leaves scalar RHS unwritten','xm_wpxp_solve.old_solution declaration shape disagrees with callers; only BiCGStab indexes it and existing lower rejects that JAX-excluded method','Positive-definite stats temporaries are source-unwritten when no adjustment occurs'],
 'matrix_adapter':known_differences,
 'manual_source_audit':['All Native routine/control/call/argument/comment/error/stat sequences','Source pressure, TA, diffusion, grid orientation, scalar, Coriolis, perturbation, LMM and clipping branches','Private coefficient and diagnosis formulas with Native grouping and boundary assignments','Array slices expanded with explicit index conversion and no mutable global state','Full Native diagnostic payload/labels, array order, max_digits10 and host-only I/O','Current source owners/lower contracts inspected; no shadow APIs/constants implemented'],
 'runtime':'No builds, installs, compilation, unit tests or model runs; root owns Serial CPU BFB gate and occasional GPU checks',
 'sha256':{name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ['advance_xm_wpxp_module.cpp','advance_xm_wpxp_module.hpp','f2k_dependencies.hpp']},
 'native_sha256':hashlib.sha256(source.encode()).hexdigest()
}
(folder/'source_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ['status','argument_counts','source_order_call_counts','all_non_solver_call_arguments_checked','stats_calls_by_routine','device_kernel_bodies_checked']},indent=2))
