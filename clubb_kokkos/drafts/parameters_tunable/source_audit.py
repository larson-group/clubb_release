"""Audit the complete supported parameter source group without executing model code."""
from pathlib import Path
import ast,hashlib,json,re,sys,textwrap,dataclasses
D=Path(__file__).resolve().parent;R=D.parents[2];N=R/'src/CLUBB_core/parameters_tunable.F90';S=N.read_text();sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
import fortran_procedures as fp,generate_procedure_bridges as gp
outline=json.loads((D/'SOURCE_OUTLINE.json').read_text());names=outline['native_inventory_and_order'];contract=json.loads((D/'parameter_contracts.json').read_text())
raw={n:re.search(r'(?ims)^  subroutine '+n+r'\b.*?^  end subroutine '+n+r'\b',S)[0] for n in names}
T=(D/'procedure_definitions.hpp.inc').read_text();C=re.sub(r'//[^\n]*','',T)
# Reuse the already source-audited expression-reader logic, not the helper module audit itself.
helper=(D.parent/'advance_helper_module/source_audit.py').read_text();ex=helper[helper.index('class Cpp:'):helper.index('\nreviews=[]')];exec(ex,globals())
def definitions(c):
 out=[]
 for m in re.finditer(r'void\s+(\w+)\s*\(([^{};]*?)\)\s*\{',c,re.S):
  end=m.end();depth=1
  while depth:depth+=(c[end]=='{')-(c[end]=='}');end+=1
  out.append({'name':m[1],'arguments':[a.strip().split()[-1] for a in split_top(re.sub(r'Kokkos::View<.*?>','HostMatrix',m[2]))],'body':c[m.end():end-1]})
 return out
Dfs=definitions(C);assert [d['name'] for d in Dfs]==names;body={d['name']:d['body'] for d in Dfs}
initializer=definitions(re.sub(r'//[^\n]*','',(D/'init_clubb_params_api.cpp.inc').read_text()));assert len(initializer)==1 and initializer[0]['name']=='init_clubb_params_api'
assert initializer[0]['arguments']==[a.name for a in read_procedure(N,'init_clubb_params_api').arguments]
body['init_clubb_params_api']=initializer[0]['body']
assert '#include \"f2k/driver_io.hpp\"' not in (D/'procedure_includes.hpp.inc').read_text()
for d in Dfs:assert d['arguments']==[a.name for a in read_procedure(N,d['name']).arguments],d['name']
# All 105 defaults retain Native order/literal values and all 210 pack/unpack mappings are exact.
reviews=[]
for n in ['set_default_parameters','pack_parameters','unpack_parameters']:
 native=[]
 for l in statements(raw[n]):
  m=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',l)
  if not m:continue
  lhs,rhs=m.groups()
  if n=='set_default_parameters':lhs+='(i)'
  native.append((fexpr(lhs),fexpr(rhs)))
 actual=[]
 for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',body[n]):
  lhs,rhs=m.groups();rhs=re.sub(r'\bparameters\b','clubb_params',rhs)
  try:actual.append((Cpp(lhs,set()).parse(),Cpp(rhs,set()).parse()))
  except Exception:raise ValueError((n,lhs,rhs))
 assert len(native)==105 and native==actual,(n,native,actual)
 reviews.append({'name':n,'arguments':len(read_procedure(N,n).arguments),'ordered_assignment_expression_trees':105})
# Source-derived immutable catalog and exact Fortran-flat bound table.
cat=(D/'constants.hpp.inc').read_text();catcode=re.sub(r'//[^\n]*','',cat)
strings=re.findall(r'"([^"]*)"',catcode);assert strings==contract['catalog_names_32'] and all(len(s)==32 for s in strings)
original=S[S.index('  params_list ='):S.index('  ! Hard configuration bounds')];assert [s.ljust(32) for s in re.findall(r'"([^"]*)"',original)]==strings
bounds=re.search(r'parameter_hard_bounds\s*=\s*\{(.*?)\};',catcode,re.S)[1];values=split_top(bounds.rstrip().rstrip(','));assert len(values)==210
for f,c in zip(contract['bounds_Fortran_flat'],values):
 f=f.replace('epsilon(one)','machine_epsilon');c=c.replace('std::numeric_limits<real_t>::epsilon()','machine_epsilon').replace('parameters_tunable_detail::','')
 assert fexpr(f)==Cpp(c,set()).parse(),(f,c)
# Derived scientific assignment ASTs preserve associations after rank/index syntax mapping.
n='calc_derived_params';native=[]
def native_normal(s):
 s=s.replace('clubb_params(:,','clubb_params(i,')
 s=re.sub(r'gr%([zk][mt]|k_(?:ub|lb)_z[mt])',r'\1',s)
 s=s.replace('gr%zm','zm').replace('gr%zt','zt').replace('gr%nzm','nzm').replace('gr%nzt','nzt')
 s=re.sub(r'real\(\s*([^,()]+),\s*kind\s*=\s*core_rknd\s*\)',r'(\1)',s)
 s=re.sub(r'\b(avg_deltaz|mult_factor_zt|mult_factor_zm)\(i\)',r'\1',s)
 return s
for l in statements(raw[n]):
 m=re.match(r'^(lmin|mixt_frac_max_mag|avg_deltaz\(i\)|mult_factor_z[mt]\(i\)|nu_vert_res_dep%\w+\(i\))\s*=\s*(.*)',l)
 if not m:continue
 lhs,rhs=m.groups()
 if lhs in {'lmin','mixt_frac_max_mag'}:lhs+='(i)'
 native.append((fexpr(native_normal(lhs)),fexpr(native_normal(rhs))))
actual=[]
for m in re.finditer(r'(?m)^\s*(minimum_length\(i\)|mixture_limit\(i\)|avg_deltaz|mult_factor_z[mt]|nu(?:1|2|6|8|9|10|_hm)\(i\))\s*=\s*([^;]*);',body[n]):
 lhs,rhs=m.groups()
 for a,b in [('minimum_length','lmin'),('mixture_limit','mixt_frac_max_mag'),('parameters','clubb_params'),('spacing','deltaz')]:lhs=re.sub(r'\b'+a+r'\b',b,lhs);rhs=re.sub(r'\b'+a+r'\b',b,rhs)
 lhs=re.sub(r'\b(nu(?:1|2|6|8|9|10|_hm))\(',r'nu_vert_res_dep_\1(',lhs)
 actual.append((Cpp(lhs,set()).parse(),Cpp(rhs,set()).parse()))
assert native==actual,(n,native,actual)
assert len(native)==24
b=body[n];assert b.index('minimum_length(i) =')<b.index('throw std::runtime_error')<b.index('mixture_limit(i) =')
assert [m[1] for m in re.finditer(r'nu_vert_res_dep\.(\w+)\s*=\s*Real1D',b)]==['nu1','nu2','nu6','nu8','nu9','nu10','nu_hm']
assert b.count('if (!nu_vert_res_dep.nu1.is_allocated())')==1
assert ['l_prescribed_avg_deltaz','grid_type == 3','grid_type == 1','grid_type == 2']==[m[1] for m in re.finditer(r'(?:else\s+)?if\s*\((l_prescribed_avg_deltaz|grid_type == [123])\)',b)]
assert all(f'const int {v} = gr.{v} - 1;' in b for v in ['k_ub_zm','k_lb_zm','k_ub_zt','k_lb_zt'])
reviews.append({'name':n,'arguments':9,'ordered_science_assignment_expression_trees':len(native),'seven_allocations_only_if_nu1_absent':True,'lmin_write_fatal_before_further_math':True,'invalid_grid_type_Native_undefinedness_retained':True})
# The four paired checks retain exact source algebra and field order.
b=body['check_parameters_api'];pairs=[('C6rt','C6thl'),('C6rtb','C6thlb'),('C6rtc','C6thlc'),('C6rt_Lscale0','C6thl_Lscale0')]
for a,z in pairs:
 expected=f'abs(clubb_params(i,i{a})-clubb_params(i,i{z})) > abs(clubb_params(i,i{a})+clubb_params(i,i{z}))/2*eps'
 cm=re.search(r'if\s*\((std::abs\(parameters\(i, i'+a+r'\).*?)\)\s*{',b,re.S)[1]
 assert fexpr(expected)==Cpp(cm.replace('std::abs','abs').replace('parameters','clubb_params'),set()).parse()
assert b.count('err_code(i) = clubb_fatal_error;')==4 and b.index('check_parameter_hard_bounds')<b.index('for (int i')
assert 'return;' not in b
b=body['check_parameter_hard_bounds'];assert b.index('for (int i')<b.index('for (int k')
assert 'has_lower && parameters(i, k) < parameter_hard_bounds[2 * k]' in b
assert 'has_upper && parameters(i, k) > parameter_hard_bounds[2 * k + 1]' in b
assert 'has_lower && has_upper' in b and b.count('err_code(i) = clubb_fatal_error;')==1
assert 'from_execution_space(err_info.err_code, err_code);' in b
for literal in [' is outside hard configuration bounds: ',' violates lower hard configuration bound: ',' violates upper hard configuration bound: ']:assert json.dumps(literal) in b
assert 'nparams_in != nparams' in body['get_parameter_hard_bounds_api']
assert 'get_parameter_hard_bounds_api: nparams mismatch' in body['get_parameter_hard_bounds_api']
# Exact default/read/pack call arguments and order; no call to a private Native helper.
b=body['init_clubb_params_api'];assert b.index('set_default_parameters(')<b.index('f2k::read_clubb_params_nl(')<b.index('parameters_tunable_detail::pack_parameters(')
def calls(s,name):
 out=[]
 for m in re.finditer(r'\b'+re.escape(name)+r'\s*\(',s):
  e=m.end();depth=1
  while depth:depth+=(s[e]=='(')-(s[e]==')');e+=1
  out.append(split_top(s[m.end():e-1]))
 return out
assert calls(b,'set_default_parameters')==[['ngrdcol']+contract['default_argument_order']]
assert calls(b,'pack_parameters')==[['ngrdcol']+[a.name for a in read_procedure(N,'pack_parameters').arguments[1:-1]]+['packed_parameters']]
assert calls(b,'read_clubb_params_nl')==[['filename','ngrdcol','iunit']+contract['namelist_member_order']+['read_status']]
assert b.count('Kokkos::subview(local_parameters, Kokkos::ALL,')==105
assert "filename.find_first_not_of(' ') != std::string::npos" in b
assert 'f2k::init_clubb_params_api' not in T and 'f2k::pack_parameters' not in T and 'f2k::unpack_parameters' not in T
# Generate I/O text only in memory with the existing owner helper, then compare isolated fragments.
cpp=[];fort=[];env=globals()|{'root':R};exec(compile(textwrap.dedent((D/'generate_driver_io.py.inc').read_text()),'staged_driver_io','exec'),env)
assert '\n'.join(cpp)+'\n'==(D/'driver_io.hpp.inc').read_text();assert '\n'.join(fort)+'\n'==(D/'driver_io.F90.inc').read_text()
f='\n'.join(fort);assert 'native_call' not in f and f.count('open(newunit=unit')==1 and f.count('read(unit,nml=clubb_params_nl')==1 and f.count('close(unit)')==1
assert not re.search(r'(?im)^\s*call\s+(set_default|pack_|unpack_|init_clubb)',f)
assert [a.name for a in env['params_proc'].arguments]==['runfile','ngrdcol','iunit']+contract['namelist_member_order']+['read_status']
# Public single/multi-column wrapper contracts; one-column rank association is documented.
a=definitions(re.sub(r'//[^\n]*','',(D/'clubb_api_module.hpp.inc').read_text()));assert len(a)==2
for d,native_name in zip(a,['calc_derived_params_api_single_col','calc_derived_params_api_multi_col']):assert d['arguments']==[x.name for x in read_procedure(R/'src/CLUBB_core/clubb_api_module.F90',native_name).arguments]
# No owning type is captured in any science kernel and no device I/O/bridge/global state exists.
kernels=0
for m in re.finditer(r'KOKKOS_LAMBDA\s*\([^)]*\)\s*{',C):
 e=m.end();depth=1
 while depth:depth+=(C[e]=='{')-(C[e]=='}');e+=1
 k=C[m.end():e-1];assert not re.search(r'\b(?:gr|nu_vert_res_dep|err_info)\s*\.|f2k::|std::(?:cerr|cout|numeric_limits)|\bfma\s*\(',k);kernels+=1
assert kernels==4 and 'thread_local' not in C and not re.search(r'\bstruct\s+nu_vertical_res_dep',C)
files=['procedure_includes.hpp.inc','constants.hpp.inc','private_declarations.hpp.inc','procedure_definitions.hpp.inc','clubb_api_module.hpp.inc','api_includes.hpp.inc','init_clubb_params_api.cpp.inc','cpp_includes.inc','generate_driver_io.py.inc','driver_io.hpp.inc','driver_io.F90.inc','parameter_contracts.json']
hashes={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in files}
record={'review_kind':'Owner full supported-source/contract/formula audit, not an independent review. Separate cleanup record.','canonical_sha256':hashlib.sha256(N.read_bytes()).hexdigest(),'hashes':hashes,'Native_supported_order':names,'exact_signatures':8,'default_assignment_ASTs':105,'pack_assignment_ASTs':105,'unpack_assignment_ASTs':105,'derived_assignment_ASTs':24,'catalog_names_with_exact32padding':105,'exact_Fortran_flat_hard_bound_values':210,'paired_equality_predicates':4,'View_only_kernels':kernels,'API_wrapper_contracts':2,'namelist_members_in_exact_source_order':105,'temporary_IO_contract_required_args':109,'reviews':reviews,'manual_review':['all column-major hard-bound validation branches, messages, NaN behavior and per-column fatal updates','four ordered C6 equality diagnostics after hard-bound check without resetting errors or early return','all7 diffusivity mappings and allocation guard; no repair of malformed partial existing state','lmin write/fatal ordering, prescribed/type3/type1/type2 average spacing, signed physical endpoints','Native invalid grid_type undefined avg_deltaz retained behind existing setup gate','Native integer exponent2/3 expanded mechanically without reassociation, FMA or speculative BFB workaround','empty filename skips I/O; required iunit retained at call-local NEWUNIT service boundary','single-column Native sequence association represented as one-column packed View','same-named out-of-line HostSpace initializer keeps driver_io after complete schema declaration, avoiding a generated conversion include cycle; generic output adapter retains caller ownership','all data case/call owned; generated nu schema reused; immutable catalogs only at module scope'],'limits':['No shared production/generator/CMake edits, build, syntax compilation or model execution. All CPU/GPU branches and BFB remain pending root integration.','Native tuner-only read_param_minmax/read_param_constraints omitted per actual JAX/root scope.','Native OPEN/READ failure termination becomes existing I/O-owner read_status followed by C++ fatal exception; compiler-specific diagnostics/spacing are not claimed identical.','set_default_parameters and private pack accept homogeneous column View types/accessibility; actual initializer uses one HostSpace owner. No claim of mixed-space arguments in a single kernel.']}
(D/'SOURCE_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n')
lexer=re.compile(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|\w+|[^\s]')
def tokens(s):return lexer.findall(re.sub(r'//[^\n]*|/\*.*?\*/','',s,flags=re.S))
cleanup=[]
for file in files:
 p=D/file
 if file.endswith(('.hpp.inc','.cpp.inc')) and file != 'driver_io.hpp.inc':
  assert all(len(l)<=100 and not l.endswith((' ','\t')) for l in p.read_text().splitlines()),file
  before=Path('/tmp/params_before_'+file)
  if before.exists():assert tokens(before.read_text())==tokens(p.read_text());cleanup.append({'file':file,'exact_executable_token_equality':True,'tokens':len(tokens(p.read_text()))})
(D/'CLEANUP_AUDIT.json').write_text(json.dumps({'hashes':hashes,'separate_owner_cleanup':True,'checks':cleanup,'Native_continuation_groups_and_intents_retained':True,'handwritten_Cpp_max_columns':100,'generated_IO_format_uses_existing_generator_owner':True,'runtime_status':'Uncompiled/untested; all comparison evidence pending root.'},indent=2)+'\n')
print('Parameters source audit passed: eight contracts,315 default/pack/unpack assignments,24 derived expressions,105 names,210 bounds,2 API overloads,4 kernels.')
