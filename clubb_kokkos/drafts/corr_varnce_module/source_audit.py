"""Audit metadata contracts/tables/expression trees without executing physics."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
reader=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();reader=reader[reader.index('class Cpp:'):reader.index('reviews=[]')];reader=reader.replace("   if name in {'max','min'}:return ('if',('bin','>' if name=='max' else '<',args[0],args[1]),args[0],args[1])", "   if name in {'max','min'}:return ('call',name,args)");exec(reader)
NATIVE=ROOT/'src/CLUBB_core/corr_varnce_module.F90';CPP=HERE/'corr_varnce_module.cpp';source=NATIVE.read_text();text=CPP.read_text();code=re.sub(r'//[^\n]*','',text);outline=json.loads((HERE/'SOURCE_OUTLINE.json').read_text());names=outline['native_inventory_and_order'];defs=[]
for m in re.finditer(r'(?:int|void)\s+(\w+)\s*\(([^{};]*)\)\s*\{',code,re.S):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 defs.append({'name':m.group(1),'args':[a.strip().split()[-1] for a in split_top(m.group(2))],'body':code[m.end():end-1]})
assert [d['name'] for d in defs]==names
aliases={n+'_host':'hm_metadata_'+n for n in ['hmp2_ip_on_hmm2_ip','hydromet_list','hydromet_tol','l_mix_rat_hm','l_frozen_hm']}
aliases.update({n+'_host':n for n in ['corr_array_n_cloud','corr_array_n_below','corr_array_n']});aliases.update({n+'_view':n for n in ['corr_array_n_cloud_def','corr_array_n_below_def']})
arrays=set(aliases)|{'corr_array_n_cloud','corr_array_n_below','corr_array_n','corr_array_n_cloud_def','corr_array_n_below_def'}
def canonical(n):
 if not isinstance(n,tuple):return n
 if n[0]=='v':return ('v',aliases.get(n[1],n[1]))
 if n[0]=='field' and n[1][0]=='v':return ('v',n[1][1]+'_'+n[2])
 if n[0]=='call':return ('call',aliases.get(n[1],n[1]).replace('std::',''),tuple(canonical(a) for a in n[2]))
 if n[0]=='unary' and n[1]=='*':return canonical(n[2])
 return tuple(canonical(x) if isinstance(x,tuple) else x for x in n)
def cpp_expr(s):return canonical(Cpp(s,arrays).parse())
def raw(name):return re.search(r'^  (?:function|subroutine) '+name+r'\(.*?^  end (?:function|subroutine) '+name+r'\b',source,re.S|re.M).group()
def assignments(name):
 expected=[]
 for l in statements(raw(name)):
  m=re.match(r'^([\w%]+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',l)
  if m:
   rhs=m.group(2).strip();expected.append((fexpr(m.group(1)),('str',rhs[1:-1]) if rhs.startswith('"') and rhs.endswith('"') else fexpr(rhs)))
 return expected

def actual_assignments(body,exclude=()):
 entries=[]
 pat=r'(?m)^\s*(?:(?:const\s+)?(?:int|bool|real_t)\s+)?([\w.]+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);'
 for m in re.finditer(pat,body):
  lhs,rhs=m.groups()
  if lhs in exclude or rhs.strip().startswith('decltype('):continue
  if lhs=='hm_metadata':assert rhs.strip()=='{}';continue
  if lhs.startswith('hydromet_list_host('):
   args=split_top(lhs[lhs.index('(')+1:-1]);assert args[0]=='character';idx=Cpp(args[1],set()).parse();assert idx[0:2]==('bin','-') and idx[3]==('n',1.)
   string=re.fullmatch(r'"([^"]*)"\[character\]',rhs.strip()).group(1);assert len(string)==10
   entries.append((m.start(),('call','hm_metadata_hydromet_list',(idx[2],)),('str',string.rstrip(' '))));continue
  entries.append((m.start(),cpp_expr(lhs),cpp_expr(rhs)))
 for m in re.finditer(r'Kokkos::deep_copy\s*\(([^;]*)\);',body,re.S):
  args=split_top(m.group(1));assert len(args)==2
  if args[0].startswith('hm_metadata.'):continue # End-of-call host mirror restoration.
  entries.append((m.start(),cpp_expr(args[0]),cpp_expr(args[1])))
 return [(a,b) for _,a,b in sorted(entries)]

reviews=[]
for d in defs:
 n=d['name'];contract=next(p for p in outline['routines'] if p['name']==n);args=[a['name'] for a in contract['arguments']];assert d['args']==args,(n,'signature',d['args'],args)
 # Exact signature names/order include Native's private argument lacking INTENT;
 # its complete body only reads metadata, documented as const in source outline.
 if n in {'def_corr_idx','set_corr_arrays_to_default','get_corr_var_index','init_pdf_hydromet_arrays_api','setup_corr_varnce_array_api'}:
  expected=assignments(n);actual=actual_assignments(d['body'])
  assert expected==actual,(n,'assignments',len(expected),len(actual),[(i,a,b) for i,(a,b) in enumerate(zip(expected,actual)) if a!=b][:3])
  count=len(expected)
 else:count='manual error/formatting review below'
 native_conditions=[]
 for l in statements(raw(n)):
  if m:=re.match(r'(?:if|elseif)\s*\((.*)\)\s*then$',l,re.I):native_conditions.append(fexpr(m.group(1)))
 actual_conditions=[]
 for m in re.finditer(r'\bif\s*\(',d['body']):
  end=m.end();depth=1
  while depth:depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  cond=d['body'][m.end():end-1];cond=re.sub(r'(\w+)\.has_value\(\)',r'present(\1)',cond)
  actual_conditions.append(cpp_expr(cond))
 if n not in {'get_corr_var_index','print_corr_matrix'}:assert native_conditions==actual_conditions,(n,'predicates',native_conditions,actual_conditions)
 elif n=='get_corr_var_index':
  cases=re.findall(r'case\(\s*"([^"]*)"\s*\)',raw(n));assert len(cases)==12
  assert [p[3][1] for p in actual_conditions]==cases
  assert "while (!trimmed.empty() && trimmed.back() == ' ')" in d['body']
  assert 'trimmed.pop_back();' in d['body']
 # Public declaration order and all names/default optional statuses remain Native.
 decls=re.sub(r'//[^\n]*','',(HERE/'procedure_declarations.hpp.inc').read_text())
 if n not in {'def_corr_idx','set_corr_arrays_to_default'}:
  declarations=re.findall(r'\b'+n+r'\s*\(([^{};]*)\)\s*[;{]',decls,re.S)
  declaration_=declarations[0]
  actualargs=[x.strip().split('=')[0].strip().split()[-1] for x in split_top(declaration_)];assert actualargs==args,(n,'header',actualargs,args)
  for a in contract['arguments']:
   if a['optional']:assert re.search(r'\b'+a['name']+r'\s*=\s*std::nullopt',declaration_)
 reviews.append({'name':n,'arguments':len(args),'ordered_assignment_trees':count,'ordered_native_if_predicates':len(native_conditions),'ordered_string_cases':12 if n=='get_corr_var_index' else 0})

# Verify the private mapping constants, source tolerance and diagnostic text.
expected_slots=re.findall(r'\b(ii_\w+)\s*=\s*(\d+)',raw('def_corr_idx'))
actual_slots=re.findall(r'constexpr int (ii_\w+)\s*=\s*(\d+)',defs[0]['body'])
assert expected_slots==actual_slots and len(actual_slots)==12
assert 'constexpr real_t tol = real_t(1.0e-6);' in defs[-2]['body']
setup_body=next(d['body'] for d in defs if d['name']=='setup_corr_varnce_array_api')
native_warnings=[''.join(re.findall(r'"([^"]*)"',l)) for l in statements(raw('setup_corr_varnce_array_api')) if l.startswith('write(fstderr,*)')]
cpp_warnings=[''.join(json.loads(q) for q in re.findall(r'"(?:\\.|[^"\\])*"',m.group(1))).rstrip('\n') for m in re.finditer(r'std::cerr\s*<<(.*?);',setup_body,re.S)]
assert native_warnings==cpp_warnings and len(cpp_warnings)==2

# Table values, source kind suffix and Fortran flat ordering are independent
# of the procedure parser and mirror alias representation.
tables=(HERE/'default_tables.hpp.inc').read_text();table_reviews=[]
for n in ['corr_array_n_cloud_def','corr_array_n_below_def']:
 native=re.search(n+r'\s*=\s*reshape\(\s*&\s*(.*?)\(/d_var_total,d_var_total/\)',source,re.S).group(1)
 expected=re.findall(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)_c',native)
 target=re.search(r'\b'+n+r'\s*=\s*\{(.*?)\};',tables,re.S).group(1);actual=re.findall(r'real_t\(([^)]*)\)',target)
 assert expected==actual and len(actual)==144 and 'inline constexpr std::array<real_t' in tables
 table_reviews.append({'name':n,'core_kind_literals_in_exact_Fortran_flat_order':144})
param=(HERE/'parameter_types.hpp.inc').read_text();type_reviews=[]
for n in ['hmp2_ip_on_hmm2_ip_ratios_type','hmp2_ip_on_hmm2_ip_slope_type','hmp2_ip_on_hmm2_ip_intrcpt_type']:
 native=re.search(r'^  type '+n+r'\b(.*?)^  end type '+n,source,re.S|re.M).group(1);expected=re.findall(r'(\w+)\s*=\s*([\d.eE+-]+)_core_rknd',native)
 target=re.search(r'struct '+n+r'\s*\{(.*?)\};',param,re.S).group(1);actual=re.findall(r'real_t\s+(\w+)\s*=\s*real_t\(([^)]*)\)',target);assert expected==actual and len(actual)==8
 type_reviews.append({'name':n,'ordered_native_core_kind_fields_and_defaults':8})
assert not re.search(r'thread_local|KOKKOS_LAMBDA|struct hm_metadata_type',code+tables+param)
assert 'using DefaultMatrix = Kokkos::View<const real_t**' in code
assert not re.search(r'(?m)^\s*(?:inline\s+)?(?:static\s+)?(?:thread_local\s+)?(?:Real\dD|HostMatrix)\s+\w+\s*[=;(]',code[:code.index('int def_corr_idx')])
# Full lower calls: default helper then cloud/below mirror, no scientific call hidden.
setup=next(d['body'] for d in defs if d['name']=='setup_corr_varnce_array_api');native_calls=[fexpr(re.sub(r'^call\s+','',l)) for l in statements(raw('setup_corr_varnce_array_api')) if re.match(r'call\s+',l,re.I)];actual_calls=[]
for m in re.finditer(r'(?:f2k::)?(?:set_corr_arrays_to_default|mirror_lower_triangular_matrix)\s*\(',setup):
 end=m.end();depth=1
 while depth:depth+=(setup[end]=='(')-(setup[end]==')');end+=1
 actual_calls.append(cpp_expr(setup[m.start():end]))
assert actual_calls==native_calls
# Native first-failure fatal semantics, two exact header choices, matrix traversal.
assertbody=next(d['body'] for d in defs if d['name']=='assert_corr_symmetric');assert assertbody.count('Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);')==2 and assertbody.count('return;')==2
assert 'to_host_space(err_info.err_header)' in assertbody and 'err_info.err_header_global' in assertbody
assert assertbody.count('corr_array_n_host(row, column)')==2
for t in ['Correlation array is non symmetric.','Correlation array is formatted incorrectly.','in subroutine assert_corr_symmetric: ']:assert t in assertbody
# Native print layout control, including the unusual LEN10/F5.2/stride6/LEN72.
printbody=next(d['body'] for d in defs if d['name']=='print_corr_matrix')
for s in ['std::fixed << std::setprecision(2) << std::setw(5)','corr_array_n_host(m - 1, n - 1)','str_array_value.resize(10, \' \');','current_character_index = current_character_index + 6;','std::cout << current_line','current_line.resize(72, \' \');']:assert s in printbody
allfiles=[CPP,HERE/'parameter_types.hpp.inc',HERE/'default_tables.hpp.inc',HERE/'procedure_declarations.hpp.inc',HERE/'f2k_dependencies.hpp']
for p in allfiles:
 t=p.read_text();assert max(map(len,t.splitlines()))<=100 and all(l.rstrip()==l and '\t' not in l for l in t.splitlines())
record={'canonical_source':str(NATIVE.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in allfiles},'routines':reviews,'tables':table_reviews,'scalar_parameter_types':type_reviews,'ordered_public_lower_calls':len(actual_calls),'metadata_owner':'Exact existing generated26-field hm_metadata_type; no duplicate;INTENT OUT reset and five local host mirrors restore caller-owned View components.','manual_review':['Full950line Native and852line JAX reads precede outline; no JAX-only compatibility metadata/helpers or reduced diagnostic interfaces imported.','Root-approved immutable default tables retain all288 core-kind values and order; no mutable/thread_local cache/global allocation. External mutation is not a new supported C++API.','Initialization retains all species index/variance/name/padding/flag/tolerance assignment order, optional slopes/intercepts and sequential PDF counts; valid Native species indices are unique and cover hydromet_dim.','Native host scalar metadata work has no GPU captures/launches/bridges; persistent fields remain Kokkos Views. Generated owner and converter untouched.','Public setup retains both-present override semantics, exact cloud/below mirroring order and warning predicates/messages.','Fatal check order/header array vs global choice/all-column status/first return/full matrix payload retained. List-directed numeric text uses local max_digits10 stream; exact compiler-dependent decimal formatting is not claimed.','Fixed-format print keeps transpose orientation,F5.2,LEN10,stride6,LEN72 and stdout; exact decimal/NaN/infinity formatting remains unverified.','Separate final cleanup reviews public contracts/optionals, Native argument groups, comments/units, declared type owner, literal kinds,100columns/whitespace and no target-only scientific helper.'],'remaining_lower':'One missing required public host bridge: mirror_lower_triangular_matrix(int nvars, rank2 real matrix InOut). Existing Native debug query remains; no new constant needed.','runtime_status':'Uncompiled/untested. No unit/model execution, CPU/GPU equality or BFB claim; no speculative arithmetic workaround.'};(HERE/'SOURCE_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'routines':reviews,'tables':table_reviews,'types':type_reviews}))
