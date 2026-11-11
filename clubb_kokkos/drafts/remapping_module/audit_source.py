"""Record full source/contract/ownership review; do not run scientific code."""
from pathlib import Path
import re,json,hashlib,sys,ast
D=Path(__file__).resolve().parent;R=D.parents[2];sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'))
from generate_derived_types import split_top
O=json.loads((D/'SOURCE_OUTLINE.json').read_text());S=(R/O['source']).read_text();C=(D/'remapping_module.cpp').read_text();H=(D/'remapping_module.hpp').read_text()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clean(s):return re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"',lambda m:re.sub(r'[^\n]',' ',m[0]) if m[0].startswith('/') else m[0],s,flags=re.S)
def balanced(s,op,l='(',r=')'):
 end=op+1;depth=1
 while depth:depth+=(s[end]==l)-(s[end]==r);end+=1
 return s[op+1:end-1],end
def defs(text):
 rows=[];T=clean(text)
 for m in re.finditer(r'\b('+ '|'.join(O['native_inventory_and_order'])+r')\s*\(',T):
  args,end=balanced(T,m.end()-1);pos=end
  while pos<len(T) and T[pos].isspace():pos+=1
  if T[pos:pos+1]!='{':continue
  body,stop=balanced(T,pos,'{','}');actual=split_top(args)
  if not all(re.search(r'\s+[*&]?\w+(?:\s*=.*)?$',a) for a in actual):continue
  rows.append((m[1],actual,body,m.start()))
 return rows
rows=defs(C);names=[r[0] for r in rows];names.remove('kmppm')
assert names==O['native_inventory_and_order'],names
contracts={p['name']:p for p in O['routines']};bodies={};contract_counts={}
for n,args,body,pos in rows:
 if n=='kmppm' and len(args)==3:
  assert [a.strip().split()[-1] for a in args]==['dm','a4','lmt'];continue
 actual=[re.sub(r'\s*=.*','',a).strip().split()[-1].lstrip('*&') for a in args]
 expected=[a['name'] for a in contracts[n]['arguments']];assert actual==expected,(n,actual,expected)
 for arg,formal in zip(contracts[n]['arguments'],args):
  assert ('*' in formal)==arg['optional'],(n,arg['name'],formal)
  if arg['optional']:assert re.search(r'=\s*nullptr',formal)
 bodies[n]=body;contract_counts[n]=len(args)
 assert '// '+actual[0]+':' in C[pos:pos+len(body)+2000],n
assert len(contract_counts)==18
assert all(not p['optional'] for p in contracts['remap_vals_to_target']['arguments'])
assert 'const grid &' in H and 'gr_source' in H and 'gr_target' in H
assert 'typename Source::memory_space' in H and 'from_execution_space(result, remapped)' in H
assert [(n,len(a)) for n,a,b,p in defs(H)]==[('remap_vals_to_target',13)]
# Exact field-name selector strings and case are authoritative.
source_name_block='\n'.join(S.splitlines()[40:80]);native_names=[]
for clause in re.findall(r'case\s*\((.*?)\)',source_name_block,re.S|re.I):native_names+=re.findall(r'"([^"]*)"',clause)
cpp_names=re.findall(r'name\s*==\s*"([^"]*)"',bodies['remap_iv_from_var_name'])
assert native_names==cpp_names
# Every Native explicit CALL is retained in order; the array limiter's scalar
# dispatch is the one additional Native-named shape/device overload.
call_counts={}
for n in bodies:
 expected=[re.match(r'call\s+(\w+)',s,re.I)[1].lower() for s in contracts[n]['ordered_calls']]
 called=set(expected);actual=[]
 for m in re.finditer(r'\b('+ '|'.join(called)+r')\s*\(',bodies[n]) if called else []:actual.append(m[1].lower())
 assert actual==expected,(n,expected,actual)
 call_counts[n]=len(actual)
# Preserve debug checks at their actual host boundary and in Native order.
assert bodies['remap_vals_to_target_helper'].index('calc_mass_over_grid_intervals(')<bodies['remap_vals_to_target_helper'].index('if (grid_remap_method == cons_ullrich_remap && R_ij == nullptr)')
assert bodies['remap_vals_to_target_helper'].index('check_conservation(')<bodies['remap_vals_to_target_helper'].index('check_monotonicity(')<bodies['remap_vals_to_target_helper'].index('check_consistency(')
assert bodies['remapping_matrix'].index('check_remap_matrix_conservation(')<bodies['remapping_matrix'].index('check_remap_matrix_consistency(')<bodies['remapping_matrix'].index('check_remap_matrix_monotonicity(')
assert bodies['ppm2m'].index('steepz(')<bodies['ppm2m'].index('if (kord >= 7)')
assert 'if (km > 32)' in bodies['ppm2m']
assert bodies['map1_ppm'].index('return;')<bodies['map1_ppm'].index('qsumk(column, k - 1) = zero')
assert 'kk <= km + 1' in bodies['map1_ppm'] and 'pe2_host(i, kk - 1)' in bodies['map1_ppm']
assert 'Kokkos::WithoutInitializing, "remap_vals_ppm::target_values_flipped"' in bodies['remap_vals_ppm']
assert re.search(r'const int .*?kord = 4',bodies['remap_vals_ppm'],re.S)
assert 'if constexpr' not in C
assert set(re.findall(r'f2k::(\w+)\(',C))=={'clubb_at_least_debug_level_api'}
assert not re.search(r'static\s+(?!assert|constexpr)',clean(C))
# Device bodies may call only the device limiter and the existing scalar
# interpolation owner. No host bridges, diagnostics or grid-object captures.
kernels=[];T=clean(C)
for m in re.finditer(r'KOKKOS_LAMBDA\s*\([^)]*\)\s*\{',T):
 body,end=balanced(T,m.end()-1,'{','}')
 assert not re.search(r'f2k::|std::(?:cerr|cout)|to_host_space|gr_source\.|gr_target\.',body)
 kernels.append({'line':C.count('\n',0,m.start())+1,'calls_scalar_kmppm':bool(re.search(r'\bkmppm\(',body))})
# Major compile-time/native-profile conventions and retained undefined paths
# were reread independently of the formula parser during the owner audit.
manual={
 'calc_mass_over_grid_intervals':['Both upper/lower spline coverage errors and full two-line payloads','Moving spline index, node equality tolerance, trapezoid sum/upper interpolation order; no added node clamp'],
 'matrix_vector_mult':['Owned result, every column/target-cell store, serial ascending source-cell fold'],
 'vertical_integral_conserve_mass':['Native SUM represented by explicit ascending fold; exact compiler association unverified'],
 'check_conservation':['Source then target integral; relative AND absolute predicates, exact warning/data order'],
 'check_monotonicity':['Global extrema and every violating value; PPM checks cells4..N-3 only'],
 'check_consistency':['Native matrix owner on ones; singleton-rank adapters retain storage lifetime'],
 'remap_vals_to_target_helper':['Cached optional operator stays read-only; mass and ordered pressures still computed before presence selection','Ullrich/PPM/unsupported method branch order; debug conservation/monotonicity/consistency per Native column'],
 'remap_vals_to_target':['Read-only grids despite missing Native INTENT','ZT uses surrounding ZM; ZM builds two exterior edges plus all interior ZT, with exact counts'],
 'check_remap_matrix_conservation':['Serial weighted column sums and first failing column; exact warning+fatal text'],
 'check_remap_matrix_consistency':['Serial row sums and first failing row; exact warning+fatal text'],
 'check_remap_matrix_monotonicity':['Any negative entry per row, then early outer exit; exact warning+fatal text'],
 'remapping_matrix':['Independent overlap cells, zero/nonzero branches and target-thickness division','Debug operator checks conservation,consistency,monotonicity order'],
 'kmppm':['All lmt0/1/2/3 branches and condition ordering','Scalar overload reads A6 only in Native-required branches; array3 returns before launches/reads'],
 'steepz':['Rat k2..km; F k2..km-1; alfa k3..km-2; updated AL k4..km-2; every shifted level checked'],
 'ppm2m':['All top/bottom iv0,-1,other constraints, km>32 sponge, AR-copy stage','Top2/interior/bottom2 limiter order; Huynh>=7 private branches retained','kord4 does not load unwritten A6 before Lin limiter initializes it; interior-only independent-cell launch fusion'],
 'map1_ppm':['Original16args, longitude/latitude physical vs scratch coordinates and k0 zero-index storage','Ascending moving search, first-column/level diagnostic, unchanged q2 before failure return','Retain undefined failed-diagnostic pe2 read when kn<km; no fabricated zeros/search clamp','Full-cell kl sum, first/last partial sums and divisions preserve Native arithmetic trees'],
 'remap_vals_ppm':['Public PPM hardcodes kord4 and singleton latitude controls','Edges/values reversed in Native order, one sequential map1 call per column, output scratch reused across columns'],
 'remap_iv_from_var_name':['Two independent case-sensitive lists, ADJUSTL/TRIM/LEN64 blank/truncation handling']}
assert set(manual)==set(O['native_inventory_and_order'])
formatting={}
for p in [D/'remapping_module.cpp',D/'remapping_module.hpp']:
 lines=p.read_text().splitlines();long=[i+1 for i,l in enumerate(lines) if len(l)>100];trail=[i+1 for i,l in enumerate(lines) if l.rstrip()!=l]
 assert not long and not trail,(p.name,long,trail)
 formatting[p.name]={'sha256':sha(p),'lines':len(lines),'max_columns':max(map(len,lines)),'trailing_whitespace':False}
record={'provenance':'Owner source/cleanup review, not an independent-agent verdict.','source_sha256':sha(R/O['source']),'jax_sha256':sha(R/O['jax']),'targets':formatting,'all_18_Native_bodies_reviewed':True,'definitions_in_Native_order':True,'Native_argument_counts':contract_counts,'total_Native_arguments':sum(contract_counts.values()),'exact_selector_names':len(cpp_names),'ordered_explicit_calls':call_counts,'device_kernel_count':len(kernels),'device_kernels':kernels,'manual_routine_review':manual,'device_only_same_kind_constant_materializations':5,'deliberate_BFB_arithmetic_changes':[],'runtime_and_BFB_status':'Pending; no live C++ caller. Serial compilation is separate evidence.','remaining_findings':[]}
(D/'source_audit.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k not in {'manual_routine_review','device_kernels'}},indent=2))
