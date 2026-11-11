"""Review pressure formulas/contracts/recurrence ownership without execution."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
r=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();r=r[r.index('class Cpp:'):r.index('reviews=[]')];r=r.replace("    assert b in [('n',2.),('n',3.)]", "    if b not in [('n',2.),('n',3.)]:return ('call','pow',(a,b))");exec(r)
NATIVE=ROOT/'src/CLUBB_core/calc_pressure.F90';TARGET=HERE/'calc_pressure.hpp';source=NATIVE.read_text();text=TARGET.read_text();code=re.sub(r'//[^\n]*','',text);defs=[]
for m in re.finditer(r'void\s+(\w+)\s*\(([^{};]*)\)\s*\{',code,re.S):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 defs.append({'name':m.group(1),'args':[a.strip().split()[-1] for a in split_top(m.group(2))],'body':code[m.end():end-1]})
assert [d['name'] for d in defs]==['init_pressure','calculate_thvm','init_pressure']
arrays={'thvm','p_sfc','p_in_pa','exner','p_in_pa_zm','exner_zm','thvm_zm','thlm','rtm','rcm','thv_ds_zt','gr_zt','gr_zm'}
def canonical(n,init=False):
 if not isinstance(n,tuple):return n
 if n==('v','thvm_epsilon'):return ('call','epsilon',(('v','thvm'),))
 if init and n[0]=='v' and n[1] in {'nzm','nzt'}:return ('v','gr_'+n[1])
 if n[0]=='call':return ('call',n[1],tuple(canonical(a,init) for a in n[2]))
 return tuple(canonical(x,init) if isinstance(x,tuple) else x for x in n)
reviews=[]
for d in defs[:2]:
 n=d['name'];p=read_procedure(NATIVE,n);assert d['args']==[a.name for a in p.arguments]
 raw=re.search(r'^  subroutine '+n+r'\(.*?^  end subroutine '+n+r'\b',source,re.S|re.M).group();expected=[]
 for l in statements(raw):
  if m:=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',l):expected.append((fexpr(m.group(1)),fexpr(m.group(2))))
 actual=[]
 for m in re.finditer(r'(?m)^\s*(?:const Real2D\s+)?(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',d['body']):
  actual.append((canonical(Cpp(m.group(1),arrays).parse(),n=='init_pressure'),canonical(Cpp(m.group(2),arrays).parse(),n=='init_pressure')))
 assert expected==actual,(n,'assignments',len(expected),len(actual),[(i,a,b) for i,(a,b) in enumerate(zip(expected,actual)) if a!=b])
 expectedifs=[fexpr(m.group(1)) for l in statements(raw) if (m:=re.match(r'if\s*\((.*)\)\s*then$',l,re.I))];actualifs=[]
 for m in re.finditer(r'\bif\s*\(',d['body']):
  end=m.end();depth=1
  while depth:depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  actualifs.append(canonical(Cpp(d['body'][m.end():end-1],arrays).parse(),n=='init_pressure'))
 assert expectedifs==actualifs,(n,'predicates')
 reviews.append({'name':n,'required_arguments':len(p.arguments),'ordered_assignment_trees':len(expected),'ordered_branch_predicates':len(expectedifs)})
init=defs[0]['body'];assert 'constexpr real_t g_ov_Cp = grav / Cp;' in init
assert 'const real_t thvm_epsilon = std::numeric_limits<real_t>::epsilon();' in init
stages=['surface','first_exner','first_pressure','thermodynamic','momentum'];positions=[init.index('"init_pressure::'+n+'"') for n in stages];assert positions==sorted(positions)
interp=init.index('zt2zm_api(');assert positions[0]<interp<positions[1]
assert 'for (int k = 2; k <= nzt; ++k)' in init
assert 'Policy(0, ngrdcol * (nzm - 1))' in init and 'const int k = cell / ngrdcol + 2;' in init
assert 'using Policy = Kokkos::RangePolicy<typename Thvm::execution_space>;' in defs[1]['body']
assert 'Thlm::rank == 2' in defs[1]['body'] and 'Thvm::rank == 2' in defs[1]['body']
# Generic initialization retains all required args, allocates shared-owner staging
# only at this boundary and restores all four original output allocations.
adapter=defs[2];p=read_procedure(NATIVE,'init_pressure');assert adapter['args']==[a.name for a in p.arguments]
for a in p.arguments:
 if a.dims:assert 'to_execution_space('+a.name+(', false' if a.intent=='out' else '')+')' in adapter['body']
call=re.search(r'init_pressure\(([^;]*)\);',adapter['body'],re.S);assert [a.strip() for a in split_top(call.group(1))]==['ngrdcol','gr','thvm_exec','p_sfc_exec','p_exec','exner_exec','p_zm_exec','exner_zm_exec']
assert re.findall(r'from_execution_space\(([^;]*)\);',adapter['body'])==['p_in_Pa, p_exec','exner, exner_exec','p_in_Pa_zm, p_zm_exec','exner_zm, exner_zm_exec']
# Captures keep only Views/scalars; the only vertical recurrence remains inside
# a column work item. No f2k, host IO or nested launch enters a device body.
kernels=[]
for m in re.finditer(r'KOKKOS_LAMBDA\([^)]*\)\s*\{',code):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 body=code[m.end():end-1];kernels.append(body);assert not re.search(r'\bgr\.|f2k::|std::|parallel_for|parallel_reduce|to_host_space',body)
assert len(kernels)==6 and not re.search(r'fma|thread_local|_exner_step',code)
assert max(map(len,text.splitlines()))<=100 and all(l.rstrip()==l and '\t' not in l for l in text.splitlines())
record={'canonical_source':str(NATIVE.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),'target_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),'routines':reviews,'kernels':len(kernels),'lower_scientific_dependency':'Existing C++zt2zm_api only; its full six-argument call checked as part of assignment trees.','manual_review':['Full265line Native and actualJAX bodies/callers read before outline. Both routines complete; no JAX-only helper or safe-branch arithmetic imported.','Surface/first-thermo/first-pressure/thermo-recurrence/momentum stages retain source order and storage indices. Thermo k dependencies remain serial within each column. Momentum cells depend only on completed thermo(k-1).','Native three epsilon predicates retain their own reference theta_v; momentum epsilon*previous-thermo value differs from the generalized JAX helper threshold and remains Native.','All output fields written, scratch owned by interpolation return, grid component aliases only. No UDT capture/global cache/bridge/stats/error path.','calculate_thvm uses output memory/execution space directly, retaining host/strided drivers and device core fields. Inputs must be accessible from that space; all current C++calls have matching residency.','Generic init adapter reuses existing array owner and restores caller outputs; Native legacy single-column sequence-association input readers remain Native external callers, not a new rank1CppAPI.','Separate cleanup reviews Native argument groups/comments/units/full hydrostatic equations, ordinary powers/logs/core constants, case-sensitive names,100columns and whitespace.'],'remaining_dependencies':'None beyond live grid/array/types/constants owners.','runtime_status':'Uncompiled/untested: CPU/GPU/model arithmetic, log/constant branches, threshold edges,initialization/host/strided caller paths all pending. No BFB claim/speculative workaround.'};(HERE/'SOURCE_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'routines':reviews,'kernels':len(kernels)}))
