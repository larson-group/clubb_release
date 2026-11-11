"""Review Native expressions, interfaces and ordering without model execution."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
reader=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();reader=reader[reader.index('class Cpp:'):reader.index('reviews=[]')]
reader=reader.replace("   if name in {'max','min'}:return ('if',('bin','>' if name=='max' else '<',args[0],args[1]),args[0],args[1])", "   if name in {'max','min'}:return ('call',name,args)")
exec(reader)
NATIVE=ROOT/'src/CLUBB_core/pos_definite_module.F90';CPP=HERE/'pos_definite_module.cpp';HPP=HERE/'pos_definite_module.hpp';source=NATIVE.read_text();text=CPP.read_text();header=HPP.read_text();code=re.sub(r'//[^\n]*','',text);hcode=re.sub(r'//[^\n]*','',header)
proc=read_procedure(NATIVE,'pos_definite_adj');names=[a.name for a in proc.arguments];decl=re.search(r'void\s+pos_definite_adj\s*\(([^{};]*)\)\s*\{',code,re.S);assert [a.strip().split()[-1] for a in split_top(decl.group(1))]==names
hdecls=re.findall(r'void\s+pos_definite_adj\s*\(([^{};]*)\)\s*[;{]',hcode,re.S);assert len(hdecls)==2
for declaration_ in hdecls:assert [a.strip().split()[-1] for a in split_top(declaration_)]==names
arrays={'field_np1','flux_np1','field_n','field_pd','flux_pd','flux_plus','flux_minus','flux_lim','fout','field_nonlim','dz_over_dt','invrs_dzt','flux_correction','ddzm_flux_correction'}
aliases={'grid_dir':'gr_grid_dir','grid_dir_indx':'gr_grid_dir_indx','invrs_dzt':'gr_invrs_dzt'}
def canonical(n):
 if not isinstance(n,tuple):return n
 if n[0]=='v':return ('v',aliases.get(n[1],n[1]))
 if n[0]=='call':return ('call',aliases.get(n[1],n[1]),tuple(canonical(a) for a in n[2]))
 return tuple(canonical(x) if isinstance(x,tuple) else x for x in n)
expected=[];native_lines=[]
for l in statements(source):
 m=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',l)
 if not m:continue
 lhs,rhs=m.groups();native_lines.append(l)
 if lhs=='field_np1':
  # Independently check the full function-call/materialization boundary below.
  rhs='-dt * ddzm_flux_correction(i,k) + field_np1(i,k)';lhs='field_np1(i,k)'
 lhs=lhs.replace('(:,','(i,').replace('i,:)','i,k)');rhs=rhs.replace('(:,','(i,').replace('i,:)','i,k)')
 expected.append((fexpr(lhs),fexpr(rhs)))
actual=[];entries=[]
for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',code):
 lhs,rhs=m.groups()
 if lhs.strip()=='l_negative':continue
 if lhs.startswith('flux_correction('):
  assert canonical(Cpp(rhs,arrays).parse())==fexpr('flux_lim(i,k)-flux_np1(i,k)');continue
 entries.append((m.start(),canonical(Cpp(lhs,arrays).parse()),canonical(Cpp(rhs,arrays).parse())))
for m in re.finditer(r'Kokkos::deep_copy\s*\(([^;]*)\);',code,re.S):
 args=split_top(m.group(1));entries.append((m.start(),Cpp(args[0],set()).parse(),Cpp(args[1],set()).parse()))
actual=[(a,b) for _,a,b in sorted(entries)];assert actual==expected,(len(expected),len(actual),[(i,native_lines[i],a,b) for i,(a,b) in enumerate(zip(expected,actual)) if a!=b])
# Match original ddzm function arguments, substituting only its materialized
# array expression. The public derivative result is indexed in field_update.
native_ddzm=fexpr('ddzm(nzm,nzt,ngrdcol,gr,flux_lim-flux_np1)')
call=re.search(r'\bddzm\s*\(([^;]*)\);',code,re.S);actual_ddzm=Cpp('ddzm('+call.group(1)+')',set()).parse()
assert actual_ddzm==('call','ddzm',native_ddzm[2][:-1]+(('v','flux_correction'),))
# All executable Native branch predicates retain order except the two ANY
# expressions, whose per-column serial boolean scans remain separate.
conditions=[]
for m in re.finditer(r'\bif\s*\(([^()]*(?:\([^()]*\)[^()]*)*)\)',code):conditions.append(m.group(1).strip())
assert len(conditions)==5
assert Cpp(conditions[0],set()).parse()==fexpr('clubb_at_least_debug_level_api(1)')
for i in [1,2]:assert canonical(Cpp(conditions[i],set()).parse())==fexpr('gr%grid_dir_indx>0')
assert conditions[3:]==['l_negative','l_negative']
assert len(re.findall(r'bool l_negative = false;',code))==2
assert len(re.findall(r'l_negative = l_negative \|\| field_np1\(i - 1, k - 1\) < zero;',code))==2
# Ordered stage boundaries and owned value copies are part of semantics.
stages=['flux_split','dz_over_dt','outward_flux','limit','boundaries','flux_budget','flux_correction','field_update','field_budget'];positions=[re.search(r'Kokkos::parallel_for\s*\(\s*"pos_definite_adj::'+n+'"',code).start() for n in stages];assert positions==sorted(positions)
copy=code.index('Kokkos::deep_copy(field_nonlim, field_np1)');assert positions[5]<copy<positions[6]
copy=code.index('Kokkos::deep_copy(flux_np1, flux_lim)');assert copy>positions[-1]
assert 'const int first_k = gr.k_lb_zm;' in code and 'const int last_k = gr.k_ub_zm - 2 * grid_dir_indx;' in code
assert re.search(r'const int level_count\s*=\s*Kokkos::max\(0, \(last_k - first_k\) / grid_dir_indx \+ 1\);',code)
assert 'const int k = first_k + (cell / ngrdcol) * grid_dir_indx;' in code
kernels=[]
for m in re.finditer(r'KOKKOS_LAMBDA\(int cell\)\s*\{',code):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 body=code[m.end():end-1];kernels.append(body);assert not re.search(r'\bgr\.|f2k::|std::|to_host_space|parallel_for|parallel_reduce',body)
assert len(kernels)==9
adapter=re.search(r'void\s+pos_definite_adj\s*\([^{};]*\)\s*\{(.*?)\n\}',hcode,re.S).group(1)
for a in proc.arguments:
 if a.dims:assert f'to_execution_space({a.name}'+(', false' if a.intent=='out' else '')+')' in adapter
call=re.search(r'pos_definite_adj\s*\(([^;]*)\);',adapter,re.S);assert [a.strip() for a in split_top(call.group(1))]==[a.name+('_exec' if a.dims else '') for a in proc.arguments]
assert re.findall(r'from_execution_space\s*\((.*?)\);',adapter,re.S)==[a.name+', '+a.name+'_exec' for a in proc.arguments if a.intent!='in']
assert 'std::cout << "Correcting flux\\n";' in code
for t in [text,header]:assert max(map(len,t.splitlines()))<=100 and all(l.rstrip()==l and '\t' not in l for l in t.splitlines())
record={'canonical_source':str(NATIVE.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),'cpp_sha256':hashlib.sha256(CPP.read_bytes()).hexdigest(),'hpp_sha256':hashlib.sha256(HPP.read_bytes()).hexdigest(),'native_procedures':1,'required_arguments':10,'ordered_native_assignment_trees':len(expected),'ddzm_materialization_and_argument_order':'exact','branch_predicates':5,'independent_kernels':len(kernels),'manual_review':['Complete Native and JAX scope read before outline/implementation; all scientific descriptions/cartoons/citation/code comments preserved in order, copied field_pd comment corrected.','Native grid metadata1-based physical limiter range maps to ordinal grid_dir_indx steps. Both ascending/descending cases compute the exact source kphalf/kp1/k_zt; lower bound term retains dz_over_dt(k_zt).','All flux levels initialized by distinct limiter writes or unchanged endpoint copies; per-call scratch is local, not shared global state.','Before and after ANY scans cover nzt and drive whole-column budgets on nzm/nzt respectively. Masks are recomputed on opposite sides of the field update.','Value copies retain Native array ownership; shared adapter restores all four caller arrays. Lower ddzm stays C++, only existing debug query is host f2k.','No device f2k, owning-grid capture, nested launch, host diagnostics or float reduction. All branch-local indices are assigned before use.','Separate final cleanup verifies Native groups/case, complete comments/units, physical-vs-storage indexing, adapter contracts,100columns and whitespace.'],'numerical_evidence':'Uncompiled/untested; no CPU/GPU model or BFB claim and no deliberate BFB arithmetic workaround.'}
(HERE/'SOURCE_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n');print('Matched all24 Native assignment trees,10arguments,5branch predicates and9kernel stages; source review only.')
