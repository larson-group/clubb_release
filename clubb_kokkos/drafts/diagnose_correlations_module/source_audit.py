"""Audit Native correlation contracts/formulas/calls without executing physics."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
r=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();r=r[r.index('class Cpp:'):r.index('reviews=[]')];r=r.replace("   if name in {'max','min'}:return ('if',('bin','>' if name=='max' else '<',args[0],args[1]),args[0],args[1])", "   if name in {'max','min'}:return ('call',name,args)");exec(r)
NATIVE=ROOT/'src/CLUBB_core/diagnose_correlations_module.F90';TARGET=HERE/'diagnose_correlations_module.hpp';source=NATIVE.read_text();text=TARGET.read_text();code=re.sub(r'//[^\n]*','',text);names=json.loads((HERE/'SOURCE_OUTLINE.json').read_text())['native_active_inventory_and_order'];defs=[]
for m in re.finditer(r'(?:void|real_t)\s+(\w+)\s*\(([^{};]*)\)\s*\{',code,re.S):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 defs.append({'name':m.group(1),'args':[a.strip().split()[-1] for a in split_top(m.group(2))],'body':code[m.end():end-1]})
assert [d['name'] for d in defs[:10]]==names
arrays={'s_1j','corr_matrix_approx','corr_matrix_prescribed','corr_mtx_approx','corr_matrix','s','corr_cholesky_mtx_t','host','corr_array_swapped','swap_array'}
def cpp_expr(s):return Cpp(s,arrays).parse()
def raw(n):return re.search(r'^  (?:function|subroutine) '+n+r'\(.*?^  end (?:function|subroutine) '+n+r'\b',source,re.S|re.M).group()
def normalize(n):
 if not isinstance(n,tuple):return n
 if n[0]=='v' and n[1]=='sqrt_sigma2_on_mu2_ip':return ('call','sqrt',(('v','sigma2_on_mu2_ip_array'),))
 if n[0]=='call':return ('call',n[1].replace('diagnose_correlations_detail::',''),tuple(normalize(a) for a in n[2]))
 return tuple(normalize(x) if isinstance(x,tuple) else x for x in n)
reviews=[];callcount=0
for d in defs[:10]:
 n=d['name'];p=read_procedure(NATIVE,n);assert d['args']==[a.name for a in p.arguments],(n,'signature')
 expected=[]
 if n in ['diagnose_corr','calc_w_corr','calc_varnce','calc_mean','calc_cholesky_corr_mtx_approx','setup_corr_cholesky_mtx']:
  for l in statements(raw(n)):
   if m:=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',l):
    lhs,rhs=m.groups()
    if n=='setup_corr_cholesky_mtx' and lhs=='corr_cholesky_mtx_t(j,i)' and rhs=='zero':lhs='corr_cholesky_mtx_t'
    expected.append((fexpr(lhs),fexpr(rhs)))
  entries=[]
  for m in re.finditer(r'(?m)^\s*(?:real_t\s+)?(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',d['body']):
   if m.group(1) not in ['f_ij','calc_w_corr'] and '(' not in m.group(1):continue
   entries.append((m.start(),cpp_expr(m.group(1)),cpp_expr(m.group(2))))
  if n in ['calc_varnce','calc_mean']:
   m=re.search(r'return\s+(.*?);',d['body'],re.S);entries.append((m.start(),('v',n),cpp_expr(m.group(1))))
  if n=='setup_corr_cholesky_mtx':
   m=re.search(r'Kokkos::deep_copy\(([^;]*)\);',d['body']);args=split_top(m.group(1));entries.append((m.start(),cpp_expr(args[0]),cpp_expr(args[1])))
  actual=[(a,b) for _,a,b in sorted(entries)];assert expected==actual,(n,'assignments',len(expected),len(actual),[(i,a,b) for i,(a,b) in enumerate(zip(expected,actual)) if a!=b])
  exprcount=len(expected)
 else:exprcount='manual array/stop/ordered-intrinsic adaptation below'
 expectedcalls=[normalize(fexpr(re.sub(r'^call\s+','',l))) for l in statements(raw(n)) if re.match(r'call\s+',l,re.I)];actualcalls=[]
 for m in re.finditer(r'(?:diagnose_correlations_detail::)?(?:rearrange_corr_array|diagnose_corr|setup_corr_cholesky_mtx|cholesky_to_corr_mtx_approx|corr_array_assertion_checks)\s*\(',d['body']):
  end=m.end();depth=1
  while depth:depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  actualcalls.append(normalize(cpp_expr(d['body'][m.start():end])))
 assert expectedcalls==actualcalls,(n,'calls',expectedcalls,actualcalls);callcount+=len(actualcalls)
 reviews.append({'name':n,'required_arguments':len(p.arguments),'ordered_assignment_trees':exprcount,'ordered_calls':len(actualcalls)})
# Scalar f_ij clipping, scalar w-correlation clipping and fatal-check predicates.
for n in ['diagnose_corr','calc_w_corr']:
 d=next(d for d in defs[:10] if d['name']==n);expected=[]
 for l in statements(raw(n)):
  if m:=re.match(r'(?:if|else if)\s*\((.*)\)\s*then$',l,re.I):
   if m.group(1).strip()!='.false.':expected.append(fexpr(m.group(1)))
 actual=[]
 for m in re.finditer(r'\bif\s*\(([^()]*)\)',d['body']):
  if 'f_ij' in m.group(1) or 'calc_w_corr' in m.group(1):actual.append(cpp_expr(m.group(1)))
 assert expected==actual,(n,'clipping')
# Native recurrences keep products serial inside each output work item.
setup=next(d['body'] for d in defs[:10] if d['name']=='setup_corr_cholesky_mtx')
for s in ['for (int i = 1; i <= j - 1; ++i)','for (int k = 1; k <= i - 1; ++k)']:assert s in setup
# The Native MATMUL adaptation contracts k in ascending order without symmetrizing
# or triangular pruning; this is an unverified compiler-arithmetic boundary.
mm=next(d['body'] for d in defs[:10] if d['name']=='cholesky_to_corr_mtx_approx')
for s in ['real_t correlation_sum = real_t(0.);','for (int k = 1; k <= n_variables; ++k)','corr_matrix_approx(i - 1, j - 1) = correlation_sum;']:assert s in mm
expr=re.search(r'correlation_sum = correlation_sum\s*\+\s*([^;]*);',mm).group(1);assert cpp_expr(expr)==fexpr('corr_cholesky_mtx_t(i,k)*corr_cholesky_mtx_t(j,k)')
# Swapping preserves six Native array-assignment stages, including the saved
# first column and reversed slice endpoint meanings. No shared/global storage.
swap=next(d['body'] for d in defs[:10] if d['name']=='rearrange_corr_array')
assert 'Kokkos::deep_copy(corr_array_swapped, corr_array);' in swap
assert 'Kokkos::deep_copy(swap_array, Kokkos::subview(corr_array_swapped, Kokkos::ALL, 0));' in swap
swap_pairs=[('corr_array_swapped(i,1)','corr_array_swapped(iiPDF_w,iiPDF_w-i+1)'),('corr_array_swapped(i,1)','corr_array_swapped(i,iiPDF_w)'),('corr_array_swapped(iiPDF_w,j)','swap_array(iiPDF_w-j+1)'),('corr_array_swapped(i,iiPDF_w)','swap_array(i)')]
actual=[(cpp_expr(m.group(1)),cpp_expr(m.group(2))) for m in re.finditer(r'(?m)^\s*(corr_array_swapped\([^;]*?\))\s*=\s*([^;]*);',swap)]
assert actual==[(fexpr(a),fexpr(b)) for a,b in swap_pairs]
assert len(re.findall(r'Policy|parallel_for',swap))==4
# Exact supported API gate and source stop strings, rather than JAX boolean checks.
parent=defs[0]['body'];assert parent.index('if (l_calc_w_corr)')<parent.index('Real2D corr_array_pre_swapped')
assert 'Native working matrix is uninitialized' in parent
assertion=next(d['body'] for d in defs[:10] if d['name']=='corr_array_assertion_checks');assert assertion.count('throw std::runtime_error')==2
expected=[m.group(1) for m in re.finditer(r'error stop "([^"]*)"',raw('corr_array_assertion_checks'))]
for message in expected:assert json.dumps(message) in assertion
assert assertion.index('clubb_at_least_debug_level_api(1)')<assertion.index('clubb_at_least_debug_level_api(2)')
assert 'constexpr real_t tol = real_t(1.e-6);' in assertion
# Kernel captures and explicit function returns; public adapters use the actual owner.
kernels=[]
for m in re.finditer(r'KOKKOS_LAMBDA\([^)]*\)\s*\{',code):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 body=code[m.end():end-1];kernels.append(body);assert not re.search(r'f2k::|std::|parallel_for|parallel_reduce|to_host_space|\.err_',body)
assert len(kernels)==14
assert 'Kokkos::fma' not in code and 'Kokkos::pow' not in code and 'thread_local' not in code
assert len(defs)==14 #10Native bodies +4thin public array adapters.
for d in defs[10:]:
 assert 'to_execution_space' in d['body'] and 'from_execution_space' in d['body']
 assert d['name'] in ['diagnose_correlations','calc_cholesky_corr_mtx_approx','setup_corr_cholesky_mtx','cholesky_to_corr_mtx_approx']
assert max(map(len,text.splitlines()))<=100 and all(l.rstrip()==l and '\t' not in l for l in text.splitlines())
record={'canonical_source':str(NATIVE.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),'target_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),'routines':reviews,'ordered_call_trees':callcount,'independent_kernels':len(kernels),'manual_review':['Full1037Native and272JAX lines read before outline. Ten active routines implemented in order; four wholly commented-out routines excluded.','Exact Native1based iiPDF_w, lower/upper physical matrix slots, reversed slices and all per-call copy ownership preserved. Each Native slice assignment is a separate queue stage.','Two per-entry scalar-product recurrences keep Native ascending order; s entries are written before every read. MATMUL is an ordered explicit contraction, with compiler arithmetic equality pending.','Native square roots and conditional clamps retained; JAX gradient floors/symmetrization/boolean-only diagnostics not imported.','Supported false path only: explicit API rejection of true, whose Native working array is uninitialized. No invented initialization.','Native debug1off-diagonal then debug2diagonal fatal exceptions preserve strings/check order/tolerance; no finite test or newerr_info/stats object.','Four generic array adapters use actual existing array owner and preserve caller output allocations. All concrete kernels use default RealView execution space; no owning UDT/global cache/device bridge.','Separate cleanup reviews signature/intent groups, complete scientific comments and disabled alternative formulas, reference citations,core-kind values,case,100columns and whitespace.'],'remaining_dependencies':'Existing host Native debug query only. No new lower scientific routine, type, accessor or physical constant.','untested':'All runtime branches/CPU/GPU/model behavior, ordered MATMUL/compiler contraction, assertions,swapped layouts and metadata boundary callers uncompiled/untested; no BFB claim/speculative workaround.'};(HERE/'SOURCE_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'routines':reviews,'calls':callcount,'kernels':len(kernels)}))
