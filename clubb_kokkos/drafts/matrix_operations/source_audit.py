"""Review Native matrix contracts, formulas and lower boundaries; no execution."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
r=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();exec(r[r.index('class Cpp:'):r.index('reviews=[]')])
NATIVE=ROOT/'src/CLUBB_core/matrix_operations.F90';TARGET=HERE/'matrix_operations.hpp';source=NATIVE.read_text();text=TARGET.read_text();code=re.sub(r'//[^\n]*','',text);names=['symm_covar_matrix_2_corr_matrix','Cholesky_factor','Symm_matrix_eigenvalues','mirror_lower_triangular_matrix'];defs=[]
for m in re.finditer(r'void\s+(\w+)\s*\(([^{};]*)\)\s*\{',code,re.S):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 defs.append({'name':m.group(1),'args':[a.strip().split()[-1] for a in split_top(m.group(2))],'body':code[m.end():end-1]})
assert [d['name'] for d in defs]==names
aliases={n+'_host':n for n in ['a_input','a_scaling','a_Cholesky','a_eigenvalues']};aliases={k.lower():v.lower() for k,v in aliases.items()}
arrays={'covar','corr','a_cholesky_host','a_scaling_host','a_input_host','a_eigenvalues','a_corr','a_scaled','a_scratch','work','matrix'}
def canonical(n):
 if not isinstance(n,tuple):return n
 if n[0]=='v':return ('v',aliases.get(n[1],n[1]))
 if n[0]=='str':return ('str',n[1].lower())
 if n[0]=='call':
  if n[1]=='real':return canonical(n[2][0])
  return ('call',aliases.get(n[1],n[1]).replace('matrix_operations_detail::',''),tuple(canonical(a) for a in n[2]))
 return tuple(canonical(x) if isinstance(x,tuple) else x for x in n)
def cpp_expr(s):return canonical(Cpp(re.sub(r"'([^']*)'",lambda m:json.dumps(m.group(1)),s),arrays).parse())
def raw(n):return re.search(r'^  subroutine '+n+r'\(.*?^  end subroutine '+n+r'\b',source,re.S|re.M).group()
reviews=[]
for d in defs:
 n=d['name'];p=read_procedure(NATIVE,n);assert d['args']==[a.name for a in p.arguments]
 expected=[]
 for l in statements(raw(n)):
  if m:=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',l):expected.append((canonical(fexpr(m.group(1))),canonical(fexpr(m.group(2)))))
  elif m:=re.match(r'if\s*\(.*\)\s+(d_smallest)\s*=\s*(.*)$',l,re.I):expected.append((fexpr(m.group(1)),fexpr(m.group(2))))
 actual=[];entries=[]
 for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',d['body']):entries.append((m.start(),cpp_expr(m.group(1)),cpp_expr(m.group(2))))
 for m in re.finditer(r'Kokkos::deep_copy\s*\(([^;]*)\);',d['body'],re.S):
  a=split_top(m.group(1));entries.append((m.start(),cpp_expr(a[0]),cpp_expr(a[1])))
 if n=='Symm_matrix_eigenvalues':
  assert 'auto a_scratch = to_host_space(a_input);' in d['body'];actual=[(('v','a_scratch'),('v','a_input'))]
 else:actual=[(a,b) for _,a,b in sorted(entries)]
 assert expected==actual,(n,'assignments',expected,actual)
 expectedcalls=[canonical(fexpr(re.sub(r'^call\s+','',l))) for l in statements(raw(n)) if re.match(r'call\s+',l,re.I)];actualcalls=[]
 for m in re.finditer(r'(?:matrix_operations_detail::)?(?:lapack_poequ|lapack_laqsy|lapack_potrf|lapack_syev|Symm_matrix_eigenvalues|symm_covar_matrix_2_corr_matrix)\s*\(',d['body']):
  end=m.end();depth=1
  while depth:depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  t=d['body'][m.start():end];t=re.sub(r'(\w+)\.data\(\)',r'\1',t);actualcalls.append(cpp_expr(t))
 assert expectedcalls==actualcalls,(n,'calls',expectedcalls,actualcalls)
 reviews.append({'name':n,'required_arguments':len(p.arguments),'ordered_assignment_trees':len(expected),'ordered_lower_calls':len(actualcalls)})
assert 'constexpr int itermax = 10;' in code and 'constexpr real_t d_coef = real_t(0.1);' in code
assert 'constexpr int lwork = 180;' in code and 'for (int iter = 1; iter <= itermax; ++iter)' in code
assert 'if (info < 0)' in defs[1]['body'] and 'throw std::runtime_error' in defs[1]['body'] and 'break;' in defs[1]['body']
assert defs[1]['body'].index('if (iter == itermax)')<defs[1]['body'].index('d_smallest =')
assert 'from_execution_space(a_scaling, a_scaling_host);' in code and 'from_execution_space(a_Cholesky, a_Cholesky_host);' in code
# Full failure/success diagnostics and native matrix/vector field payloads.
for label in ['a_factored (worked)=','Cholesky_factor: leading minor of order ',' is not positive definite.','factorization failed.','a_input=','a_Cholesky=','a_eigenvalues=','a_correlations=','iteration =','itermax =','Fatal error in Cholesky_factor','Attempting to modify matrix to allow factorization.','a_modified eigenvalues=','Symm_matrix_eigenvalues: Algorithm failed to converge.']:assert label in text,label
assert code.count('KOKKOS_LAMBDA')==2
for m in re.finditer(r'KOKKOS_LAMBDA\(int cell\)\s*\{',code):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 assert not re.search(r'f2k::|std::cerr|to_host_space|lapack_',code[m.end():end-1])
assert 'using Policy = Kokkos::RangePolicy<typename Correlation::execution_space>;' in code
assert 'using Policy = Kokkos::RangePolicy<typename Matrix::execution_space>;' in code
assert 'if (nvars > 1)' in code and 'if (col >= 2 && row < col)' in code
for p in [TARGET,HERE/'lapack_dependencies.hpp']:
 t=p.read_text();assert max(map(len,t.splitlines()))<=100 and all(l.rstrip()==l and '\t' not in l for l in t.splitlines())
rec={'canonical_source':str(NATIVE.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [TARGET,HERE/'lapack_dependencies.hpp']},'routines':reviews,'excluded_unused_public_JAX_routines':['row_mult_lower_tri_matrix','set_lower_triangular_matrix','get_lower_triangular_matrix','print_lower_triangular_matrix'],'manual_review':['Full565line Native and173line actualJAX scope read before outline; root approved four-routine supported scope.','Native scaling/character flags/info branches, ten attempts total, retry restore/minimum scan/tau association and final-failed-iteration modification retained.','Complete debug1/2 matrix/eigen/correlation diagnostics with field-specific payloads; no JAX retry-count/trace omissions imported.','Private eigen helper preserves fixed lwork180 and returns after diagnostics for either LAPACK status failure; source-owned scratch copies remain per call.','Raw buffers cross only public lapack_interfaces boundary; transient ABI extension is pending in existing shim owner, no raw Fortran symbols in science.','Two simple kernels retain matrix input/output memory space, unchanged lower/diagonal ownership and explicit void enclosing types. No mutable cache/global state/device bridge.','Separate cleanup covers Native arguments/order/case, core-kind constants, scientific comments, complete payloads,100columns/whitespace.'],'remaining_lower':'Four typed LAPACK entry extensions in existing core/shim owners: potrf,poequ,laqsy,syev; existing host Native debug query. No new physical constant.','limits':'Diagnostic numbers use local width/precision streams; exact compiler-specific G/E formatting is unverified. Uncompiled/untested; no runtime/BFB claim or speculative arithmetic workaround.'};(HERE/'SOURCE_AUDIT.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(reviews))
