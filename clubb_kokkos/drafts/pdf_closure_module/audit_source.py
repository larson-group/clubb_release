from pathlib import Path
import ast,copy,hashlib,json,re,sys
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import find_procedure,statements
from generate_derived_types import split_top
names=['calc_wp4_pdf','calc_wp2xp_pdf','calc_wpxp2_pdf','calc_wpxpyp_pdf']
source=(root/'src/CLUBB_core/pdf_closure_module.F90').read_text()
target=(folder/'pdf_closure_module.cpp').read_text()
clean=lambda s:re.sub(r'//[^\n]*|/\*.*?\*/','',s,flags=re.S)
cpp=clean(target)
ms=list(re.finditer(r'void\s+(\w+)\s*\(([^)]*)\)\s*\{',cpp));assert [m[1] for m in ms]==names
class Powers(ast.NodeTransformer):
 def visit_BinOp(self,node):
  node=self.generic_visit(node)
  if isinstance(node.op,ast.Pow):
   assert isinstance(node.right,ast.Constant) and node.right.value in (2,4)
   square=ast.BinOp(left=copy.deepcopy(node.left),op=ast.Mult(),right=copy.deepcopy(node.left))
   return square if node.right.value==2 else ast.BinOp(left=copy.deepcopy(square),op=ast.Mult(),right=copy.deepcopy(square))
  return node
counts={}
for m in ms:
 name=m[1];native=re.search(r'^  subroutine '+name+r'\(.*?^  end subroutine '+name+r'\b',source,re.M|re.S).group()
 args=[a.strip().split()[-1] for a in split_top(m[2])];assert args==[a.name for a in find_procedure(name).arguments]
 output=args[-1]
 ns=next(s for s in statements(native) if re.match(r'\s*'+output+r'\s*\(',s) and '=' in s).split('=',1)[1]
 cs=re.search(r'\b'+output+r'\s*\(\s*i\s*,\s*k\s*\)\s*=(.*?);',cpp,re.S)[1].replace('Kokkos::sqrt','sqrt')
 a=Powers().visit(ast.parse(ns.strip(),mode='eval'));b=ast.parse(re.sub(r'\s+',' ',cs).strip(),mode='eval')
 assert ast.dump(a)==ast.dump(b),name
 counts[name]=len(args)
for n in ['pdf_closure_module.cpp','pdf_closure_module.hpp']:
 s=(folder/n).read_text();assert all(len(l)<=100 and l.rstrip()==l and '\t' not in l for l in s.splitlines())
 code=re.sub(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'','',clean(s))
 for op,cl in [('(',')'),('{','}'),('[',']')]:
  d=0
  for ch in code:d+=(ch==op)-(ch==cl);assert d>=0,(n,op)
  assert d==0,(n,op,d)
assert target.count('KOKKOS_LAMBDA')==4
assert 'f2k::' not in target and 'Kokkos::fma' not in target and 'pow(' not in cpp
record={'status':'source-audited draft, NOT compiled/integrated/runtime tested','source':'src/CLUBB_core/pdf_closure_module.F90:1828–2347 (four selected public routines)','native_order':names,'argument_counts':counts,'all_signatures_and_complete_assignment_trees_match':True,'all_integer_powers_expanded_with_surrounding_association_intact':True,'device_kernels':4,'scientific_lower_calls':[],'manual_review':['Full four Native procedures read before drafting; preimplementation outline retained.','Complete equation/physical-parameter descriptions, references blocks and all input/output units retained.','Each formula matches both component means/variances/correlations and mixture factors in Native order, with no added floors, limits, aliases, temporaries or arithmetic fixes.','Native final arrays overwritten in every cell; no branch/stats/error behavior exists in these four routines.','Native wp4/wpxpyp missing INTENT is not used to infer optional or omitted outputs; complete bodies assign them.','Flat (column,level) loops, explicit scalar indices and only required View captures.','Public definitions must be interleaved with private anonymous scopes in production to preserve full Native order; documented in integration notes.','Separate cleanup simplified redundant power-expression parentheses and argument-comment padding; formatter and syntax-token checks pass.'],'numerical_limits':'Ordinary C++ expressions only; no BFB or GPU claim until maintained model integration.','hashes':{n:hashlib.sha256((folder/n).read_bytes()).hexdigest() for n in ['pdf_closure_module.cpp','pdf_closure_module.hpp']},'native_sha256':hashlib.sha256((root/'src/CLUBB_core/pdf_closure_module.F90').read_bytes()).hexdigest()}
(folder/'source_audit.json').write_text(json.dumps(record,indent=2)+'\n');print({k:record[k] for k in ['status','argument_counts','device_kernels']})
