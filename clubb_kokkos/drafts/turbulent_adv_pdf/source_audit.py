"""Compare Native contracts and expression trees without executing physics."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements,declaration
from generate_derived_types import split_top
reader=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();reader=reader[reader.index('class Cpp:'):reader.index('reviews=[]')]
reader=reader.replace("   if name in {'max','min'}:return ('if',('bin','>' if name=='max' else '<',args[0],args[1]),args[0],args[1])", "   if name in {'max','min'}:return ('call',name,args)")
exec(reader)
NATIVE=ROOT/'src/CLUBB_core/turbulent_adv_pdf.F90';CPP=HERE/'turbulent_adv_pdf.cpp';HPP=HERE/'turbulent_adv_pdf.hpp'
source=NATIVE.read_text();text=CPP.read_text();header=HPP.read_text();code=re.sub(r'//[^\n]*','',text);hcode=re.sub(r'//[^\n]*','',header)
outline=json.loads((HERE/'SOURCE_OUTLINE.json').read_text());names=outline['native_inventory_and_order'];defs=[]
for m in re.finditer(r'void\s+(\w+)\s*\(([^{};]*)\)\s*\{',code,re.S):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 defs.append({'name':m.group(1),'args':[a.strip().split()[-1] for a in split_top(m.group(2))],'body':code[m.end():end-1]})
assert [d['name'] for d in defs]==names
aliases={'invrs_dzm':'gr_invrs_dzm','invrs_dzt':'gr_invrs_dzt','weights_zm2zt':'gr_weights_zm2zt','grid_dir':'gr_grid_dir'}
def canonical(n):
 if not isinstance(n,tuple):return n
 if n[0]=='v':return ('v',aliases.get(n[1],n[1]))
 if n[0]=='call':return ('call',aliases.get(n[1],n[1]),tuple(canonical(a) for a in n[2]))
 return tuple(canonical(x) if isinstance(x,tuple) else x for x in n)
reviews=[]
for d in defs:
 name=d['name'];proc=read_procedure(NATIVE,name);args=[a.name for a in proc.arguments];assert d['args']==args
 raw=re.search(r'^  subroutine '+name+r'\(.*?^  end subroutine '+name+r'\b',source,re.S|re.M).group();decls={}
 for l in statements(raw):decls.update(declaration(l))
 arrays={a.name for a in proc.arguments if a.dims}|{n for n,a in decls.items() if a.dims}|{'invrs_dzm','invrs_dzt','weights_zm2zt'}
 expected=[]
 for l in statements(raw):
  m=re.match(r'^(\w+\([^()]*\))\s*=(?!=)\s*(.*)$',l)
  if m:expected.append((fexpr(m.group(1)),fexpr(m.group(2))))
 actual=[]
 for m in re.finditer(r'(?m)^\s*(\w+\([^;{}]*?\))\s*=(?!=)\s*([^;]*);',d['body']):actual.append((canonical(Cpp(m.group(1),arrays).parse()),canonical(Cpp(m.group(2),arrays).parse())))
 assert expected==actual,(name,sorted(arrays),len(expected),len(actual),next((i for i,(a,b) in enumerate(zip(expected,actual)) if a!=b),None))
 expectedifs=[fexpr(m.group(1)) for l in statements(raw) if (m:=re.match(r'if\s*\((.*)\)\s*then$',l,re.I))];actualifs=[]
 for m in re.finditer(r'\bif\s*\(',d['body']):
  end=m.end();depth=1
  while depth:depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  actualifs.append(canonical(Cpp(d['body'][m.end():end-1],arrays).parse()))
 assert expectedifs==actualifs,(name,'branches')
 declarations=re.findall(r'void\s+'+name+r'\s*\(([^{};]*)\)\s*[;{]',hcode,re.S);assert len(declarations)==2
 for declaration_ in declarations:assert [a.strip().split()[-1] for a in split_top(declaration_)]==args
 # Generic adapter must stage every input, skip Out value reads, call Native order,
 # and restore exactly the original caller-owned output.
 match=re.search(r'void\s+'+name+r'\s*\([^{};]*\)\s*\{(.*?)\n\}',hcode,re.S);body=match.group(1)
 expectedstage=[a for a in proc.arguments if a.dims]
 for a in expectedstage:assert f'to_execution_space({a.name}'+(', false' if a.intent=='out' else '')+')' in body
 call=re.search(r'\b'+name+r'\s*\(([^;]*)\);',body,re.S)
 assert [a.strip() for a in split_top(call.group(1))]==[a.name+('_exec' if a.dims else '') for a in proc.arguments]
 output=proc.arguments[-1].name;assert re.findall(r'from_execution_space\s*\((.*?)\);',body,re.S)==[output+', '+output+'_exec']
 # All flattened kernels assign the Native independent column/storage indices.
 kernels=[]
 for m in re.finditer(r'KOKKOS_LAMBDA\(int cell\)\s*\{',d['body']):
  end=m.end();depth=1
  while depth:depth+=(d['body'][end]=='{')-(d['body'][end]=='}');end+=1
  body_=d['body'][m.end():end-1];kernels.append(body_)
  assert not re.search(r'\bgr\.|f2k::|to_host_space|parallel_for|parallel_reduce',body_)
  assert 'const int i = cell' in body_
 assert len(kernels)==len(re.findall(r'!\$acc parallel loop',raw))
 reviews.append({'name':name,'arguments':len(args),'ordered_indexed_assignment_trees':len(expected),'ordered_if_predicates':len(expectedifs),'independent_kernels':len(kernels),'header_contracts_and_generic_call_order':'exact'})
for t in [text,header]:assert max(map(len,t.splitlines()))<=100 and all(l.rstrip()==l and '\t' not in l for l in t.splitlines())
assert 'from_level' not in text and 'Optional Arguements' not in text
record={'canonical_source':str(NATIVE.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),'cpp_sha256':hashlib.sha256(CPP.read_bytes()).hexdigest(),'hpp_sha256':hashlib.sha256(HPP.read_bytes()).hexdigest(),'routines':reviews,'separate_manual_review':['Full four Native routines and JAX scope read; Native scientific descriptions and executable comments retained in place except documented stale/incorrect comments.','Native band-first operator output, super/main/sub categories and all two endpoints checked; no transpose.','Grid weight categories1:2 are upper/lower slots, not physical offsets. View conversion is the only index-base change.','Flattened interior k=2+cell/ngrdcol,i=1+cell%ngrdcol covers Native2:nzm-1 exactly for either grid direction. Boundary flattening covers all three bands at both storage endpoints.','Upwind >zero branch, zero else, Godunov min/max ordering and thermodynamic sign profile are unchanged.','Concrete execution-space Views pair with execution-space grid components; generic adapter uses existing array owner and retains caller allocations. No UDT/device bridge/host IO capture.','Separate cleanup: public/definition continuation groups, argument units/shapes, scientific descriptions, Native block comments/dividers, case-sensitive calls, semicolon placement,100columns and whitespace.'],'numerical_evidence':'Pending: no compilation, tests, model execution or BFB claim. No deliberate BFB-specific arithmetic adaptation.'}
(HERE/'SOURCE_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(reviews))
