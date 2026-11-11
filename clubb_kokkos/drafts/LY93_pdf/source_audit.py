"""Compare LY93 source contracts and formula trees; do not execute physics."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,read_procedure,declaration
from generate_derived_types import split_top
p=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();part=p[p.index('class Cpp:'):p.index('reviews=[]')]
part=part.replace("    assert b in [('n',2.),('n',3.)]", "    if b not in [('n',2.),('n',3.)]:return ('call','pow',(a,b))")
part=part.replace("   if name in {'max','min'}:return ('if',('bin','>' if name=='max' else '<',args[0],args[1]),args[0],args[1])", "   if name in {'max','min'}:return ('call',name,args)")
exec(part)
native=ROOT/'src/CLUBB_core/LY93_pdf.F90';target=HERE/'LY93_pdf.hpp';source=native.read_text();text=target.read_text();code=re.sub(r'//[^\n]*','',text);outline=json.loads((HERE/'SOURCE_OUTLINE.json').read_text());defs=[]
for name in outline['native_inventory_and_order']:
 m=re.search(r'\b'+name+r'\s*\(([^{};]*?)\)\s*\{',code,re.S);assert m,name;end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 defs.append({'name':name,'args':[a.strip().split()[-1] for a in split_top(m.group(1))],'body':code[m.end():end-1],'position':m.start()})
assert [d['name'] for d in sorted(defs,key=lambda d:d['position'])]==outline['native_inventory_and_order']
def scalarize(n,fields):
 if not isinstance(n,tuple):return n
 if n[0]=='v' and n[1] in fields:return ('call',n[1],(('v','k'),))
 return tuple(scalarize(x,fields) if isinstance(x,tuple) else x for x in n)
def canonical(n):
 if not isinstance(n,tuple):return n
 if n[0]=='call':
  args=tuple(canonical(a) for a in n[2])
  if n[1]=='max' and args[0][0:2]==('call','max'):return ('call','max',args[0][2]+args[1:])
  return ('call',n[1],args)
 return tuple(canonical(x) if isinstance(x,tuple) else x for x in n)
reviews=[]
for d in defs:
 name=d['name'];p=read_procedure(native,name);assert d['args']==[a.name for a in p.arguments]
 raw=re.search(r'(?ms)^  (?:subroutine|function) '+name+r'\(.*?^  end (?:subroutine|function) '+name+r'\b',source).group();decls={}
 for l in statements(raw):decls.update(declaration(l))
 arrays={n for n,a in decls.items() if a.dims};scalar_fields=arrays-{'sgn_skx','b_x'};expected=[]
 for line in statements(raw):
  m=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',line)
  if not m:continue
  lhs,rhs=m.groups();expected.append((canonical(scalarize(fexpr(lhs),scalar_fields)),canonical(scalarize(fexpr(rhs),scalar_fields))))
 actual=[];entries=[]
 for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',d['body']):entries.append((m.start(),m.group(1),m.group(2)))
 for m in re.finditer(r'Kokkos::deep_copy\s*\(',d['body']):
  end=m.end();depth=1
  while depth:depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  args=split_top(d['body'][m.end():end-1]);entries.append((m.start(),*args))
 for _,lhs,rhs in sorted(entries):
  # Whole-profile result assignment retains an owned returned View and caller output.
  lhsnode=Cpp(lhs,set()).parse();rhsnode=Cpp(rhs,set()).parse()
  if name=='LY93_driver' and lhsnode==('v','mixt_frac'):lhsnode=scalarize(lhsnode,scalar_fields)
  actual.append((canonical(lhsnode),canonical(rhsnode)))
 # The returned mixture API consumes the whole Sk_max profile rather than one cell.
 if name=='LY93_driver':
  expected[1]=(expected[1][0],('call','calc_mixt_frac_ly93',(('v','nz'),('v','sk_max'))))
 assert expected==actual,(name,'assignments',[(i,a,b) for i,(a,b) in enumerate(zip(expected,actual)) if a!=b],len(expected),len(actual))
 reviews.append({'name':name,'arguments':len(p.arguments),'ordered_assignment_trees':len(expected)})
# All three Native private-free public calls retain exact field order.
d=defs[0];expected=[fexpr(re.sub(r'^call\s+','',l,flags=re.I)) for l in statements(re.search(r'(?ms)^  subroutine LY93_driver\(.*?^  end subroutine LY93_driver\b',source).group()) if re.match(r'call\s+',l,re.I)];actual=[]
for m in re.finditer(r'calc_params_LY93\s*\(',d['body']):
 end=m.end();depth=1
 while depth:depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
 actual.append(Cpp(d['body'][m.start():end],set()).parse())
assert expected==actual
for m in re.finditer(r'KOKKOS_LAMBDA\([^)]*\)\s*\{',code):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 assert not re.search(r'f2k::|to_host_space|std::cerr|parallel_for|parallel_reduce',code[m.end():end-1])
assert max(map(len,text.splitlines()))<=100 and all(l.rstrip()==l and '\t' not in l for l in text.splitlines())
a={'canonical_source':str(native.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(native.read_bytes()).hexdigest(),'target_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'routines':reviews,'ordered_parameter_calls':len(actual),'manual_review':'Bisection bounds/uncapped termination/tolerance, WHERE sign branch inclsignedzero, profile0-based execution mapping, Nativearray scratch lowered to work-item scalars, return ownership/memoryspace, Nativegrouping/comments/citations/units; no scopes/options omitted.','runtime':'No builds/tests/numerical evidence.'}
(HERE/'expression_audit.json').write_text(json.dumps(a,indent=2)+'\n');print('3signatures,',sum(r['ordered_assignment_trees'] for r in reviews),'assignmenttrees,3field-specificcalls,3safe kernels;source only')
