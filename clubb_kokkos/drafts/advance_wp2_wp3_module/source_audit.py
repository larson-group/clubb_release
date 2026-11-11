"""Audit Native/draft source structure, not compiled or numerical behavior."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,declaration
from generate_derived_types import split_top
utility=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text()
part=utility[utility.index('class Cpp:'):utility.index('reviews=[]')]
part=part.replace("   if isinstance(n.op,ast.Pow):\n    assert b in [('n',2.),('n',3.)]", "   if isinstance(n.op,ast.Pow):\n    if b==('n',4.):return ('bin','*',('bin','*',a,a),('bin','*',a,a))\n    if b not in [('n',2.),('n',3.)]:return ('call','pow',(a,b))")
part=part.replace("   if name in {'max','min'}:return ('if',('bin','>' if name=='max' else '<',args[0],args[1]),args[0],args[1])", "   if name in {'max','min'}:\n    out=args[0]\n    for value in args[1:]:out=('if',('bin','>' if name=='max' else '<',out,value),out,value)\n    return out")
exec(part)
NATIVE=ROOT/'src/CLUBB_core/advance_wp2_wp3_module.F90';TARGET=HERE/'advance_wp2_wp3_module.cpp'
source=NATIVE.read_text();text=TARGET.read_text();code=re.sub(r'//[^\n]*','',text);outline=json.loads((HERE/'SOURCE_OUTLINE.json').read_text());defs=[]
for m in re.finditer(r'\bvoid\s+(\w+)\s*\(([^{};]*?)\)\s*\{',code,re.S):
 d=1;end=m.end()
 while d:d+=(code[end]=='{')-(code[end]=='}');end+=1
 defs.append({'name':m.group(1),'args':[a.strip().split()[-1].lstrip('&*') for a in split_top(m.group(2))],'body':code[m.end():end-1]})
assert [d['name'] for d in defs]==outline['native_inventory_and_order']
reviews=[]
for d,p in zip(defs,outline['routines']):
 name=d['name'];assert d['args']==[a['name'] for a in p['arguments']]
 raw=re.search(r'(?ms)^  subroutine '+name+r'\(.*?^  end subroutine '+name+r'\b',source).group()
 if name in {'advance_wp2_wp3','wp23_solve','wp23_rhs'}:
  reviews.append({'name':name,'arguments':len(d['args']),'review':'Manual full parent/solver/RHS staging, call/field/error/stats review; adapters tracked separately'});continue
 decls={}
 for line in statements(raw):decls.update(declaration(line))
 arrays={k for k,a in decls.items() if a.dims};arrays.update(re.findall(r'\b(gr_\w+)\s*\(',d['body']))
 expected=[]
 for line in statements(raw):
  m=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',line)
  if not m:continue
  lhs,rhs=m.groups()
  if ':' in lhs:lhs=lhs.replace(':','band+1')
  expected.append((fexpr(lhs),fexpr(rhs)))
 body=d['body']
 for name2 in arrays:body=re.sub(r'\b'+name2+r'\(band,',name2+'((band+1)-1,',body,flags=re.I)
 actual=[]
 for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',body):
  actual.append((Cpp(m.group(1),arrays).parse(),Cpp(m.group(2),arrays).parse()))
 # Whole-band Native fills become explicit band loops in a device work item.
 for m in re.finditer(r'for\s*\(int band=0;[^;]*;\+\+band\)\s*(\w+\([^;]*?\))\s*=\s*([^;]*);',body):
  actual.insert(0,(Cpp(m.group(1),arrays).parse(),Cpp(m.group(2),arrays).parse()))
 if expected!=actual:
  diffs=[(i,a,b) for i,(a,b) in enumerate(zip(expected,actual)) if a!=b]
  raise ValueError((name,'assignment trees',len(expected),len(actual),diffs[:1]))
 native_conditions=[]
 for line in statements(raw):
  m=re.match(r'(?:if|elseif|else if)\s*\((.*)\)\s*then$',line,re.I)
  if m:native_conditions.append(fexpr(m.group(1)))
 actual_conditions=[]
 for m in re.finditer(r'\bif\s*\(',body):
  pos=m.end();end=pos;depth=1
  while depth:depth+=(body[end]=='(')-(body[end]==')');end+=1
  actual_conditions.append(Cpp(body[pos:end-1],arrays).parse())
 assert native_conditions==actual_conditions,(name,'branches')
 reviews.append({'name':name,'arguments':len(d['args']),'ordered_science_assignment_trees':len(expected),'ordered_branch_predicates':len(native_conditions)})
kernels=[]
for m in re.finditer(r'KOKKOS_LAMBDA\s*\([^)]*\)\s*\{',text):
 pos=m.end();end=pos;depth=1
 while depth:depth+=(text[end]=='{')-(text[end]=='}');end+=1
 body=text[pos:end-1]
 assert not re.search(r'\b(?:gr|stats|pdf_implicit_coefs_terms|nu_vert_res_dep|err_info)\s*(?:\.|->)|f2k::',body),(text[:m.start()].count('\n')+1,body)
 # Mutable work-item scalars must be assigned before their first read. Values
 # selected on the host must be captured, not redeclared as zero locals.
 body_without_comments=re.sub(r'//[^\n]*','',body)
 for local in re.finditer(r'\b(?:int|real_t|bool)\s+(\w+)\{\};',body_without_comments):
  name=local.group(1);tail=body_without_comments[local.end():]
  first=re.search(r'\b'+re.escape(name)+r'\b',tail)
  assert first and re.match(r'\s*=(?!=)',tail[first.end():]),(name,'kernel local read before write')
 kernels.append(text[:m.start()].count('\n')+1)
native_names={name.lower():name for name in outline['native_inventory_and_order']}
for name in re.findall(r'\b(\w+)\s*\(',code):
 if name.lower() in native_names:assert name==native_names[name.lower()],(name,'C++ symbol case')
a={'canonical_source':str(NATIVE.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),'target':str(TARGET.relative_to(ROOT)),'target_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),'routines':reviews,'kernels_without_owning_udt_or_f2k_capture':len(kernels),'kernel_local_first_use_and_exact_internal_symbol_case_checked':True,'status':'Source audit only. Compilation and maintained CPU-vs-CPU model comparison pending.'}
(HERE/'expression_audit.json').write_text(json.dumps(a,indent=2)+'\n')
print('Checked23 Native signatures/order and',sum(r.get('ordered_science_assignment_trees',0) for r in reviews),'ordered assignment trees;',len(kernels),'device bodies')
