"""Source-only contract/formula audit; no compiled or numerical evidence."""
from pathlib import Path
import ast,re,json,hashlib,sys
from dataclasses import asdict
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements,declaration
from generate_derived_types import split_top
# Reuse the independent expression reader, not any source preparation script.
utility=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text()
part=utility[utility.index('class Cpp:'):utility.index('reviews=[]')]
part=part.replace("   if isinstance(n.op,ast.Pow):\n    assert b in [('n',2.),('n',3.)]", "   if isinstance(n.op,ast.Pow):\n    if b not in [('n',2.),('n',3.)]:return ('call','pow',(a,b))")
part=part.replace("   if name in {'max','min'}:return ('if',('bin','>' if name=='max' else '<',args[0],args[1]),args[0],args[1])", "   if name in {'max','min'}:\n    out=args[0]\n    for value in args[1:]:out=('if',('bin','>' if name=='max' else '<',out,value),out,value)\n    return out")
exec(part)
native=ROOT/'src/CLUBB_core/sfc_varnce_module.F90';target=HERE/'sfc_varnce_module.cpp';source=native.read_text();text=target.read_text();code=re.sub(r'//[^\n]*','',text)
p=read_procedure(native,'calc_sfc_varnce');signature=re.search(r'void\s+calc_sfc_varnce\s*\(([^{};]*)\)\s*\{',code,re.S)
assert [a.strip().split()[-1].lstrip('&*') for a in split_top(signature.group(1))]==[a.name for a in p.arguments]
assert 'inline void calc_sfc_varnce' not in text
body=code[signature.end():];decls={}
for line in statements(source):decls.update(declaration(line))
arrays={n for n,a in decls.items() if a.dims};arrays.add('gr_zm')
expected=[]
for line in statements(source):
 m=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',line)
 if m:expected.append((fexpr(m.group(1)),fexpr(m.group(2))))
actual=[]
for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',body):
 if re.search(r'_sf_(?:begin|final)\(',m.group(1)):continue
 actual.append((Cpp(m.group(1),arrays).parse(),Cpp(m.group(2),arrays).parse()))
assert expected==actual,[(i,a,b) for i,(a,b) in enumerate(zip(expected,actual)) if a!=b][:1]
# Compare every science branch, including the dormant Andre1978 body. The final
# diagnostic ANY/status branches are explicit host-language adaptations.
native_conditions=[]
for line in statements(source):
 m=re.match(r'(?:if|elseif|else if)\s*\((.*)\)\s*then$',line,re.I)
 if m:
  value=m.group(1)
  if 'err_info%' not in value:native_conditions.append(fexpr(value))
actual_conditions=[]
for m in re.finditer(r'\bif\s*\(',body):
 end=m.end();pos=end;depth=1
 while depth:depth+=(body[end]=='(')-(body[end]==')');end+=1
 value=body[pos:end-1]
 if 'err_info.' in value:continue
 value=value.replace('stats.flags().l_sample','stats_l_sample')
 actual_conditions.append(Cpp(value,arrays).parse())
assert native_conditions==actual_conditions
# Stats expression materialization must be exact and preserve six begin/six
# finalize call order/field-specific mapping, rather than merely label counts.
science_fields=['thlp2','rtp2','rtpthlp','up2','vp2','wp2']
for method,phase in [('stats_begin_budget','begin'),('stats_finalize_budget','final')]:
 found=re.findall(r'f2k::'+method+r'\("(\w+)",\s*(\w+),\s*stats\)',body)
 assert found==[(n+'_sf',n+'_sf_'+phase) for n in science_fields]
 for name in science_fields:
  assert re.search(name+'_sf_'+phase+r'\(column, level\)\s*=\s*'+name+r'\(column, level\)\s*/\s*dt;',body)
# Device ownership is checked independently from the Native arithmetic trees.
kernels=[]
for m in re.finditer(r'KOKKOS_LAMBDA\s*\([^)]*\)\s*\{',text):
 pos=m.end();end=pos;depth=1
 while depth:depth+=(text[end]=='{')-(text[end]=='}');end+=1
 kernel=text[pos:end-1]
 assert not re.search(r'\b(?:gr|stats|sclr_idx|err_info)\s*(?:\.|->)|f2k::',kernel)
 kernels.append(text[:m.start()].count('\n')+1)
contract=read_procedure(ROOT/'src/CLUBB_core/numerical_check.F90','sfc_varnce_check')
record={'canonical_source':str(native.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(native.read_bytes()).hexdigest(),'target':str(target.relative_to(ROOT)),'target_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'native_signature_arguments':len(p.arguments),'scientific_assignment_trees':len(expected),'ordered_branch_predicates':len(native_conditions),'ordered_sampled_budget_calls':12,'device_bodies_without_owning_udt_or_f2k_capture':len(kernels),'required_public_native_checker_contract':asdict(contract),'status':'Source-only audit. No build or physics execution; CPU-vs-CPU BFB evidence pending.'}
(HERE/'expression_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print('Surface source audit:',len(p.arguments),'arguments;',len(expected),'assignment trees;',len(native_conditions),'branch predicates;',len(kernels),'device bodies.')
