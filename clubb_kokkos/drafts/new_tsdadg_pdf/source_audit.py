"""Audit Native/TSDADG contracts, expressions, calls and diagnostics without running physics."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,read_procedure
from generate_derived_types import split_top
util=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();part=util[util.index('class Cpp:'):util.index('reviews=[]')]
part=part.replace("   if name in {'max','min'}:return ('if',('bin','>' if name=='max' else '<',args[0],args[1]),args[0],args[1])", "   if name in {'max','min'}:return ('call',name,args)");exec(part)
native=ROOT/'src/CLUBB_core/new_tsdadg_pdf.F90';target=HERE/'new_tsdadg_pdf.hpp';source=native.read_text();text=target.read_text();code=re.sub(r'//[^\n]*',lambda m:' '*len(m.group()),text);outline=json.loads((HERE/'SOURCE_OUTLINE.json').read_text());defs=[]
for name in outline['native_inventory_and_order']:
 m=re.search(r'\bvoid '+name+r'\s*\(([^{};]*?)\)\s*\{',code,re.S);assert m,name;end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 defs.append({'name':name,'args':[a.strip().split()[-1].lstrip('&*') for a in split_top(m.group(1))],'body':code[m.end():end-1],'start':m.start(),'end':end})
assert [d['name'] for d in sorted(defs,key=lambda d:d['start'])]==outline['native_inventory_and_order'];mapping={n.lower():n for n in outline['native_inventory_and_order']}
for name in re.findall(r'\b(\w+)\s*\(',code):
 if name.lower() in mapping:assert name==mapping[name.lower()],(name,'C++case')
local_arrays={n.lower() for n in ['big_L_w_1','big_L_w_2','big_L_rt_1','big_L_rt_2','big_L_thl_1','big_L_thl_2','sgn_wprtp','sgn_wpthlp','sgn_wp2','coef_sigma_w_1_sqd','coef_sigma_w_2_sqd','coef_sigma_rt_1_sqd','coef_sigma_rt_2_sqd','coef_sigma_thl_1_sqd','coef_sigma_thl_2_sqd']}
def lower(n):
 if not isinstance(n,tuple):return n
 if n[0]=='call' and n[1] in local_arrays:return ('v',n[1])
 return tuple(lower(a) if isinstance(a,tuple) else a for a in n)
def exprs(body):return [(m,Cpp(m.group(1),set()).parse(),Cpp(m.group(2),set()).parse()) for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;{}]*?\))?)\s*=(?!=)\s*([^;]*);',body)]
reviews=[]
for d in defs:
 name=d['name'];p=read_procedure(native,name);assert d['args']==[a.name for a in p.arguments]
 raw=re.search(r'(?ms)^  subroutine '+name+r'\(.*?^  end subroutine '+name+r'\b',source).group();expected=[]
 for line in statements(raw):
  m=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',line)
  if m:
   lhs,rhs=m.groups()
   if lhs.lower().startswith('small_l_'):continue
   expected.append((lower(fexpr(lhs)),lower(fexpr(rhs))))
 actual=[(a,b) for m,a,b in exprs(d['body']) if not m.group(1).startswith('negative_variances')]
 assert expected==actual,(name,'assignments',[(i,a,b) for i,(a,b) in enumerate(zip(expected,actual)) if a!=b],len(expected),len(actual))
 native_conditions=[]
 for line in statements(raw):
  m=re.match(r'(?:if|elseif)\s*\((.*)\)\s*then$',line,re.I)
  if m:native_conditions.append(lower(fexpr(m.group(1))))
 actual_conditions=[]
 body=d['body']
 if name=='tsdadg_pdf_driver':
  # The driver adds2 WHERE sign branches and a host warning-count/6print gate.
  kernel_start=body.index('KOKKOS_LAMBDA');begin=body.index('{',kernel_start)+1;end=begin;depth=1
  while depth:depth+=(body[end]=='{')-(body[end]=='}');end+=1
  body=body[begin:end-1]
 for m in re.finditer(r'\bif\s*\(',body):
  end=m.end();depth=1
  while depth:depth+=(body[end]=='(')-(body[end]==')');end+=1
  actual_conditions.append(Cpp(body[m.end():end-1],set()).parse())
 if name=='tsdadg_pdf_driver':actual_conditions=actual_conditions[2:]
 assert native_conditions==actual_conditions,(name,'predicates',native_conditions,actual_conditions)
 reviews.append({'name':name,'arguments':len(p.arguments),'ordered_assignment_trees':len(expected),'ordered_Native_IF_predicates':len(native_conditions)})
# Exact immutable Native small_l setup moved outside the work-item kernel.
for field in ('w','rt','thl'):
 assert re.search(r'const real_t small_l_'+field+r'_1\s*=\s*real_t\(\.75\),\s*small_l_'+field+r'_2\s*=\s*real_t\(\.5\)',code)
# Every scalar helper call keeps Native field order after local-array lowering.
d=defs[0];raw=re.search(r'(?ms)^  subroutine tsdadg_pdf_driver\(.*?^  end subroutine tsdadg_pdf_driver\b',source).group();expected=[lower(fexpr(re.sub(r'^call\s+','',l,flags=re.I))) for l in statements(raw) if re.match(r'call\s+',l,re.I)];actual=[]
for m in re.finditer(r'\b(calc_L_x_Skx_fnc|calc_setter_parameters|calc_respnder_parameters)\s*\(',d['body']):
 end=m.end();depth=1
 while depth:depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
 actual.append(Cpp(d['body'][m.start():end],set()).parse())
assert expected==actual
helper_calls=len(actual);helper_arguments=sum(len(c[2]) for c in actual)
# Native warning messages concatenate literals. C++ may split long literals for
# formatting; compare their exact concatenated message/data labels in order.
expected_labels=[]
for l in statements(raw):
 if re.match(r'write\s*\(',l,re.I):expected_labels.append(''.join(re.findall(r'"([^"]*)"',l)))
actual_labels=[]
for m in re.finditer(r'(?m)^\s*diagnostic\s*<<\s*([^;]*);',d['body']):
 if 'std::setprecision' in m.group(1):continue
 actual_labels.append(''.join(json.loads('"'+a+'"') for a in re.findall(r'"((?:\\.|[^"\\])*)"',m.group(1))))
assert expected_labels==actual_labels,(expected_labels,actual_labels)
fields=['sigma_w_1_sqd','sigma_w_2_sqd','sigma_rt_1_sqd','sigma_rt_2_sqd','sigma_thl_1_sqd','sigma_thl_2_sqd']
for idx,name in enumerate(fields):
 assert re.search(r'negative_variances\(k, '+str(idx)+r'\)\s*=\s*'+name+r'\(k\)',d['body'])
 assert re.search(r'if\s*\(warnings_host\(k, '+str(idx)+r'\) < zero\)',d['body'])
 assert re.search(r'<< warnings_host\(k, '+str(idx)+r'\)',d['body'])
 assert not re.search(r'\b'+re.escape('coef_'+name)+r'\(k\)',d['body'])
# Verify both complete Native mathematical-method sections remain in their
# original position, word/equation order and spelling (ignoring whitespace).
for before,after,index in [('tsdadg_pdf_driver','calc_setter_parameters',1),('calc_L_x_Skx_fnc','calc_respnder_parameters',3)]:
 begin=source.index('  !=============================================================================',source.index('  end subroutine '+before));end=source.index('  subroutine '+after,begin)
 expected=''.join(l.split('!',1)[1] if '!' in l else '' for l in source[begin:end].splitlines())
 region=text[defs[index-1]['end']:defs[index]['start']];actual=''.join(re.findall(r'//([^\n]*)',region))
 assert re.sub(r'\s+','',expected) in re.sub(r'\s+','',actual),(before,after,'method comments')
assert not re.search(r'f2k::|std::cerr|to_host_space|parallel_for|parallel_reduce',code[code.index('KOKKOS_LAMBDA'):code.index('},\n      warning_count')])
for d in defs[1:]:assert not re.search(r'f2k::|std::|Kokkos::View|parallel_for|parallel_reduce',d['body'])
assert max(map(len,text.splitlines()))<=100 and all(l.rstrip()==l and '\t' not in l for l in text.splitlines())
a={'canonical_source':str(native.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(native.read_bytes()).hexdigest(),'target_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'routines':reviews,'immutable_small_l_values':6,'ordered_helper_calls':helper_calls,'ordered_helper_arguments':helper_arguments,'exact_warning_and_data_labels':len(actual_labels),'warning_payload_fields':6,'complete_source_method_comment_regions':2,'runtime':'Source only; compilation/numerical evidence pending.'}
(HERE/'expression_audit.json').write_text(json.dumps(a,indent=2)+'\n');print('4signatures,',sum(r['ordered_assignment_trees'] for r in reviews),'assignmenttrees,14predicates,12calls/129args,12diagnosticlabels; source only')
