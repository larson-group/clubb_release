"""Compare scalar metadata contracts and decision trees; do not execute physics."""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
reader=(ROOT/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();exec(reader[reader.index('class Cpp:'):reader.index('reviews=[]')])
NATIVE=ROOT/'src/CLUBB_core/index_mapping.F90';CPP=HERE/'index_mapping.cpp';HPP=HERE/'index_mapping.hpp';source=NATIVE.read_text();text=CPP.read_text();header=HPP.read_text();code=re.sub(r'//[^\n]*','',text);hcode=re.sub(r'//[^\n]*','',header);names=json.loads((HERE/'SOURCE_OUTLINE.json').read_text())['native_inventory_and_order'];defs=[]
for m in re.finditer(r'(int|real_t)\s+(\w+)\s*\(([^{};]*)\)\s*\{',code,re.S):
 end=m.end();depth=1
 while depth:depth+=(code[end]=='{')-(code[end]=='}');end+=1
 defs.append({'name':m.group(2),'type':m.group(1),'args':[a.strip().split()[-1] for a in split_top(m.group(3))],'body':code[m.end():end-1]})
assert [d['name'] for d in defs]==names

def canonical(n):
 if not isinstance(n,tuple):return n
 if n[0]=='field' and n[1]==('v','hm_metadata'):return ('v','hm_metadata_'+n[2])
 return tuple(canonical(x) if isinstance(x,tuple) else x for x in n)
reviews=[]
for d in defs:
 name=d['name'];p=read_procedure(NATIVE,name);assert d['args']==[a.name for a in p.arguments]
 h=re.search(r'\b'+name+r'\s*\(([^;]*)\);',hcode);assert [a.strip().split()[-1] for a in split_top(h.group(1))]==d['args']
 raw=re.search(r'^  function '+name+r'\(.*?^  end function '+name+r'\b',source,re.S|re.M).group();expected=[];pred=[]
 for l in statements(raw):
  if m:=re.match(r'^(\w+)\s*=(?!=)\s*(.*)$',l):expected.append((fexpr(m.group(1)),fexpr(m.group(2))))
  if m:=re.match(r'(?:if|elseif)\s*\((.*)\)\s*then$',l,re.I):pred.append(fexpr(m.group(1)))
 actual=[(canonical(Cpp(m.group(1),set()).parse()),canonical(Cpp(m.group(2),set()).parse())) for m in re.finditer(r'(?m)^\s*(\w+)\s*=(?!=)\s*([^;]*);',d['body'])]
 actualpred=[canonical(Cpp(m.group(1),set()).parse()) for m in re.finditer(r'\bif\s*\(([^()]*)\)',d['body'])]
 assert expected==actual and pred==actualpred,(name,'assignment or predicate')
 assert re.search(r'\breturn\s+'+p.result.name+r'\s*;',d['body'])
 # Case-sensitive member names must be present in the actual generated owner.
 owner=(ROOT/'clubb_kokkos/src/CLUBB_core/corr_varnce_module.hpp').read_text()
 for member in re.findall(r'hm_metadata\.(\w+)',d['body']):assert re.search(r'\b'+member+r'\b',owner),member
 reviews.append({'name':name,'required_arguments':len(p.arguments),'ordered_assignments':len(actual),'ordered_branch_predicates':len(pred),'result_name_and_type':'exact','case_sensitive_generated_members':'verified'})
assert not re.search(r'f2k::|KOKKOS_LAMBDA|std::cout|std::cerr|thread_local|static ',code)
for t in [text,header]:assert max(map(len,t.splitlines()))<=100 and all(l.rstrip()==l and '\t' not in l for l in t.splitlines())
rec={'canonical_source':str(NATIVE.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),'cpp_sha256':hashlib.sha256(CPP.read_bytes()).hexdigest(),'hpp_sha256':hashlib.sha256(HPP.read_bytes()).hexdigest(),'routines':reviews,'manual_review':['Full Native/JAX scope; no JAX-only helper, guarded sentinel or0based conversion imported.','Field-specific branch order including inverse mapping order differences and scalar radius OR predicates retained.','Existing generated26-field metadata owner reused; functions read scalar metadata only, no array staging/global owner/device captures.','All five full procedure descriptions and species comments retained. Stale module ownership/initialization comments corrected to existing hm_metadata_type/init_pdf_hydromet_arrays_api.','Separate cleanup:100columns, whitespace, exact case-sensitive names, same argument order, no target-only helper or optional.'],'remaining_dependencies':'None beyond existing generated metadata type and Native-valued constants_clubb maximum radii/zero.','numerical_evidence':'No build, unit/model test, runtime or BFB claim.'};(HERE/'SOURCE_AUDIT.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(reviews))
