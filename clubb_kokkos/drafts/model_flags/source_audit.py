"""Audit Native flag contracts/order/printing mappings; execute no model code."""
from pathlib import Path
import json,re,sys,hashlib
D=Path(__file__).resolve().parent;R=D.parents[2];P=R/'src/CLUBB_core/model_flags.F90';S=P.read_text();T=(D/'model_flags.cpp').read_text();H=(D/'procedure_declarations.hpp.inc').read_text();sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'));from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
pattern=r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"'
def clean(s):return re.sub(pattern,lambda m:re.sub(r'[^\n]',' ',m[0]) if m[0].startswith('/') else m[0],s,flags=re.S)
C=clean(T);names=json.loads((D/'SOURCE_OUTLINE.json').read_text())['native_inventory_and_order'];defs=[]
for m in re.finditer(r'void\s+('+ '|'.join(names)+r')\s*\((.*?)\)\s*\{',C,re.S):
 e=m.end();depth=1
 while depth:depth+=(C[e]=='{')-(C[e]=='}');e+=1
 args=[a.strip().split()[-1] for a in split_top(m[2])];assert args==[a.name for a in read_procedure(P,m[1]).arguments];defs.append({'name':m[1],'args':m[2],'body':C[m.end():e-1]})
assert [d['name'] for d in defs]==names+[names[-1]]
raw={n:re.search(r'^[ \t]*subroutine\s+'+n+r'\b.*?^[ \t]*end subroutine\s+'+n+r'\b',S,re.M|re.S).group() for n in names};reviews=[]
for d in defs[:2]:
 native=[(m[1].replace('%','.'),m[2].strip().replace('%','.').replace('.true.','true').replace('.false.','false')) for line in statements(raw[d['name']]) if (m:=re.match(r'([\w%]+)\s*=(?!=)\s*(.*)',line))];cpp=[(m[1],m[2].strip()) for m in re.finditer(r'(?m)^\s*([\w.]+)\s*=\s*([^;]+);',d['body'])];assert native==cpp and len(native)==66,(d['name'],native,cpp)
 reviews.append({'name':d['name'],'arguments':len(read_procedure(P,d['name']).arguments),'ordered_assignments':len(native)})
b=defs[1]['body'];native_pred=[re.sub(r'\s+','',m[1].replace('.or.','||')) for line in statements(raw[names[1]]) if (m:=re.match(r'if\s*\((.*?)\)\s*then',line))];cpp_pred=[re.sub(r'\s+','',m[1]) for m in re.finditer(r'if\s*\((.*?)\)\s*{',b)];assert native_pred==cpp_pred and len(native_pred)==2
native_errors=re.findall(r"error stop '(.*)'",raw[names[1]]);cpp_errors=[json.loads('"'+m+'"') for m in re.findall(r'throw std::runtime_error\(\s*"([^"\\]*(?:\\.[^"\\]*)*)"\s*\);',b)];assert native_errors==cpp_errors and len(native_errors)==2
assert b.index('throw std::runtime_error')<b.index('clubb_config_flags.iiPDF_type')
assert 'int iunit' in defs[2]['args'] and 'std::ostream& iunit' in defs[3]['args']
assert 'f2k::print_clubb_config_flags_api(iunit, clubb_config_flags);' in defs[2]['body']
writes=[]
for line in statements(raw[names[-1]]):
 if m:=re.match(r'write\(iunit,\*\)\s*"([^"]*)",\s*clubb_config_flags%(\w+)',line,re.I):writes.append(m.groups())
printed=re.findall(r'output\s*<<\s*"([^"]*)"\s*<<\s*\(?clubb_config_flags\.(\w+)',defs[3]['body']);assert writes==printed and len(writes)==64
assert ('l_mono_flux_lim_um = ','l_mono_flux_lim_vm') in writes and ('l_mono_flux_lim_vm = ','l_mono_flux_lim_um') in writes
assert not any(label.startswith(('saturation_formula','fill_holes_type')) for label,f in writes)
reviews.append({'name':names[-1],'required_arguments':2,'integer_unit_public_Native_IO_boundary':True,'stream_adapter_ordered_row_field_pairs':64,'Native_omitted_rows_and_wind_swap_preserved':True})
for n in names:
 for m in re.finditer(r'void\s+'+n+r'\s*\((.*?)\)\s*;',clean(H),re.S):assert [a.strip().split()[-1] for a in split_top(m[1])]==[a.name for a in read_procedure(P,n).arguments]
assert not re.search(r'\bstruct\s+clubb_config_flags_type|\bthread_local\b|KOKKOS_LAMBDA|\bfma\s*\(',C)
result={'scope':'All three Native behavior APIs plus same-named stream adapter justified by actual C++ case-owned output sinks. Existing mirror type reused.','canonical_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [D/'model_flags.cpp',D/'procedure_declarations.hpp.inc',D/'procedure_includes.hpp.inc',D/'f2k_manifest.py.inc',D/'generator_constants.py.inc']},'reviews':reviews,'ownership':'Native modules have no mutable module instances. All scalar outputs/configuration fields belong to existing caller/ClubbDriver. Immutable source constants remain generated in their original owner.','errors':'Two exact Native selector fatal checks precede all field assignments; unsupported invalid IntentOut state is not initialized artificially.','formatting_boundary':'int-iunit overload keeps exact Native I/O through public bridge. ostream overload preserves64row/value mappings; compiler-specific list-directed spacing unverified. No C++ logical-unit registry.','review_kind':'Owner source review; separate cleanup/independent review tracked separately.','runtime_status':'Draft only. No shared edits/builds/tests/CPU/GPU/BFB claim. Missing two selected constants/public I/O bridge are explicit integration dependencies.'};(D/'SOURCE_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print('Flags source audit passed: 3 Native procedures, 132 assignments, 2 fatal checks, 64 print mappings.')
