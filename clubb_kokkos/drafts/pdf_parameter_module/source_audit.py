from pathlib import Path
import re,sys,json,hashlib,ast
R=Path(__file__).resolve().parents[3];sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'));from fortran_procedures import statements,read_procedure
from generate_derived_types import split_top

def clean(t):return re.sub(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*.*?\*/',lambda m:m.group() if m.group()[0] in '\"\'' else '',t,flags=re.S)
def balanced(t,start,op='(',cl=')'):
 d=1;i=start;quote=''
 while d:
  c=t[i]
  if quote:
   if c==quote and t[i-1]!='\\':quote=''
  elif c in '\"\'':quote=c
  else:d+=(c==op)-(c==cl)
  i+=1
 return i

def raw(S,n):return re.search(r'^\s*subroutine '+n+r'\s*\(.*?^\s*end subroutine '+n+r'\b',S,re.S|re.M).group()
def audit(module):
 D=R/'clubb_kokkos/drafts'/module;P=R/('src/CLUBB_core/'+module+'.F90');S=P.read_text();O=json.loads((D/'SOURCE_OUTLINE.json').read_text());T=(D/(module+'.cpp')).read_text();C=clean(T);names=O['native_inventory_and_order'];defs=[]
 for m in re.finditer(r'void\s+(\w+)\s*\(',C):
  if m.group(1) not in names:continue
  e=balanced(C,m.end());tail=e
  while C[tail].isspace():tail+=1
  if C[tail]!='{':continue
  b=tail+1;f=balanced(C,b,'{','}');defs.append({'name':m.group(1),'args':[x.strip().split()[-1] for x in split_top(C[m.end():e-1])],'body':C[b:f-1]})
 assert [x['name'] for x in defs]==names
 for d,n in zip(defs,names):assert d['args']==[a.name for a in read_procedure(P,n).arguments]
 body={d['name']:re.sub(r'\s+\)',')',re.sub(r'\(\s+','(',re.sub(r'\s+',' ',d['body']))) for d in defs};reviews=[];details=[]
 if module=='pdf_parameter_module':
  zeropdf=raw(S,'zero_pdf_params_api');pdf=re.findall(r'pdf_params%(\w+)\(i,k\)\s*=\s*zero',zeropdf);assert len(pdf)==47
  zeroimp=raw(S,'zero_pdf_implicit_coefs_terms_api');imp2=re.findall(r'pdf_implicit_coefs_terms%(\w+)\(i,k\)\s*=\s*zero',zeroimp);imp3=re.findall(r'pdf_implicit_coefs_terms%(\w+)\(i,k,sclr\)\s*=\s*zero',zeroimp);assert len(imp2)==19 and len(imp3)==8
  for n,owner,ff,ff3 in [('init_pdf_params','pdf_params',pdf,[]),('init_pdf_implicit_coefs_terms_api','pdf_implicit_coefs_terms',imp2,imp3)]:
   expected=re.findall(owner+r'%(\w+)\(\s*(?:1:)?ngrdcol\s*,\s*(?:1:)?nz(?:\s*,\s*1:sclr_dim)?\s*\)',raw(S,n));actual=re.findall(owner+r'\.(\w+)\s*=\s*decltype\('+owner+r'\.\1\)',body[n]);assert expected==actual==ff+ff3,(n,expected,actual)
   assert body[n].index(owner+' = {};')<body[n].index(owner+'.'+ff[0]);assert re.search(owner+r'\.ngrdcol = ngrdcol;',body[n]);assert re.search(owner+r'\.nz = nz;',body[n])
   if ff3:assert 'if (sclr_dim > 0)' in body[n] and owner+'.sclr_dim = sclr_dim;' in body[n]
   call='zero_pdf_params_api' if owner=='pdf_params' else 'zero_pdf_implicit_coefs_terms_api';assert call+'('+owner+');' in body[n];reviews.append({'name':n,'allocation_field_order':len(actual),'metadata_and_nested_reset':True})
  for n,owner,ff,ff3 in [('zero_pdf_params_api','pdf_params',pdf,[]),('zero_pdf_implicit_coefs_terms_api','pdf_implicit_coefs_terms',imp2,imp3)]:
   actual=re.findall(owner+r'_(\w+)\(i, k(?:, sclr)?\) = zero;',body[n]);assert actual==ff+ff3
   for f in ff+ff3:assert 'const auto '+owner+'_'+f+' = '+owner+'.'+f+';' in re.sub(r'\s+',' ',body[n])
   if ff3:assert 'if ('+owner+'.sclr_dim > 0)' in body[n]
   assert not re.search(owner+r'\.[a-z0-9_]+\s*=',body[n]);reviews.append({'name':n,'ordered_zero_field_assignments':len(actual),'allocation_preserved':True})
  for n in ['pack_pdf_params_api','unpack_pdf_params_api']:
   native=raw(S,n)
   if n.startswith('pack'):
    expected=[(int(slot),f) for slot,f in re.findall(r'r_param_array\(:,(\d+)\)\s*=\s*pdf_params%(\w+)\(1,k_start:k_end\)',native)];actual=[(int(slot)+1,f) for slot,f in re.findall(r'r_param_array\(level, (\d+)\)\s*=\s*pdf_params_(\w+)\(0, k\);',body[n])]
   else:
    expected=[(int(slot),f) for f,slot in re.findall(r'pdf_params%(\w+)\(1,k_start:k_end\)\s*=\s*r_param_array\(:,(\d+)\)',native)];actual=[(int(slot)+1,f) for f,slot in re.findall(r'pdf_params_(\w+)\(0, k\)\s*=\s*r_param_array\(level, (\d+)\);',body[n])]
   assert expected==actual==list(enumerate(pdf,1)),n
   assert 'if (k_start_in.has_value() && k_end_in.has_value())' in body[n]
   assert 'k_start = *k_start_in;' in body[n] and 'k_end = *k_end_in;' in body[n]
   assert 'k_start = 1;' in body[n] and 'k_end = nz;' in body[n]
   assert 'const int k = k_start - 1 + level;' in body[n] and 'k_end - k_start + 1 != nz' in body[n]
   assert 'pdf_params = {}' not in body[n] and not re.search(r'pdf_params\.(?:nz|ngrdcol)\s*=',body[n]);reviews.append({'name':n,'ordered_slot_mappings':47,'Native_first_column_and_optional_range':True})
  pp=body['print_pdf_params'];native_fields=re.findall(r'pdf_params%(\w+)\(i,:\)',raw(S,'print_pdf_params'));actual_labels=re.findall(r'output << "(\w+)\(i,:\) = "',pp);assert native_fields==actual_labels==pdf
  for f in pdf:assert 'to_host_space(pdf_params.'+f+')' in pp and f+'_host(i - 1, k)' in pp
  reviews.insert(4,{'name':'print_pdf_params','ordered_field_labels_and_values':47,'per_column_output_order':True})
  decl=(D/'procedure_declarations.hpp.inc').read_text();assert 'inline constexpr int num_pdf_params = 47;' in decl
  for n in ['pack_pdf_params_api','unpack_pdf_params_api']:
   assert n in decl and 'buffer_exec' in decl
  assert 'to_execution_space(r_param_array, false)' in decl and 'from_execution_space(r_param_array, buffer_exec);' in decl
  wrapper=re.sub(r'\(\s+','(',re.sub(r'\s+',' ',(D/'clubb_api_module.hpp.inc').read_text()));assert 'init_pdf_params_api' in wrapper and 'init_pdf_params(nz, ngrdcol,' in wrapper
  details=['All 47 PDF allocation/zero/packing/printing field orders and19+8 implicit-coefficient fields independently read from Native and matched against target statements.','Native first-column47slot mapping and both-present optional range retained; target explicitly rejects nonconformable slice length outside the defined Native array-assignment contract.','Unpack changes neither other columns nor metadata and does not import JAX zeroed reconstruction. Scalar implicit fields remain unallocated/untouched when sclr_dim<=0.','Existing Native API forwarder init_pdf_params_api staged separately for its actual two driver callers; no duplicate mirror schema or copy helper added.','The two unused Native copy helpers have no current repository CALL sites/JAX bodies; their Native exports remain unchanged. C++ printing has complete labels/values/order but unverified list-directed spacing.']
 elif module=='hydromet_pdf_parameter_module':
  n='init_hydromet_pdf_params';ff=re.findall(r'hydromet_pdf_params%(\w+)\s*=\s*zero',raw(S,n));actual=re.findall(r'Kokkos::deep_copy\(hydromet_pdf_params\.(\w+), zero\);|hydromet_pdf_params\.(\w+) = zero;',body[n]);actual=[a or b for a,b in actual];assert ff==actual and len(ff)==18
  assert 'hydromet_pdf_parameter_collection::max_hydromet_dim' in body[n]
  for f in ff[:14]:assert 'if (!hydromet_pdf_params.'+f+'.is_allocated())' in body[n]
  reviews.append({'name':n,'ordered_zero_fields':18,'fixed_view_contracts':14})
  n='init_precip_fracs_api';ff=['precip_frac','precip_frac_1','precip_frac_2'];actual=re.findall(r'precip_fracs\.(\w+)\s*=\s*decltype\(precip_fracs\.\1\)',body[n]);assert actual==ff
  assert 'precip_fracs = {};' in body[n] and 'precip_fracs.ngrdcol = ngrdcol;' in body[n] and 'precip_fracs.nzt = nzt;' in body[n] and 'zero_precip_fracs_api(precip_fracs);' in body[n];reviews.append({'name':n,'ordered_allocation_fields':3,'metadata_and_reset':True})
  n='zero_precip_fracs_api';assert re.findall(r'precip_fracs_(\w+)\(i, k\) = zero;',body[n])==ff;reviews.append({'name':n,'ordered_zero_fields':3,'allocation_preserved':True})
  decl=(D/'procedure_declarations.hpp.inc').read_text();assert 'hydromet_pdf_params.zero();' in decl
  details=['Native fixed-size 14 array fields reuse existing mirror Views and collection max-8 constant; allocate missing fixed storage once, reject malformed fixed shapes, and reset18 fields in source order.','Precipitation IntentOut releases/recreates exactly three Views, sets original metadata, then calls its Native-named reset. Reset mutates existing allocations.','Whole-grid hydrometeor initialization is a thin overload of the existing generated collection::zero owner, justified by actual Native per-cell initialization and existing C++ collection.zero calls. No new schema/global cache.']
 else:
  assert 'const int mpi_rank_dummy = 1;' in body['init_default_err_info_api'] and 'const int chunk_idx_dummy = 1;' in body['init_default_err_info_api']
  assert re.search(r'init_err_info_api\(ngrdcol, mpi_rank_dummy, chunk_idx_dummy,',body['init_default_err_info_api'])
  b=body['init_err_info_api'];assert re.findall(r'err_info\.(\w+) = decltype\(err_info\.\1\)',b)==['host_lat','host_lon','err_code','err_header'];assert 'Kokkos::deep_copy(err_info.err_code, clubb_no_error);' in b
  b=body['set_err_info_values_api'];assert re.search(r'if \(!err_info.err_header.is_allocated\(\) \|\| !err_info.host_lat.is_allocated\(\) \|\|\s*!err_info.err_header.is_allocated\(\)\)',b)
  stages=['init_default_err_info_api(','if (mpi_rank_in.has_value())','if (chunk_idx_in.has_value())','if (lat_in.has_value())','if (lon_in.has_value())','update_headers(ngrdcol, err_info);'];loc=[b.index(x) for x in stages];assert loc==sorted(loc)
  for coord,limit in [('lat',90),('lon',180)]:
   assert 'if ('+coord+'_values(i) > real_t('+str(limit)+'.))' in b
   assert 'tmp = '+coord+'_values(i) - real_t('+str(limit)+'.);' in b
   assert 'else ' in b and 'tmp = '+coord+'_values(i);' in b
  formatter=body['update_headers'];assert 'auto formatted = err_info;' in formatter and 'f2k::format_err_info_headers(ngrdcol, formatted);' in formatter
  assert 'Kokkos::deep_copy(err_info.err_header, formatted.err_header);' in formatter and 'err_info.err_header_global = formatted.err_header_global;' in formatter
  assert re.findall(r'err_info\.(\w+) = \{\};',body['cleanup_err_info_api'])==['err_header','err_code','host_lat','host_lon']
  manifest=(D/'f2k/procedure_manifest.py.inc').read_text();assert '"format_err_info_headers": ("set_err_info_values_api", set()),' in manifest
  assert not re.search(r'f2k::(?:update_headers|f2k_get|f2k_set)',C)
  # Prove public absent-optional path is formatting-only from actual Native.
  native_set=raw(S,'set_err_info_values_api');assert 'if ( present(lat_in) ) then' in native_set and 'if ( present(lon_in) ) then' in native_set and 'call update_headers(ngrdcol, err_info)' in native_set
  reviews=[{'name':d['name'],'arguments':len(d['args']),'contract_and_source_order_reviewed':True} for d in defs]
  details=['Default MPI/chunk values both1 and actual Native dummy argument order retained; required initializer allocation order/no-error fill then setter call preserved.','Native setter guard visibly repeats err_header and omits host_lon/err_code; optional MPI/chunk/lat/lon processing order and exact >90/-90,>180/-180 formulas retained.','Explicit optional View adapters preserve independent absence and caller memory-space staging. Actual driver lat/lon-only keyword call needs two nullopt integer placeholders at promotion.','Private update_headers delegates exact fixed-width I/O through an alias of PUBLIC Native setter with all optionals absent; actual state/allocation/normalization is C++. The alias neither exposes the Native private helper nor adds permanent private accessors.','Only formatted header bytes/string return to the original case from a call-local snapshot. No C++ mutable global state.','Cleanup preserves Native component release order and scalar context. Kokkos refcount release has no Native DEALLOCATE STAT status; no fabricated failure tests/warnings. Already allocated InOut initialization throws instead of silently reallocating.']
 # Captures: owner objects stay on host; explicit View aliases enter kernels.
 kernels=[]
 for m in re.finditer(r'KOKKOS_LAMBDA\([^)]*\)\s*\{',C):
  e=balanced(C,m.end(),'{','}');b=C[m.end():e-1];assert not re.search(r'\b(?:pdf_params|pdf_implicit_coefs_terms|precip_fracs|err_info)\.|f2k::|std::cerr|to_host_space|parallel_for|parallel_reduce',b),b;kernels.append(b)
 assert not re.search(r'\bstruct\s+(?:pdf_parameter|implicit_coefs_terms|hydromet_pdf_parameter|precipitation_fractions|err_info_type)\s*\{',C)
 assert not re.search(r'\bthread_local\b|\bfma\s*\(',C)
 files=[D/(module+'.cpp'),D/'procedure_declarations.hpp.inc',D/'procedure_includes.hpp.inc']
 if module=='pdf_parameter_module':files.append(D/'clubb_api_module.hpp.inc')
 if module=='err_info_type_module':files.append(D/'f2k/procedure_manifest.py.inc')
 result={'canonical_source':str(P.relative_to(R)),'canonical_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'hashes':{str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},'routines':reviews,'kernels_reviewed':len(kernels),'manual_review':details,'review_kind':'Owner full Native/JAXread, source contract/field/order review; separate cleanup and independent reviewer tracked separately.','runtime_status':'Uncompiled and untested. No CPU/GPU/BFB claim. All lifecycle, optional bounds, errors, printer and formatting backend behavior remain pending integration.'}
 (D/'SOURCE_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__ == '__main__':
 result = audit(Path(__file__).resolve().parent.name)
 print('Owner source audit passed:', len(result['routines']), 'routines;', result['kernels_reviewed'], 'kernels. No compilation or numerical claim.')
