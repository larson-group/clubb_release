from pathlib import Path
import sys,re,json,hashlib
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import find_procedure
from generate_derived_types import split_top
names=['ADG1_pdf_driver','ADG2_pdf_driver','Luhar_3D_pdf_driver','ADG1_w_closure','calc_Luhar_params','close_Luhar_pdf','ADG1_ADG2_responder_params','backsolve_Luhar_params','max_cubic_root']
source=(root/'src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90').read_text();target=(folder/'adg1_adg2_3d_luhar_pdf.cpp').read_text()
clean=lambda s:re.sub(r'//[^\n]*|/\*.*?\*/','',s,flags=re.S)
cpp=clean(target)
match=list(re.finditer(r'(void|real_t)\s+(\w+)\s*\(([^)]*)\)\s*\{',cpp))
assert [m[2] for m in match]==names
blocks={}
for m in match:
 n=m[2];args=[a.strip().split()[-1] for a in split_top(m[3])]
 assert args==[a.name for a in find_procedure(n).arguments],n
 depth=1;i=m.end()
 while depth:
  depth+=(cpp[i]=='{')-(cpp[i]=='}');i+=1
 blocks[n]=cpp[m.end():i-1]
def calls(s,is_cpp):
 s=clean(s) if is_cpp else '\n'.join(l.split('!')[0] for l in s.splitlines())
 pat=(r'\b('+'|'.join(names)+r')\s*\(') if is_cpp else (r'call\s+('+'|'.join(names)+r')\s*\(')
 result=[]
 for m in re.finditer(pat,s,re.I):
  start=m.end();i=start;d=1
  while d:
   d+=(s[i]=='(')-(s[i]==')');i+=1
  args=s[start:i-1]
  if is_cpp and re.search(r'\b(?:int|Real\dD|HostReal1D|real_t)\s+\w',args):continue
  # Source function, not subroutine, is recorded separately below.
  if m[1]=='max_cubic_root':continue
  result.append((m[1],split_top(args.replace('&',''))))
 return result
def norm(a,is_cpp):
 if is_cpp:
  a=re.sub(r'Kokkos::subview\(\s*(\w+)\s*,\s*std::make_pair\(\s*k\s*,\s*k\s*\+\s*1\s*\)\s*\)',r'\1(k)',a)
  a=re.sub(r'Kokkos::subview\(\s*(\w+)\s*,\s*i\s*,\s*Kokkos::ALL\s*\)',r'\1(i,:)',a)
  a=re.sub(r'Kokkos::subview\(\s*(\w+)\s*,\s*Kokkos::ALL\s*,\s*Kokkos::ALL\s*,\s*sclr\s*\)',r'\1(:,:,sclr)',a)
  a=a.replace('sclr_tol_host','sclr_tol').replace('thl_tol * thl_tol','thl_tol**2').replace('rt_tol * rt_tol','rt_tol**2')
 return re.sub(r'[\s()]','',a).lower()
counts={}
for n in names:
 native=re.search(r'^  (?:subroutine|function) '+n+r'\(.*?^  end (?:subroutine|function) '+n+r'\b',source,re.M|re.S).group()
 nc=calls(native,False);cc=calls(blocks[n],True)
 assert [a for a,b in nc]==[a for a,b in cc],(n,nc,cc)
 for (nn,na),(cn,ca) in zip(nc,cc):
  assert [norm(a,False) for a in na]==[norm(a,True) for a in ca],(n,nn,na,ca)
 counts[n]=len(nc)
# All public/native-only profile adapters retain exact positional contracts.
head=clean((folder/'adg1_adg2_3d_luhar_pdf.hpp').read_text())
for n in names[:6]:
 signatures=re.findall(r'void\s+'+n+r'\s*\(([^)]*)\)',head)
 assert len(signatures)==(2 if n in {'Luhar_3D_pdf_driver','calc_Luhar_params','close_Luhar_pdf'} else 1)
 for args in signatures:assert [a.strip().split()[-1] for a in split_top(args)]==[a.name for a in find_procedure(n).arguments]
for name in ['adg1_adg2_3d_luhar_pdf.cpp','adg1_adg2_3d_luhar_pdf.hpp','f2k_dependencies.hpp']:
 text=(folder/name).read_text();code=clean(text);code=re.sub(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'','',code)
 for a,b in [('(',')'),('{','}'),('[',']')]:
  depth=0
  for c in code:
   depth+=(c==a)-(c==b);assert depth>=0,(name,a)
  assert depth==0,(name,a,depth)
 assert all(len(l)<=100 and l.rstrip()==l and '\t' not in l for l in text.splitlines())
assert 'f2k::ADG1_w_closure' not in target
assert target.count('small_m_root = max_cubic_root(')==1
assert 'std::numeric_limits<real_t>::max()' in target
assert 'constexpr bool l_use_cubic_backsolve = true' in target
assert 'No Native BFB claim' in target
record={'state':'source-complete draft; NOT compiled or runtime validated','authority':'src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90 (all1537lines read)','routine_order':names,'argument_counts':{n:len(find_procedure(n).arguments) for n in names},'all_required_no_optionals':True,'call_counts_in_native_order':counts,'every_native_subroutine_call_argument':True,'public_api_and_profile_adapter_order':True,'token_balance_and_100cols':True,'definedness_deviation':'close_Luhar_pdf initializes covariance sign for nonvarying normalized outputs, matching JAX; native leaves it undefined','pending_lower_dependencies':['f2k::cubic_solve returns host complex(nz,3)','f2k::quadratic_solve returns host complex(nz,2)'],'manual_source_audit':['all nine routines source order/signatures/intents/calls','ADG1 debug all-column errors, parent fatal gate, field mappings/native wind thl_tol','ADG2 source column loop, sigma overwrite and scalar flags/tolerances','Luhar three setter branches/tie order and all covariance/tolerance/sequence associations','ADG1_w threshold, mixture clamp, signed normalized means, equal variances','calc_Luhar varying/nonvarying formulas and real power1/3','close_Luhar two branches, definedness boundary, unchanged physical means/variances','backsolve near-half/zero/ordinary cubic and inactive quadratic branches, Native HUGE','max-root coefficient thresholds, real/complex choice, quadratic and unguarded linear fallback','scientific descriptions/citations, inactive branch comments, local state ownership','profile host boundary/adapters and single-evaluation root/candidate temporaries classified'],'runtime':'No builds/tests run; root owns hookup and model comparisons','sha256':{p:hashlib.sha256((folder/p).read_bytes()).hexdigest() for p in ['adg1_adg2_3d_luhar_pdf.cpp','adg1_adg2_3d_luhar_pdf.hpp','f2k_dependencies.hpp']}}
(folder/'source_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k in ['state','routine_order','argument_counts','call_counts_in_native_order','pending_lower_dependencies']}))
