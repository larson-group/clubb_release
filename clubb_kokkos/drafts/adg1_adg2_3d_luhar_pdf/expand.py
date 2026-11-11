from pathlib import Path
import sys,re
root=Path(__file__).resolve().parents[3]; folder=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import find_procedure,statements,declaration
from generate_derived_types import split_top
src=(root/'src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90').read_text()
names=['ADG1_pdf_driver','ADG2_pdf_driver','Luhar_3D_pdf_driver','ADG1_w_closure','calc_Luhar_params','close_Luhar_pdf','ADG1_ADG2_responder_params','backsolve_Luhar_params','max_cubic_root']
public=names[:6]; profiles={'Luhar_3D_pdf_driver','calc_Luhar_params','close_Luhar_pdf'}
blocks={n:re.search(r'^  (?:subroutine|function) '+n+r'\(.*?^  end (?:subroutine|function) '+n+r'\b',src,re.M|re.S).group() for n in names}
procs={n:find_procedure(n) for n in names}
entities={n:{a.name.lower():a for st in statements(b) for a in declaration(st).values()} for n,b in blocks.items()}
def ctype(a,n):
 if a.kind=='integer':return 'int'
 if a.kind=='logical':return 'bool'
 if a.kind=='err_info_type':return 'err_info_type&'
 if a.kind=='real':
  if a.dims:return ('HostReal1D' if len(a.dims)==1 and (n in profiles or n=='max_cubic_root') else 'Real'+str(len(a.dims))+'D')
  return 'real_t'+('&' if a.intent!='in' else '')
 raise ValueError(a)
def sig(n):
 args={a.name.lower():a for a in procs[n].arguments};head='\n'.join(l.split('!')[0] for l in blocks[n].splitlines());head=head[:head.index(')')];groups=[]
 for l in head.splitlines():
  code=re.sub(r'^\s*(?:subroutine|function)\s+\w+\(','',l.replace('&',''));parts=[p.strip().lower() for p in code.split(',') if p.strip()]
  if parts:groups.append((', '.join(ctype(args[x],n)+' '+args[x].name for x in parts), args[parts[-1]].intent))
 return ('real_t' if n=='max_cubic_root' else 'void')+' '+n+'('+''.join(g+(', // '+{'in':'In','out':'Out','inout':'InOut'}[intent]+'\n    ' if i<len(groups)-1 else ')') for i,(g,intent) in enumerate(groups))
def description(n):
 b=blocks[n];i=b.find('    ! Description:')
 if i<0:i=b.find('    ! Calculates')
 if i<0:return '  // Native '+n+' procedure.\n'
 end=b.index('    use ',i)
 return '\n'.join('  //'+l.strip()[1:] if l.strip().startswith('!') else '' for l in b[i:end].splitlines())+'\n'
def power(s):
 while '**' in s:
  pos=s.index('**');end=pos
  while s[end-1].isspace():end-=1
  start=end-1
  if s[start]==')':
   lev=1;start-=1
   while lev:
    lev+=(s[start]==')')-(s[start]=='(');start-=1
   start+=1;j=start-1
   while j>=0 and (s[j].isalnum() or s[j] in '_.%'):j-=1
   if j<start-1:start=j+1
  else:
   while start>=0 and (s[start].isalnum() or s[start] in '_.%'):start-=1
   start+=1
  m=re.match(r'\s*(\w+(?:\.\d+)?)',s[pos+2:]);assert m,s
  base=s[start:end];exponent=m[1];endexp=pos+2+m.end()
  if exponent in {'2','3','4'}:repl='('+' * '.join(base for i in range(int(exponent)))+')'
  else:repl='Kokkos::pow('+base+', '+exponent+')'
  s=s[:start]+repl+s[endexp:]
 return s

def expr(s,n,elementwise=False):
 s=power(s)
 s=re.sub(r'((?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)_core_rknd',r'real_t(\1)',s)
 s=s.replace('&',' ').replace('%','.')
 for a,b in [('.and.','&&'),('.or.','||'),('.not.','!'),('.true.','true'),('.false.','false'),('/=','!=')]:s=re.sub(re.escape(a),b,s,flags=re.I)
 s=re.sub(r'huge\([^)]*\)','std::numeric_limits<real_t>::max()',s,flags=re.I)
 for f in ['sqrt','abs','exp']:s=re.sub(r'(?<![\w:])'+f+r'\s*\(','Kokkos::'+f+'(',s,flags=re.I)
 s=re.sub(r'aimag\s*\(','std::imag(',s,flags=re.I)
 s=re.sub(r'real\(\s*(.*?),\s*kind\s*=\s*core_rknd\s*\)',r'std::real(\1)',s,flags=re.I)
 for f,op in [('min','<'),('max','>')]:
  while m:=re.search(r'(?<![\w:])'+f+r'\s*\(',s,re.I):
   pos=m.end();lev=1
   while lev:lev+=(s[pos]=='(')-(s[pos]==')');pos+=1
   args=split_top(s[m.end():pos-1]);value=args[0]
   for a in args[1:]:value='(('+value+') '+op+' ('+a+') ? ('+value+') : ('+a+'))'
   s=s[:m.start()]+value+s[pos:]
 for low,a in entities[n].items():s=re.sub(r'\b'+re.escape(low)+r'\b',a.name,s,flags=re.I)
 # Native fixed singleton coefficient/root indices retain explicit zero-based offsets.
 s=re.sub(r'((?:[abcd]_coef(?:_in)?))\(1\)',r'\1(0)',s)
 s=re.sub(r'((?:cubic|quadratic)_roots)\(1,([123])\)',lambda m:m[1]+'(0,'+str(int(m[2])-1)+')',s)
 if elementwise:
  for low,a in entities[n].items():
   if a.dims and len(a.dims)==1:s=re.sub(r'\b'+re.escape(a.name)+r'\b(?!\s*\()',a.name+'(k)',s)
 return re.sub(r'\s+',' ',s).strip()

def call(st,n):
 m=re.match(r'call\s+(\w+)\s*\((.*)\)',st,re.S|re.I);callee=m[1];values=split_top(m[2]);proc=procs[callee]
 out=[]
 for val,a in zip(values,proc.arguments):
  if a.dims:
   mprofile=re.fullmatch(r'(\w+)\(k\)',val.strip())
   mcolumn=re.fullmatch(r'(\w+)\(i,:\)',val.strip())
   msclr=re.fullmatch(r'(\w+)\(:,:,sclr\)',val.strip())
   if mprofile:val='Kokkos::subview('+mprofile[1]+', std::make_pair(k,k+1))'
   elif mcolumn:val='Kokkos::subview('+mcolumn[1]+', i, Kokkos::ALL)'
   elif msclr:val='Kokkos::subview('+msclr[1]+', Kokkos::ALL, Kokkos::ALL, sclr)'
  val=expr(val,n)
  if n=='ADG2_pdf_driver':val=val.replace('sclr_tol(sclr)','sclr_tol_host(sclr)')
  out.append(val)
 return callee+'('+',\n    '.join(out)+');'

def translated(sec,n,elementwise=False):
 out=[];pending=''
 for line in sec.splitlines():
  if line.strip().startswith('!$acc'):continue
  if line.strip().startswith('!'):out.append('  //'+line.strip()[1:]);continue
  code=line.split('!')[0].strip()
  if not code:continue
  pending+=' '+code.strip('& ')
  if code.endswith('&'):continue
  st=pending.strip();pending=''
  if st.lower()=='return':continue
  if re.match(r'end\s*(?:do|if|where)',st,re.I):out.append('  }');continue
  if re.match(r'elsewhere|else\s*$',st,re.I):out.append('  } else {');continue
  m=re.match(r'(else\s*if|elseif|if|where)\s*\((.*)\)\s*(?:then)?$',st,re.I)
  if m:
   out.append(('  } else if' if m[1].lower().replace(' ','') in {'elseif'} else '  if')+' ('+expr(m[2],n,elementwise)+') {');continue
  m=re.match(r'do\s+(\w+)\s*=\s*1\s*,\s*(\w+)(?:\s*,\s*1)?$',st,re.I)
  if m:
   idx,dim=m.groups();out.append('  for (int '+idx+' = 0; '+idx+' < '+dim+'; ++'+idx+') {');continue
  if st.lower().startswith('call '):out.append('  '+call(st,n));continue
  if '=' in st:out.append('  '+expr(st,n,elementwise)+';');continue
  raise ValueError((n,st))
 assert not pending
 return '\n'.join(out)+'\n'
def local_views(n):
 b=blocks[n];m=re.search(r'![- ]*Local [Vv]ariable[s]?[- ]*',b)
 if n=='Luhar_3D_pdf_driver':start=b.index('    real( kind = core_rknd ), dimension(nz) ::')
 elif m:start=m.end()
 else:return ''
 sec=b[start:];sec=sec[:sec.find('    !----- Begin') if '    !----- Begin' in sec else len(sec)]
 args={a.name for a in procs[n].arguments};out=[]
 for st in statements(sec):
  for a in declaration(st).values():
   if a.name in args or a.kind!='real' or not a.dims:continue
   t='HostReal1D' if len(a.dims)==1 else 'Real'+str(len(a.dims))+'D'
   out.append('  '+t+' '+a.name+'("'+n+'::'+a.name+'", '+','.join(a.dims)+');')
 return '\n'.join(out)+'\n'
# Preserve compact original definitions unchanged for incremental root hookup.
compact=(folder/'incremental_adg1/adg1_adg2_3d_luhar_pdf.cpp').read_text()
def cpp_definition(n):
 start=compact.index('void '+n+'(',compact.index('} // namespace') if n=='ADG1_pdf_driver' else compact.index('namespace {',compact.index('void ADG1_pdf_driver(')))
 # Find the definition rather than its preceding prototype.
 pos=compact.index('{',start);depth=1;i=pos+1
 while depth:
  depth+=(compact[i]=='{')-(compact[i]=='}');i+=1
 return compact[start:i]
# First definition has no comment braces before function body; private comments
# include braces, so use the final known namespace suffix to delimit it.
adg1=compact[compact.index('void ADG1_pdf_driver('):compact.index('\nnamespace {',compact.index('void ADG1_pdf_driver('))].strip().replace('f2k::ADG1_w_closure','ADG1_w_closure')
resp_start=compact.index('void ADG1_ADG2_responder_params(',compact.index('\nnamespace {',compact.index('void ADG1_pdf_driver(')))
responder=compact[resp_start:compact.rindex('\n} // namespace\n')].strip()

bodies={}
n='ADG2_pdf_driver';b=blocks[n];sec=b[b.index('    ! Calculate the mixture fraction'):b.rindex('    return')]
# Profile whole-array overwrite is one independent device loop per source column.
sec=sec.replace('sigma_sqd_w(i,:) = min( one / ( one + small_m_w(i,:)**2 ), 0.99_core_rknd )','DRAFT_SIGMA_OVERWRITE')
translated_adg2=translated(sec.replace('DRAFT_SIGMA_OVERWRITE','sigma_overwrite = 0'),'ADG2_pdf_driver')
translated_adg2=translated_adg2.replace('sigma_overwrite = 0;', '''Kokkos::parallel_for("ADG2_pdf_driver::sigma_overwrite", nz, KOKKOS_LAMBDA(const int k) {
      const real_t width = one / (one + small_m_w(i,k)*small_m_w(i,k));
      sigma_sqd_w(i,k) = width < real_t(0.99) ? width : real_t(0.99);
    });''')
bodies[n]=description(n)+local_views(n)+'  const auto sclr_tol_host = to_host_space(sclr_tol);\n'+translated_adg2
n='Luhar_3D_pdf_driver';b=blocks[n];sec=b[b.index('    do k = 1, nz, 1'):b.rindex('    return')]
bodies[n]=description(n)+'''  // Host profile orchestration is an explicit temporary boundary while the
  // complex calc_roots backend remains native. Public adapters own staging.
'''+local_views(n)+translated(sec,n)
n='ADG1_w_closure';b=blocks[n];sec=b[b.index('    !$acc parallel loop'):b.index('    !$acc end parallel loop')]
# Reuse source scalar formulas inside the two-dimensional execution-space loop.
inner=sec[sec.index('        ! Width'):];inner=inner[:inner.rindex('      end do')]
bodies[n]=description(n)+'''  if (f2k::clubb_at_least_debug_level_api(2)) {
    int below_tolerance = 0;
    Kokkos::parallel_reduce("ADG1_w_closure::variance_check", ngrdcol*nz,
      KOKKOS_LAMBDA(const int cell, int& count) {
        const int i=cell%ngrdcol; const int k=cell/ngrdcol;
        if (wp2(i,k) < w_tol_sqd) ++count;
      }, below_tolerance);
    if (below_tolerance) {
      std::cerr << err_info.err_header_global << '\\n';
      std::cerr << "Error in ADG1_w_closure: wp2 must not be less than the w tolerance value.\\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  Kokkos::parallel_for("ADG1_w_closure::component_parameters", ngrdcol*nz,
    KOKKOS_LAMBDA(const int cell) {
      const int i=cell%ngrdcol; const int k=cell/ngrdcol;
'''+translated(inner,n)+'    });\n'
n='calc_Luhar_params';b=blocks[n];sec=b[b.index('    where ( xp2 >'):b.rindex('    return')]
sec=sec.replace('small_m = max( two_thirds * abs( Skx )**one_third, 0.05_core_rknd )','small_m_candidate = two_thirds * abs( Skx )**one_third\n       small_m = max( small_m_candidate, 0.05_core_rknd )')
bodies[n]=description(n)+local_views(n)+'  for (int k=0; k<nz; ++k) {\n    real_t small_m_candidate; // Evaluate the native MAX operand once.\n'+translated(sec,n,True)+'  }\n'
n='close_Luhar_pdf';b=blocks[n];sec=b[b.index('    where ( xp2 >'):b.rindex('    return')]
# Native reads this sign uninitialized in its nonvarying branch; choose the
# actual JAX definition, explicitly at this exact boundary, without new floors.
sec=sec.replace('       where ( wpxp >= zero )\n          sgn_wpxp = one\n       elsewhere ! <w\'x\'> < 0\n          sgn_wpxp = -one\n       endwhere ! <w\'x\'> >= 0','       ! Sign is initialized for both branches immediately above.')
bodies[n]=description(n)+local_views(n)+'''  for (int k=0; k<nz; ++k) {
    // Definedness deviation: Native initializes sign only for xp2>tolerance
    // then reads it below tolerance too. Match JAX's defined covariance sign;
    // physical means/variances retain Native formulas. No Native BFB claim
    // applies to its undefined nonvarying normalized means (Native1022–1060).
    sgn_wpxp(k) = wpxp(k) >= zero ? one : -one;
'''+translated(sec,n,True)+'  }\n'
n='backsolve_Luhar_params';b=blocks[n];sec=b[b.index('    if ( l_use_cubic_backsolve )'):b.rindex('    return')]
sec=re.sub(r'small_m_x\s*&\s*= sqrt\( max\( max_cubic_root\( a_coef_in, b_coef_in, &\s*c_coef_in, d_coef_in \), &\s*0\.05_core_rknd\*\*2 \) \)', 'small_m_root = max_cubic_root( a_coef_in, b_coef_in, c_coef_in, d_coef_in )\n          small_m_x = sqrt( max( small_m_root, 0.05_core_rknd**2 ) )', sec)
bodies[n]=description(n)+'''  real_t alpha, a, b, c, alpha_upr, alpha_low, discrim;
  real_t small_m_root; // Own the native function-result temporary: call once.
  HostReal1D a_coef_in("backsolve_Luhar_params::a_coef_in",1);
  HostReal1D b_coef_in("backsolve_Luhar_params::b_coef_in",1);
  HostReal1D c_coef_in("backsolve_Luhar_params::c_coef_in",1);
  HostReal1D d_coef_in("backsolve_Luhar_params::d_coef_in",1);
  constexpr bool l_use_cubic_backsolve = true;
'''+translated(sec,n)
n='max_cubic_root';b=blocks[n];sec=b[b.index('    ! Calculate a minimum threshold'):b.rindex('    return')]
bodies[n]=description(n)+'''  f2k::HostComplex2D cubic_roots;
  f2k::HostComplex2D quadratic_roots;
  real_t max_root, a_coef_thresh, b_coef_thresh;
'''+translated(sec,n).replace('cubic_solve(', 'f2k::cubic_solve(').replace('quadratic_solve(', 'f2k::quadratic_solve(')+'  return max_root;\n'

# Profile API overloads stage arbitrary caller Views once at their explicit
# host/native-root boundary. Native names/argument order/ownership are retained.
def profile_adapter(n):
 p=procs[n];arrays=[a for a in p.arguments if a.dims];mapping={a.name:'Array'+str(i) for i,a in enumerate(arrays)}
 args=[]
 for a in p.arguments:args.append((mapping[a.name] if a.dims else ctype(a,n))+' '+a.name)
 out='template <'+', '.join('class '+mapping[a.name] for a in arrays)+'>\nvoid '+n+'('+',\n    '.join(args)+') {\n'
 out+='  // Stage profiles at the temporary host/complex-root boundary.\n'
 for a in arrays:
  out+=('  auto '+a.name+'_host = to_host_space('+a.name+');\n' if a.intent=='in' else '  HostReal1D '+a.name+'_host("'+n+'::'+a.name+'",nz);\n')
 out+='  '+n+'('+', '.join(a.name+'_host' if a.dims else a.name for a in p.arguments)+');\n'
 for a in arrays:
  if a.intent!='in':out+='  from_execution_space('+a.name+', '+a.name+'_host);\n'
 return out+'}\n'
head='''// ADG1, ADG2 and 3-D Luhar component PDF closures.
//
// ADG1_pdf_driver and ADG2_pdf_driver close two-dimensional column/level fields;
// Luhar_3D_pdf_driver closes one profile. ADG1_w_closure supplies the vertical-
// velocity mixture, while calc_Luhar_params and close_Luhar_pdf supply Luhar
// mixture/response parameters. Means retain field units, variances squared units,
// and normalized means/mixture/width factors are dimensionless.
// ADG1 kernels use column-first LayoutLeft execution-space Views. The three
// Luhar profile APIs currently stage arbitrary Views to host C++ calculations
// while complex calc_roots remains a transient native dependency. No stats or
// persistent/global model storage is owned here. See DRAFT_STATUS.md for scope,
// missing lower bridges, definedness deviation and absent runtime validation.
#pragma once
#include "array_types.hpp"
#include "err_info_type_module.hpp"
namespace clubb_kokkos {
using HostReal1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
'''
head+='\n\n'.join(sig(n)+';' for n in public)+'\n\n'+'\n'.join(profile_adapter(n) for n in public if n in profiles)+'\n} // namespace clubb_kokkos\n'
(folder/'adg1_adg2_3d_luhar_pdf.hpp').write_text(head)
cpp='''// Component PDF closures; see the public header and draft status for boundaries.
#include "adg1_adg2_3d_luhar_pdf.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include "f2k_dependencies.hpp"
#include <iostream>
#include <limits>
#include <utility>
namespace clubb_kokkos {
namespace {
'''
cpp+='\n\n'.join(sig(n)+';' for n in names[6:])+'\n} // namespace\n\n'
for n in names:
 if n=='ADG1_ADG2_responder_params':cpp+='namespace {\n\n'
 cpp+=(adg1 if n=='ADG1_pdf_driver' else responder if n=='ADG1_ADG2_responder_params' else sig(n)+' {\n'+bodies[n]+'}')+'\n\n'
cpp+='} // namespace\n} // namespace clubb_kokkos\n'
(folder/'adg1_adg2_3d_luhar_pdf.cpp').write_text(cpp)
(folder/'f2k_dependencies.hpp').write_text('''// Draft-only declarations for the unported public complex root backend.
// Root integration must supply these transient adapters; there is no fallback.
#pragma once
#include "adg1_adg2_3d_luhar_pdf.hpp"
#include <complex>
namespace f2k {
using HostComplex2D = Kokkos::View<std::complex<clubb_kokkos::real_t>**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
HostComplex2D cubic_solve(int nz, clubb_kokkos::HostReal1D a_coef,
  clubb_kokkos::HostReal1D b_coef, clubb_kokkos::HostReal1D c_coef, clubb_kokkos::HostReal1D d_coef);
HostComplex2D quadratic_solve(int nz, clubb_kokkos::HostReal1D a_coef,
  clubb_kokkos::HostReal1D b_coef, clubb_kokkos::HostReal1D c_coef);
} // namespace f2k
''')
