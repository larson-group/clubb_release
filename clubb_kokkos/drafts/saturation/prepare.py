from pathlib import Path
import re,sys
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from generate_derived_types import split_top
source=(root/'src/CLUBB_core/saturation.F90').read_text()
names=['sat_mixrat_liq_k','sat_mixrat_liq_2D','sat_vapor_press_liq','sat_vapor_press_liq_flatau','sat_vapor_press_liq_bolton','sat_vapor_press_liq_gfdl','sat_vapor_press_liq_lookup','sat_mixrat_ice_k','sat_mixrat_ice_2D','rcm_sat_adj']
blocks={n:re.search(r'^  (?:elemental )?(?:subroutine|function) '+n+r'\(.*?^  end (?:subroutine|function) '+n+r'\b',source,re.M|re.S).group() for n in names}
def powexpr(s):
 while '**' in s:
  pos=s.rindex('**');end=pos
  while s[end-1].isspace():end-=1
  start=end-1
  if s[start]==')':
   lev=1;start-=1
   while lev:lev+=(s[start]==')')-(s[start]=='(');start-=1
   start+=1;j=start-1
   while j>=0 and (s[j].isalnum() or s[j] in '_.'):j-=1
   if j<start-1:start=j+1
  else:
   while start>=0 and (s[start].isalnum() or s[start] in '_.'):start-=1
   start+=1
  i=pos+2
  while s[i].isspace():i+=1
  expstart=i
  if s[i]=='(':
   lev=1;i+=1
   while lev:lev+=(s[i]=='(')-(s[i]==')');i+=1
  else:
   while i<len(s) and (s[i].isalnum() or s[i] in '_.'):i+=1
  exponent=s[expstart:i];base=s[start:end]
  result='('+base+' * '+base+')' if exponent=='2' else 'Kokkos::pow('+base+', '+exponent+')'
  s=s[:start]+result+s[i:]
 return s

def expr(s):
 s=s.replace('&',' ')
 # Powers are expanded before kind-suffix conversion, so Native real bases
 # and parenthesized exponents retain their complete source operands.
 s=powexpr(s)
 s=re.sub(r'((?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)_core_rknd',lambda m:'real_t('+m[1]+')',s)
 for a,b in [('.and.','&&'),('.or.','||'),('.not.','!'),('.true.','true'),('.false.','false'),('/=','!=')]:s=re.sub(re.escape(a),b,s,flags=re.I)
 for f in ['sqrt','exp','log10','log','abs']:s=re.sub(r'(?<![\w:])'+f+r'\s*\(','Kokkos::'+f+'(',s,flags=re.I)
 for f,op in [('min','<'),('max','>')]:
  while m:=re.search(r'(?<![\w:])'+f+r'\s*\(',s,re.I):
   pos=m.end();d=1
   while d:d+=(s[pos]=='(')-(s[pos]==')');pos+=1
   a=split_top(s[m.end():pos-1]);assert len(a)==2,s
   value='(('+a[0]+') '+op+' ('+a[1]+') ? ('+a[0]+') : ('+a[1]+'))'
   s=s[:m.start()]+value+s[pos:]
 s=re.sub(r'int\(\s*anint\(\s*T_in_K\s*\)\s*\)','static_cast<int>(Kokkos::round(T_in_K))',s,flags=re.I)
 return re.sub(r'\s+',' ',s).strip()
def translate(sec):
 out=[];pending=''
 for raw in sec.splitlines():
  if raw.strip().startswith('!$acc'):continue
  if raw.strip().startswith('!'):out.append('  //'+raw.strip()[1:]);continue
  code=raw.split('!')[0].strip()
  if not code:continue
  pending+=' '+code.strip('& ')
  if code.endswith('&'):continue
  st=pending.strip();pending=''
  if st.lower()=='return':continue
  if re.match(r'select case',st,re.I):out.append('  switch (saturation_formula) {');continue
  if re.match(r'end select',st,re.I):out.append('  break;\n  }');continue
  m=re.match(r'case\s*\(\s*(\w+)\s*\)',st,re.I)
  if m:
   if any('case ' in a for a in out):out.append('  break;')
   out.append('  case '+m[1]+':');continue
  if st.lower()=='case default':out.append('  break;\n  default:');continue
  if re.match(r'end\s*if',st,re.I):out.append('  }');continue
  if re.match(r'else\s*$',st,re.I):out.append('  } else {');continue
  m=re.match(r'(else\s*if|elseif|if)\s*\((.*)\)\s*then',st,re.I)
  if m:out.append(('  } else if' if m[1].lower()!='if' else '  if')+' ('+expr(m[2])+') {');continue
  if st.lower().startswith('call '):out.append('  '+expr(st[5:])+';');continue
  if '=' in st:out.append('  '+expr(st)+';');continue
  raise ValueError(st)
 assert not pending
 return '\n'.join(out)+'\n'
def bodypart(n,start):
 b=blocks[n];return b[b.index(start):b.rindex('    return')] if '    return' in b else b[b.index(start):b.index('  end function')]
def desc(n):
 b=blocks[n];m=re.search(r'^\s*! Description:',b,re.M)
 if not m:return ''
 end=re.search(r'^\s*(?:use|implicit none)\b',b[m.end():],re.M).start()+m.end()
 return '\n'.join('  //'+l.strip()[1:] if l.strip().startswith('!') else '' for l in b[m.start():end].splitlines())+'\n'
# Exact immutable Native table, emitted as a device getter to avoid host-only
# constexpr array addresses or mutable global initialization in scalar kernels.
rawtable=source[source.index('    (/ 0.049560547'):source.index(' /)',source.index('    (/ 0.049560547'))]
values=re.findall(r'([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?\d+)?)_core_rknd',rawtable)
assert len(values)==156,len(values)
lookup='KOKKOS_INLINE_FUNCTION real_t svp_liq_lookup_table(int T_in_K_int) {\n  // Return vapor pressure [Pa] from the immutable Native 188:343 K table.\n  // Switch literals avoid host-only table addresses in device scalar calls.\n  switch (T_in_K_int) {\n'
for i,val in enumerate(values):lookup+='  case '+str(i+188)+': return real_t('+val+');\n'
lookup+='  default: Kokkos::abort("saturation lookup index out of range"); return real_t(0);\n  }\n}\n'
# Native liquid mixing-ratio scalar and array SVP bodies are identical.
liquid=blocks['sat_mixrat_liq_k'];start=liquid.index('    ! Calculate the SVP');end=liquid.index('    ! If esat exceeds',start)
flatau_parameters='\n'.join('  //'+line.strip()[1:] for line in liquid[liquid.index('    ! Relative error norm expansion'):liquid.index('    real( kind = core_rknd ), parameter ::')].splitlines() if line.strip().startswith('!'))+'\n'
liquid_esat='''KOKKOS_INLINE_FUNCTION real_t liquid_esat_for_mixrat(real_t T_in_K, int saturation_formula) {
  // Shared device-cell extraction of the Native scalar/2D liquid SVP blocks.
  // Temperature [K] and the formula selector give vapor pressure [Pa].
  // Its GFDL 173.15 K floor intentionally differs from the public SVP routine.
'''+flatau_parameters+'''
  real_t T_in_C, T_in_C_sqd, T_in_K_clipped, esat;
  int T_in_K_int;
  constexpr real_t min_T_in_C = real_t(-85.0);
  constexpr real_t min_T_in_K = real_t(173.15);
'''+translate(liquid[start:end])+'  return esat;\n}\n'
head='''// Saturation vapor pressures, liquid/ice mixing ratios and cloud-water adjustment.
//
// sat_mixrat_liq_api and sat_mixrat_ice accept scalar pressure [Pa], temperature [K]
// and a Native formula selector, or column-first Kokkos View fields. The liquid
// array API additionally accepts grid metadata and an optional Native 1-based
// start index. sat_vapor_press_liq returns vapor pressure [Pa]; rcm_sat_adj returns
// cloud-water mixing ratio [kg/kg]. Scalar forms are device callable. Array outputs
// retain the pressure input's memory space; unselected liquid levels are unwritten.
// Native selectors/errors and both specific-humidity branches are retained.
// References: Emanuel 1994, Eqs. 4.4.14–15; Flatau et al. 1992, JAM 31, 1507–1513;
// Bolton 1980; Goff/Gratch 1946. See DRAFT_STATUS.md before integration.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "constants_clubb.hpp"
#include "model_flags.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <optional>
#include <type_traits>
namespace clubb_kokkos {
namespace saturation_detail {
'''+lookup+'\n'+liquid_esat+'''
KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq_flatau(real_t T_in_K, real_t& esat);
KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq_bolton(real_t T_in_K, real_t& esat);
KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq_gfdl(real_t T_in_K, real_t& esat);
KOKKOS_INLINE_FUNCTION real_t sat_vapor_press_liq_lookup(real_t T_in_K);
} // namespace saturation_detail

'''
head+='''KOKKOS_INLINE_FUNCTION real_t sat_mixrat_liq_k(real_t p_in_Pa, real_t T_in_K, int saturation_formula) {
'''+desc('sat_mixrat_liq_k')+'''  // Input pressure [Pa], temperature [K] and selector; return mixing ratio [kg/kg].
  // -------------------- Begin Code --------------------
  const real_t esat = saturation_detail::liquid_esat_for_mixrat(T_in_K, saturation_formula);
  real_t result;
  if (p_in_Pa-esat < real_t(1)) {
    // If esat exceeds the air pressure, then assume esat~=0.5*pressure
    // and set rsat = ep = 0.622.
    result = ep;
  } else {
    if (I_sat_sphum) {
      // GFDL uses specific humidity. Retain the Native scalar parentheses.
      result = ep * (esat / (p_in_Pa - (real_t(1)-ep)*esat));
    } else {
      // Saturation mixing ratio: rs = epsilon*esat/(p-esat), epsilon = Rd/Rv.
      result = ep * esat / (p_in_Pa-esat);
    }
  }
  return result;
}
KOKKOS_INLINE_FUNCTION real_t sat_mixrat_liq_api(real_t p_in_Pa, real_t T_in_K, int saturation_formula) {
  // Return liquid saturation mixing ratio [kg/kg] from pressure [Pa],
  // temperature [K] and selector through the Native scalar generic.
  return sat_mixrat_liq_k(p_in_Pa,T_in_K,saturation_formula);
}

'''
# All array calculations are independent cells; fuse SVP and conversion without
# changing arithmetic inside either source expression or dropping selectors.
head+='''template <class Pressure, class Temperature>
auto sat_mixrat_liq_2D(int nz, int ngrdcol, const grid& gr, Pressure p_in_Pa, Temperature T_in_K,
                       int saturation_formula, std::optional<int> start_index_in=std::nullopt) {
'''+desc('sat_mixrat_liq_2D')+'''  // Input column-first pressure [Pa] and temperature [K] fields on gr;
  // return mixing ratio [kg/kg]. The optional start index remains 1-based.
  // Levels outside the selected range are unwritten in Native.
  // -------------------- Begin Code --------------------
  // start_index is optional and used for choosing the sub-arrays.
  int start_index=0, stop_index=nz;
  if (start_index_in) {
    if (gr.grid_dir_indx>0) {
      // Ascending grid.
      start_index=*start_index_in-1;
    } else {
      // Descending grid.
      stop_index=*start_index_in;
    }
  }
  auto pressure=to_execution_space(p_in_Pa);
  auto temperature=to_execution_space(T_in_K);
  Real2D result(Kokkos::view_alloc(Kokkos::WithoutInitializing,"sat_mixrat_liq_2D"),ngrdcol,nz);
  // Fuse independent SVP and conversion loops without an intermediate esat View.
  // An empty Native level range launches no work.
  const int cells=start_index<stop_index?ngrdcol*(stop_index-start_index):0;
  Kokkos::parallel_for("sat_mixrat_liq_2D",cells,
    KOKKOS_LAMBDA(const int cell) {
      const int i=cell%ngrdcol; const int k=start_index+cell/ngrdcol;
      const real_t esat=saturation_detail::liquid_esat_for_mixrat(temperature(i,k),saturation_formula);
      if (pressure(i,k)-esat < real_t(1)) result(i,k)=ep;
      else {
        if (I_sat_sphum) {
          // GFDL uses specific humidity. Retain the source 2D denominator,
          // which differs from the scalar expression (Native lines 529–530).
          result(i,k)=ep*(esat/(pressure(i,k)-(real_t(1)-ep*esat)));
        } else result(i,k)=ep*esat/(pressure(i,k)-esat);
      }
    });
  if constexpr (std::is_same_v<typename Pressure::memory_space,typename Real2D::memory_space>) return result;
  else {
    Kokkos::View<real_t**,Kokkos::LayoutLeft,typename Pressure::memory_space> output(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,"sat_mixrat_liq_2D"),ngrdcol,nz);
    from_execution_space(output,result);
    return output;
  }
}
template <class Pressure, class Temperature>
auto sat_mixrat_liq_api(int nz,int ngrdcol,const grid& gr,Pressure p_in_Pa,Temperature T_in_K,
                       int saturation_formula,std::optional<int> start_index_in=std::nullopt) {
  // Return the column-first liquid saturation field [kg/kg] through the Native
  // array generic, retaining optional 1-based start_index_in and grid direction.
  return sat_mixrat_liq_2D(nz,ngrdcol,gr,p_in_Pa,T_in_K,saturation_formula,start_index_in);
}

KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq(real_t T_in_K,int saturation_formula,real_t& esat) {
'''+desc('sat_vapor_press_liq')+'''  // Input temperature [K] and selector; write vapor pressure [Pa].
'''+translate(blocks['sat_vapor_press_liq'][blocks['sat_vapor_press_liq'].index('    ! Saturation Vapor Pressure'):blocks['sat_vapor_press_liq'].rindex('    return')]).replace('sat_vapor_press_liq_bolton(', 'saturation_detail::sat_vapor_press_liq_bolton(').replace('sat_vapor_press_liq_flatau(', 'saturation_detail::sat_vapor_press_liq_flatau(').replace('sat_vapor_press_liq_gfdl(', 'saturation_detail::sat_vapor_press_liq_gfdl(').replace('sat_vapor_press_liq_lookup(', 'saturation_detail::sat_vapor_press_liq_lookup(')+'''
}

namespace saturation_detail {
'''
for n in ['sat_vapor_press_liq_flatau','sat_vapor_press_liq_bolton','sat_vapor_press_liq_gfdl']:
 startmarker='    ! ---------------------- Begin Code' if n.endswith('flatau') else '    ! --------------------------- Begin Code'
 if n.endswith('bolton'):startmarker='    ! --------------------------- Begin Code'
 sec=bodypart(n,startmarker)
 const='  constexpr real_t min_T_in_C=real_t(-85.0);\n  real_t T_in_C,T_in_C_sqd;\n' if n.endswith('flatau') else '  constexpr real_t min_T_in_K=real_t(203.15);\n  real_t T_in_K_clipped;\n' if n.endswith('gfdl') else ''
 head+='KOKKOS_INLINE_FUNCTION void '+n+'(real_t T_in_K,real_t& esat) {\n'+desc(n)+'  // Input temperature [K]; write vapor pressure [Pa].\n'+(flatau_parameters if n.endswith('flatau') else '')+const+translate(sec)+'}\n\n'
n='sat_vapor_press_liq_lookup';sec=bodypart(n,'    ! ---- Begin Code')
head+='KOKKOS_INLINE_FUNCTION real_t '+n+'(real_t T_in_K) {\n'+desc(n)+'  // Input temperature [K]; return tabulated vapor pressure [Pa].\n  int T_in_K_int;\n  real_t esat;\n'+translate(sec)+'  return esat;\n}\n} // namespace saturation_detail\n\n'
# Scalar ice directly extracts the source2D cell calculation. No allocation
# inside a device parcel call; array generic reuses exactly the same cell body.
ice=blocks['sat_mixrat_ice_2D'];sec=ice[ice.index('    ! Determine the SVP'):ice.index('    !$acc parallel loop gang vector collapse(2) default(present)\n    do k = 1, nz',ice.index('    end select'))]
sec='\n'.join(l for l in sec.splitlines() if not re.match(r'\s*(?:do\s+[ik]\s*=|end\s*do)',l,re.I))
sec=sec.replace('T_in_K(i,k)','T_in_K').replace('esat_ice(i,k)','esat_ice')
# Native a coefficients are multiplied by100 at core precision.
a=[x for x in re.findall(r'([\d.]+(?:[eE][+-]?\d+)?)_core_rknd',ice[ice.index('100._core_rknd * (/'):ice.index('/)',ice.index('100._core_rknd * (/'))])][1:]
assert len(a)==9,a
sec='\n'.join(line if line.strip().startswith('!') else re.sub(r'a\((\d)\)',lambda m:'(real_t(100)*real_t('+a[int(m[1])-1]+'))',line) for line in sec.splitlines())
head+='''KOKKOS_INLINE_FUNCTION real_t sat_mixrat_ice_k(real_t p_in_Pa,real_t T_in_K,int saturation_formula) {
'''+desc('sat_mixrat_ice_k')+'''  // Input pressure [Pa], temperature [K] and selector; return mixing ratio [kg/kg].
  // Device adaptation: inline the Native 1x1 array call's cell body instead
  // of allocating arrays in each parcel. Array overload uses this same body.
  // Flatau ice coefficients use the relative-error expansion (-90 to 0 deg_C),
  // Table 4, p. 1511 of Flatau et al. 1992; factor 100 converts mb to Pa.
  constexpr real_t min_T_in_C=real_t(-90.0), min_T_in_K=real_t(173.15);
  real_t T_in_C,T_in_K_clipped,esat_ice,result;
'''+translate(sec)+'''  if (p_in_Pa-esat_ice<real_t(1)) result=ep;
  else {
    if (I_sat_sphum) result=ep*(esat_ice/(p_in_Pa-(real_t(1)-ep)*esat_ice));
    else result=ep*esat_ice/(p_in_Pa-esat_ice);
  }
  return result;
}
KOKKOS_INLINE_FUNCTION real_t sat_mixrat_ice(real_t p_in_Pa,real_t T_in_K,int saturation_formula) {
  // Return ice saturation mixing ratio [kg/kg] from pressure [Pa],
  // temperature [K] and selector through the Native scalar generic.
  return sat_mixrat_ice_k(p_in_Pa,T_in_K,saturation_formula);
}

template <class Pressure,class Temperature>
auto sat_mixrat_ice_2D(int nz,int ngrdcol,Pressure p_in_Pa,Temperature T_in_K,int saturation_formula) {
'''+desc('sat_mixrat_ice_2D')+'''  // Input column-first pressure [Pa] and temperature [K] fields;
  // return mixing ratio [kg/kg] in the pressure field's memory space.
  // ------------------------ Begin Code ------------------------
  auto pressure=to_execution_space(p_in_Pa);
  auto temperature=to_execution_space(T_in_K);
  Real2D result(Kokkos::view_alloc(Kokkos::WithoutInitializing,"sat_mixrat_ice_2D"),ngrdcol,nz);
  Kokkos::parallel_for("sat_mixrat_ice_2D",ngrdcol*nz,KOKKOS_LAMBDA(const int cell) {
    const int i=cell%ngrdcol;const int k=cell/ngrdcol;
    result(i,k)=sat_mixrat_ice_k(pressure(i,k),temperature(i,k),saturation_formula);
  });
  if constexpr (std::is_same_v<typename Pressure::memory_space,typename Real2D::memory_space>) return result;
  else {
    Kokkos::View<real_t**,Kokkos::LayoutLeft,typename Pressure::memory_space> output(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,"sat_mixrat_ice_2D"),ngrdcol,nz);
    from_execution_space(output,result);return output;
  }
}
template <class Pressure,class Temperature>
auto sat_mixrat_ice(int nz,int ngrdcol,Pressure p_in_Pa,Temperature T_in_K,int saturation_formula) {
  // Return the column-first ice saturation field [kg/kg] through the Native
  // array generic for pressure [Pa], temperature [K] and the formula selector.
  return sat_mixrat_ice_2D(nz,ngrdcol,p_in_Pa,T_in_K,saturation_formula);
}

KOKKOS_INLINE_FUNCTION real_t rcm_sat_adj(real_t thlm,real_t rtm,real_t p_in_Pa,real_t exner,
                                        int saturation_formula) {
'''+desc('rcm_sat_adj')+'''  // Input liquid-water potential temperature thlm [K], total-water mixing
  // ratio rtm [kg/kg], pressure [Pa], Exner function [-] and formula selector.
  // Return cloud-water mixing ratio rcm [kg/kg].
  constexpr real_t tolerance=real_t(0.001); // Theta calculation tolerance [K].
  constexpr int itermax=1000000;
  real_t rcm,theta=thlm,answer,too_low=real_t(0),too_high=real_t(0);
  int iteration;
  // ----- Begin Code -----
  // Default theta/bracket initialization is retained above.
  for (iteration=1;iteration<=itermax;++iteration) {
    const real_t excess=rtm-sat_mixrat_liq_api(p_in_Pa,theta*exner,saturation_formula);
    answer=theta-(Lv/(Cp*exner))*(excess>zero_threshold?excess:zero_threshold);
    if (Kokkos::abs(answer-thlm)<=tolerance) break;
    else if (answer-thlm>tolerance) too_high=theta;
    else if (thlm-answer>tolerance) too_low=theta;
    // First iteration must set an upper bound that is way too high.
    if (iteration==1) too_high=theta+real_t(20);
    theta=(too_low+too_high)/real_t(2);
  }
  // Preserve the exact Native loop-counter check, including its last-iteration
  // exit versus natural-exhaustion distinction; this is not JAX's fixed scan.
  if (iteration==itermax) {
    rcm=real_t(0);
    Kokkos::abort("Error in rcm_sat_adj: could not determine rcm");
  } else {
    const real_t excess=rtm-sat_mixrat_liq_api(p_in_Pa,theta*exner,saturation_formula);
    rcm=excess>zero_threshold?excess:zero_threshold;
  }
  return rcm;
}

} // namespace clubb_kokkos
'''
(folder/'saturation.hpp').write_text(head)
(folder/'lookup_tokens.txt').write_text('\n'.join(values)+'\n')
