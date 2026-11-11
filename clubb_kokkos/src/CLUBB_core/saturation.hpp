// Saturation vapor pressures, liquid/ice mixing ratios and cloud-water adjustment.
//
// sat_mixrat_liq_api and sat_mixrat_ice accept scalar pressure [Pa], temperature [K]
// and a Native formula selector, or column-first Kokkos View fields. The liquid
// array API additionally accepts grid metadata and an optional zero-based
// start index. sat_vapor_press_liq returns vapor pressure [Pa]; rcm_sat_adj returns
// cloud-water mixing ratio [kg/kg]. Scalar forms are device callable. Array outputs
// retain the pressure input's memory space; unselected liquid levels are unwritten.
// Native selectors/errors and both specific-humidity branches are retained.
// References: Emanuel 1994, Eqs. 4.4.14–15; Flatau et al. 1992, JAM 31, 1507–1513;
// Bolton 1980; Goff/Gratch 1946.
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
KOKKOS_INLINE_FUNCTION real_t svp_liq_lookup_table(int T_in_K_int) {
  // Return vapor pressure [Pa] from the immutable Native 188:343 K table.
  // Switch literals avoid host-only table addresses in device scalar calls.
  switch (T_in_K_int) {
  case 188:
    return real_t(0.049560547);
  case 189:
    return real_t(0.059753418);
  case 190:
    return real_t(0.070129395);
  case 191:
    return real_t(0.083618164);
  case 192:
    return real_t(0.09814453);
  case 193:
    return real_t(0.11444092);
  case 194:
    return real_t(0.13446045);
  case 195:
    return real_t(0.15686035);
  case 196:
    return real_t(0.18218994);
  case 197:
    return real_t(0.21240234);
  case 198:
    return real_t(0.24725342);
  case 199:
    return real_t(0.28668213);
  case 200:
    return real_t(0.33184814);
  case 201:
    return real_t(0.3826294);
  case 202:
    return real_t(0.4416504);
  case 203:
    return real_t(0.50775146);
  case 204:
    return real_t(0.58343506);
  case 205:
    return real_t(0.6694946);
  case 206:
    return real_t(0.7668457);
  case 207:
    return real_t(0.87750244);
  case 208:
    return real_t(1.0023804);
  case 209:
    return real_t(1.1434937);
  case 210:
    return real_t(1.3028564);
  case 211:
    return real_t(1.482544);
  case 212:
    return real_t(1.6847534);
  case 213:
    return real_t(1.9118042);
  case 214:
    return real_t(2.1671143);
  case 215:
    return real_t(2.4535522);
  case 216:
    return real_t(2.774231);
  case 217:
    return real_t(3.1330566);
  case 218:
    return real_t(3.5343628);
  case 219:
    return real_t(3.9819336);
  case 220:
    return real_t(4.480713);
  case 221:
    return real_t(5.036072);
  case 222:
    return real_t(5.6540527);
  case 223:
    return real_t(6.340088);
  case 224:
    return real_t(7.1015015);
  case 225:
    return real_t(7.9450684);
  case 226:
    return real_t(8.8793335);
  case 227:
    return real_t(9.91217);
  case 228:
    return real_t(11.053528);
  case 229:
    return real_t(12.313049);
  case 230:
    return real_t(13.70166);
  case 231:
    return real_t(15.231018);
  case 232:
    return real_t(16.91394);
  case 233:
    return real_t(18.764038);
  case 234:
    return real_t(20.795898);
  case 235:
    return real_t(23.025574);
  case 236:
    return real_t(25.470093);
  case 237:
    return real_t(28.147766);
  case 238:
    return real_t(31.078003);
  case 239:
    return real_t(34.282043);
  case 240:
    return real_t(37.782593);
  case 241:
    return real_t(41.60382);
  case 242:
    return real_t(45.771606);
  case 243:
    return real_t(50.31366);
  case 244:
    return real_t(55.259644);
  case 245:
    return real_t(60.641174);
  case 246:
    return real_t(66.492004);
  case 247:
    return real_t(72.84802);
  case 248:
    return real_t(79.74756);
  case 249:
    return real_t(87.23126);
  case 250:
    return real_t(95.34259);
  case 251:
    return real_t(104.12747);
  case 252:
    return real_t(113.634796);
  case 253:
    return real_t(123.91641);
  case 254:
    return real_t(135.02725);
  case 255:
    return real_t(147.02563);
  case 256:
    return real_t(159.97308);
  case 257:
    return real_t(173.93488);
  case 258:
    return real_t(188.97995);
  case 259:
    return real_t(205.18109);
  case 260:
    return real_t(222.61517);
  case 261:
    return real_t(241.36334);
  case 262:
    return real_t(261.51108);
  case 263:
    return real_t(283.14853);
  case 264:
    return real_t(306.37054);
  case 265:
    return real_t(331.27698);
  case 266:
    return real_t(357.97278);
  case 267:
    return real_t(386.56842);
  case 268:
    return real_t(417.17978);
  case 269:
    return real_t(449.9286);
  case 270:
    return real_t(484.94254);
  case 271:
    return real_t(522.3556);
  case 272:
    return real_t(562.30804);
  case 273:
    return real_t(604.947);
  case 274:
    return real_t(650.42645);
  case 275:
    return real_t(698.9074);
  case 276:
    return real_t(750.55835);
  case 277:
    return real_t(805.55554);
  case 278:
    return real_t(864.0828);
  case 279:
    return real_t(926.3325);
  case 280:
    return real_t(992.5052);
  case 281:
    return real_t(1062.8102);
  case 282:
    return real_t(1137.4657);
  case 283:
    return real_t(1216.6995);
  case 284:
    return real_t(1300.7483);
  case 285:
    return real_t(1389.8594);
  case 286:
    return real_t(1484.2896);
  case 287:
    return real_t(1584.3064);
  case 288:
    return real_t(1690.1881);
  case 289:
    return real_t(1802.224);
  case 290:
    return real_t(1920.7146);
  case 291:
    return real_t(2045.9724);
  case 292:
    return real_t(2178.3218);
  case 293:
    return real_t(2318.099);
  case 294:
    return real_t(2465.654);
  case 295:
    return real_t(2621.3489);
  case 296:
    return real_t(2785.5596);
  case 297:
    return real_t(2958.6758);
  case 298:
    return real_t(3141.101);
  case 299:
    return real_t(3333.2534);
  case 300:
    return real_t(3535.5657);
  case 301:
    return real_t(3748.4863);
  case 302:
    return real_t(3972.4792);
  case 303:
    return real_t(4208.024);
  case 304:
    return real_t(4455.616);
  case 305:
    return real_t(4715.7686);
  case 306:
    return real_t(4989.0127);
  case 307:
    return real_t(5275.8945);
  case 308:
    return real_t(5576.9795);
  case 309:
    return real_t(5892.8535);
  case 310:
    return real_t(6224.116);
  case 311:
    return real_t(6571.3926);
  case 312:
    return real_t(6935.3213);
  case 313:
    return real_t(7316.5674);
  case 314:
    return real_t(7715.8105);
  case 315:
    return real_t(8133.755);
  case 316:
    return real_t(8571.125);
  case 317:
    return real_t(9028.667);
  case 318:
    return real_t(9507.15);
  case 319:
    return real_t(10007.367);
  case 320:
    return real_t(10530.132);
  case 321:
    return real_t(11076.282);
  case 322:
    return real_t(11646.683);
  case 323:
    return real_t(12242.221);
  case 324:
    return real_t(12863.808);
  case 325:
    return real_t(13512.384);
  case 326:
    return real_t(14188.913);
  case 327:
    return real_t(14894.385);
  case 328:
    return real_t(15629.823);
  case 329:
    return real_t(16396.268);
  case 330:
    return real_t(17194.799);
  case 331:
    return real_t(18026.516);
  case 332:
    return real_t(18892.55);
  case 333:
    return real_t(19794.07);
  case 334:
    return real_t(20732.262);
  case 335:
    return real_t(21708.352);
  case 336:
    return real_t(22723.592);
  case 337:
    return real_t(23779.273);
  case 338:
    return real_t(24876.709);
  case 339:
    return real_t(26017.258);
  case 340:
    return real_t(27202.3);
  case 341:
    return real_t(28433.256);
  case 342:
    return real_t(29711.578);
  case 343:
    return real_t(31038.766);
  default:
    Kokkos::abort("saturation lookup index out of range");
    return real_t(0);
  }
}

KOKKOS_INLINE_FUNCTION real_t liquid_esat_for_mixrat(real_t T_in_K, int saturation_formula) {
  // Shared device-cell extraction of the Native scalar/2D liquid SVP blocks.
  // Temperature [K] and the formula selector give vapor pressure [Pa].
  // Its GFDL 173.15 K floor intentionally differs from the public SVP routine.
  // Relative error norm expansion (-50 to 50 deg_C) from
  // Table 3 of pp. 1510 of Flatau et al. 1992 (Water Vapor)
  // (The 100 coefficient converts from mb to Pa)
  //   real, dimension(7), parameter :: a = &
  //   100.* (/ 6.11176750,      0.443986062,     0.143053301E-01, &
  //            0.265027242E-03, 0.302246994E-05, 0.203886313E-07, &
  //            0.638780966E-10 /)
  // Relative error norm expansion (-85 to 70 deg_C) from
  // Table 4 of pp. 1511 of Flatau et al.
  //real( kind = core_rknd ), dimension(9), parameter :: a = &
  //100._core_rknd * &
  //  Commented out because the form has been redone, causing these number to no longer be needed,
  //  leaving them in for now for reference.
  //         (/ 6.11583699_core_rknd,      0.444606896_core_rknd,     0.143177157E-01_core_rknd, &
  //         0.264224321E-03_core_rknd, 0.299291081E-05_core_rknd, 0.203154182E-07_core_rknd, &
  //         0.702620698E-10_core_rknd, 0.379534310E-13_core_rknd,-0.321582393E-15_core_rknd /)

  real_t T_in_C, T_in_C_sqd, T_in_K_clipped, esat;
  int T_in_K_int;
  constexpr real_t min_T_in_C = real_t(-85.0);
  constexpr real_t min_T_in_K = real_t(173.15);
  // Calculate the SVP for water vapor.
  switch (saturation_formula) {
  case saturation_flatau:
    // Using the Flatau, et al. polynomial approximation for SVP over vapor
    // Determine deg K - 273.15
    T_in_C = T_in_K - T_freeze_K;
    // Since this approximation is only good out to -85 degrees Celsius we
    // truncate the result here (Flatau, et al. 1992)
    T_in_C = ((T_in_C) > (min_T_in_C) ? (T_in_C) : (min_T_in_C));
    // Polynomial approx. (Flatau, et al. 1992)
    // This is the generalized formula but is not computationally efficient.
    // Based on Wexler's expressions(2.1)-(2.4) (See Flatau et al. p 1508)
    // e_{sat} = a_1 + a_2 ( T - T_0 ) + ... + a_{n+1} ( T - T_0 )^n
    // esat = a(1)
    // do i = 2, size( a ) , 1
    //   esat = esat + a(i) * ( T_in_C )**(i-1)
    // end do
    // The 8th order polynomial fit.  When running deep
    // convective cases I noticed that absolute temperature often dips below
    // -50 deg_C at higher altitudes, where the 6th order approximation is
    // not accurate.  -dschanen 20 Nov 2008
    //esat = a(1) + T_in_C*( a(2) + T_in_C*( a(3) + T_in_C*( a(4) + T_in_C &
    //*( a(5) + T_in_C*( a(6) + T_in_C*( a(7) + T_in_C*( a(8) + T_in_C*( a(9) ) ) ) ) ) ) ) )
    // Factoring the polynomial above and changing it into this form allows the cpu
    // to complete the calculations out of order. This is because modern cpus can complete
    // multiple instructions at once if they do not depend on eachother, in the above case
    // each instruction relies on the result of the last. In this version however, the terms
    // in the parentheses could potentially be calculated in parallel by different execution
    // units in the cpu, then only when those terms are being multiplied together do the
    // instructions need to be done one at a time. See clubb issue 834 for more info.
    //   - Gunther Huebler, Aug 2018
    T_in_C_sqd = (T_in_C * T_in_C);
    esat       = -real_t(3.21582393e-14) * (T_in_C - real_t(646.5835252598777)) *
                 (T_in_C + real_t(90.72381630364440)) *
                 (T_in_C_sqd + real_t(111.0976961559954) * T_in_C + real_t(6459.629194243118)) *
                 (T_in_C_sqd + real_t(152.3131930092453) * T_in_C + real_t(6499.774954705265)) *
                 (T_in_C_sqd + real_t(174.4279584934021) * T_in_C + real_t(7721.679732114084));
    break;
  case saturation_bolton:
    // Using the Bolton 1980 approximations for SVP over vapor
    // Generally this more computationally expensive than the Flatau polnomial expansion
    esat = real_t(611.2) *
           Kokkos::exp((real_t(17.67) * (T_in_K - T_freeze_K)) / (T_in_K - real_t(29.65)));
    // ---> h1g
    break;
  case saturation_gfdl:
    // Using GFDL polynomial approximation for SVP with respect to liquid
    // Since the Goff-Gratch approximation is valid only down to -70 degrees Celsius,
    //   we threshold the temperature.  This will yield a minimal saturation at
    //   cold temperatures.
    T_in_K_clipped = ((min_T_in_K) > (T_in_K) ? (min_T_in_K) : (T_in_K));
    // Goff Gratch equation, uncertain below -70 C
    esat =
        Kokkos::pow(
            real_t(10.),
            (-real_t(7.90298) * (real_t(373.16) / T_in_K_clipped - real_t(1.)) +
             real_t(5.02808) * Kokkos::log10(real_t(373.16) / T_in_K_clipped) -
             real_t(1.3816e-7) *
                 (Kokkos::pow(real_t(10.),
                              (real_t(11.344) * (real_t(1.) - T_in_K_clipped / real_t(373.16)))) -
                  real_t(1.)) +
             real_t(8.1328e-3) *
                 (Kokkos::pow(real_t(10.),
                              (-real_t(3.49149) * (real_t(373.16) / T_in_K_clipped - real_t(1.)))) -
                  real_t(1.)) +
             Kokkos::log10(real_t(1013.246)))) *
        real_t(100.);
    // <--- h1g
    break;
  case saturation_lookup:
    T_in_K_int = static_cast<int>(Kokkos::round(T_in_K));
    // Since this approximation is only good out to -85 degrees Celsius we
    // truncate the result here
    T_in_K_int = ((((T_in_K_int) > (188) ? (T_in_K_int) : (188))) < (343)
                      ? (((T_in_K_int) > (188) ? (T_in_K_int) : (188)))
                      : (343));
    // Use the lookup table to determine the saturation vapor pressure.
    esat = svp_liq_lookup_table(T_in_K_int);
    break;
  default:
    // Undefined approximation
    esat = -real_t(99999.999);
    break;
  }
  return esat;
}

KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq_flatau(real_t T_in_K, real_t& esat);
KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq_bolton(real_t T_in_K, real_t& esat);
KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq_gfdl(real_t T_in_K, real_t& esat);
KOKKOS_INLINE_FUNCTION real_t sat_vapor_press_liq_lookup(real_t T_in_K);
} // namespace saturation_detail

KOKKOS_INLINE_FUNCTION real_t sat_mixrat_liq_k(real_t p_in_Pa, real_t T_in_K,
                                               int saturation_formula) {
  // Description:
  //   Used to compute the saturation mixing ratio of liquid water.

  // References:
  //   Formula from Emanuel 1994, 4.4.14
  //-------------------------------------------------------------------------
  // Input pressure [Pa], temperature [K] and selector; return mixing ratio [kg/kg].
  // -------------------- Begin Code --------------------
  const real_t esat = saturation_detail::liquid_esat_for_mixrat(T_in_K, saturation_formula);
  real_t result;
  if (p_in_Pa - esat < real_t(1)) {
    // If esat exceeds the air pressure, then assume esat~=0.5*pressure
    // and set rsat = ep = 0.622.
    result = ep;
  } else {
    if (I_sat_sphum) {
      // GFDL uses specific humidity. Retain the Native scalar parentheses.
      result = ep * (esat / (p_in_Pa - (real_t(1) - ep) * esat));
    } else {
      // Saturation mixing ratio: rs = epsilon*esat/(p-esat), epsilon = Rd/Rv.
      result = ep * esat / (p_in_Pa - esat);
    }
  }
  return result;
}
KOKKOS_INLINE_FUNCTION real_t sat_mixrat_liq_api(real_t p_in_Pa, real_t T_in_K,
                                                 int saturation_formula) {
  // Return liquid saturation mixing ratio [kg/kg] from pressure [Pa],
  // temperature [K] and selector through the Native scalar generic.
  return sat_mixrat_liq_k(p_in_Pa, T_in_K, saturation_formula);
}

template <class Pressure, class Temperature>
void sat_mixrat_liq_2D(int nz, int ngrdcol, const grid& gr, Pressure p_in_Pa, Temperature T_in_K,
                       int saturation_formula, Real2D result,
                       std::optional<int> start_index_in = std::nullopt) {
  // Description:
  //   Used to compute the saturation mixing ratio of liquid water.
  //
  // References:
  //   Formula from Emanuel 1994, 4.4.14
  //-------------------------------------------------------------------------
  // Input column-first pressure [Pa] and temperature [K] fields on gr;
  // return mixing ratio [kg/kg]. The optional start index is a zero-based inclusive slot.
  // Levels outside the selected range are unwritten in Native.
  // -------------------- Begin Code --------------------
  // start_index is optional and used for choosing the sub-arrays.
  int start_index = 0, stop_index = nz;
  if (start_index_in) {
    if (gr.grid_dir_indx > 0) {
      // Ascending grid.
      start_index = *start_index_in;
    } else {
      // Descending grid.
      stop_index = *start_index_in + 1;
    }
  }
  auto pressure    = to_execution_space(p_in_Pa);
  auto temperature = to_execution_space(T_in_K);
  // Fuse independent SVP and conversion loops without an intermediate esat View.
  // An empty Native level range launches no work.
  const int cells = start_index < stop_index ? ngrdcol * (stop_index - start_index) : 0;
  Kokkos::parallel_for(
      "sat_mixrat_liq_2D", cells, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = start_index + cell / ngrdcol;
        const real_t esat =
            saturation_detail::liquid_esat_for_mixrat(temperature(i, k), saturation_formula);
        if (pressure(i, k) - esat < real_t(1))
          result(i, k) = ep;
        else {
          if (I_sat_sphum) {
            // GFDL uses specific humidity. Retain the source 2D denominator,
            // which differs from the scalar expression (Native lines 529–530).
            result(i, k) = ep * (esat / (pressure(i, k) - (real_t(1) - ep * esat)));
          } else
            result(i, k) = ep * esat / (pressure(i, k) - esat);
        }
      });
  return;
}

template <class Pressure, class Temperature>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Pressure::memory_space>
sat_mixrat_liq_2D(int nz, int ngrdcol, const grid& gr, Pressure p_in_Pa, Temperature T_in_K,
                  int saturation_formula, std::optional<int> start_index_in = std::nullopt) {
  // Preserve the owned-result/memory-space contract through the same 2D kernel.
  Real2D result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sat_mixrat_liq_2D"), ngrdcol, nz);
  sat_mixrat_liq_2D(nz, ngrdcol, gr, p_in_Pa, T_in_K, saturation_formula, result, start_index_in);
  if constexpr (std::is_same_v<typename Pressure::memory_space, typename Real2D::memory_space>)
    return result;
  else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Pressure::memory_space> output(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "sat_mixrat_liq_2D"), ngrdcol, nz);
    from_execution_space(output, result);
    return output;
  }
}
template <class Pressure, class Temperature>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Pressure::memory_space>
sat_mixrat_liq_api(int nz, int ngrdcol, const grid& gr, Pressure p_in_Pa, Temperature T_in_K,
                   int saturation_formula, std::optional<int> start_index_in = std::nullopt) {
  // Return the column-first liquid saturation field [kg/kg] through the Native
  // array generic, retaining optional zero-based start_index_in and grid direction.
  return sat_mixrat_liq_2D(nz, ngrdcol, gr, p_in_Pa, T_in_K, saturation_formula, start_index_in);
}

template <class Pressure, class Temperature>
inline void sat_mixrat_liq_api(int nz, int ngrdcol, const grid& gr, Pressure p_in_Pa,
                               Temperature T_in_K, int saturation_formula, Real2D result,
                               std::optional<int> start_index_in = std::nullopt) {
  // Fill selected execution-space cells; retain zero-based optional selection.
  sat_mixrat_liq_2D(nz, ngrdcol, gr, p_in_Pa, T_in_K, saturation_formula, result, start_index_in);
}

KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq(real_t T_in_K, int saturation_formula,
                                                real_t& esat) {

  // Description:
  //   Computes SVP for water vapor. Calls one of the other functions
  //   that calculate an approximation to SVP.

  // References:
  //   None
  //-----------------------------------------------------------------
  // Input temperature [K] and selector; write vapor pressure [Pa].
  // Saturation Vapor Pressure over Water [Pa]
  // ------------------------ Being Code ------------------------
  // Saturation Vapor Pressure, esat, can be found to be approximated
  // in many different ways.
  switch (saturation_formula) {
  case saturation_bolton:
    // Using the Bolton 1980 approximations for SVP over vapor
    saturation_detail::sat_vapor_press_liq_bolton(T_in_K, esat);
    //Earthworks case
    break;
  case saturation_flatau:
    // Using the Flatau, et al. polynomial approximation for SVP over vapor
    saturation_detail::sat_vapor_press_liq_flatau(T_in_K, esat);
    // ---> h1g
    break;
  case saturation_gfdl:
    // Using GFDL polynomial approximation for SVP with respect to liquid
    saturation_detail::sat_vapor_press_liq_gfdl(T_in_K, esat);
    // <--- h1g
    break;
  case saturation_lookup:
    // Use the lookup table to determine the saturation vapor pressure.
    esat = saturation_detail::sat_vapor_press_liq_lookup(T_in_K);
    break;
  default:
    // Undefined approximation
    esat = -real_t(99999.999);
    break;
  }
}

namespace saturation_detail {
KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq_flatau(real_t T_in_K, real_t& esat) {

  // Description:
  //   Computes SVP for water vapor.

  // References:
  //   ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
  //     and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
  //     pp. 1507--1513
  //------------------------------------------------------------------------
  // Input temperature [K]; write vapor pressure [Pa].
  // Relative error norm expansion (-50 to 50 deg_C) from
  // Table 3 of pp. 1510 of Flatau et al. 1992 (Water Vapor)
  // (The 100 coefficient converts from mb to Pa)
  //   real, dimension(7), parameter :: a = &
  //   100.* (/ 6.11176750,      0.443986062,     0.143053301E-01, &
  //            0.265027242E-03, 0.302246994E-05, 0.203886313E-07, &
  //            0.638780966E-10 /)
  // Relative error norm expansion (-85 to 70 deg_C) from
  // Table 4 of pp. 1511 of Flatau et al.
  //real( kind = core_rknd ), dimension(9), parameter :: a = &
  //100._core_rknd * &
  //  Commented out because the form has been redone, causing these number to no longer be needed,
  //  leaving them in for now for reference.
  //         (/ 6.11583699_core_rknd,      0.444606896_core_rknd,     0.143177157E-01_core_rknd, &
  //         0.264224321E-03_core_rknd, 0.299291081E-05_core_rknd, 0.203154182E-07_core_rknd, &
  //         0.702620698E-10_core_rknd, 0.379534310E-13_core_rknd,-0.321582393E-15_core_rknd /)
  constexpr real_t min_T_in_C = real_t(-85.0);
  real_t T_in_C, T_in_C_sqd;
  // ---------------------- Begin Code ----------------------
  // Determine deg K - 273.15
  T_in_C = T_in_K - T_freeze_K;
  // Since this approximation is only good out to -85 degrees Celsius we
  // truncate the result here (Flatau, et al. 1992)
  T_in_C = ((T_in_C) > (min_T_in_C) ? (T_in_C) : (min_T_in_C));
  // Polynomial approx. (Flatau, et al. 1992)
  // This is the generalized formula but is not computationally efficient.
  // Based on Wexler's expressions(2.1)-(2.4) (See Flatau et al. p 1508)
  // e_{sat} = a_1 + a_2 ( T - T_0 ) + ... + a_{n+1} ( T - T_0 )^n
  // esat = a(1)
  // do i = 2, size( a ) , 1
  //   esat = esat + a(i) * ( T_in_C )**(i-1)
  // end do
  // The 8th order polynomial fit.  When running deep
  // convective cases I noticed that absolute temperature often dips below
  // -50 deg_C at higher altitudes, where the 6th order approximation is
  // not accurate.  -dschanen 20 Nov 2008
  //esat = a(1) + T_in_C*( a(2) + T_in_C*( a(3) + T_in_C*( a(4) + T_in_C &
  //*( a(5) + T_in_C*( a(6) + T_in_C*( a(7) + T_in_C*( a(8) + T_in_C*( a(9) ) ) ) ) ) ) ) )
  // Factoring the polynomial above and changing it into this form allows the cpu
  // to complete the calculations out of order. This is because modern cpus can complete
  // multiple instructions at once if they do not depend on eachother, in the above case
  // each instruction relies on the result of the last. In this version however, the terms
  // in the parentheses could potentially be calculated in parallel by different execution
  // units in the cpu, then only when those terms are being multiplied together do the
  // instructions need to be done one at a time. See clubb issue 834 for more info.
  //   - Gunther Huebler, Aug 2018
  T_in_C_sqd = (T_in_C * T_in_C);
  esat       = -real_t(3.21582393e-14) * (T_in_C - real_t(646.5835252598777)) *
               (T_in_C + real_t(90.72381630364440)) *
               (T_in_C_sqd + real_t(111.0976961559954) * T_in_C + real_t(6459.629194243118)) *
               (T_in_C_sqd + real_t(152.3131930092453) * T_in_C + real_t(6499.774954705265)) *
               (T_in_C_sqd + real_t(174.4279584934021) * T_in_C + real_t(7721.679732114084));
}

KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq_bolton(real_t T_in_K, real_t& esat) {
  // Description:
  //   Computes SVP for water vapor.
  // References:
  //   Bolton 1980
  //------------------------------------------------------------------------
  // Input temperature [K]; write vapor pressure [Pa].
  // --------------------------- Begin Code ---------------------------
  // (Bolton 1980) approx.
  // Generally this more computationally expensive than the Flatau polnomial expansion
  esat = real_t(611.2) *
         Kokkos::exp((real_t(17.67) * (T_in_K - T_freeze_K)) / (T_in_K - real_t(29.65)));
}

KOKKOS_INLINE_FUNCTION void sat_vapor_press_liq_gfdl(real_t T_in_K, real_t& esat) {
  // Description:
  // copy from "GFDL polysvp.F90"
  //  Compute saturation vapor pressure with respect to liquid  by using
  // function from Goff and Gratch (1946)

  //  Polysvp returned in units of pa.
  //  T_in_K  is input in units of K.
  //------------------------------------------------------------------------
  // Input temperature [K]; write vapor pressure [Pa].
  constexpr real_t min_T_in_K = real_t(203.15);
  real_t T_in_K_clipped;
  // --------------------------- Begin Code ---------------------------
  // Since the Goff-Gratch approximation is valid only down to -70 degrees Celsius,
  //   we threshold the temperature.  This will yield a minimal saturation at
  //   cold temperatures.
  T_in_K_clipped = ((min_T_in_K) > (T_in_K) ? (min_T_in_K) : (T_in_K));
  // Goff Gratch equation, uncertain below -70 C
  esat =
      Kokkos::pow(real_t(10.),
                  (-real_t(7.90298) * (real_t(373.16) / T_in_K_clipped - real_t(1.)) +
                   real_t(5.02808) * Kokkos::log10(real_t(373.16) / T_in_K_clipped) -
                   real_t(1.3816e-7) *
                       (Kokkos::pow(real_t(10.), (real_t(11.344) *
                                                  (real_t(1.) - T_in_K_clipped / real_t(373.16)))) -
                        real_t(1.)) +
                   real_t(8.1328e-3) *
                       (Kokkos::pow(real_t(10.), (-real_t(3.49149) *
                                                  (real_t(373.16) / T_in_K_clipped - real_t(1.)))) -
                        real_t(1.)) +
                   Kokkos::log10(real_t(1013.246)))) *
      real_t(100.);
}

KOKKOS_INLINE_FUNCTION real_t sat_vapor_press_liq_lookup(real_t T_in_K) {

  // Description:
  //   Computes SVP for water vapor, using a lookup table.
  //
  //   The lookup table was constructed using the Flatau approximation.

  // References:
  //   ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
  //     and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
  //     pp. 1507--1513
  //------------------------------------------------------------------------
  // Input temperature [K]; return tabulated vapor pressure [Pa].
  int T_in_K_int;
  real_t esat;
  // ---- Begin Code ----
  T_in_K_int = static_cast<int>(Kokkos::round(T_in_K));
  // Since this approximation is only good out to -85 degrees Celsius we
  // truncate the result here
  T_in_K_int = ((((T_in_K_int) > (188) ? (T_in_K_int) : (188))) < (343)
                    ? (((T_in_K_int) > (188) ? (T_in_K_int) : (188)))
                    : (343));
  // Use the lookup table to determine the saturation vapor pressure.
  esat = svp_liq_lookup_table(T_in_K_int);
  return esat;
}
} // namespace saturation_detail

KOKKOS_INLINE_FUNCTION real_t sat_mixrat_ice_k(real_t p_in_Pa, real_t T_in_K,
                                               int saturation_formula) {

  // Description:
  //   Used to compute the saturation mixing ratio of ice.

  // References:
  //   Formula from Emanuel 1994, 4.4.15
  //-------------------------------------------------------------------------
  // Input pressure [Pa], temperature [K] and selector; return mixing ratio [kg/kg].
  // Device adaptation: inline the Native 1x1 array call's cell body instead
  // of allocating arrays in each parcel. Array overload uses this same body.
  // Flatau ice coefficients use the relative-error expansion (-90 to 0 deg_C),
  // Table 4, p. 1511 of Flatau et al. 1992; factor 100 converts mb to Pa.
  constexpr real_t min_T_in_C = real_t(-90.0), min_T_in_K = real_t(173.15);
  real_t T_in_C, T_in_K_clipped, esat_ice, result;
  // Determine the SVP for the given temperature
  switch (saturation_formula) {
  case saturation_bolton:
    // Using the Bolton 1980 approximations for SVP over ice
    // Exponential approx.
    esat_ice = real_t(100.0) * Kokkos::exp(real_t(23.33086) - (real_t(6111.72784) / T_in_K) +
                                           (real_t(0.15215) * Kokkos::log(T_in_K)));
    break;
  case saturation_flatau:
    // Using the Flatau, et al. polynomial approximation for SVP over ice
    // Determine deg K - 273.15
    T_in_C = T_in_K - T_freeze_K;
    // Since this approximation is only good out to -90 degrees Celsius we
    // truncate the result here (Flatau, et al. 1992)
    T_in_C = ((T_in_C) > (min_T_in_C) ? (T_in_C) : (min_T_in_C));
    // Polynomial approx. (Flatau, et al. 1992)
    //   esati = a(1)
    //   do i = 2, size( a ), 1
    //     esati = esati + a(i) * ( T_in_C )**(i-1)
    //   end do
    esat_ice =
        (real_t(100) * real_t(6.09868993)) +
        T_in_C *
            ((real_t(100) * real_t(0.499320233)) +
             T_in_C *
                 ((real_t(100) * real_t(0.184672631E-01)) +
                  T_in_C *
                      ((real_t(100) * real_t(0.402737184E-03)) +
                       T_in_C *
                           ((real_t(100) * real_t(0.565392987E-05)) +
                            T_in_C * ((real_t(100) * real_t(0.521693933E-07)) +
                                      T_in_C * ((real_t(100) * real_t(0.307839583E-09)) +
                                                T_in_C * ((real_t(100) * real_t(0.105785160E-11)) +
                                                          T_in_C * (real_t(100) *
                                                                    real_t(0.161444444E-14)))))))));
    // ---> h1g, 2010-06-16
    break;
  case saturation_gfdl:
    // Using GFDL polynomial approximation for SVP with respect to ice
    // Since the Goff-Gratch ice approximation is valid only down to -100 degrees Celsius,
    //   we threshold the temperature.  This will yield a minimal saturation at
    //   cold temperatures.
    T_in_K_clipped = ((min_T_in_K) > (T_in_K) ? (min_T_in_K) : (T_in_K));
    // Goff Gratch equation (good down to -100 C)
    esat_ice = Kokkos::pow(real_t(10.),
                           (-real_t(9.09718) * (real_t(273.16) / T_in_K_clipped - real_t(1.)) -
                            real_t(3.56654) * Kokkos::log10(real_t(273.16) / T_in_K_clipped) +
                            real_t(0.876793) * (real_t(1.) - T_in_K_clipped / real_t(273.16)) +
                            Kokkos::log10(real_t(6.1071)))) *
               real_t(100.);
    // <--- h1g, 2010-06-16
    break;
  default:
    // Undefined approximation
    esat_ice = -real_t(99999.999);
    break;
  }
  if (p_in_Pa - esat_ice < real_t(1))
    result = ep;
  else {
    if (I_sat_sphum)
      result = ep * (esat_ice / (p_in_Pa - (real_t(1) - ep) * esat_ice));
    else
      result = ep * esat_ice / (p_in_Pa - esat_ice);
  }
  return result;
}
KOKKOS_INLINE_FUNCTION real_t sat_mixrat_ice(real_t p_in_Pa, real_t T_in_K,
                                             int saturation_formula) {
  // Return ice saturation mixing ratio [kg/kg] from pressure [Pa],
  // temperature [K] and selector through the Native scalar generic.
  return sat_mixrat_ice_k(p_in_Pa, T_in_K, saturation_formula);
}

template <class Pressure, class Temperature>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Pressure::memory_space>
sat_mixrat_ice_2D(int nz, int ngrdcol, Pressure p_in_Pa, Temperature T_in_K,
                  int saturation_formula) {

  // Description:
  //   Used to compute the saturation mixing ratio of ice.

  // References:
  //   Formula from Emanuel 1994, 4.4.15
  //-------------------------------------------------------------------------
  // Input column-first pressure [Pa] and temperature [K] fields;
  // return mixing ratio [kg/kg] in the pressure field's memory space.
  // ------------------------ Begin Code ------------------------
  auto pressure    = to_execution_space(p_in_Pa);
  auto temperature = to_execution_space(T_in_K);
  Real2D result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sat_mixrat_ice_2D"), ngrdcol, nz);
  Kokkos::parallel_for(
      "sat_mixrat_ice_2D", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i  = cell % ngrdcol;
        const int k  = cell / ngrdcol;
        result(i, k) = sat_mixrat_ice_k(pressure(i, k), temperature(i, k), saturation_formula);
      });
  if constexpr (std::is_same_v<typename Pressure::memory_space, typename Real2D::memory_space>)
    return result;
  else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Pressure::memory_space> output(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "sat_mixrat_ice_2D"), ngrdcol, nz);
    from_execution_space(output, result);
    return output;
  }
}
template <class Pressure, class Temperature>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Pressure::memory_space>
sat_mixrat_ice(int nz, int ngrdcol, Pressure p_in_Pa, Temperature T_in_K, int saturation_formula) {
  // Return the column-first ice saturation field [kg/kg] through the Native
  // array generic for pressure [Pa], temperature [K] and the formula selector.
  return sat_mixrat_ice_2D(nz, ngrdcol, p_in_Pa, T_in_K, saturation_formula);
}

KOKKOS_INLINE_FUNCTION real_t rcm_sat_adj(real_t thlm, real_t rtm, real_t p_in_Pa, real_t exner,
                                          int saturation_formula) {

  // Description:
  //
  //   This function uses an iterative method to find the value of rcm
  //   from an initial profile that has saturation at some point.
  //
  // References:
  //   None
  //-------------------------------------------------------------------------
  // Input liquid-water potential temperature thlm [K], total-water mixing
  // ratio rtm [kg/kg], pressure [Pa], Exner function [-] and formula selector.
  // Return cloud-water mixing ratio rcm [kg/kg].
  constexpr real_t tolerance = real_t(0.001); // Theta calculation tolerance [K].
  constexpr int itermax      = 1000000;
  real_t rcm, theta = thlm, answer, too_low = real_t(0), too_high = real_t(0);
  int iteration;
  // ----- Begin Code -----
  // Default theta/bracket initialization is retained above.
  for (iteration = 1; iteration <= itermax; ++iteration) {
    const real_t excess = rtm - sat_mixrat_liq_api(p_in_Pa, theta * exner, saturation_formula);
    answer = theta - (Lv / (Cp * exner)) * (excess > zero_threshold ? excess : zero_threshold);
    if (Kokkos::abs(answer - thlm) <= tolerance)
      break;
    else if (answer - thlm > tolerance)
      too_high = theta;
    else if (thlm - answer > tolerance)
      too_low = theta;
    // First iteration must set an upper bound that is way too high.
    if (iteration == 1)
      too_high = theta + real_t(20);
    theta = (too_low + too_high) / real_t(2);
  }
  // Preserve the exact Native loop-counter check, including its last-iteration
  // exit versus natural-exhaustion distinction; this is not JAX's fixed scan.
  if (iteration == itermax) {
    rcm = real_t(0);
    Kokkos::abort("Error in rcm_sat_adj: could not determine rcm");
  } else {
    const real_t excess = rtm - sat_mixrat_liq_api(p_in_Pa, theta * exner, saturation_formula);
    rcm                 = excess > zero_threshold ? excess : zero_threshold;
  }
  return rcm;
}

} // namespace clubb_kokkos
