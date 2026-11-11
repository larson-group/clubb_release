// Cloudy updraft and downdraft diagnostics for aerosol activation.
//
// calc_w_up_in_cloud integrates each PDF component's positive and negative
// vertical velocities, then weights them by its cloud and mixture fractions.
// Fields use (column, level) LayoutLeft Views; all scalar scratch is cell local.
// The Native caller owns PDF selection and statistics gates. See INTEGRATION.md
// before merging this source-only fragment into the permanent module owner.
// Reference: https://www.overleaf.com/project/614a136d47846639af22ae34
#include "pdf_closure_module.hpp"
#include "constants_clubb.hpp"
#include <Kokkos_MathematicalFunctions.hpp>

namespace clubb_kokkos {

//=============================================================================
void calc_w_up_in_cloud(int nz, int ngrdcol,                                          // In
                        Real2D mixt_frac, Real2D cloud_frac_1, Real2D cloud_frac_2,   // In
                        Real2D w_1, Real2D w_2, Real2D varnce_w_1, Real2D varnce_w_2, // In
                        Real2D w_up_in_cloud, Real2D w_down_in_cloud,                 // Out
                        Real2D cloudy_updraft_frac, Real2D cloudy_downdraft_frac) {   // Out
  // Description:
  // Subroutine that computes the mean cloudy updraft (and also calculates
  // the mean cloudy downdraft).
  //
  // In order to activate aerosol, we'd like to feed the activation scheme
  // a vertical velocity that's representative of cloudy updrafts. For skewed
  // layers, like cumulus layers, this might be an improvement over the square
  // root of wp2 that's currently used. At the same time, it would be simpler
  // and less expensive than feeding SILHS samples into the aerosol code
  // (see larson-group/e3sm#19 and larson-group/e3sm#26).
  //
  // The formulas are only valid for certain PDFs in CLUBB (ADG1, ADG2,
  // new hybrid), hence we omit calculation if another PDF type is used.
  //
  // References: https://www.overleaf.com/project/614a136d47846639af22ae34
  //----------------------------------------------------------------------

  // nz: Number of vertical levels
  // ngrdcol: Number of grid columns
  //----------- Input Variables -----------
  // mixt_frac: mixture fraction [-]
  // cloud_frac_1: cloud fraction (1st PDF component) [-]
  // cloud_frac_2: cloud fraction (2nd PDF component) [-]
  // w_1: upward velocity (1st PDF component) [m/s]
  // w_2: upward velocity (2nd PDF component) [m/s]
  // varnce_w_1: variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: variance of w (2nd PDF component) [m^2/s^2]

  //----------- Output Variables -----------
  // w_up_in_cloud: mean cloudy updraft speed [m/s]
  // w_down_in_cloud: mean cloudy downdraft speed [m/s]
  // cloudy_updraft_frac: cloudy updraft fraction [-]
  // cloudy_downdraft_frac: cloudy downdraft fraction [-]

  Kokkos::parallel_for(
      "calc_w_up_in_cloud", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;

        //----------- Local Variables -----------
        real_t w_up_1, w_up_2;        // integral of w and Heaviside fnc, where w > 0
        real_t w_down_1, w_down_2;    // integral of w and Heaviside fnc, where w < 0
        real_t stdev_w_1, stdev_w_2;  // Standard deviation of w
        real_t ratio_w_1;             // mu_w_1 / (sqrt(2) * sigma_w_1)
        real_t ratio_w_2;             // mu_w_2 / (sqrt(2) * sigma_w_2)
        real_t erf_ratio_w_1;         // erf(ratio_w_1)
        real_t erf_ratio_w_2;         // erf(ratio_w_2)
        real_t exp_neg_ratio_w_1_sqd; // exp(-ratio_w_1^2)
        real_t exp_neg_ratio_w_2_sqd; // exp(-ratio_w_2^2)
        real_t updraft_frac_1;        // Fraction of 1st PDF comp. where w > 0
        real_t updraft_frac_2;        // Fraction of 2nd PDF comp. where w > 0
        real_t downdraft_frac_1;      // Fraction of 1st PDF comp. where w < 0
        real_t downdraft_frac_2;      // Fraction of 2nd PDF comp. where w < 0

        stdev_w_1 = Kokkos::sqrt(varnce_w_1(i, k));
        stdev_w_2 = Kokkos::sqrt(varnce_w_2(i, k));

        // Calculate quantities in the 1st PDF component.
        if (w_1(i, k) > max_num_stdevs * stdev_w_1) {

          // The mean of w in the 1st PDF component is more than
          // max_num_stdevs standard deviations above 0.
          // The entire 1st PDF component is found in an updraft (w > 0).
          w_up_1           = w_1(i, k);
          updraft_frac_1   = one;
          w_down_1         = zero;
          downdraft_frac_1 = zero;

        } else if (w_1(i, k) < -max_num_stdevs * stdev_w_1) {

          // The mean of w in the 1st PDF component is more than
          // max_num_stdevs standard deviations below 0.
          // The entire 1st PDF component is found in a downdraft (w < 0).
          w_up_1           = zero;
          updraft_frac_1   = zero;
          w_down_1         = w_1(i, k);
          downdraft_frac_1 = one;

        } else {

          // The 1st PDF component contains both updraft and downdraft.
          ratio_w_1             = w_1(i, k) / (sqrt_2 * (eps > stdev_w_1 ? eps : stdev_w_1));
          erf_ratio_w_1         = Kokkos::erf(ratio_w_1);
          exp_neg_ratio_w_1_sqd = Kokkos::exp(-(ratio_w_1 * ratio_w_1));

          w_up_1 = one_half * w_1(i, k) * (one + erf_ratio_w_1) +
                   (stdev_w_1 / sqrt_2pi) * exp_neg_ratio_w_1_sqd;

          updraft_frac_1 = one_half * (one + erf_ratio_w_1);

          w_down_1 = one_half * w_1(i, k) * (one - erf_ratio_w_1) -
                     (stdev_w_1 / sqrt_2pi) * exp_neg_ratio_w_1_sqd;

          // downdraft_frac_1 = one_half * (one - erf_ratio_w_1)
          downdraft_frac_1 = one - updraft_frac_1;
        }

        // Calculate quantities in the 2nd PDF component.
        if (w_2(i, k) > max_num_stdevs * stdev_w_2) {

          // The mean of w in the 2nd PDF component is more than
          // max_num_stdevs standard deviations above 0.
          // The entire 2nd PDF component is found in an updraft (w > 0).
          w_up_2           = w_2(i, k);
          updraft_frac_2   = one;
          w_down_2         = zero;
          downdraft_frac_2 = zero;

        } else if (w_2(i, k) < -max_num_stdevs * stdev_w_2) {

          // The mean of w in the 2nd PDF component is more than
          // max_num_stdevs standard deviations below 0.
          // The entire 2nd PDF component is found in a downdraft (w < 0).
          w_up_2           = zero;
          updraft_frac_2   = zero;
          w_down_2         = w_2(i, k);
          downdraft_frac_2 = one;

        } else {

          // The 2nd PDF component contains both updraft and downdraft.
          ratio_w_2             = w_2(i, k) / (sqrt_2 * (eps > stdev_w_2 ? eps : stdev_w_2));
          erf_ratio_w_2         = Kokkos::erf(ratio_w_2);
          exp_neg_ratio_w_2_sqd = Kokkos::exp(-(ratio_w_2 * ratio_w_2));

          w_up_2 = one_half * w_2(i, k) * (one + erf_ratio_w_2) +
                   (stdev_w_2 / sqrt_2pi) * exp_neg_ratio_w_2_sqd;

          updraft_frac_2 = one_half * (one + erf_ratio_w_2);

          w_down_2 = one_half * w_2(i, k) * (one - erf_ratio_w_2) -
                     (stdev_w_2 / sqrt_2pi) * exp_neg_ratio_w_2_sqd;

          // downdraft_frac_2 = one_half * (one - erf_ratio_w_2)
          downdraft_frac_2 = one - updraft_frac_2;
        }

        // Calculate the total cloudy updraft fraction.
        cloudy_updraft_frac(i, k) = mixt_frac(i, k) * cloud_frac_1(i, k) * updraft_frac_1 +
                                    (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * updraft_frac_2;

        // Calculate the total cloudy downdraft fraction.
        cloudy_downdraft_frac(i, k) =
            mixt_frac(i, k) * cloud_frac_1(i, k) * downdraft_frac_1 +
            (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * downdraft_frac_2;

        // Calculate the mean vertical velocity found in a cloudy updraft.
        w_up_in_cloud(i, k) = (mixt_frac(i, k) * cloud_frac_1(i, k) * w_up_1 +
                               (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * w_up_2) /
                              (eps > cloudy_updraft_frac(i, k) ? eps : cloudy_updraft_frac(i, k));

        // Calculate the mean vertical velocity found in a cloudy downdraft.
        w_down_in_cloud(i, k) =
            (mixt_frac(i, k) * cloud_frac_1(i, k) * w_down_1 +
             (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * w_down_2) /
            (eps > cloudy_downdraft_frac(i, k) ? eps : cloudy_downdraft_frac(i, k));
      });

  return;
}

} // namespace clubb_kokkos
