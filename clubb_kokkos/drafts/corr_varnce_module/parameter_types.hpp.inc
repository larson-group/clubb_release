// Merge inside clubb_kokkos, before the existing generated hm_metadata_type.
// These are the three Native scalar parameter types, not metadata replacements.

struct hmp2_ip_on_hmm2_ip_ratios_type {
  // Prescribed parameters for hydrometeor values of <hm|_ip'^2> / <hm|_ip>^2,
  // where  <hm|_ip'^2>  is the in-precip. variance of the hydrometeor and
  //        <hm|_ip>     is the in-precip. mean of the hydrometeor
  //
  // These values are dependent on the horizontal grid spacing of the run, and are calculated
  // using a slope and intercept corresponding to each hydrometer.
  real_t rr = real_t(1.0); // For rain water mixing ratio [-]
  real_t Nr = real_t(1.0); // For rain drop concentration [-]
  real_t ri = real_t(1.0); // For ice mixing ratio        [-]
  real_t Ni = real_t(1.0); // For ice concentration       [-]
  real_t rs = real_t(1.0); // For snow mixing ratio       [-]
  real_t Ns = real_t(1.0); // For snow concentration      [-]
  real_t rg = real_t(1.0); // For graupel mixing ratio    [-]
  real_t Ng = real_t(1.0); // For graupel concentration   [-]
};

// These slopes and intercepts below are used to calculate the hmp2_ip_on_hmm2_ip_ratios_type
// values that are defined above. This functionality is described by equations 8, 10, & 11
// in ``Parameterization of the Spatial Variability of Rain for Large-Scale Models and
// Remote Sensing`` Lebo, et al, October 2015
// https://journals.ametsoc.org/doi/pdf/10.1175/JAMC-D-15-0066.1
// see clubb:ticket:830 for more detail
//
// hmp2_ip_on_hmm2_ip(iirr) = hmp2_ip_on_hmm2_ip_intrcpt%rr + &
//                            hmp2_ip_on_hmm2_ip_slope%rr * max( host_dx, host_dy )
//
// In Lebo et al. the suggested values were
//     slope = 2.12e-5 [1/m]
//     intercept = 0.54 [-]
//
// In CLUBB standalone, these parameters can be set based on the value for a
// given case in the CASE_model.in file.

struct hmp2_ip_on_hmm2_ip_slope_type {
  real_t rr = real_t(2.12e-5); // For rain water mixing ratio [1/m]
  real_t Nr = real_t(2.12e-5); // For rain drop concentration [1/m]
  real_t ri = real_t(2.12e-5); // For ice mixing ratio        [1/m]
  real_t Ni = real_t(2.12e-5); // For ice concentration       [1/m]
  real_t rs = real_t(2.12e-5); // For snow mixing ratio       [1/m]
  real_t Ns = real_t(2.12e-5); // For snow concentration      [1/m]
  real_t rg = real_t(2.12e-5); // For graupel mixing ratio    [1/m]
  real_t Ng = real_t(2.12e-5); // For graupel concentration   [1/m]
};

struct hmp2_ip_on_hmm2_ip_intrcpt_type {
  real_t rr = real_t(0.54); // For rain water mixing ratio [-]
  real_t Nr = real_t(0.54); // For rain drop concentration [-]
  real_t ri = real_t(0.54); // For ice mixing ratio        [-]
  real_t Ni = real_t(0.54); // For ice concentration       [-]
  real_t rs = real_t(0.54); // For snow mixing ratio       [-]
  real_t Ns = real_t(0.54); // For snow concentration      [-]
  real_t rg = real_t(0.54); // For graupel mixing ratio    [-]
  real_t Ng = real_t(0.54); // For graupel concentration   [-]
};
