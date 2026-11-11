// Merge inside clubb_kokkos after the generated hm_metadata_type; include <array>.
// Repository standalone/JAX code only reads these default tables. They are
// immutable C++ constants; external mutation of the public Native arrays is
// not a supported C++ API. Fortran-flat order retains its RESHAPE semantics.
inline constexpr int d_var_total = 12; // Size of the default correlation arrays.

inline constexpr std::array<real_t, d_var_total * d_var_total> corr_array_n_cloud_def = {
    real_t(1.),   real_t(-.6),  real_t(.09),   real_t(.09),
    real_t(.788), real_t(.675), real_t(.240),  real_t(.222),
    real_t(.240), real_t(.222), real_t(.240),  real_t(.222), // chi Native column
    real_t(0.),   real_t(1.),   real_t(.027),  real_t(.027),
    real_t(.114), real_t(.115), real_t(-.029), real_t(.093),
    real_t(.022), real_t(.013), real_t(0.),    real_t(0.), // eta Native column
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.34),
    real_t(.315), real_t(.270), real_t(.120),  real_t(.167),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.), // w Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.),
    real_t(0.),   real_t(0.),   real_t(.464),  real_t(.320),
    real_t(.168), real_t(.232), real_t(0.),    real_t(0.), // Ncn Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(1.),   real_t(.821), real_t(0.),    real_t(0.),
    real_t(.173), real_t(.164), real_t(.319),  real_t(.308), // rr Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(1.),   real_t(.152),  real_t(.143),
    real_t(0.),   real_t(0.),   real_t(.285),  real_t(.273), // Nr Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.758),
    real_t(.585), real_t(.571), real_t(.379),  real_t(.363), // ri Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.),
    real_t(.571), real_t(.550), real_t(.363),  real_t(.345), // Ni Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(1.),   real_t(.758), real_t(.485),  real_t(.470), // rs Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(1.),   real_t(.470),  real_t(.450), // Ns Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.758), // rg Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.), // Ng Native column
};

inline constexpr std::array<real_t, d_var_total * d_var_total> corr_array_n_below_def = {
    real_t(1.),   real_t(.3),   real_t(.09),   real_t(.09),
    real_t(.788), real_t(.675), real_t(.240),  real_t(.222),
    real_t(.240), real_t(.222), real_t(.240),  real_t(.222), // chi Native column
    real_t(0.),   real_t(1.),   real_t(.027),  real_t(.027),
    real_t(.114), real_t(.115), real_t(-.029), real_t(.093),
    real_t(.022), real_t(.013), real_t(0.),    real_t(0.), // eta Native column
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.34),
    real_t(.315), real_t(.270), real_t(.120),  real_t(.167),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.), // w Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.),
    real_t(0.),   real_t(0.),   real_t(.464),  real_t(.320),
    real_t(.168), real_t(.232), real_t(0.),    real_t(0.), // Ncn Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(1.),   real_t(.821), real_t(0.),    real_t(0.),
    real_t(.173), real_t(.164), real_t(.319),  real_t(.308), // rr Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(1.),   real_t(.152),  real_t(.143),
    real_t(0.),   real_t(0.),   real_t(.285),  real_t(.273), // Nr Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.758),
    real_t(.585), real_t(.571), real_t(.379),  real_t(.363), // ri Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.),
    real_t(.571), real_t(.550), real_t(.363),  real_t(.345), // Ni Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(1.),   real_t(.758), real_t(.485),  real_t(.470), // rs Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(1.),   real_t(.470),  real_t(.450), // Ns Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.758), // rg Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.), // Ng Native column
};
