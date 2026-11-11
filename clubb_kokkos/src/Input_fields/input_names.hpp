// Column-header strings used by the Native-compatible case input readers.
#pragma once
namespace clubb_kokkos {
inline constexpr char z_name[] = "z[m]";
inline constexpr char pressure_name[] = "Press[Pa]";
inline constexpr char press_mb_name[] = "Press[mb]";
inline constexpr char temperature_name[] = "T[K]";
inline constexpr char theta_name[] = "thm[K]";
inline constexpr char thetal_name[] = "thlm[K]";
inline constexpr char temperature_f_name[] = "T_f[K/s]";
inline constexpr char thetal_f_name[] = "thlm_f[K/s]";
inline constexpr char theta_f_name[] = "thm_f[K/s]";
inline constexpr char rt_name[] = "rt[kg/kg]";
inline constexpr char sp_humidity_name[] = "sp_hmdty[kg/kg]";
inline constexpr char rt_f_name[] = "rtm_f[kg/kg/s]";
inline constexpr char sp_humidity_f_name[] = "sp_hmdty_f[kg/kg/s]";
inline constexpr char um_name[] = "u[m/s]";
inline constexpr char vm_name[] = "v[m/s]";
inline constexpr char ug_name[] = "ug[m/s]";
inline constexpr char vg_name[] = "vg[m/s]";
inline constexpr char um_ref_name[] = "um_ref[m/s]";
inline constexpr char vm_ref_name[] = "vm_ref[m/s]";
inline constexpr char um_f_name[] = "um_f[m/s^2]";
inline constexpr char vm_f_name[] = "vm_f[m/s^2]";
inline constexpr char wm_name[] = "w[m/s]";
inline constexpr char omega_name[] = "omega[Pa/s]";
inline constexpr char omega_mb_hr_name[] = "omega[mb/hr]";
inline constexpr char CO2_name[] = "CO2[ppmv]";
inline constexpr char CO2_umol_name[] = "CO2[umol/m^2/s]";
inline constexpr char ozone_name[] = "o3[kg/kg]";
inline constexpr char time_name[] = "Time[s]";
inline constexpr char latent_ht_name[] = "latent_ht[W/m^2]";
inline constexpr char sens_ht_name[] = "sens_ht[W/m^2]";
inline constexpr char upwp_sfc_name[] = "upwp_sfc[(m/s)^2]";
inline constexpr char vpwp_sfc_name[] = "vpwp_sfc[(m/s)^2]";
inline constexpr char T_sfc_name[] = "T_sfc[K]";
inline constexpr char wpthlp_sfc_name[] = "wpthlp_sfc[mK/s]";
inline constexpr char wpqtp_sfc_name[] = "wpqtp_sfc[(kg/kg)m/s]";
} // namespace clubb_kokkos
