// Read SCAM IOP observations and calculate their sounding/surface/forcing fields.
//
// The four core readers retain Native real4 working precision and array order.
// Dimension/type selection, surface extrapolation, thermodynamics, forcing
// precedence and data reshaping are C++; raw reads use the NetCDF-C service.
// The three NetCDF helpers preserve
// required/optional diagnostics. calcdate advances an integer YYMMDD date.
// State belongs to one readiopdata_state. Storage indices are zero-based;
// Native NetCDF record ordinals are supplied only at the file-service boundary.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include <string>
namespace clubb_kokkos {
struct readiopdata_state;
using IopReal1D = Kokkos::View<float*, Kokkos::LayoutLeft>;
using IopReal2D = Kokkos::View<float**, Kokkos::LayoutLeft>;
using IopInt1D  = Kokkos::View<int*, Kokkos::LayoutLeft>;
void readiopdata_core_dims(ScratchWorkspace& workspace, bool masterproc,
                           const std::string& iopfilepath,
                           bool scamiop_from_global_cam,                                // In
                           int& ntime, int& nlev, int& nlat, int& nlon, IopInt1D& tsec, // Out
                           int& bdate, IopReal1D& dplevs, IopReal1D& lat_in,
                           IopReal1D& lon_in,       // Out
                           readiopdata_state& iop); // Native module InOut
void readiopdata_core_snd(
    ScratchWorkspace& workspace,                                                   // InOut
    bool masterproc, const std::string& iopfilepath, bool scamiop_from_global_cam, // In
    int ntime, int nlev, IopReal1D dplevs, float rgas, float cp, float fac_cond,
    float fac_sub,                                                                       // In
    IopReal2D& zsnd, IopReal2D& usnd, IopReal2D& vsnd, IopReal2D& tsnd,                  // Out
    IopReal2D& qsnd, IopReal2D& psnd, IopReal2D& ugls, IopReal2D& vgls, IopReal2D& wgls, // Out
    bool* get_add_surface_data, bool* have_omega, bool* wgls_holds_omega, // Optional Out
    bool* have_w, bool* wgls_holds_w, std::string& subs_type, // UWM_MISC Optional/required Out
    readiopdata_state& iop);                                  // Native module InOut
void readiopdata_core_sfc(ScratchWorkspace& workspace, bool masterproc,
                          const std::string& iopfilepath,
                          bool scamiop_from_global_cam,                           // In
                          bool SFC_FLX_FXD, bool get_add_surface_data, int ntime, // In
                          IopReal1D& sstsfc, IopReal1D& shsfc, IopReal1D& lhsfc,
                          IopReal1D& tausfc,       // Out
                          readiopdata_state& iop); // Native module InOut
void readiopdata_core_frc(ScratchWorkspace& workspace, bool masterproc,
                          const std::string& iopfilepath, bool get_add_surface_data,
                          bool have_omega,                              // In
                          IopReal1D dplevs, int nlev, int ntime,        // In
                          IopReal2D& dtls, IopReal2D& dqls,             // Out
                          readiopdata_state& iop);                      // Native module InOut
void get_netcdf_dimlength(int NCID, const std::string& dimName,         // In
                          int& dimlength, int& status,                  // Out
                          bool required, const readiopdata_state& iop); // In, Native module In
void get_netcdf_var1d_real(ScratchWorkspace& workspace, int NCID, const std::string& varName,
                           IopReal1D var,
                           const bool& use_NF90_REAL, // In, In, InOut, In
                           int& status, bool required,
                           const readiopdata_state& iop); // Out, In, Native module In
void get_netcdf_var2d_real(ScratchWorkspace& workspace,   // InOut
                           int NCID, const std::string& varName, int ntime, int nlev, // In
                           IopReal2D var, const bool& use_NF90_REAL,                  // InOut, In
                           int& status, bool required,
                           const readiopdata_state& iop);    // Out, In, Native module In
void calcdate(int inDate, int inSecs,                        // In
              int& outDate, int& outSecs, float& outCalday); // Out
[[noreturn]] void task_abort(); // Native external fatal service, permanent C++.
} // namespace clubb_kokkos
