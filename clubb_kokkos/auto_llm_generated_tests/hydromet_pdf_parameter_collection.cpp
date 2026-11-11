#include "CLUBB_core/hydromet_pdf_parameter_module.hpp"
#include "f2k/derived_type_bridge.hpp"

#include <Kokkos_Core.hpp>

#include <stdexcept>

int main() {
  // Check the device field layout, first-index contiguity, and batch reset.
  Kokkos::initialize();
  {
    clubb_kokkos::hydromet_pdf_parameter_collection params;
    params.allocate(3, 2);
    if (params.hm_1.extent(0) != 3 || params.hm_1.extent(1) != 2 ||
        params.hm_1.extent(2) != params.max_hydromet_dim ||
        params.corr_hmx_hmy_1.extent(3) != params.max_hydromet_dim) {
      throw std::runtime_error("hydrometeor PDF collection shape");
    }

    auto hm_1 = params.hm_1;
    auto corr_hmx_hmy_1 = params.corr_hmx_hmy_1;
    auto mu_Ncn_1 = params.mu_Ncn_1;
    Kokkos::parallel_for("fill_hydromet_pdf_params", Kokkos::RangePolicy<>(0, 3 * 2 * 8 * 8),
                         KOKKOS_LAMBDA(const int linear) {
      const int column = linear % 3;
      const int level = (linear / 3) % 2;
      const int species = (linear / 6) % 8;
      const int other_species = linear / 48;
      corr_hmx_hmy_1(column, level, species, other_species) =
          static_cast<clubb_kokkos::real_t>(linear);
      if (other_species == 0) {
        hm_1(column, level, species) = static_cast<clubb_kokkos::real_t>(linear);
      }
      if (species == 0 && other_species == 0) {
        mu_Ncn_1(column, level) = static_cast<clubb_kokkos::real_t>(linear);
      }
    });

    f2k::roundtrip(params);
    auto host_hm = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, params.hm_1);
    auto host_corr =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, params.corr_hmx_hmy_1);
    if (&host_hm(1, 0, 0) - &host_hm(0, 0, 0) != 1 ||
        &host_corr(1, 0, 0, 0) - &host_corr(0, 0, 0, 0) != 1 ||
        host_hm(2, 1, 7) != 47 || host_corr(2, 1, 7, 7) != 383) {
      throw std::runtime_error("hydrometeor PDF column-first indexing");
    }

    params.zero();
    host_hm = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, params.hm_1);
    host_corr = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, params.corr_hmx_hmy_1);
    auto host_mu = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, params.mu_Ncn_1);
    if (host_hm(2, 1, 7) != 0 || host_corr(2, 1, 7, 7) != 0 || host_mu(2, 1) != 0) {
      throw std::runtime_error("hydrometeor PDF batch reset");
    }
  }
  Kokkos::finalize();
}
