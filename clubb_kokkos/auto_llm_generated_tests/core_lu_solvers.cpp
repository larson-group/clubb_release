// Compare all five LU solver routines with compiled native Fortran.
// Column-varying diagonally dominant systems exercise minimal and longer row
// counts, one/many columns, one/many/zero RHSs, and the native band-index mapping.
// Solutions must match native bits; residuals independently check A*x=rhs.
#include "scratch_workspace_owner.hpp"
#include "CLUBB_core/tridiag_lu_solver.hpp"
#include "CLUBB_core/penta_lu_solver.hpp"

#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>

using clubb_kokkos::real_t;
extern "C" void f2k_lu_reference_c(int, int, int, int, const real_t*, real_t*, real_t*, real_t*,
                                   real_t*);

int main(int argc, char* argv[]) {
  // Solve each shape through real C++ device and native overloads, preserving
  // matrix/RHS inputs and checking exact output and independent row residuals.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    clubb_kokkos::test::ScratchWorkspaceOwner scratch_owner;
    auto& workspace = scratch_owner.workspace;
    using namespace clubb_kokkos;
    using Host2D = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    using Host3D = Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;

    for (int nbands : {3, 5})
      for (int ndim : {2, 4, 7, 33}) {
        if (nbands == 5 && ndim < 4)
          continue;

        for (int nc : {1, 5, 129})
          for (int nrhs : {0, 1, 3}) {
            Host3D lhs("lhs", nc, nbands, ndim), native_lhs("native_lhs", nbands, nc, ndim);
            Host3D rhs("rhs", nc, ndim, nrhs), native_soln("native_soln", nc, ndim, nrhs);
            Host2D native_single("native_single", nc, ndim),
                native_legacy("native_legacy", nc, ndim);
            const int diagonal = nbands / 2;
            for (int k = 0; k < ndim; ++k)
              for (int i = 0; i < nc; ++i) {
                for (int b = 0; b < nbands; ++b) {
                  const int col       = k - (b - diagonal);
                  lhs(i, b, k)        = (b == diagonal) ? real_t(4) + real_t(0.013) * real_t(i + k)
                                        : (col >= 0 && col < ndim)
                                            ? real_t(0.0017) * real_t(1 + (i + b + k) % 19) *
                                                  ((i + b + k) % 2 ? real_t(-1) : real_t(1))
                                            : real_t(9876.5);
                  native_lhs(b, i, k) = lhs(i, b, k);
                }
                for (int j = 0; j < nrhs; ++j)
                  rhs(i, k, j) =
                      real_t(0.173) * real_t(1 + i + k + j) * (j % 2 ? real_t(-1) : real_t(1));
              }
            auto lhs_d = to_execution_space(lhs), rhs_d = to_execution_space(rhs);
            Real3D soln("soln", nc, ndim, nrhs);
            f2k_lu_reference_c(ndim, nrhs, nc, nbands, rhs.data(), native_lhs.data(),
                               native_soln.data(), native_single.data(), native_legacy.data());
            if (nbands == 3)
              tridiag_lu_solve(workspace, ndim, nrhs, nc, lhs_d, rhs_d, soln);
            else
              penta_lu_solve(workspace, ndim, nrhs, nc, lhs_d, rhs_d, soln);

            auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, soln);
            if (host.extent(0) != std::size_t(nc) || host.extent(1) != std::size_t(ndim) ||
                host.extent(2) != std::size_t(nrhs))
              throw std::runtime_error("LU result shape mismatch");
            for (int j = 0; j < nrhs; ++j)
              for (int k = 0; k < ndim; ++k)
                for (int i = 0; i < nc; ++i) {
                  if (std::memcmp(&host(i, k, j), &native_soln(i, k, j), sizeof(real_t))) {
                    std::cerr << "LU differs: bands=" << nbands << " rows=" << ndim << " col=" << i
                              << " row=" << k << " RHS=" << j << std::hexfloat
                              << " C++=" << host(i, k, j) << " native=" << native_soln(i, k, j)
                              << '\n';
                    throw std::runtime_error("LU solution differs from compiled native bits");
                  }
                  real_t row = real_t(0), scale = std::abs(rhs(i, k, j));
                  for (int b = 0; b < nbands; ++b) {
                    const int col = k - (b - diagonal);
                    if (col >= 0 && col < ndim) {
                      const real_t product = lhs(i, b, k) * host(i, col, j);
                      row += product;
                      scale += std::abs(product);
                    }
                  }
                  if (std::abs(row - rhs(i, k, j)) >
                      real_t(500) * std::numeric_limits<real_t>::epsilon() * scale)
                    throw std::runtime_error("LU residual failed");
                }
            if (nrhs) {
              auto rhs_one = Kokkos::subview(rhs_d, Kokkos::ALL, Kokkos::ALL, 0);
              Real2D single("single", nc, ndim);
              if (nbands == 3)
                tridiag_lu_solve(workspace, ndim, nc, lhs_d, rhs_one, single);
              else
                penta_lu_solve(workspace, ndim, nc, lhs_d, rhs_one, single);
              auto single_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, single);
              if (std::memcmp(single_h.data(), native_single.data(), nc * ndim * sizeof(real_t)))
                throw std::runtime_error("single-RHS LU differs from compiled native bits");

              if (nbands == 3) {
                // The one-matrix interface has no column axis; pack a selected
                // matrix/profile on the device and compare the native overload.
                Real2D matrix("single_matrix", nbands, ndim);
                Real1D vector("single_vector", ndim), legacy("legacy", ndim);
                const int column = nc - 1;
                Kokkos::deep_copy(matrix, Kokkos::subview(lhs_d, column, Kokkos::ALL, Kokkos::ALL));
                Kokkos::deep_copy(vector, Kokkos::subview(rhs_d, column, Kokkos::ALL, 0));
                tridiag_lu_solve(workspace, ndim, matrix, vector, legacy);
                auto legacy_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, legacy);
                for (int k = 0; k < ndim; ++k)
                  if (std::memcmp(&legacy_h(k), &native_legacy(column, k), sizeof(real_t)))
                    throw std::runtime_error("single-matrix LU differs from compiled native bits");
              }
            }
            auto matrix_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, lhs_d);
            auto rhs_h    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, rhs_d);
            if (std::memcmp(matrix_h.data(), lhs.data(), lhs.size() * sizeof(real_t)) ||
                (rhs.size() && std::memcmp(rhs_h.data(), rhs.data(), rhs.size() * sizeof(real_t))))
              throw std::runtime_error("LU modified caller matrix/RHS inputs");
            for (int k = 0; k < ndim; ++k)
              for (int i = 0; i < nc; ++i)
                for (int b = 0; b < nbands; ++b)
                  if (native_lhs(b, i, k) != lhs(i, b, k))
                    throw std::runtime_error("native LU matrix contract differs");
          }
      }
    std::cout << "All five LU overloads match native bits and satisfy A*x=rhs\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
