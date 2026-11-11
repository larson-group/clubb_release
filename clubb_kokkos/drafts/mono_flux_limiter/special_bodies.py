"""Explicit Native spike stage and actual public matrix-layout boundary."""
spike_body=r'''
// Ensure there are no spikes at the top of the domain. Each integral and
// subsequent whole-column scaling retains Native sequential level order.
Kokkos::View<int*,Kokkos::LayoutLeft> spike_warning(
  Kokkos::view_alloc(Kokkos::WithoutInitializing,"monotonic_turbulent_flux_limit::spike_warning"),ngrdcol);
Kokkos::parallel_for("monotonic_turbulent_flux_limit::spike",ngrdcol,KOKKOS_LAMBDA(int column) {
  const int i=column+1;
  spike_warning(column)=0;
  if (Kokkos::abs(xm(i-1,gr_k_ub_zt-1)-xm_enter_mfl(i-1,gr_k_ub_zt-1)) > real_t(10.)*xm_tol) {
    const real_t dz=gr_zm(i-1,gr_k_ub_zm-1)-gr_zm(i-1,gr_k_ub_zm-gr_grid_dir_indx-1);
    const real_t xm_density_weighted=rho_ds_zt(i-1,gr_k_ub_zt-1)
      *(xm(i-1,gr_k_ub_zt-1)-xm_enter_mfl(i-1,gr_k_ub_zt-1))*dz;
    real_t xm_vert_integral=real_t(0.0);
    for (int k=gr_k_lb_zt;gr_grid_dir_indx*(k-(gr_k_ub_zt-gr_grid_dir_indx))<=0;k+=gr_grid_dir_indx)
      xm_vert_integral=xm_vert_integral+rho_ds_zt(i-1,k-1)*xm(i-1,k-1)*gr_grid_dir*gr_dzt(i-1,k-1);
    // Check the integral before dividing, without a new numerical denominator guard.
    if (Kokkos::abs(xm_vert_integral)<eps) {
      spike_warning(column)=1;
      xm(i-1,gr_k_ub_zt-1)=xm_enter_mfl(i-1,gr_k_ub_zt-1);
    } else {
      real_t xm_adj_coef=xm_density_weighted/xm_vert_integral;
      // The source bound is -0.99, rather than the theoretical -1.
      if (xm_adj_coef < -real_t(0.99)) {
        spike_warning(column)=2;
        xm_adj_coef=-real_t(0.99);
      }
      for (int k=1;k<=nzt;++k) xm(i-1,k-1)=xm(i-1,k-1)*(real_t(1.)+xm_adj_coef);
      xm(i-1,gr_k_ub_zt-1)=xm_enter_mfl(i-1,gr_k_ub_zt-1);
    }
  }
});
// Native CPU warnings are preserved at a host boundary; device kernels never print.
auto spike_warning_host=to_host_space(spike_warning);
for (int i=1;i<=ngrdcol;++i) {
  if (spike_warning_host(i-1)==1)
    std::cerr << "Vertical integral of xm is zero;"
                 "mfl will remove spike at top of domain,"
                 "but it will not conserve xm.\n";
  if (spike_warning_host(i-1)==2)
    std::cerr << "xm_adj_coef in mfl less than -0.99, mx_adj_coef set to -0.99\n";
}
'''
solver_body=r'''
switch(solve_type) {
case mono_flux_rtm:solve_type_str="rtm";break;
case mono_flux_thlm:solve_type_str="thlm";break;
default:solve_type_str="scalars";break;
}
// Matrix solves differ by grid direction. Native forces both vertical rows
// and bands to descending order when requested, before solving.
const bool reverse=l_force_descending_solves && gr.grid_dir_indx>0;
Real3D solver_lhs(Kokkos::view_alloc(Kokkos::WithoutInitializing,"mfl_xm_solve::solver_lhs"),ngrdcol,ndiags3,nzt);
Real2D solver_rhs(Kokkos::view_alloc(Kokkos::WithoutInitializing,"mfl_xm_solve::solver_rhs"),ngrdcol,nzt);
Real2D solver_xm(Kokkos::view_alloc(Kokkos::WithoutInitializing,"mfl_xm_solve::solver_xm"),ngrdcol,nzt);
// This is the actual matrix API boundary: Native bands-first -> column-first.
Kokkos::parallel_for("mfl_xm_solve::pack",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
  const int i=cell%ngrdcol,k=cell/ngrdcol;
  const int row=reverse ? nzt-1-k:k;
  solver_rhs(i,k)=rhs(i,row);
  for (int b=0;b<ndiags3;++b) solver_lhs(i,b,k)=lhs(reverse ? ndiags3-1-b:b,i,row);
});
tridiag_solve(solve_type_str,tridiag_solve_method,
               ngrdcol,nzt,l_implemented, // In
               solver_lhs,solver_rhs,err_info, // InOut
               solver_xm); // Out
// Restore Native InOut arrays in their post-flip/post-solver order. They do
// not flip back. Restore the still-unflipped solution before the fatal return.
Kokkos::parallel_for("mfl_xm_solve::restore",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
  const int i=cell%ngrdcol,k=cell/ngrdcol;
  rhs(i,k)=solver_rhs(i,k);
  xm(i,k)=solver_xm(i,k);
  for (int b=0;b<ndiags3;++b) lhs(b,i,k)=solver_lhs(i,b,k);
});
if (f2k::clubb_at_least_debug_level_api(0)) {
  if (any_equal(err_info.err_code,clubb_fatal_error)) {
    std::cerr << err_info.err_header_global << "\n"
              << "calling tridiag_solve in mfl_xm_solve\n";
    return;
  }
}
// Only a successful Native return reaches the final solution direction reversal.
if (reverse) {
  Kokkos::parallel_for("mfl_xm_solve::flip_solution",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
    const int i=cell%ngrdcol,k=cell/ngrdcol;
    xm(i,k)=solver_xm(i,nzt-1-k);
  });
}
return;
'''
# The public source declarations and description are emitted by prepare.py;
# this specialized body owns its Native scalar/string locals directly.
