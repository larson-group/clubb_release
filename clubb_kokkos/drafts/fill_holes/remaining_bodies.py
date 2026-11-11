"""Remaining Native-order bodies, consumed only by the owned prepare script."""

# The scalar-level device overload keeps the exact public three-argument contract.
# Host diagnostics are attached by its concrete wrapper or the batched caller.
device_one_level=r'''
// Device form of the same Native one-level API, used by the batched frozen fill.
// It preserves the ordered species sums; host wrappers own warning I/O.
template<class Input,class Output>
KOKKOS_INLINE_FUNCTION void hole_filling_hm_one_lev(
    int num_hm_fill,Input hm_one_lev, // In
    Output hm_one_lev_filled)        // Out
{
  real_t total_hole=real_t(0.0),total_mass=real_t(0.0);
  int num_neg_hm=0;
  // Initialize and determine missing mass and available same-phase material.
  for (int i=1;i<=num_hm_fill;++i) hm_one_lev_filled(i-1)=real_t(0.0);
  for (int i=1;i<=num_hm_fill;++i) {
    if (hm_one_lev(i-1)<zero) {
      total_hole=total_hole+hm_one_lev(i-1);
      num_neg_hm=num_neg_hm+1;
    } else {
      total_mass=total_mass+hm_one_lev(i-1);
    }
  }
  // No material: leave the field unchanged.
  if (Kokkos::abs(total_mass)<eps) {
    for (int i=1;i<=num_hm_fill;++i) hm_one_lev_filled(i-1)=hm_one_lev(i-1);
    return;
  }
  // Insufficient mass fills holes partly; otherwise remove them and scale donors.
  for (int i=1;i<=num_hm_fill;++i) {
    if (Kokkos::abs(total_hole)>total_mass)
      hm_one_lev_filled(i-1)=(hm_one_lev(i-1)<zero ? hm_one_lev(i-1):zero)*(one+total_mass/total_hole);
    else
      hm_one_lev_filled(i-1)=(hm_one_lev(i-1)>zero ? hm_one_lev(i-1):zero)*(one+total_hole/total_mass);
  }
}
'''
header_top=header_top.replace('#include <limits>','#include <limits>\n#include <Kokkos_MathematicalFunctions.hpp>')
old_write=write
def write():
    old_write()
    h=(folder/'fill_holes.hpp').read_text();position=h.index('// rtm:') if '// rtm:' in h else h.rindex('} // namespace clubb_kokkos')
    # The device form is visible before C++ definitions, with the same Native name.
    h=h[:position]+device_one_level+'\n'+h[position:]
    adapters=r'''
// Moment callers pass strided scalar/column sections. Adapt once at the public
// boundary using the existing array owner; the Native axis meanings do not change.
template<class Dz,class Rho,class Field>
void fill_holes_vertical_api(int nz,int ngrdcol,real_t threshold,
                             int lower_hf_level,int upper_hf_level,
                             Dz dz,Rho rho_ds,int grid_dir_indx,
                             int fill_holes_type,Field field) {
  auto dz_exec=to_execution_space(dz);
  auto rho_exec=to_execution_space(rho_ds);
  auto field_exec=to_execution_space(field);
  fill_holes_vertical_api(nz,ngrdcol,threshold,
                           lower_hf_level,upper_hf_level,
                           dz_exec,rho_exec,grid_dir_indx,
                           fill_holes_type,field_exec);
  from_execution_space(field,field_exec);
}
template<class W,class U,class V>
void fill_holes_wp2_from_horz_tke(int nz,int ngrdcol,real_t threshold,
                                  int lower_hf_level,int upper_hf_level,
                                  W wp2,U up2,V vp2) {
  auto wp2_exec=to_execution_space(wp2);
  auto up2_exec=to_execution_space(up2);
  auto vp2_exec=to_execution_space(vp2);
  fill_holes_wp2_from_horz_tke(nz,ngrdcol,threshold,
                                lower_hf_level,upper_hf_level,
                                wp2_exec,up2_exec,vp2_exec);
  from_execution_space(wp2,wp2_exec);
  from_execution_space(up2,up2_exec);
  from_execution_space(vp2,vp2_exec);
}
'''
    pos=h.rindex('} // namespace clubb_kokkos');h=h[:pos]+adapters+'\n'+h[pos:]
    (folder/'fill_holes.hpp').write_text(h)
    c=(folder/'fill_holes.cpp').read_text();pos=c.index('constexpr bool l_debug')
    helper=r'''
// Shared host rendering of the Native one-level warnings. It changes no fields.
template<class Input,class Output>
void print_hm_fill_diagnostics(int num_hm_fill,Input hm_one_lev,Output hm_one_lev_filled) {
  std::cout << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  real_t total_hole=zero,total_mass=zero;
  for (int i=1;i<=num_hm_fill;++i) {
    if (hm_one_lev(i-1)<zero) total_hole=total_hole+hm_one_lev(i-1);
    else total_mass=total_mass+hm_one_lev(i-1);
  }
  if (Kokkos::abs(total_mass)<eps) {
    std::cout << "Warning: One-level hole filling was not successful! total_mass ~= 0\n";
    return;
  }
  if (Kokkos::abs(total_hole)>total_mass) {
    for (int i=1;i<=num_hm_fill;++i)
      std::cout << "Warning: One-level hole filling was not able to fill holes completely!"
                   " The holes were filled partially. |total_hole| > total_mass\n";
  }
  real_t sum_before=zero,sum_after=zero;
  bool l_negative=false;
  for (int i=1;i<=num_hm_fill;++i) {
    sum_before=sum_before+hm_one_lev(i-1);
    sum_after=sum_after+hm_one_lev_filled(i-1);
    if (hm_one_lev_filled(i-1)<zero) l_negative=true;
  }
  if (Kokkos::abs(sum_before-sum_after)>Kokkos::abs(sum_before+sum_after)*eps/2)
    std::cout << "Warning: Hole filling was not conservative!\n";
  if (l_negative) std::cout << "Warning: Hole filling failed! A hole could not be filled.\n";
}
'''
    c=c[:pos]+helper+c[pos:]
    # Diagnostics preserve values, with C++ max_digits10 rendering on the host.
    c=c.replace('if (l_clubb_at_least_debug_level_3) {\n    // Mirror',
      'if (l_clubb_at_least_debug_level_3) {\n    std::cerr << std::setprecision(std::numeric_limits<real_t>::max_digits10);\n    // Mirror')
    (folder/'fill_holes.cpp').write_text(c)

bodies['hole_filling_hm_one_lev']=r'''
  // The concrete entry keeps diagnostics on the host. A batched caller invokes
  // the identical device overload on its per-level sections without kernel launches.
  Kokkos::parallel_for("hole_filling_hm_one_lev::level",1,KOKKOS_LAMBDA(int) {
    auto input=Kokkos::subview(hm_one_lev,Kokkos::ALL);
    auto output=Kokkos::subview(hm_one_lev_filled,Kokkos::ALL);
    hole_filling_hm_one_lev(num_hm_fill,input,output);
  });
  if (f2k::clubb_at_least_debug_level_api(2)) {
    auto input=to_host_space(hm_one_lev);
    auto output=to_host_space(hm_one_lev_filled);
    print_hm_fill_diagnostics(num_hm_fill,input,output);
  }
'''
# A full-range rank1 subview can be the exact Real1D type on Serial; explicitly
# select the device template rather than recursively resolving the host entry.
bodies['hole_filling_hm_one_lev']=bodies['hole_filling_hm_one_lev'].replace('hole_filling_hm_one_lev(num_hm_fill,input,output);','hole_filling_hm_one_lev<decltype(input),decltype(output)>(num_hm_fill,input,output);')

bodies['fill_holes_hydromet_api']=r'''
  auto frozen_flags=to_host_space(l_frozen_hm);
  auto mixing_flags=to_host_space(l_mix_rat_hm);
  int num_frozen_hm=0;
  // Determine the number of frozen mixing ratios, preserving species order.
  for (int i=1;i<=hydromet_dim;++i)
    if (frozen_flags(i-1) && mixing_flags(i-1)) num_frozen_hm=num_frozen_hm+1;
  Real2D hydromet_frozen(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_hydromet_api::hydromet_frozen"),nzt,num_frozen_hm);
  Real2D hydromet_frozen_filled(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_hydromet_api::hydromet_frozen_filled"),nzt,num_frozen_hm);
  // Pack the exact Native frozen subset; zero padding would alter sum association.
  int j=1;
  for (int i=1;i<=hydromet_dim;++i) {
    if (frozen_flags(i-1) && mixing_flags(i-1)) {
      Kokkos::deep_copy(Kokkos::subview(hydromet_frozen,Kokkos::ALL,j-1),
                        Kokkos::subview(hydromet,Kokkos::ALL,i-1));
      j=j+1;
    }
  }
  // Same-phase redistribution is sequential in species, independent in height.
  Kokkos::parallel_for("fill_holes_hydromet_api::levels",nzt,KOKKOS_LAMBDA(int level) {
    bool l_holes=false;
    for (int j=1;j<=num_frozen_hm;++j) if (hydromet_frozen(level,j-1)<zero) l_holes=true;
    auto input=Kokkos::subview(hydromet_frozen,level,Kokkos::ALL);
    auto output=Kokkos::subview(hydromet_frozen_filled,level,Kokkos::ALL);
    if (l_holes) {
      hole_filling_hm_one_lev<decltype(input),decltype(output)>(num_frozen_hm,input,output);
    } else {
      for (int j=1;j<=num_frozen_hm;++j) output(j-1)=input(j-1);
    }
  });
  if (f2k::clubb_at_least_debug_level_api(2)) {
    auto input=to_host_space(hydromet_frozen);
    auto output=to_host_space(hydromet_frozen_filled);
    for (int i=1;i<=nzt;++i) {
      bool l_holes=false;
      for (int j=1;j<=num_frozen_hm;++j) if (input(i-1,j-1)<zero) l_holes=true;
      if (l_holes) print_hm_fill_diagnostics(num_frozen_hm,
                    Kokkos::subview(input,i-1,Kokkos::ALL),
                    Kokkos::subview(output,i-1,Kokkos::ALL));
    }
  }
  // Reassemble the filled array; liquid species and concentrations are unchanged.
  j=1;
  for (int i=1;i<=hydromet_dim;++i) {
    if (frozen_flags(i-1) && mixing_flags(i-1)) {
      Kokkos::deep_copy(Kokkos::subview(hydromet_filled,Kokkos::ALL,i-1),
                        Kokkos::subview(hydromet_frozen_filled,Kokkos::ALL,j-1));
      j=j+1;
    } else {
      Kokkos::deep_copy(Kokkos::subview(hydromet_filled,Kokkos::ALL,i-1),
                        Kokkos::subview(hydromet,Kokkos::ALL,i-1));
    }
  }
  // Native liquid-phase same-level hole filling remains unimplemented.
'''
bodies['fill_holes_wv']=r'''
  // An unknown species terminates only upon reaching an actual hole, after the
  // Native vapor-tendency update. Keep this exceptional path ordered on the host.
  std::string name=hydromet_name;
  while (!name.empty() && name.back()==' ') name.pop_back();
  const bool l_rain=name=="rrm";
  const bool l_frozen=name=="rim" || name=="rsm" || name=="rgm";
  if (!l_rain && !l_frozen) {
    auto hm=to_host_space(hydromet);
    for (int k=1;k<=nzt;++k) {
      if (hm(k-1)<zero_threshold) {
        auto vapor=to_host_space(rvm_mc);
        const real_t rvm_clip_tndcy=hm(k-1)/dt;
        vapor(k-1)=vapor(k-1)+rvm_clip_tndcy;
        from_execution_space(rvm_mc,to_execution_space(vapor));
        std::cerr << "Fatal error in microphys_driver\n";
        std::exit(EXIT_FAILURE); // Native ERROR STOP, after the vapor update.
      }
    }
    return;
  }
  Kokkos::parallel_for("fill_holes_wv::levels",nzt,KOKKOS_LAMBDA(int level) {
    const int k=level+1;
    if (hydromet(k-1)<zero_threshold) {
      // Tendency transferred from the hydrometeor to vapor.
      const real_t rvm_clip_tndcy=hydromet(k-1)/dt;
      rvm_mc(k-1)=rvm_mc(k-1)+rvm_clip_tndcy;
      // Evaporation and sublimation have different latent heat compensation.
      if (l_rain) thlm_mc(k-1)=thlm_mc(k-1)-rvm_clip_tndcy*(Lv/(Cp*exner(k-1)));
      else thlm_mc(k-1)=thlm_mc(k-1)-rvm_clip_tndcy*(Ls/(Cp*exner(k-1)));
      hydromet(k-1)=zero_threshold;
    }
  });
'''

bodies['fill_holes_driver_api']=r'''
  Real3D hydromet_filled(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_driver_api::hydromet_filled"),ngrdcol,nzt,hydromet_dim);
  Real3D hydromet_clipped(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_driver_api::hydromet_clipped"),ngrdcol,nzt,hydromet_dim);
  real_t max_velocity;
  std::string name_bt,name_hf,name_wvhf,name_cl,name_mc,hydromet_name;
  constexpr bool l_hole_fill=true;
  auto hydromet_list_host=to_host_space(hm_metadata.hydromet_list);
  auto l_mix_rat_hm_host=to_host_space(hm_metadata.l_mix_rat_hm);
  auto l_frozen_hm_host=to_host_space(hm_metadata.l_frozen_hm);
  auto hydromet_tol_host=to_host_space(hm_metadata.hydromet_tol);
  // Stats payloads are staged only at the existing host stats API boundary.
  auto budget=[&](const std::string& name,int ihm,bool begin) -> void {
    Real2D values(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_driver_api::stats_values"),ngrdcol,nzt);
    Kokkos::parallel_for("fill_holes_driver_api::stats_values",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
      const int i=cell%ngrdcol,k=cell/ngrdcol;
      values(i,k)=hydromet(i,k,ihm-1)/dt;
    });
    if (begin) f2k::stats_begin_budget(name,values,stats);
    else f2k::stats_finalize_budget(name,values,stats);
  };
  // Start _hf budgets before any same-phase or vertical filling.
  if (stats.flags().l_sample) {
    for (int ihm=1;ihm<=hydromet_dim;++ihm) {
      setup_stats_names(ihm,hydromet_dim,hm_metadata.hydromet_list,
                        max_velocity,&name_bt,&name_hf,&name_wvhf,&name_cl,&name_mc);
      budget(name_hf,ihm,true);
    }
  }
  if (any_below(hydromet,zero_threshold) && l_fill_holes_hm) {
    // Preserve single-column public APIs at their actual boundaries. Packing a
    // strided row reuses the established array owner; no Fortran scalar call.
    for (int i=1;i<=ngrdcol;++i) {
      auto input=to_execution_space(Kokkos::subview(hydromet,i-1,Kokkos::ALL,Kokkos::ALL));
      auto output=to_execution_space(Kokkos::subview(hydromet_filled,i-1,Kokkos::ALL,Kokkos::ALL),false);
      fill_holes_hydromet_api(nzt,hydromet_dim,input,
                              hm_metadata.l_frozen_hm,hm_metadata.l_mix_rat_hm,
                              output);
      from_execution_space(Kokkos::subview(hydromet_filled,i-1,Kokkos::ALL,Kokkos::ALL),output);
    }
    Kokkos::deep_copy(hydromet,hydromet_filled);
  }
  for (int ihm=1;ihm<=hydromet_dim;++ihm) {
    setup_stats_names(ihm,hydromet_dim,hm_metadata.hydromet_list,
                      max_velocity,&name_bt,&name_hf,&name_wvhf,&name_cl,&name_mc);
    hydromet_name.assign(hydromet_list_host.data()+(ihm-1)*hydromet_list_host.extent(0),hydromet_list_host.extent(0));
    while (!hydromet_name.empty() && hydromet_name.back()==' ') hydromet_name.pop_back();
    if (f2k::clubb_at_least_debug_level_api(1)) {
      auto values=to_host_space(hydromet);
      for (int k=1;k<=nzt;++k) for (int i=1;i<=ngrdcol;++i)
        if (values(i-1,k-1,ihm-1)<zero_threshold)
          std::cerr << hydromet_name << " < " << zero_threshold
                    << " in fill_holes_driver_api at i, k= " << i << " " << k << "\n";
    }
    auto species=to_execution_space(Kokkos::subview(hydromet,Kokkos::ALL,Kokkos::ALL,ihm-1));
    if (any_below(species,zero_threshold) && !hydromet_name.empty() && hydromet_name[0]=='r' && l_hole_fill) {
      for (int i=1;i<=ngrdcol;++i) {
        auto column=to_execution_space(Kokkos::subview(hydromet,std::make_pair(i-1,i),Kokkos::ALL,ihm-1));
        if (any_below(column,zero_threshold)) {
          auto dz=to_execution_space(Kokkos::subview(gr.dzt,std::make_pair(i-1,i),Kokkos::ALL));
          auto rho=to_execution_space(Kokkos::subview(rho_ds_zt,std::make_pair(i-1,i),Kokkos::ALL));
          // Native deliberately passes1..nzt and positive direction for these rows.
          fill_holes_vertical_api(nzt,1,zero_threshold,
                                   1,nzt,dz,rho,1,
                                   fill_holes_type,column);
          from_execution_space(Kokkos::subview(hydromet,std::make_pair(i-1,i),Kokkos::ALL,ihm-1),column);
        }
      }
    }
    if (stats.flags().l_sample) {
      budget(name_hf,ihm,false);
      budget(name_wvhf,ihm,true);
    }
    if (!hydromet_name.empty() && hydromet_name[0]=='r' && l_hole_fill) {
      for (int i=1;i<=ngrdcol;++i) {
        auto hm=to_execution_space(Kokkos::subview(hydromet,i-1,Kokkos::ALL,ihm-1));
        if (any_below(hm,zero_threshold)) {
          // Residual holes can require vapor, notably ASEX A209 with LHS sampling.
          auto exner_column=to_execution_space(Kokkos::subview(exner,i-1,Kokkos::ALL));
          auto vapor=to_execution_space(Kokkos::subview(rvm_mc,i-1,Kokkos::ALL));
          auto temperature=to_execution_space(Kokkos::subview(thlm_mc,i-1,Kokkos::ALL));
          fill_holes_wv(nzt,dt,exner_column,hydromet_name,
                         vapor,temperature,hm);
          from_execution_space(Kokkos::subview(rvm_mc,i-1,Kokkos::ALL),vapor);
          from_execution_space(Kokkos::subview(thlm_mc,i-1,Kokkos::ALL),temperature);
          from_execution_space(Kokkos::subview(hydromet,i-1,Kokkos::ALL,ihm-1),hm);
        }
      }
    }
    if (stats.flags().l_sample) budget(name_wvhf,ihm,false);
    if (l_mix_rat_hm_host(ihm-1)) {
      if (stats.flags().l_sample) budget(name_cl,ihm,true);
      const bool l_frozen=l_frozen_hm_host(ihm-1);
      const real_t hydromet_tol=hydromet_tol_host(ihm-1);
      Kokkos::parallel_for("fill_holes_driver_api::mixing_ratio_clip",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
        const int i=cell%ngrdcol,k=cell/ngrdcol;
        // Clip negative values, then transfer very small positive mass to vapor.
        if (hydromet(i,k,ihm-1)<zero_threshold) hydromet(i,k,ihm-1)=zero_threshold;
        if (hydromet(i,k,ihm-1)<=hydromet_tol) {
          rvm_mc(i,k)=rvm_mc(i,k)+hydromet(i,k,ihm-1)/dt;
          if (!l_frozen)
            thlm_mc(i,k)=thlm_mc(i,k)-(Lv/(Cp*exner(i,k)))*(hydromet(i,k,ihm-1)/dt);
          else
            thlm_mc(i,k)=thlm_mc(i,k)-(Ls/(Cp*exner(i,k)))*(hydromet(i,k,ihm-1)/dt);
          hydromet(i,k,ihm-1)=zero;
        }
      });
      if (stats.flags().l_sample) budget(name_cl,ihm,false);
    }
  }
  // Concentration clipping does not change the mixing-ratio fields.
  for (int i=1;i<=ngrdcol;++i) {
    auto input=to_execution_space(Kokkos::subview(hydromet,i-1,Kokkos::ALL,Kokkos::ALL));
    auto output=to_execution_space(Kokkos::subview(hydromet_clipped,i-1,Kokkos::ALL,Kokkos::ALL),false);
    clip_hydromet_conc_mvr(nzt,hydromet_dim,hm_metadata,input,output);
    from_execution_space(Kokkos::subview(hydromet_clipped,i-1,Kokkos::ALL,Kokkos::ALL),output);
  }
  for (int ihm=1;ihm<=hydromet_dim;++ihm) {
    if (!l_mix_rat_hm_host(ihm-1)) {
      setup_stats_names(ihm,hydromet_dim,hm_metadata.hydromet_list,
                        max_velocity,&name_bt,&name_hf,&name_wvhf,&name_cl,&name_mc);
      if (stats.flags().l_sample) budget(name_cl,ihm,true);
      Kokkos::deep_copy(Kokkos::subview(hydromet,Kokkos::ALL,Kokkos::ALL,ihm-1),
                        Kokkos::subview(hydromet_clipped,Kokkos::ALL,Kokkos::ALL,ihm-1));
      if (stats.flags().l_sample) budget(name_cl,ihm,false);
    }
  }
'''

bodies['clip_hydromet_conc_mvr']=r'''
  auto mix=to_host_space(hm_metadata.l_mix_rat_hm);
  auto frozen=to_host_space(hm_metadata.l_frozen_hm);
  for (int idx=1;idx<=hydromet_dim;++idx) {
    if (!mix(idx-1)) {
      // Precompute public Native metadata functions on the host once per species.
      // No scalar bridge enters the independent level kernel.
      real_t Nxm_min_coef;
      const real_t mvr_hydromet_max=f2k::mvr_hm_max(idx,hm_metadata);
      if (!frozen(idx-1)) {
        // <Nr>=<rr>/((4/3)*pi*rho_lw*mvr_rain_max^3).
        Nxm_min_coef=one/(four_thirds*pi*rho_lw*((mvr_hydromet_max*mvr_hydromet_max)*mvr_hydromet_max));
      } else {
        // Frozen species use rho_ice and their species-specific maximum radius.
        Nxm_min_coef=one/(four_thirds*pi*rho_ice*((mvr_hydromet_max*mvr_hydromet_max)*mvr_hydromet_max));
      }
      const int rx_idx=f2k::Nx2rx_hm_idx(idx,hm_metadata);
      Kokkos::parallel_for("clip_hydromet_conc_mvr::levels",nzt,KOKKOS_LAMBDA(int level) {
        const int k=level+1;
        if (hydromet(k-1,rx_idx-1)>zero) {
          const real_t minimum=Nxm_min_coef*hydromet(k-1,rx_idx-1);
          hydromet_clipped(k-1,idx-1)=hydromet(k-1,idx-1)>minimum ? hydromet(k-1,idx-1):minimum;
        } else {
          hydromet_clipped(k-1,idx-1)=zero;
        }
      });
    } else {
      Kokkos::deep_copy(Kokkos::subview(hydromet_clipped,Kokkos::ALL,idx-1),
                        Kokkos::subview(hydromet,Kokkos::ALL,idx-1));
    }
  }
'''

# Keep all nine case bodies and their literal strings/speeds directly from Native.
parts=[]
for s in statements(blocks['setup_stats_names']):
    if s.lower().startswith('select case'):break
    if re.match(r'(?:max_velocity|name_\w+_local)\s*=',s):parts.append(s)
setup='\n'.join(expr(s,'setup_stats_names')+';' for s in parts)
setup+='\nauto hydromet_list_host=to_host_space(hydromet_list);\nstd::string species(hydromet_list_host.data()+(ihm-1)*hydromet_list_host.extent(0),hydromet_list_host.extent(0));\nwhile (!species.empty() && species.back()==\' \') species.pop_back();\n'
cases=re.search(r'select case.*?end select',blocks['setup_stats_names'],re.S|re.I)[0]
current=None;first=True
for token in tokens(cases):
    s=token.text
    if token.kind=='comment':setup+='  //'+s+'\n';continue
    if s.lower().startswith('select'):continue
    if s.lower().startswith('case'):
        if current is not None:setup+='}\n'
        if 'default' in s.lower():setup+='else {\n'
        else:
            label=s[s.index('(')+1:s.rindex(')')];setup+=('if' if first else 'else if')+' (species=='+label.strip()+') {\n';first=False
        current=s
    elif s.lower().startswith('end select'):setup+='}\n'
    else:setup+='  '+expr(s,'setup_stats_names')+';\n'
for field in ('bt','hf','wvhf','cl','mc'):setup+='if (name_'+field+' != nullptr) *name_'+field+'=name_'+field+'_local;\n'
bodies['setup_stats_names']='  std::string name_bt_local,name_hf_local,name_wvhf_local,name_cl_local,name_mc_local;\n'+'\n'.join('  '+s for s in setup.splitlines())
