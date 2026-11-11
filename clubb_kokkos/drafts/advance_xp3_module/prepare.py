"""Assemble the isolated third-moment draft with canonical scientific comments."""
from pathlib import Path
import re,sys,textwrap
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,declaration,statements
path=root/'src/CLUBB_core/advance_xp3_module.F90';source=path.read_text()
names=['advance_xp3','diagnose_xp3','advance_xp3_simplified','term_tp_rhs','term_ac_rhs']
blocks={n:re.search(r'^  (?:subroutine|function) '+n+r'\(.*?^  end (?:subroutine|function) '+n+r'\b',source,re.M|re.S).group() for n in names}
contracts={n:read_procedure(path,n) for n in names}
metadata={};docs={}
for n,b in blocks.items():
 info={}
 for st in statements(b):
  if '::' in st:info.update(declaration(st))
 metadata[n]=info;docs[n]={}
 last=b.rfind('::');end=b.find('\n\n',last);prefix=b[:end if end>=0 else len(b)]
 for key,a in info.items():
  for raw in prefix.splitlines():
   if '!' not in raw or raw.strip().startswith('!'):continue
   code,comment=raw.split('!',1)
   note=comment.strip()
   if re.fullmatch(r'(?:intent\(\w+\)|in|out|in/?out)',note,re.I):continue
   if re.search(r'\b'+re.escape(a.name)+r'\s*(?:,\s*&|&)?\s*$',code):docs[n][key]=note;break

def description(n):
 b=blocks[n];begin=b.index('! Description');end=re.search(r'^\s*use\b',b[begin:],re.M).start()+begin
 return '\n'.join('//'+line.strip()[1:] for line in b[begin:end].splitlines() if line.strip().startswith('!'))+'\n'

def signature(n,prototype=False):
 p=contracts[n];array_args=[a for a in p.arguments if a.dims];parameters=[]
 for a in p.arguments:
  if a.dims:typ=('Array'+str(array_args.index(a))) if n=='advance_xp3_simplified' else 'Real'+str(len(a.dims))+'D'
  elif a.kind=='grid':typ='const grid&'
  elif a.kind=='stats_type':typ='f2k::stats_type&'
  else:typ={'real':'real_t','integer':'int','logical':'bool'}[a.kind]
  parameters.append(typ+' '+a.name)
 prefix='template<'+','.join('class Array'+str(j) for j in range(len(array_args)))+'>\n' if n=='advance_xp3_simplified' else 'KOKKOS_INLINE_FUNCTION ' if n.startswith('term_') else ''
 result='real_t' if p.result else 'void'
 if prototype:return prefix+result+' '+n+'('+','.join(parameters)+')'
 lines=[]
 for j,(a,decl) in enumerate(zip(p.arguments,parameters)):
  note=a.intent.title()+(': '+docs[n][a.name.lower()] if a.name.lower() in docs[n] else '')
  if len(decl+note)>83 or j==len(parameters)-1:
   lines.extend('// '+line for line in textwrap.wrap(note,width=83));lines.append(decl+(',' if j<len(parameters)-1 else ')'))
  else:lines.append(decl+', // '+note)
 return prefix+result+' '+n+'(\n'+'\n'.join(lines)

def locals(n):
 args={a.name.lower() for a in contracts[n].arguments};result=contracts[n].result
 out=['//----------------------- Local Variables -----------------------']
 for key,a in metadata[n].items():
  if key in args or (result and key==result.name.lower()):continue
  if a.dims:
   returned=bool(re.search(r'\b'+re.escape(a.name)+r'\s*=\s*(?:zm2zt_api|zt2zm_api|ddzm)\(',{'advance_xp3':advance,'diagnose_xp3':diagnose,'advance_xp3_simplified':simplified}[n]))
   if returned or a.name=='xp3_old':out.append('Real'+str(len(a.dims))+'D '+a.name+';'+(' // '+docs[n][key] if key in docs[n] else ''))
   else:out.append('Real'+str(len(a.dims))+'D '+a.name+'(Kokkos::view_alloc(Kokkos::WithoutInitializing,"'+n+'::'+a.name+'"),'+','.join(a.dims)+');'+(' // '+docs[n][key] if key in docs[n] else ''))
 # Native loop counters become kernel-local/host-loop scoped; no mutable capture.
 return '\n'.join(out)+'\n'

advance='''wp2_zt = zm2zt_api(nzm,nzt,ngrdcol,gr,wp2,w_tol_sqd);

// Advance <rt'^3> one model timestep or calculate <rt'^3> using a
// steady-state approximation.
advance_xp3_simplified(nzm,nzt,ngrdcol,gr,xp3_rtp3,dt,
 rtm,rtp2,wprtp,wprtp2,rho_ds_zm,invrs_rho_ds_zt,
 invrs_tau_zt,tau_max_zt,rt_tol,l_lmm_stepping,stats,rtp3);

// Advance <thl'^3> one model timestep or calculate <thl'^3> using a
// steady-state approximation.
advance_xp3_simplified(nzm,nzt,ngrdcol,gr,xp3_thlp3,dt,
 thlm,thlp2,wpthlp,wpthlp2,rho_ds_zm,invrs_rho_ds_zt,
 invrs_tau_zt,tau_max_zt,thl_tol,l_lmm_stepping,stats,thlp3);

// Advance <sclr'^3> one model timestep or calculate <sclr'^3> using a
// steady-state approximation. Stage scalar thresholds at the typed scalar
// boundary once; the source's per-scalar lower-call order is retained.
auto sclr_tol_host=to_host_space(sclr_tol);
for(int sclr=0;sclr<sclr_dim;++sclr) {
 advance_xp3_simplified(nzm,nzt,ngrdcol,gr,xp3_sclrp3,dt,
  Kokkos::subview(sclrm,Kokkos::ALL,Kokkos::ALL,sclr),
  Kokkos::subview(sclrp2,Kokkos::ALL,Kokkos::ALL,sclr),
  Kokkos::subview(wpsclrp,Kokkos::ALL,Kokkos::ALL,sclr),
  Kokkos::subview(wpsclrp2,Kokkos::ALL,Kokkos::ALL,sclr),rho_ds_zm,
  invrs_rho_ds_zt,invrs_tau_zt,tau_max_zt,sclr_tol_host(sclr),
  l_lmm_stepping,stats,Kokkos::subview(sclrp3,Kokkos::ALL,Kokkos::ALL,sclr));
}

// Use a modified form of the Larson and Golaz (2005) ansatz for the
// ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
Skx_func(nzt,ngrdcol,wp2_zt,wp3,w_tol,clubb_params,Skw_zt);
wpthlp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,wpthlp);
wprtp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,wprtp);
// Positive def. quantity.
thlp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,thlp2,thl_tol*thl_tol);
// Positive def. quantity.
rtp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,rtp2,rt_tol*rt_tol);
upwp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,upwp);
vpwp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,vpwp);
// Positive def. quantity.
up2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,up2,w_tol_sqd);
// Positive def. quantity.
vp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,vp2,w_tol_sqd);
thvm_zm=zt2zm_api(nzm,nzt,ngrdcol,gr,thvm,zero_threshold);
ddzm_thvm_zm=ddzm(nzm,nzt,ngrdcol,gr,thvm_zm);
Kokkos::parallel_for("advance_xp3::buoyancy_frequency",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
 const int i=cell%ngrdcol,k=cell/ngrdcol;
 const real_t stability=(grav/thvm(i,k))*ddzm_thvm_zm(i,k);
 brunt_vaisala_freq_sqd_zt(i,k)=stability>zero?stability:zero;
});

// The xp3_coef_fnc provides some extra tunability to the simple xp3
// equation.
// When xp3_coef_fnc goes to 0, the value of Skx goes to the smallest
// magnitude permitted by the function. When xp3_coef_fnc goes to 1, the
// magnitude of Skx becomes huge.
Kokkos::parallel_for("advance_xp3::xp3_coef_fnc",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
 const int i=cell%ngrdcol,k=cell/ngrdcol;
 xp3_coef_fnc(i,k)=clubb_params(i,ixp3_coef_base)
  +(one-clubb_params(i,ixp3_coef_slope))
   *(one-Kokkos::exp(brunt_vaisala_freq_sqd_zt(i,k)/clubb_params(i,ixp3_coef_slope)));
});
xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,upwp_zt,wp2_zt,
 up2_zt,xp3_coef_fnc,clubb_params,w_tol,up3);
xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,vpwp_zt,wp2_zt,
 vp2_zt,xp3_coef_fnc,clubb_params,w_tol,vp3);

if(stats.flags().l_sample) {
 f2k::stats_update("thlp2_zt",thlp2_zt,stats);
 f2k::stats_update("wpthlp_zt",wpthlp_zt,stats);
 f2k::stats_update("wprtp_zt",wprtp_zt,stats);
 f2k::stats_update("rtp2_zt",rtp2_zt,stats);
 f2k::stats_update("up2_zt",up2_zt,stats);
 f2k::stats_update("vp2_zt",vp2_zt,stats);
 f2k::stats_update("upwp_zt",upwp_zt,stats);
 f2k::stats_update("vpwp_zt",vpwp_zt,stats);
}
'''

diagnose='''//----------------------------- Begin Code ------------------------------
wp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,wp2,w_tol_sqd);
// The ADG1 PDF must use this option.
Skx_func(nzt,ngrdcol,wp2_zt,wp3,w_tol,clubb_params,Skw_zt);
wpthlp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,wpthlp);
wprtp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,wprtp);
// Positive def. quantity.
thlp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,thlp2,thl_tol*thl_tol);
// Positive def. quantity.
rtp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,rtp2,rt_tol*rt_tol);
upwp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,upwp);
vpwp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,vpwp);
// Positive def. quantity.
up2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,up2,w_tol_sqd);
// Positive def. quantity.
vp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,vp2,w_tol_sqd);
// The typed ansatz API receives one threshold scalar per Native scalar call.
auto sclr_tol_host=to_host_space(sclr_tol);

if(iiPDF_type==iiPDF_ADG1) {
 // Use the Larson and Golaz (2005) ansatz for the ADG1 PDF to
 // calculate <rt'^3>, <thl'^3>, <u'^3>, <v'^3>, and <sclr'^3>.
 sigma_sqd_w_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,sigma_sqd_w,zero_threshold);
 xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,wpthlp_zt,wp2_zt,
  thlp2_zt,sigma_sqd_w_zt,clubb_params,thl_tol,thlp3);
 xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,wprtp_zt,wp2_zt,
  rtp2_zt,sigma_sqd_w_zt,clubb_params,rt_tol,rtp3);
 xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,upwp_zt,wp2_zt,
  up2_zt,sigma_sqd_w_zt,clubb_params,w_tol,up3);
 xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,vpwp_zt,wp2_zt,
  vp2_zt,sigma_sqd_w_zt,clubb_params,w_tol,vp3);
 for(int sclr=0;sclr<sclr_dim;++sclr) {
  // Keep the Native ADG1 scalar-flux minimum as well as the variance minimum.
  wpsclrp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,
   Kokkos::subview(wpsclrp,Kokkos::ALL,Kokkos::ALL,sclr),
   sclr_tol_host(sclr)*sclr_tol_host(sclr));
  sclrp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,
   Kokkos::subview(sclrp2,Kokkos::ALL,Kokkos::ALL,sclr),
   sclr_tol_host(sclr)*sclr_tol_host(sclr));
  // Existing array owner handles the typed ansatz output boundary. A compatible
  // contiguous scalar section aliases its allocation; otherwise restore it.
  auto sclrp3_section=Kokkos::subview(sclrp3,Kokkos::ALL,Kokkos::ALL,sclr);
  auto sclrp3_field=to_execution_space(sclrp3_section,false);
  xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,wpsclrp_zt,wp2_zt,
   sclrp2_zt,sigma_sqd_w_zt,clubb_params,sclr_tol_host(sclr),sclrp3_field);
  from_execution_space(sclrp3_section,sclrp3_field);
 }
} else {
 // Use a modified form of the Larson and Golaz (2005) ansatz for the
 // ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
 thvm_zm=zt2zm_api(nzm,nzt,ngrdcol,gr,thvm,zero_threshold);
 ddzm_thvm_zm=ddzm(nzm,nzt,ngrdcol,gr,thvm_zm);
 Kokkos::parallel_for("diagnose_xp3::buoyancy_frequency",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
  const int i=cell%ngrdcol,k=cell/ngrdcol;
  const real_t stability=(grav/thvm(i,k))*ddzm_thvm_zm(i,k);
  brunt_vaisala_freq_sqd_zt(i,k)=stability>zero?stability:zero;
 });
 // Initialize sigma_sqd_w_zt to zero so we don't break output.
 sigma_sqd_w_zt=Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing,"diagnose_xp3::sigma_sqd_w_zt"),ngrdcol,nzt);
 Kokkos::deep_copy(sigma_sqd_w_zt,zero);

 // The xp3_coef_fnc is used in place of sigma_sqd_w_zt when the
 // ADG1 PDF is not being used. The xp3_coef_fnc provides some extra
 // tunability to the simple xp3 equation.
 // When xp3_coef_fnc goes to 0, the value of Skx goes to the smallest
 // magnitude permitted by the function. When xp3_coef_fnc goes to 1,
 // the magnitude of Skx becomes huge.
 // The value of Skx becomes large near cloud top, where there is a
 // higher degree of static stability. The exp{ } portion of the
 // xp3_coef_fnc allows the xp3_coef_fnc to become larger in regions
 // of high static stability, producing larger magnitude values of Skx.
 Kokkos::parallel_for("diagnose_xp3::xp3_coef_fnc",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
  const int i=cell%ngrdcol,k=cell/ngrdcol;
  xp3_coef_fnc(i,k)=clubb_params(i,ixp3_coef_base)
   +(one-clubb_params(i,ixp3_coef_slope))
    *(one-Kokkos::exp(brunt_vaisala_freq_sqd_zt(i,k)/clubb_params(i,ixp3_coef_slope)));
 });
 xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,wpthlp_zt,wp2_zt,
  thlp2_zt,xp3_coef_fnc,clubb_params,thl_tol,thlp3);
 xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,wprtp_zt,wp2_zt,
  rtp2_zt,xp3_coef_fnc,clubb_params,rt_tol,rtp3);
 xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,upwp_zt,wp2_zt,
  up2_zt,xp3_coef_fnc,clubb_params,w_tol,up3);
 xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,vpwp_zt,wp2_zt,
  vp2_zt,xp3_coef_fnc,clubb_params,w_tol,vp3);
 for(int sclr=0;sclr<sclr_dim;++sclr) {
  // Native non-ADG1 scalar-flux interpolation has no minimum.
  wpsclrp_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,
   Kokkos::subview(wpsclrp,Kokkos::ALL,Kokkos::ALL,sclr));
  sclrp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,
   Kokkos::subview(sclrp2,Kokkos::ALL,Kokkos::ALL,sclr),
   sclr_tol_host(sclr)*sclr_tol_host(sclr));
  auto sclrp3_section=Kokkos::subview(sclrp3,Kokkos::ALL,Kokkos::ALL,sclr);
  auto sclrp3_field=to_execution_space(sclrp3_section,false);
  xp3_LG_2005_ansatz(nzt,ngrdcol,Skw_zt,wpsclrp_zt,wp2_zt,
   sclrp2_zt,xp3_coef_fnc,clubb_params,sclr_tol_host(sclr),sclrp3_field);
  from_execution_space(sclrp3_section,sclrp3_field);
 }
}

if(stats.flags().l_sample) {
 f2k::stats_update("thlp2_zt",thlp2_zt,stats);
 f2k::stats_update("wpthlp_zt",wpthlp_zt,stats);
 f2k::stats_update("wprtp_zt",wprtp_zt,stats);
 f2k::stats_update("rtp2_zt",rtp2_zt,stats);
 f2k::stats_update("up2_zt",up2_zt,stats);
 f2k::stats_update("vp2_zt",vp2_zt,stats);
 f2k::stats_update("upwp_zt",upwp_zt,stats);
 f2k::stats_update("vpwp_zt",vpwp_zt,stats);
}
'''

simplified='''// Coefficient in the <x'^3> turbulent dissipation term [-].
constexpr real_t C_xp3_dissipation=real_t(1.0);
// Flag to either predict <x'^3> or use steady-state approximation.
constexpr bool l_predict_xp3=false;
std::string name_bt,name_tp,name_ac,name_dp;
// The saved field is read only by the Native predictive/LMM branch.
if(l_predict_xp3 && l_lmm_stepping)
 xp3_old=Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing,"advance_xp3_simplified::xp3_old"),ngrdcol,nzt);
// Capture the grid component rather than its owning UDT in device kernels.
const auto invrs_dzt=gr.invrs_dzt;
//----------------------- Begin Code -----------------------
switch(solve_type) {
 case xp3_rtp3:
  // Budget stats for rtp3.
  name_bt="rtp3_bt";name_tp="rtp3_tp";name_ac="rtp3_ac";name_dp="rtp3_dp";
  break;
 case xp3_thlp3:
  // Budget stats for thlp3.
  name_bt="thlp3_bt";name_tp="thlp3_tp";name_ac="thlp3_ac";name_dp="thlp3_dp";
  break;
 default:
  // Budgets aren't setup for the passive scalars.
  name_bt="";name_tp="";name_ac="";name_dp="";
  break;
}
if(l_predict_xp3) {
 if(stats.flags().l_sample) {
  Real2D stats_bt_before(Kokkos::view_alloc(Kokkos::WithoutInitializing,"advance_xp3_simplified::bt_before"),ngrdcol,nzt);
  Kokkos::parallel_for("advance_xp3_simplified::bt_before",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
   const int i=cell%ngrdcol,k=cell/ngrdcol;
   stats_bt_before(i,k)=xp3(i,k)/dt;
  });
  f2k::stats_begin_budget(name_bt,stats_bt_before,stats);
 }
}
// Initialize variables.
Kokkos::deep_copy(term_tp,zero);
Kokkos::deep_copy(term_ac,zero);
// Interpolate <x> to momentum levels.
xm_zm=zt2zm_api(nzm,nzt,ngrdcol,gr,xm,zero_threshold);
// Interpolate <x'^2> to thermodynamic levels. Positive definite quantity.
xp2_zt=zm2zt_api(nzm,nzt,ngrdcol,gr,xp2,x_tol*x_tol);

// Independent cells retain Native k=1:nzt-1, with no vertical recurrence.
Kokkos::parallel_for("advance_xp3_simplified::terms",ngrdcol*(nzt>0?nzt-1:0),KOKKOS_LAMBDA(int cell) {
 const int i=cell%ngrdcol,k=cell/ngrdcol;
 // Define the km1 index.
 const int kp1=k+1<nzt-1?k+1:nzt-1;
 // Calculate the <x'^3> turbulent production (tp) term.
 term_tp(i,k)=term_tp_rhs(xp2_zt(i,k),wpxp(i,kp1),wpxp(i,k),
  rho_ds_zm(i,kp1),rho_ds_zm(i,k),invrs_rho_ds_zt(i,k),invrs_dzt(i,k));
 // Calculate the <x'^3> accumulation (ac) term.
 term_ac(i,k)=term_ac_rhs(xm_zm(i,kp1),xm_zm(i,k),wpxp2(i,k),invrs_dzt(i,k));
 if(l_predict_xp3) {
  if(l_lmm_stepping)xp3_old(i,k)=xp3(i,k);
  // Advance <x'^3> one time step.
  xp3(i,k)=((xp3(i,k)/dt)+term_tp(i,k)+term_ac(i,k))
   /((one/dt)+(C_xp3_dissipation*invrs_tau_zt(i,k)));
  if(l_lmm_stepping)xp3(i,k)=one_half*(xp3_old(i,k)+xp3(i,k));
 } else {
  // Calculate <x'^3> using the steady-state approximation.
  const real_t reciprocal_tau=one/invrs_tau_zt(i,k);
  xp3(i,k)=(reciprocal_tau<tau_max_zt(i,k)?reciprocal_tau:tau_max_zt(i,k))
   *one/C_xp3_dissipation*(term_tp(i,k)+term_ac(i,k));
 }
});
// Set Upper Boundary Condition: retain Native storage level nzt in both
// grid directions, rather than substituting a grid boundary index.
Kokkos::parallel_for("advance_xp3_simplified::upper_boundary",ngrdcol,KOKKOS_LAMBDA(int i) {
 xp3(i,nzt-1)=zero;
});
if(stats.flags().l_sample) {
 f2k::stats_update(name_tp,term_tp,stats);
 f2k::stats_update(name_ac,term_ac,stats);
 Real2D stats_dp(Kokkos::view_alloc(Kokkos::WithoutInitializing,"advance_xp3_simplified::dissipation"),ngrdcol,nzt);
 Kokkos::parallel_for("advance_xp3_simplified::dissipation",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
  const int i=cell%ngrdcol,k=cell/ngrdcol;
  stats_dp(i,k)=-(C_xp3_dissipation*invrs_tau_zt(i,k))*xp3(i,k);
 });
 f2k::stats_update(name_dp,stats_dp,stats);
 if(l_predict_xp3) {
  Real2D stats_bt_after(Kokkos::view_alloc(Kokkos::WithoutInitializing,"advance_xp3_simplified::bt_after"),ngrdcol,nzt);
  Kokkos::parallel_for("advance_xp3_simplified::bt_after",ngrdcol*nzt,KOKKOS_LAMBDA(int cell) {
   const int i=cell%ngrdcol,k=cell/ngrdcol;
   stats_bt_after(i,k)=xp3(i,k)/dt;
  });
  f2k::stats_finalize_budget(name_bt,stats_bt_after,stats);
 }
}
'''

bodies=[advance,diagnose,simplified,
'''// The <x'^3> turbulent production term.
return +three*xp2_zt*invrs_rho_ds_zt*invrs_dzt
 *(rho_ds_zmp1*wpxpp1-rho_ds_zm*wpxp);
''',
'''// The <x'^3> accumulation term.
return -three*wpxp2*invrs_dzt*(xm_zmp1-xm_zm);
''']
file_description='''// Predict and diagnose third-order scalar and horizontal-velocity moments.
//
// advance_xp3 advances rt/thl/passive-scalar third moments with the Native
// simplified equation, then diagnoses u/v third moments. diagnose_xp3 applies
// the Larson–Golaz ansatz to all those fields with the selected PDF width.
// Inputs/outputs are execution-space LayoutLeft Views indexed (column,level)
// or (column,level,scalar). Native units are listed in the interfaces.
// Private production/accumulation expressions execute within independent cell
// kernels. Stats retain the existing Fortran backend; no state is kept here.
// Native predictive/LMM and steady-state branches are retained. The last
// storage level is the Native upper boundary in either grid direction.
// Reference: Larson and Golaz (2005) for the diagnostic ansatz.
'''
head=file_description+'''#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "f2k/stats_bridge.hpp"
namespace clubb_kokkos {
'''
for name in names[:2]:head+=signature(name)+';\n\n'
head+='} // namespace clubb_kokkos\n'
cpp=file_description+'''#include "advance_xp3_module.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "model_flags.hpp"
#include "Skx_module.hpp"
#include "f2k/procedure_bridge.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <string>
namespace clubb_kokkos {
namespace {
constexpr int xp3_rtp3=1; // Named constant for solving rtp3.
constexpr int xp3_thlp3=2; // Named constant for solving thlp3.
constexpr int xp3_sclrp3=3; // Named constant for solving sclrp3.
'''
for name in names[2:]:cpp+=signature(name,True)+';\n'
cpp+='} // namespace\n\n'
for name,body in zip(names,bodies):
 if name==names[2]:cpp+='namespace {\n'
 cpp+=signature(name)+' {\n'+description(name)
 if not name.startswith('term_'):cpp+=locals(name)
 cpp+=body+'}\n\n'
cpp+='} // namespace\n} // namespace clubb_kokkos\n'
(folder/'advance_xp3_module.hpp').write_text(head)
(folder/'advance_xp3_module.cpp').write_text(cpp)
