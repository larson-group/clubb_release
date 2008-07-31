!-----------------------------------------------------------------------
! $Id: stats_zt.F90,v 1.7 2008-07-31 16:10:44 faschinj Exp $
 
module stats_zt

implicit none

private ! Default Scope

public :: stats_init_zt

contains

!-----------------------------------------------------------------------
subroutine stats_init_zt( vars_zt, lerror )

!     Description:
!     Initializes array indices for zt 
!-----------------------------------------------------------------------

use stats_variables, only: & 
    ithlm,  & ! Variable(s)
    iT_in_K, & 
    ithvm, & 
    irtm, & 
    ircm, & 
    ium, & 
    ivm, & 
    iwmt, & 
    iug, & 
    ivg, & 
    icf, & 
    ip, & 
    iexner, & 
    iLscale, & 
    iwp3, & 
    iwpthlp2, & 
    iwp2thlp, & 
    iwprtp2, & 
    iwp2rtp, & 
    iLscale_up, & 
    iLscale_down, & 
    itaut, & 
    iKht, & 
    iwp2thvp, & 
    iwp2rcp, & 
    iwprtpthlp, & 
    isct, & 
    irho, & 
    iNcm, & 
    iNcnm, & 
    iNim, & 
    isnowslope, & 
    iNsnowm, & 
    ised_rcm, & 
    irsat, & 
    irrainm, & 
    iNrm, & 
    irain_rate, & 
    iAKm, & 
    iAKm_est, & 
    iradht, & 
    iradht_LW, & 
    iradht_SW, & 
    idiam, & 
    imass_ice_cryst, & 
    ircm_icedfs, & 
    iu_T_cm, & 
    imean_vol_rad_rain, & 
    imean_vol_rad_cloud, & 
    irsnowm, & 
    irgraupelm, & 
    iricem, & 
    irtm_bt, & 
    irtm_ma, & 
    irtm_ta, & 
    irtm_forcing, & 
    irtm_mc, & 
    irtm_cl, & 
    irtm_pd, & 
    ithlm_bt, & 
    ithlm_ma, & 
    ithlm_ta, & 
    ithlm_forcing, & 
    ithlm_mc, & 
    ithlm_cl,       & 
    iwp3_bt, & 
    iwp3_ma, & 
    iwp3_ta, & 
    iwp3_tp, & 
    iwp3_ac, & 
    iwp3_bp, & 
    iwp3_pr1, & 
    iwp3_pr2, & 
    iwp3_dp1, & 
    iwp3_cl, & 
    irrainm_bt, & 
    irrainm_ma, & 
    irrainm_sd, & 
    irrainm_dff, & 
    irrainm_cond, & 
    irrainm_auto, & 
    irrainm_accr, & 
    irrainm_cond_adj, & 
    irrainm_mc, & 
    irrainm_cl, & 
    iNrm_bt, & 
    iNrm_ma, & 
    iNrm_sd, & 
    iNrm_dff, & 
    iNrm_cond, & 
    iNrm_auto, & 
    iNrm_cond_adj, & 
    iNrm_mc, & 
    iNrm_cl, & 
    irsnowm_bt, & 
    irsnowm_ma, & 
    irsnowm_sd, & 
    irsnowm_dff

use stats_variables, only: & 
    irsnowm_mc, & 
    irsnowm_cl, & 
    irgraupelm_bt, & 
    irgraupelm_ma, & 
    irgraupelm_sd, & 
    irgraupelm_dff, & 
    irgraupelm_mc, & 
    irgraupelm_cl, & 
    iricem_bt, & 
    iricem_ma, & 
    iricem_sd, & 
    iricem_dff, & 
    iricem_mc, & 
    iricem_cl, & 
    ivm_bt, & 
    ivm_ma, & 
    ivm_gf, & 
    ivm_cf, & 
    ivm_ta, & 
    ium_bt, & 
    ium_ma, & 
    ium_gf, & 
    ium_cf, & 
    ium_ta, & 
    ia, & 
    iw1, & 
    iw2, & 
    isw1, & 
    isw2, & 
    ithl1, & 
    ithl2, & 
    isthl1, & 
    isthl2, & 
    irt1, & 
    irt2, & 
    isrt1, & 
    isrt2, & 
    irc1, & 
    irc2, & 
    irsl1, & 
    irsl2, & 
    iR1, & 
    iR2, & 
    is1, & 
    is2, & 
    iss1, & 
    iss2, & 
    irrtthl, & 
    iwp2_zt, & 
    ithlp2_zt, & 
    iwpthlp_zt, & 
    iwprtp_zt, & 
    irtp2_zt, & 
    irtpthlp_zt, & 
    zt, & 
    isclram, & 
    isclram_f, & 
    isclrbm, & 
    isclrbm_f, & 
    iedsclram, & 
    iedsclrbm

use stats_type, only: & 
    stat_assign ! Procedure

implicit none

integer, parameter :: nvarmax = 250

!Input Variable
character(len= * ), dimension(nvarmax), intent(in) :: vars_zt

!Output Variable	
logical, intent(inout) :: lerror

!Local Varables
integer :: i, k

! Default initialization for array indices for zt

ithlm         = 0
iT_in_K    = 0
ithvm         = 0
irtm          = 0
ircm          = 0
ium           = 0
ivm           = 0
iwmt          = 0
iug           = 0
ivg           = 0
icf           = 0
ip            = 0
iexner        = 0
iLscale       = 0
iwp3          = 0
iwpthlp2      = 0
iwp2thlp      = 0
iwprtp2       = 0
iwp2rtp       = 0
iLscale_up          = 0
iLscale_down        = 0
itaut         = 0
iKht          = 0
iwp2thvp      = 0
iwp2rcp       = 0
iwprtpthlp    = 0
isct          = 0
irho         = 0
iNcm          = 0  ! Brian
iNcnm         = 0
iNim          = 0
isnowslope    = 0  ! Adam Smith, 22 April 2008
iNsnowm       = 0  ! Adam Smith, 22 April 2008
ised_rcm      = 0  ! Brian
irsat          = 0  ! Brian
irrainm          = 0  ! Brian
iNrm          = 0  ! Brian
irain_rate    = 0  ! Brian
iAKm          = 0  ! analytic Kessler.  Vince Larson 22 May 2005
iAKm_est      = 0  ! LH Kessler.  Vince Larson 22 May 2005
iradht        = 0
iradht_LW     = 0
iradht_SW     = 0

idiam           = 0
imass_ice_cryst = 0
ircm_icedfs     = 0
iu_T_cm         = 0

imean_vol_rad_rain  = 0  ! Brian
imean_vol_rad_cloud = 0

irsnowm       = 0
irgraupelm    = 0
iricem        = 0

irtm_bt       = 0
irtm_ma       = 0
irtm_ta       = 0
irtm_forcing  = 0
irtm_mc       = 0
irtm_cl       = 0 ! Josh
irtm_pd       = 0
ithlm_bt      = 0
ithlm_ma      = 0
ithlm_ta      = 0
ithlm_forcing = 0
ithlm_mc      = 0
ithlm_cl      = 0 ! Josh

iwp3_bt       = 0
iwp3_ma       = 0
iwp3_ta       = 0
iwp3_tp       = 0
iwp3_ac       = 0
iwp3_bp       = 0
iwp3_pr1      = 0
iwp3_pr2      = 0
iwp3_dp1      = 0
iwp3_cl       = 0

irrainm_bt       = 0
irrainm_ma       = 0
irrainm_sd       = 0
irrainm_dff      = 0
irrainm_cond     = 0
irrainm_auto     = 0
irrainm_accr     = 0
irrainm_cond_adj = 0
irrainm_mc       = 0
irrainm_cl       = 0

iNrm_bt       = 0
iNrm_ma       = 0
iNrm_sd       = 0
iNrm_dff      = 0
iNrm_cond     = 0
iNrm_auto     = 0
iNrm_cond_adj = 0
iNrm_mc       = 0
iNrm_cl       = 0

irsnowm_bt    = 0
irsnowm_ma    = 0
irsnowm_sd    = 0
irsnowm_dff   = 0
irsnowm_mc    = 0
irsnowm_cl    = 0

irgraupelm_bt = 0
irgraupelm_ma = 0
irgraupelm_sd = 0
irgraupelm_dff= 0
irgraupelm_mc = 0
irgraupelm_cl = 0

iricem_bt     = 0
iricem_ma     = 0
iricem_sd     = 0
iricem_dff    = 0
iricem_mc     = 0
iricem_cl     = 0

ivm_bt = 0
ivm_ma = 0
ivm_gf = 0
ivm_cf = 0
ivm_ta = 0

ium_bt = 0
ium_ma = 0
ium_gf = 0
ium_cf = 0
ium_ta = 0

ia            = 0
iw1           = 0
iw2           = 0
isw1          = 0
isw2          = 0
ithl1         = 0
ithl2         = 0
isthl1        = 0
isthl2        = 0
irt1          = 0
irt2          = 0
isrt1         = 0
isrt2         = 0
irc1          = 0
irc2          = 0
irsl1         = 0
irsl2         = 0
iR1           = 0
iR2           = 0
is1           = 0
is2           = 0
iss1          = 0
iss2          = 0
irrtthl       = 0

iwp2_zt     = 0
ithlp2_zt   = 0
iwpthlp_zt  = 0
iwprtp_zt   = 0
irtp2_zt    = 0
irtpthlp_zt = 0

isclram     = 0
isclram_f   = 0
isclrbm     = 0
isclrbm_f   = 0

iedsclram   = 0
iedsclrbm   = 0

lerror = .false.

!     Assign pointers for statistics variables zt

k = 1
do i=1,zt%nn

  select case ( trim(vars_zt(i)) )
  case ('thlm')
    ithlm = k
    call stat_assign( ithlm, "thlm",  & 
          "thetal (K)", "K", zt)
    k = k + 1

  case ('T_in_K')
    iT_in_K = k
    call stat_assign( iT_in_K, "T_in_K",  & 
          "absolute temperature in K", "K", zt )
    k = k + 1

  case ('thvm')
    ithvm = k
    call stat_assign( ithvm, "thvm", & 
          "virtual potential temperature (K)","K",zt)
    k = k + 1

  case ('rtm')
    irtm = k

    call stat_assign( irtm, "rtm", & 
          "total water mixing ratio (kg/kg)","kg/kg",zt)
    
    !zt%f%var(irtm)%ptr => zt%x(:,k)
    !zt%f%var(irtm)%name = "rtm"
    !zt%f%var(irtm)%description 
    != "total water mixing ratio (kg/kg)"
    !zt%f%var(irtm)%units = "kg/kg"
    
    k = k + 1

  case ('rcm')
    ircm = k
    call stat_assign( ircm, "rcm", & 
          "liquid water mixing ratio (kg/kg)","kg/kg",zt)
    k = k + 1
  case ('um')
    ium = k
    call stat_assign( ium, "um", & 
          "u wind (m/s)","m/s",zt)
    k = k + 1
  case ('vm')
    ivm = k
    call stat_assign( ivm,  "vm", & 
          "v wind (m/s)","m/s",zt)
    k = k + 1
  case ('wmt')
    iwmt = k
    call stat_assign( iwmt, "wm", & 
          "w wind (m/s)","m/s",zt)
    k = k + 1
  case ('ug')
    iug = k
    call stat_assign( iug,"ug", & 
         "u geostrophic wind (m/s)", "m/s", zt)
    k = k + 1
  case ('vg')
    ivg = k
    call stat_assign(ivg,"vg", & 
         "v geostrophic wind (m/s)", "m/s",zt)
    k = k + 1
  case ('cf')
    icf = k
    call stat_assign(icf,"cf", & 
         "cloud fraction", "count",zt)
    k = k + 1
  case ('p')
    ip = k
    call stat_assign(ip,"p", & 
         "pressure (Pa)","Pa",zt)
    k = k + 1
  case ('exner')
    iexner = k
    call stat_assign(iexner,"exner", & 
         "Exner","count",zt)
    k = k + 1
  case ('Lscale')
    iLscale = k
    call stat_assign(iLscale,"Lscale", & 
         "Mixing length","m",zt)
    k = k + 1
  case ('thlm_forcing')
    ithlm_forcing = k
    call stat_assign(ithlm_forcing,"thlm_f", & 
         "thetal forcing", "K/s",zt)
    k = k + 1
  case ('thlm_mc')
    ithlm_mc = k
    call stat_assign(ithlm_mc,"thlm_mc", & 
         "thetal micro (not in budget)", "K/s",zt)
    k = k + 1
  case ('rtm_forcing')
    irtm_forcing = k
    call stat_assign(irtm_forcing,"rtm_f", & 
         "rt forcing", "kg/(kg s)",zt)
    k = k + 1
  case ('rtm_mc')
    irtm_mc = k
    call stat_assign(irtm_mc,"rtm_mc", & 
         "rt micro (not in budget)", "kg/(kg s)",zt)
    k = k + 1

  case ('wp3')
    iwp3 = k
    call stat_assign(iwp3,"wp3", & 
         "w third order moment", "(m^3)/(s^3)",zt)
    k = k + 1

  case ('wpthlp2')
    iwpthlp2 = k
    call stat_assign(iwpthlp2,"wpthlp2", & 
         "wpthlp2 covariance", "(m K^2)/s",zt)
    k = k + 1

  case ('wp2thlp')
    iwp2thlp = k
    call stat_assign(iwp2thlp,"wp2thlp", & 
         "wp2thlp covariance", "(m^2 K)/(s^2)",zt)
    k = k + 1

  case ('wprtp2')
    iwprtp2 = k
    call stat_assign(iwprtp2,"wprtp2", & 
         "wprtp2 covariance", "(m kg^2)/(s kg^2)",zt)
    k = k + 1

  case ('wp2rtp')
    iwp2rtp = k
    call stat_assign(iwp2rtp,"wp2rtp", & 
         "wp2rtp covariance", "(m2 kg)/(s2 kg)",zt)
    k = k + 1

  case ('Lscale_up')
    iLscale_up = k
    call stat_assign(iLscale_up,"Lscale_up", & 
         "Upward mixing length","m",zt)
    k = k + 1

  case ('Lscale_down')
    iLscale_down = k
    call stat_assign(iLscale_down,"Lscale_down", & 
         "Downward mixing length","m",zt)
    k = k + 1

  case ('taut')
    itaut = k
    call stat_assign(itaut,"taut", & 
         "Dissipation time","s",zt)
    k = k + 1

  case ('kht')
    iKht = k
    call stat_assign(iKht,"Kht", & 
         "Eddy diffusivity","m^2/s",zt)
    k = k + 1

  case ('wp2thvp')
    iwp2thvp = k
    call stat_assign(iwp2thvp,"wp2thvp", & 
         "wp2thvp","(m2 K)/s2",zt)
    k = k + 1

  case ('wp2rcp')
    iwp2rcp = k
    call stat_assign(iwp2rcp,"wp2rcp", & 
         "wp2rcp","(m2 kg)/(s2 kg)",zt)
    k = k + 1

  case ('wprtpthlp')
    iwprtpthlp = k
    call stat_assign(iwprtpthlp,"wprtpthlp", & 
         "wprtpthlp","(m kg K)/(s kg)",zt)
    k = k + 1

  case ('sc')
    isct = k
    call stat_assign(isct,"Sc", & 
         "Sc","count",zt)
    k = k + 1

  case ('rho')
    irho = k
    call stat_assign(irho,"rho", & 
         "density","kg/m^3",zt)
    k = k + 1

  case ('Ncm')           ! Brian
    iNcm = k
    call stat_assign(iNcm,"Ncm", & 
         "Cloud droplet number concentration (num/kg)", & 
         "count/kg",zt)
    k = k + 1

  case ('Ncnm')
    iNcnm = k
    call stat_assign(iNcnm,"Ncnm", & 
         "Cloud nuclei number concentration (num/m^3)", & 
         "count/m^3",zt)
    k = k + 1

  case ('Nim')           ! Brian
    iNim = k
    call stat_assign(iNim,"Nim", & 
         "Ice crystal number concentration (num/m^3)", & 
         "count/m^3",zt)
    k = k + 1

  case ('snowslope')     ! Adam Smith, 22 April 2008
    isnowslope = k
    call stat_assign(isnowslope,"snowslope", & 
         "COAMPS snow slope parameter (1/m)", & 
         "1/m",zt)
    k = k + 1

  case ('Nsnowm')        ! Adam Smith, 22 April 2008
    iNsnowm = k
    call stat_assign(iNsnowm,"Nsnowm", & 
         "Snow particle number concentration (num/m^3)", & 
         "count/m^3",zt)
    k = k + 1

  case ('sed_rcm')       ! Brian
    ised_rcm = k
    call stat_assign(ised_rcm,"sed_rcm", & 
         "d(rcm)/dt due to cloud sedimentation (kg / [m^2 s])", & 
         "kg/(m^2 s)",zt)
    k = k + 1

  case ('rsat')           ! Brian
    irsat = k
    call stat_assign(irsat,"rsat", & 
         "Saturation mixing ratio (kg/kg)","kg/kg",zt)
    k = k + 1

  case ('rrainm')           ! Brian
    irrainm = k
    call stat_assign(irrainm,"rrainm", & 
         "Rain water mixing ratio (kg/kg)","kg/kg",zt)
    k = k + 1

  case ('rsnowm')
    irsnowm = k
    call stat_assign(irsnowm,"rsnowm", & 
         "Snow water mixing ratio (kg/kg)","kg/kg",zt)
    k = k + 1

  case ('ricem')
    iricem = k
    call stat_assign(iricem,"ricem", & 
         "Pristine ice water mixing ratio (kg/kg)","kg/kg",zt)
    k = k + 1

  case ('rgraupelm')
    irgraupelm = k
    call stat_assign(irgraupelm,"rgraupelm", & 
         "Graupel water mixing ratio (kg/kg)","kg/kg",zt)
    k = k + 1

  case ('Nrm')           ! Brian
    iNrm = k
    call stat_assign(iNrm,"Nrm", & 
         "Rain drop number concentration (num/kg)", & 
         "count/kg",zt)
    k = k + 1

  case ('mean_vol_rad_rain')  ! Brian
    imean_vol_rad_rain = k
    call stat_assign(imean_vol_rad_rain,"mvrr", & 
         "Rain drop mean volume radius (m)","m",zt)
    k = k + 1

  case ('mean_vol_rad_cloud')
    imean_vol_rad_cloud = k

    call stat_assign(imean_vol_rad_cloud,"mvrc", & 
         "Cloud drop mean volume radius (m)","m",zt)
    k = k + 1

  case ('rain_rate')     ! Brian
    irain_rate = k

    call stat_assign(irain_rate,"rain_rate", & 
         "Rain rate (mm/day)","mm/day",zt)
    k = k + 1
 
  case ('AKm')           ! Vince Larson 22 May 2005
    iAKm = k
    call stat_assign(iAKm,"AKm", & 
         "Analytic Kessler ac [kg/kg]","kg/kg",zt)
    k = k + 1
 
  case ('AKm_est')       ! Vince Larson 22 May 2005
    iAKm_est = k

    call stat_assign(iAKm_est,"AKm_est", & 
         "LH Kessler estimate [kg/kg]","kg/kg",zt)
    k = k + 1

  case ('radht')
    iradht = k

    call stat_assign(iradht,"radht", & 
         "Heating rate","K/s",zt)
    k = k + 1

  case ('radht_LW')
    iradht_LW = k

    call stat_assign(iradht_LW,"radht_LW", & 
         "radht_LW, Long-wave heating rate","K/s",zt)

    k = k + 1

  case ('radht_SW')
    iradht_SW = k
    call stat_assign(iradht_SW,"radht_SW", & 
         "radht_SW, Short-wave heating rate","K/s",zt)
    k = k + 1

  case ('diam')
    idiam = k

    call stat_assign(idiam,"diam", & 
         "Ice crystal diameter","m",zt)
    k = k + 1

  case ('mass_ice_cryst')
    imass_ice_cryst = k
    call stat_assign(imass_ice_cryst,"mass_ice_cryst", & 
         "Ice crystal mass","kg",zt)
    k = k + 1

  case ('rcm_icedfs')

    ircm_icedfs = k
    call stat_assign(ircm_icedfs,"rcm_icedfs", & 
         "Change in liquid due to ice","kg/kg/s",zt)
    k = k + 1

  case ('u_T_cm')
    iu_T_cm = k
    call stat_assign(iu_T_cm,"u_T_cm", & 
         "Ice crystal fallspeed","cm/s",zt)
    k = k + 1

  case ('rtm_bt')
    irtm_bt = k

    call stat_assign(irtm_bt,"rtm_bt", & 
         "rtm budget","kg/kg/s",zt)
    k = k + 1

  case ('rtm_ma')
    irtm_ma = k

    call stat_assign(irtm_ma,"rtm_ma", & 
         "rtm ma","kg/kg/s",zt)
    k = k + 1

  case ('rtm_ta')
    irtm_ta = k

    call stat_assign(irtm_ta,"rtm_ta", & 
         "rtm ta","kg/kg/s",zt)
    k = k + 1
  case ('rtm_cl')
    irtm_cl = k

    call stat_assign(irtm_cl, "rtm_cl", & 
         "rtm clipping", "kg/kg/s", zt)

    k = k + 1
  case ('rtm_pd')
    irtm_pd = k

    call stat_assign(irtm_pd, "rtm_pd", & 
         "rtm positive definite adj.", "kg/kg/s", zt)

    k = k + 1

  case ('thlm_bt')
    ithlm_bt = k

    call stat_assign(ithlm_bt,"thlm_bt", & 
         "thlm bt","kg/kg/s",zt)
    k = k + 1

  case ('thlm_ma')
    ithlm_ma = k

    call stat_assign(ithlm_ma,"thlm_ma", & 
         "thlm ma","kg/kg/s",zt)
    k = k + 1

  case ('thlm_ta')
    ithlm_ta = k

    call stat_assign(ithlm_ta,"thlm_ta", & 
         "thlm ta","kg/kg/s",zt)
    k = k + 1
  case ('thlm_cl')
    ithlm_cl = k

    call stat_assign(ithlm_cl,"thlm_cl", & 
          "thlm_cl", "kg/kg/s",zt)
    k = k + 1 
  case ('wp3_bt')
    iwp3_bt = k

    call stat_assign(iwp3_bt,"wp3_bt", & 
         "wp3 budget","(m^3)/(s^4)",zt)
    k = k + 1
 
  case ('wp3_ma')
    iwp3_ma = k

    call stat_assign(iwp3_ma,"wp3_ma", & 
         "wp3 mean advection","(m^3)/(s^4)",zt)
    k = k + 1
 
  case ('wp3_ta')
    iwp3_ta = k

    call stat_assign(iwp3_ta,"wp3_ta", & 
         "wp3 turbulent advection","(m^3)/(s^4)",zt)

    k = k + 1
 
  case ('wp3_tp')
    iwp3_tp = k
    call stat_assign(iwp3_tp,"wp3_tp", & 
         "wp3 turbulent transport","(m^3)/(s^4)",zt)
    k = k + 1
 
  case ('wp3_ac')
    iwp3_ac = k
    call stat_assign(iwp3_ac,"wp3_ac", & 
         "wp3 accumulation term","(m^3)/(s^4)",zt)
    k = k + 1
 
  case ('wp3_bp')
    iwp3_bp = k
    call stat_assign(iwp3_bp,"wp3_bp", & 
         "wp3 buoyancy production","(m^3)/(s^4)",zt)
    k = k + 1
 
  case ('wp3_pr1')
    iwp3_pr1 = k
    call stat_assign(iwp3_pr1,"wp3_pr1", & 
         "wp3 pressure term 1","(m^3)/(s^4)",zt)
    k = k + 1
 
  case ('wp3_pr2')
    iwp3_pr2 = k
    call stat_assign(iwp3_pr2,"wp3_pr2", & 
         "wp3 pressure term 2","(m^3)/(s^4)",zt)

    k = k + 1

  case ('wp3_dp1')
    iwp3_dp1 = k
    call stat_assign(iwp3_dp1,"wp3_dp1", & 
         "wp3 dissipation term 1","(m^3)/(s^4)",zt)
    k = k + 1
 
  case ('wp3_cl')
    iwp3_cl = k
    call stat_assign(iwp3_cl,"wp3_cl", & 
         "wp3 clipping term","(m^3)/(s^4)",zt)
    k = k + 1
 
  case ('rrainm_bt')
    irrainm_bt = k
    call stat_assign(irrainm_bt,"rrainm_bt", & 
         "rrainm budget","(kg/kg)/(s)",zt)
    k = k + 1
 
  case ('rrainm_ma')
    irrainm_ma = k

    call stat_assign(irrainm_ma,"rrainm_ma", & 
         "rrainm mean advection","(kg/kg)/(s)",zt)
    k = k + 1
 
  case ('rrainm_sd')
    irrainm_sd = k

    call stat_assign(irrainm_sd,"rrainm_sd", & 
         "rrainm sedimentation","(kg/kg)/(s)",zt)
    k = k + 1
 
  case ('rrainm_dff')
    irrainm_dff = k

    call stat_assign(irrainm_dff,"rrainm_dff", & 
         "rrainm diffusion","(kg/kg)/(s)",zt)
    k = k + 1
 
  case ('rrainm_cond')
    irrainm_cond = k

    call stat_assign(irrainm_cond,"rrainm_cond", & 
         "rrainm cond/evap","(kg/kg)/(s)",zt)
    k = k + 1
 
  case ('rrainm_auto')
    irrainm_auto = k

    call stat_assign(irrainm_auto,"rrainm_auto", & 
         "rrainm autoconversion","(kg/kg)/(s)",zt)
    k = k + 1
 
  case ('rrainm_accr')
    irrainm_accr = k
    call stat_assign(irrainm_accr,"rrainm_accr", & 
         "rrainm accretion","(kg/kg)/(s)",zt)
    k = k + 1

  case ('rrainm_cond_adj')
    irrainm_cond_adj = k

    call stat_assign(irrainm_cond_adj,"rrainm_cond_adj", & 
         "rrainm cond/evap adjustment due to over-evaporation", & 
         "(kg/kg)/(s)",zt)
    k = k + 1
 
  case ('rrainm_mc')
    irrainm_mc = k

    call stat_assign(irrainm_mc,"rrainm_mc", & 
         "rrainm total microphysical tendency","(kg/kg)/(s)",zt)

    k = k + 1

  case ('rrainm_cl')
    irrainm_cl = k
    call stat_assign(irrainm_cl,"rrainm_cl", & 
         "rrainm clipping term","(kg/kg)/(s)",zt)

    k = k + 1
 
  case ('Nrm_bt')
    iNrm_bt = k
    call stat_assign(iNrm_bt,"Nrm_bt", & 
         "Nrm budget","(count/kg)/s",zt)

    k = k + 1
 
  case ('Nrm_ma')
    iNrm_ma = k

    call stat_assign(iNrm_ma,"Nrm_ma", & 
         "Nrm mean advection","(count/kg)/s",zt)
    k = k + 1
 
  case ('Nrm_sd')
    iNrm_sd = k

    call stat_assign(iNrm_sd,"Nrm_sd", & 
         "Nrm sedimentation","(count/kg)/s",zt)

    k = k + 1
 
  case ('Nrm_dff')
    iNrm_dff = k
    call stat_assign(iNrm_dff,"Nrm_dff", & 
         "Nrm diffusion","(count/kg)/s",zt)

    k = k + 1
 
  case ('Nrm_cond')
    iNrm_cond = k

    call stat_assign(iNrm_cond,"Nrm_cond", & 
         "Nrm cond/evap","(count/kg)/s",zt)
    k = k + 1
 
  case ('Nrm_auto')
    iNrm_auto = k

    call stat_assign(iNrm_auto,"Nrm_auto", & 
         "Nrm autoconversion","(count/kg)/s",zt)

    k = k + 1

  case ('Nrm_cond_adj')
    iNrm_cond_adj = k

    call stat_assign(iNrm_cond_adj,"Nrm_cond_adj", & 
         "Nrm cond/evap adjustment due to over-evaporation", & 
         "(count/kg)/s",zt)
    k = k + 1
 
  case ('Nrm_mc')
    iNrm_mc = k
    call stat_assign(iNrm_mc,"Nrm_mc", & 
         "Nrm micro","(count/kg)/s",zt)

    k = k + 1

  case ('Nrm_cl')
    iNrm_cl = k
    call stat_assign(iNrm_cl,"Nrm_cl", & 
         "Nrm clipping term","(count/kg)/s",zt)
    k = k + 1

  case ('rsnowm_bt')
    irsnowm_bt = k
    call stat_assign(irsnowm_bt,"rsnowm_bt", & 
         "rsnowm budget","(kg/kg)/s",zt)

    k = k + 1
 
  case ('rsnowm_ma')
    irsnowm_ma = k

    call stat_assign(irsnowm_ma,"rsnowm_ma", & 
         "rsnowm mean advection","(kg/kg)/s",zt)
    k = k + 1
 
  case ('rsnowm_sd')
    irsnowm_sd = k
    call stat_assign(irsnowm_sd,"rsnowm_sd", & 
         "rsnowm sedimentation","(kg/kg)/s",zt)
    k = k + 1
 
  case ('rsnowm_dff')
    irsnowm_dff = k

    call stat_assign(irsnowm_dff,"rsnowm_dff", & 
         "rsnowm diffusion","(kg/kg)/s",zt)
    k = k + 1

  case ('rsnowm_mc')
    irsnowm_mc = k

    call stat_assign(irsnowm_mc,"rsnowm_mc", & 
         "rsnowm micro","(kg/kg)/s",zt)
    k = k + 1

  case ('rsnowm_cl')
    irsnowm_cl = k

    call stat_assign(irsnowm_cl,"rsnowm_cl", & 
         "rsnowm clipping term","(kg/kg)/s",zt)
    k = k + 1

 case ('ricem_bt')
    iricem_bt = k

    call stat_assign(iricem_bt,"ricem_bt", & 
         "ricem budget","(kg/kg)/s",zt)
 
    k = k + 1
 
  case ('ricem_ma')
    iricem_ma = k

    call stat_assign(iricem_ma,"ricem_ma", & 
         "ricem mean advection","(kg/kg)/s",zt)
    k = k + 1
 
  case ('ricem_sd')
    iricem_sd = k

    call stat_assign(iricem_sd,"ricem_sd", & 
         "ricem sedimentation","(kg/kg)/s",zt)
    k = k + 1
 
  case ('ricem_dff')
    iricem_dff = k

    call stat_assign(iricem_dff,"ricem_dff", & 
         "ricem diffusion","(kg/kg)/s",zt)
    k = k + 1

  case ('ricem_mc')
    iricem_mc = k

    call stat_assign(iricem_mc,"ricem_mc", & 
         "ricem micro","(kg/kg)/s",zt)
    k = k + 1

  case ('ricem_cl')
    iricem_cl = k

    call stat_assign(iricem_cl,"ricem_cl", & 
         "ricem clipping term","(kg/kg)/s",zt)
    k = k + 1

  case ('rgraupelm_bt')
    irgraupelm_bt = k

    call stat_assign(irgraupelm_bt,"rgraupelm_bt", & 
         "rgraupelm budget","(kg/kg)/s",zt)
    k = k + 1
 
  case ('rgraupelm_ma')
    irgraupelm_ma = k

    call stat_assign(irgraupelm_ma,"rgraupelm_ma", & 
         "rgraupelm mean advection","(kg/kg)/s",zt)
    k = k + 1
 
  case ('rgraupelm_sd')
    irgraupelm_sd = k

    call stat_assign(irgraupelm_sd,"rgraupelm_sd", & 
         "rgraupelm sedimentation","(kg/kg)/s",zt)
    k = k + 1
 
  case ('rgraupelm_dff')
    irgraupelm_dff = k

    call stat_assign(irgraupelm_dff,"rgraupelm_dff", & 
         "rgraupelm diffusion","(kg/kg)/s",zt)
    k = k + 1

  case ('rgraupelm_mc')
    irgraupelm_mc = k

    call stat_assign(irgraupelm_mc,"rgraupelm_mc", & 
         "rgraupelm micro","(kg/kg)/s",zt)
    k = k + 1

  case ('rgraupelm_cl')
    irgraupelm_cl = k

    call stat_assign(irgraupelm_cl,"rgraupelm_cl", & 
         "rgraupelm clipping term","(kg/kg)/s",zt)
    k = k + 1

  case ('vm_bt')
    ivm_bt = k

    call stat_assign(ivm_bt,"vm_bt", & 
         "vm time tendency(note: does not balance at zt=2)", & 
         "m/s",zt)
    k = k + 1

  case ('vm_ma')
    ivm_ma = k
    call stat_assign(ivm_ma,"vm_ma", & 
         "vm mean advection","m/s",zt)
    k = k + 1

  case ('vm_gf')
    ivm_gf = k

    call stat_assign(ivm_gf,"vm_gf", & 
         "vm geostrophic forcing","m/s",zt)
    k = k + 1

  case ('vm_cf')
    ivm_cf = k

    call stat_assign(ivm_cf,"vm_cf", & 
         "vm coriolis forcing","m/s",zt)
    k = k + 1

  case ('vm_ta')
    ivm_ta = k

    call stat_assign(ivm_ta,"vm_ta", & 
         "vm turbulent transport","m/s",zt)
    k = k + 1

  case ('um_bt')
    ium_bt = k

    call stat_assign(ium_bt,"um_bt", & 
         "um time tendency (note: does not balance at zt=2)", & 
         "m/s",zt)
    k = k + 1

  case ('um_ma')
    ium_ma = k

    call stat_assign(ium_ma,"um_ma", & 
         "um mean advection","m/s",zt)
    k = k + 1

  case ('um_gf')
    ium_gf = k
    call stat_assign(ium_gf,"um_gf", & 
         "um geostrophic forcing","m/s",zt)
    k = k + 1

  case ('um_cf')
    ium_cf = k
    call stat_assign(ium_cf,"um_cf", & 
         "um coriolis forcing","m/s",zt)
    k = k + 1

  case ('um_ta')
    ium_ta = k
    call stat_assign(ium_ta,"um_ta", & 
         "um turbulent transport","m/s",zt)
    k = k + 1
 
  case ('a')
    ia = k
    call stat_assign(ia,"a", & 
         "pdf parameter a","count",zt)
    k = k + 1
 
  case ('w1')
    iw1 = k
    call stat_assign(iw1,"w1", & 
         "pdf parameter w1","m/s",zt)

    k = k + 1
 
  case ('w2')
    iw2 = k

    call stat_assign(iw2,"w2", & 
         "pdf parameter w2","m/s",zt)
    k = k + 1
 
  case ('sw1')
    isw1 = k
    call stat_assign(isw1,"sw1", & 
         "pdf parameter sw1","m^2/s^2",zt)

    k = k + 1
 
  case ('sw2')
    isw2 = k

    call stat_assign(isw2,"sw2", & 
         "pdf parameter sw2","m^2/s^2",zt)
    k = k + 1
 
  case ('thl1')
    ithl1 = k

    call stat_assign(ithl1,"thl1", & 
         "pdf parameter thl1","K",zt)

    k = k + 1
 
  case ('thl2')
    ithl2 = k

    call stat_assign(ithl2,"thl2", & 
         "pdf parameter thl2","K",zt)
    k = k + 1
 
  case ('sthl1')
    isthl1 = k

    call stat_assign(isthl1,"sthl1", & 
         "pdf parameter sthl1","K^2",zt)

    k = k + 1
 
  case ('sthl2')
    isthl2 = k
    call stat_assign(isthl2,"sthl2", & 
         "pdf parameter sthl2","K^2",zt)

    k = k + 1
 
  case ('rt1')
    irt1 = k
    call stat_assign(irt1,"rt1", & 
         "pdf parameter rt1","kg/kg",zt)

    k = k + 1
 
  case ('rt2')
    irt2 = k

    call stat_assign(irt2,"rt2", & 
         "pdf parameter rt2","kg/kg",zt)
    k = k + 1
 
  case ('srt1')
    isrt1 = k
    call stat_assign(isrt1,"srt1", & 
         "pdf parameter srt1","(kg^2)/(kg^2)",zt)
    k = k + 1
 
  case ('srt2')
    isrt2 = k

    call stat_assign(isrt2,"srt2", & 
         "pdf parameter srt2","(kg^2)/(kg^2)",zt)
    k = k + 1
 
  case ('rc1')
    irc1 = k

    call stat_assign(irc1,"rc1", & 
         "pdf parameter rc1","kg/kg",zt)
    k = k + 1
 
  case ('rc2')
    irc2 = k

    call stat_assign(irc2,"rc2", & 
         "pdf parameter rc2","kg/kg",zt)
    k = k + 1
 
  case ('rsl1')
    irsl1 = k

    call stat_assign(irsl1,"rsl1", & 
         "pdf parameter rsl1","kg/kg",zt)
    k = k + 1
 
  case ('rsl2')
    irsl2 = k

    call stat_assign(irsl2,"rsl2", & 
         "pdf parameter rsl2","kg/kg",zt)
    k = k + 1
 
  case ('R1')
    iR1 = k
    call stat_assign(iR1,"R1", & 
         "pdf parameter R1","count",zt)
    k = k + 1
 
  case ('R2')
    iR2 = k

    call stat_assign(iR2,"R2", & 
         "pdf parameter R2","count",zt)
    k = k + 1
 
  case ('s1')
    is1 = k

    call stat_assign(is1,"s1", & 
         "pdf parameter s1","kg/kg",zt)
    k = k + 1
 
  case ('s2')
    is2 = k

    call stat_assign(is2,"s2", & 
         "pdf parameter s2","kg/kg",zt)
    k = k + 1
 
  case ('ss1')
    iss1 = k

    call stat_assign(iss1,"ss1", & 
         "pdf parameter ss1","kg/kg",zt)
    k = k + 1
 
  case ('ss2')
    iss2 = k

    call stat_assign(iss2,"ss2", & 
         "pdf parameter ss2","kg/kg",zt)
    k = k + 1
 
  case ('rrtthl')
    irrtthl = k

    call stat_assign(irrtthl,"rrtthl", & 
         "pdf parameter rrtthl","count",zt)
    k = k + 1

  case('wp2_zt')
    iwp2_zt = k

    call stat_assign(iwp2_zt,"wp2_zt", & 
         "wp2_zt","m^2/s^2",zt)
    k = k + 1

  case('thlp2_zt')
    ithlp2_zt = k

    call stat_assign(ithlp2_zt,"thlp2_zt", & 
         "thlp2_zt","K^2",zt)
    k = k + 1

  case('wpthlp_zt')
    iwpthlp_zt = k

    call stat_assign(iwpthlp_zt,"wpthlp_zt", & 
         "wpthlp_zt","(m K)/s",zt)
    k = k + 1

  case('wprtp_zt')   
    iwprtp_zt = k

    call stat_assign(iwprtp_zt,"wprtp_zt", & 
         "wprtp_zt","(m kg)/(s kg)",zt)
    k = k + 1

  case('rtp2_zt')   
    irtp2_zt = k

    call stat_assign(irtp2_zt,"rtp2_zt", & 
         "rtp2_zt","kg/kg",zt)
    k = k + 1

  case('rtpthlp_zt')   
    irtpthlp_zt = k

    call stat_assign(irtpthlp_zt,"rtpthlp_zt", & 
         "rtpthlp_zt","(kg K)/kg",zt)
    k = k + 1

  case ('sclram')
    isclram = k

    call stat_assign(isclram,"sclram", & 
         "passive scalar (currently thlm)","K",zt)
    k = k + 1

  case ('sclram_f')
    isclram_f = k

    call stat_assign(isclram_f,"sclram_f", & 
         "passive scalar (currently thlm_forcing)","K/s",zt)
    k = k + 1

  case ('sclrbm')
    isclrbm = k

    call stat_assign(isclrbm,"sclrbm", & 
         "passive scalar (currently rtm)","kg/kg",zt)
    k = k + 1

  case ('sclrbm_f')
    isclrbm_f = k
    call stat_assign(isclrbm_f,"sclrbm_f", & 
         "passive scalar (currently rtm_forcing)","kg/kg",zt)
    k = k + 1

  case ('edsclram')
    iedsclram = k

    call stat_assign(iedsclram,"edsclram", & 
         "eddy-diff scalar (currently thlm)","K",zt)
    k = k + 1

  case ('edsclrbm')
    iedsclrbm = k

    call stat_assign(iedsclrbm,"edsclrbm", & 
         "eddy-diff scalar (currently rt)","kg/kg",zt)
    k = k + 1

  case default
    write(0,*) 'Error: unrecognized variable in vars_zt: ', & 
       trim(vars_zt(i))
    lerror = .true.

  end select

end do

return

end subroutine stats_init_zt

end module stats_zt
