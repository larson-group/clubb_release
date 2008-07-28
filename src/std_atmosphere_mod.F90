!----------------------------------------------------------------------------
! $Id: std_atmosphere_mod.F90,v 1.2 2008-07-28 19:45:10 faschinj Exp $
module std_atmosphere_mod
  
  implicit none
  
  private ! Default Scope

  public :: std_alt, std_T_in_K, std_sp_hmdty, std_pinmb,  & 
            std_o3l, std_atmos_dim, std_atmosphere
  ! Size of U.S. Standard Atmo. data
  integer, parameter :: std_atmos_dim = 50
        

  ! Parameters from U.S. Standard Atmosphere, 1976; 
  ! Starting at 1Km altitude

  ! Altitude in meters
  double precision, parameter, dimension(std_atmos_dim) ::  & 
   std_alt = (/   & 
     1000.00,  2000.00,  3000.00,  4000.00,  5000.00,   & 
     6000.00,  7000.00,  8000.00,  9000.00,  10000.0,   & 
     11000.0,  12000.0,  13000.0,  14000.0,  15000.0,   & 
     16000.0,  17000.0,  18000.0,  19000.0,  20000.0,   & 
     21000.0,  22000.0,  23000.0,  24000.0,  25000.0,   & 
     26000.0,  27000.0,  28000.0,  29000.0,  30000.0,   & 
     31000.0,  32000.0,  33000.0,  34000.0,  35000.0,   & 
     36000.0,  37000.0,  38000.0,  39000.0,  40000.0,   & 
     41000.0,  42000.0,  43000.0,  44000.0,  45000.0,   & 
     46000.0,  47000.0,  48000.0,  49000.0,  50000.0   /)

  ! Temperature in degrees Kelvin
  double precision, parameter, dimension(std_atmos_dim) ::  & 
   std_T_in_K = (/   & 
     281.600,  275.100,  268.700,  262.200,  255.700,   & 
     249.200,  242.700,  236.200,  229.700,  223.200,   & 
     216.800,  216.600,  216.600,  216.600,  216.600,   & 
     216.600,  216.600,  216.600,  216.600,  216.600,   & 
     217.600,  218.600,  219.600,  220.600,  221.600,   & 
     222.580,  223.560,  224.540,  225.520,  226.500,   & 
     228.500,  230.500,  232.500,  234.500,  236.500,   & 
     239.280,  242.060,  244.840,  247.620,  250.400,   & 
     253.160,  255.920,  258.680,  261.440,  264.200,   & 
     265.480,  266.760,  268.040,  269.320,  270.600   /)
  
  ! Specific Humidity ( Water Vapor / Density )
  double precision, parameter, dimension(std_atmos_dim) ::  & 
   std_sp_hmdty = (/  & 
     0.378038E-02,  0.287984E-02,  0.197954E-02,  0.134261E-02, & 
     0.869093E-03,  0.575670E-03,  0.355932E-03,  0.228224E-03,   & 
     0.984800E-04,  0.435308E-04,  0.224781E-04,  0.118628E-04,   & 
     0.675169E-05,  0.368583E-05,  0.369610E-05,  0.366366E-05,   & 
     0.365425E-05,  0.361842E-05,  0.423077E-05,  0.494882E-05,  & 
     0.633914E-05,  0.806077E-05,  0.103636E-04,  0.129953E-04,   & 
     0.164671E-04,  0.173018E-04,  0.181366E-04,  0.189714E-04,   & 
     0.198062E-04,  0.206410E-04,  0.202939E-04,  0.199469E-04,   & 
     0.195999E-04,  0.192529E-04,  0.189058E-04,  0.184780E-04,   & 
     0.180502E-04,  0.176224E-04,  0.171946E-04,  0.167668E-04,  & 
     0.166688E-04,  0.165707E-04,  0.164727E-04,  0.163747E-04,   & 
     0.162767E-04,  0.153583E-04,  0.144398E-04,  0.135214E-04,   & 
     0.126030E-04,  0.116845E-04  /)

  ! Pressure in millibars
  double precision, parameter, dimension(std_atmos_dim) :: & 
   std_pinmb = (/   & 
     898.600,  795.000,  701.200,  616.600,  540.500,   & 
     472.200,  411.100,  356.500,  308.000,  265.000,   & 
     227.000,  194.000,  165.800,  141.700,  121.100,   & 
     103.500,  88.5000,  75.6500,  64.6700,  55.2900,   & 
     47.2900,  40.4700,  34.6700,  29.7200,  25.4600,   & 
     22.7620,  20.0640,  17.3660,  14.6680,  11.9700,   & 
     10.7252,  9.48040,  8.23560,  6.99080,  5.74600,   & 
     5.17100,  4.59600,  4.02100,  3.44600,  2.87100,   & 
     2.59500,  2.31900,  2.04300,  1.76700,  1.49100,   & 
     1.35236,  1.21372,  1.07508,  0.936440, 0.797800   /)



  ! Ozone ( O_3 / Density )
  double precision, parameter, dimension(std_atmos_dim) ::   & 
    std_o3l= (/   & 
      0.486049E-07,  0.536246E-07,  0.549874E-07,  0.561455E-07,   & 
      0.611081E-07,  0.681715E-07,  0.813559E-07,  0.988969E-07,   & 
      0.152002E-06,  0.217654E-06,  0.356360E-06,  0.512985E-06,   & 
      0.637659E-06,  0.833699E-06,  0.107803E-05,  0.138138E-05,   & 
      0.196767E-05,  0.263158E-05,  0.336538E-05,  0.427398E-05,  & 
      0.501849E-05,  0.604557E-05,  0.690909E-05,  0.766937E-05,   & 
      0.848303E-05,  0.895916E-05,  0.943528E-05,  0.991141E-05,   & 
      0.103875E-04,  0.108637E-04,  0.112905E-04,  0.117173E-04,   & 
      0.121441E-04,  0.125709E-04,  0.129978E-04,  0.128507E-04,   & 
      0.127036E-04,  0.125565E-04,  0.124094E-04,  0.122623E-04,  & 
      0.115392E-04,  0.108162E-04,  0.100931E-04,  0.937005E-05,   & 
      0.864700E-05,  0.769657E-05,  0.674614E-05,  0.579570E-05,   & 
      0.484527E-05,  0.389484E-05  /)
  
  contains
 
!-----------------------------------------------------------------------
  subroutine std_atmosphere( alt, theta, rtm )
!
!       Description: 
!       Given a specific altitude this subroutine will return interpolated 
!       values for theta and rtm from U.S. Standard Atmosphere data.

!       References:
!       McClatchey, et al., (1972) _Environmental Research Papers_, 
!       No. 411, p.94
!-----------------------------------------------------------------------
  use constants, only:  & 
      p0,  & ! Variable(s) 
      kappa

  use interpolation, only:  & 
      linint,  & ! Procedure(s) 
      binary_search
  
  implicit none
  
  intrinsic :: real

  ! Input Variable
  real,intent(in) :: alt      ! Altitude                      [m]

  ! Output Variables
  real,intent(out) ::  & 
  theta,                       & ! Potential Temperature         [K]        
  rtm                         ! Total Water Mixing Ratio      [kg/kg]
  
  ! Local Variables
  real ::  & 
  exner,               & ! Exner function                [-]
  press,               & ! Pressure                      [hPa]
  sp_humidity,         & ! Specific humidity             [kg/kg]
  tabs0               ! Temperature                   [K]
  
  ! These variables are used to make the calls to linint cleaner    
  real, dimension(std_atmos_dim) :: & 
  T_in_K, & 
  pinmb, & 
  sp_hmdty, & 
  height

  integer :: varindex
  
  T_in_K    = real( std_T_in_K )
  pinmb    = real( std_pinmb )
  sp_hmdty = real( std_sp_hmdty )
  height   = real( std_alt )
  
  varindex = binary_search( std_atmos_dim, height, alt )

  if( varindex < 0 ) then 
     stop "Cannot find altitude in Standard Atmosphere"
  endif

  ! Compute thlm from Standard Atmosphere

  tabs0 = linint( alt, height(varindex), height(varindex-1),  & 
                  T_in_K(varindex), T_in_K(varindex-1) )
  
  press = 100. *  & 
          linint( alt, height(varindex), height(varindex-1), & 
                  pinmb(varindex), pinmb(varindex-1) )

  exner = (press/p0)**kappa

  theta = tabs0/exner

  ! Compute rtm

  sp_humidity = linint( alt, height(varindex), height(varindex-1), & 
                sp_hmdty(varindex), sp_hmdty(varindex-1))

  rtm = sp_humidity/( 1. - sp_humidity)


  return

  end subroutine std_atmosphere

end module std_atmosphere_mod
