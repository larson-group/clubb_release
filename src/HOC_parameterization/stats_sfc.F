!-----------------------------------------------------------------------
! $Id: stats_sfc.F,v 1.10 2008-06-24 18:13:58 faschinj Exp $
      
      module stats_sfc
#ifdef STATS
 
      implicit none
      
      private ! Set Default Scope

      public :: stats_init_sfc
      
      contains 

!-----------------------------------------------------------------------
      subroutine stats_init_sfc( vars_sfc, lerror )

!     Description:
!     Initializes array indices for sfc
!-----------------------------------------------------------------------

      use stats_variables, only:
     .    sfc, ! Variables
     .    iustar,
     .    ilh,
     .    ish,
     .    icc,
     .    ilwp,
     .    izb,
     .    izi,
     .    irain,
     .    ipflux,
     .    irrm_sfc,
     .    iwp23_cn,
     .    irtm_cn,
     .    ithlm_cn,
     .    irtp2_cn,
     .    ithlp2_cn,
     .    irtpthlp_cn,
     .    iup2_cn,
     .    ivp2_cn
      use stats_type, only:
     .    stat_assign ! Procedure

      implicit none

      integer, parameter :: nvarmax = 250  ! Max variables
      !Input Variable
      character(len= * ), dimension(nvarmax), intent(in) :: vars_sfc

      !Output Variable	
      logical, intent(inout) :: lerror
      
      !Local Varables
      integer :: i, k

! Default initialization for array indices for sfc

      iustar   = 0
      ilh      = 0
      ish      = 0
      icc      = 0
      ilwp     = 0
      izb      = 0
      izi      = 0
      irain    = 0   ! Brian
      ipflux   = 0   ! Brian
      irrm_sfc = 0   ! Brian

      ! These are estimates of the condition number on each implicit
      ! matrices, and not located at the surface of the domain.
      iwp23_cn    = 0
      irtm_cn     = 0
      ithlm_cn    = 0
      irtp2_cn    = 0
      ithlp2_cn   = 0
      irtpthlp_cn = 0
      iup2_cn     = 0
      ivp2_cn     = 0

      ! Assign pointers for statistics variables sfc

      k = 1
      do i=1,sfc%nn

        select case ( trim(vars_sfc(i)) )

        case ('ustar')
          iustar = k

          call stat_assign(iustar,"ustar",
     .         "ustar [m/s]","m/s",sfc)
          k = k + 1

        case ('lh')
          ilh = k
          call stat_assign(ilh,"lh",
     .         "Surface latent heating [W/m^2]","W/m2",sfc)
          k = k + 1

        case ('sh')
          ish = k
          call stat_assign(ish,"sh",
     .         "Surface sensible heating [W/m^2]","W/m2",sfc)
          k = k + 1

        case ('cc')
          icc = k
          call stat_assign(icc,"cc",
     .         "Cloud cover","count",sfc)
          k = k + 1

        case ('lwp')
          ilwp = k
          call stat_assign(ilwp,"lwp",
     .         "Liquid water path [kg/m^2]","kg/m2",sfc)
          k = k + 1

        case ('zb')
          izb = k
          call stat_assign(izb,"zb",
     .         "zb","m",sfc)
          k = k + 1

        case ('zi')
          izi = k
          call stat_assign(izi,"zi",
     .         "zi","m",sfc)
          k = k + 1

        case ('rain')          ! Brian
          irain = k
          call stat_assign(irain,"rain_rate",
     .         "Surface rainfall rate [mm/day]","mm/day",sfc)
          k = k + 1

        case ('pflux')         ! Brian
          ipflux = k

          call stat_assign(ipflux,"prec_flux",
     .         "Surface precipitation flux [W/m^2]","W/m^2",sfc)
          k = k + 1

        case ('rrm_sfc')       ! Brian
          irrm_sfc = k

          call stat_assign(irrm_sfc,"rrm_sfc",
     .         "Surface rain water mixing ratio","kg/kg",sfc)
          k = k + 1

        case ('wp23_cn')
          iwp23_cn = k
          call stat_assign(iwp23_cn,"wp23_cn",
     .         "Estimate of the condition number for wp2/3","count",sfc)
          k = k + 1

        case ('thlm_cn')
          ithlm_cn = k
          call stat_assign(ithlm_cn,"thlm_cn",
     .         "Estimate of the condition number for thlm/wpthlp",
     .         "count",sfc)
          k = k + 1

        case ('rtm_cn')
          irtm_cn = k

          call stat_assign(irtm_cn,"rtm_cn",
     .         "Estimate of the condition number for rtm/wprtp",
     .         "count",sfc)
          k = k + 1

        case ('thlp2_cn')
          ithlp2_cn = k

          call stat_assign(ithlp2_cn,"thlp2_cn",
     .         "Estimate of the condition number for thlp2",
     .         "count",sfc)
          k = k + 1

        case ('rtp2_cn')
          irtp2_cn = k
          call stat_assign(irtp2_cn,"rtp2_cn",
     .         "Estimate of the condition number for rtp2",
     .         "count",sfc)
          k = k + 1

        case ('rtpthlp_cn')
          irtpthlp_cn = k
           call stat_assign(irtpthlp_cn,"rtpthlp_cn",
     .         "Estimate of the condition number for rtpthlp",
     .         "count",sfc)
          k = k + 1

        case ('up2_cn')
          iup2_cn = k
          call stat_assign(iup2_cn,"up2_cn",
     .         "Estimate of the condition number for up2","count",sfc)
          k = k + 1

        case ('vp2_cn')
          ivp2_cn = k  
          call stat_assign(ivp2_cn,"vp2_cn",
     .         "Estimate of the condition number for vp2","count",sfc)

          k = k + 1

        case default
          write(0,*) 'Error: unrecognized variable in vars_sfc: ',
     .          trim(vars_sfc(i))
          lerror = .true.

        end select

      end do
      
      return
      
      end subroutine stats_init_sfc
      
#endif /*STATS*/
      end module stats_sfc

