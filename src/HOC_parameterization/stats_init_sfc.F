!-----------------------------------------------------------------------
! $Id: stats_init_sfc.F,v 1.2 2007-10-05 15:49:27 faschinj Exp $
#ifdef STATS
!-----------------------------------------------------------------------
      subroutine stats_init_sfc( vars_sfc, lerror )

!     Description:
!     Initializes array indices for sfc
!-----------------------------------------------------------------------
      use hoc_stats
      
      implicit none

      integer, parameter :: nvarmax = 250  ! Max variables
      !Input Variable
      character(len= * ), dimension(nvarmax), intent(in) :: vars_sfc

      !Output Variable	
      logical, intent(out) :: lerror
      
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
     .         "ustar","m/s",sfc)
          k = k + 1

        case ('lh')
          ilh = k
          call stat_assign(ilh,"lh",
     .         "lh","W/m2",sfc)
          k = k + 1

        case ('sh')
          ish = k
          call stat_assign(ish,"sh",
     .         "sh","W/m2",sfc)
          k = k + 1

        case ('cc')
          icc = k
          call stat_assign(icc,"cc",
     .         "cc","count",sfc)
          k = k + 1

        case ('lwp')
          ilwp = k
          call stat_assign(ilwp,"lwp",
     .         "lwp","kg/m2",sfc)
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
     .         "rainfall rate (mm/day)","mm/day",sfc)
          k = k + 1

        case ('pflux')         ! Brian
          ipflux = k

          call stat_assign(ipflux,"prec_flux",
     .         "Precipitation flux (W/m^2)","W/m^2",sfc)
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
