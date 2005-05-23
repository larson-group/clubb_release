!------------------------------------------------------------------------
! Definition of a grid class and associated functions
!
! The grid specification is as follows:
!
!   +              ================== zm(nnzp) =========GP=======
!   |
!   |
! dzt(nnzp)  +     ------------------ zt(nnzp) ---------GP-------
!   |        |
!   |        |
!   +  dzm(nnzp-1) ================== zm(nnzp-1) ================
!            |
!            |
!            +     ------------------ zt(nnzp-1) ----------------
! 
!                                       .
!                                       .
!                                       .
!                                       .
! 
!   +              ================== zm(2) =====================
!   |
!   |
! dzt(2)     +     ------------------ zt(2) ---------------------
!   |        |
!   |        |
!   +      dzm(1)  ================== zm(1) ============GP=======  surface
!            |     //////////////////////////////////////////////
!            |
!            +     ------------------ zt(1) ------------GP-------
!
!
! GP indicates ghost points. Variables located at those levels are
! not prognosed, but only used for boundary conditions.
!
! Chris Golaz, 7/17/99
! modified 9/10/99
!
! modified 2/28/01 to get rid of UNDEFINED.
! Experience has shown that the 'UNDEFINED' idea was not shuch a good idea
! after all.
!
        module grid_class

           include 'grid.h'

           type grid
              integer :: nnzp
              real :: zm(1:NMAX),zt(1:NMAX),dzm(1:NMAX-1),dzt(1:NMAX)
           end type grid

           type (grid) gr       ! The grid is defined here so that it
                                ! is common throughout the module
                                ! The implication is that only one grid
                                ! can be defined !

           interface zt2zm
              module procedure interpolated_azm, interpolated_azmk
           end interface

           interface zm2zt
              module procedure interpolated_azt, interpolated_aztk
           end interface

           interface ddzm
              module procedure gradzm
           end interface

           interface ddzt
              module procedure gradzt
           end interface

           contains

!          --------------------------------------------------------------
!           Grid Constructor
!
!           This subroutine sets up a uniformly spaced vertical grid similar
!           to RAMS grid.
!
!          Input parameters:
!
!             integer nnzp        number of vertical levels in grid
!             real deltaz        vertical grid spacing

           subroutine gridsetup( nnzp, deltaz, runtype )
           implicit none

!           Input

           integer, intent(in) :: nnzp
           real, intent(in)    :: deltaz
           integer, intent(in) :: runtype

!          Internal

           integer i

!           Define grid size

           if ( nnzp > NMAX ) then
              write(*,*) 'ERROR: nnzp > NMAX'
              write(*,*) 'nnzp = ',nnzp
              write(*,*) 'NMAX = ',NMAX
              write(*,*) 'Check file grid.h'
              stop 'gridsetup'
           end if
           gr%nnzp = nnzp
!          Modified code from Brian. Turn on if desired

!           if (runtype == 1) then
!
!              gr%zm(1) =     0.0000000
!              gr%zm(2) =      33.01319
!              gr%zm(3) =      66.16915
!              gr%zm(4) =      99.46445
!              gr%zm(5) =      132.8614
!              gr%zm(6) =      166.3603
!              gr%zm(7) =      199.9620
!              gr%zm(8) =      233.6674
!              gr%zm(9) =      267.4767
!              gr%zm(10) =     301.3913
!              gr%zm(11) =     335.4112
!              gr%zm(12) =     369.5378
!              gr%zm(13) =     403.7718
!              gr%zm(14) =     438.1139
!              gr%zm(15) =     472.5654
!              gr%zm(16) =     507.1265
!              gr%zm(17) =     541.7983
!              gr%zm(18) =     576.5813
!              gr%zm(19) =     611.4764
!              gr%zm(20) =     646.4847
!              gr%zm(21) =     681.6051
!              gr%zm(22) =     716.8403
!              gr%zm(23) =     752.1945
!              gr%zm(24) =     787.6716
!              gr%zm(25) =     823.2727
!              gr%zm(26) =     858.9983
!              gr%zm(27) =     894.8497
!              gr%zm(28) =     930.8275
!              gr%zm(29) =     966.9328
!              gr%zm(30) =     1003.166
!              gr%zm(31) =     1039.529
!              gr%zm(32) =     1076.023
!              gr%zm(33) =     1112.648
!              gr%zm(34) =     1149.405
!              gr%zm(35) =     1186.295
!              gr%zm(36) =     1223.320
!              gr%zm(37) =     1260.480
!              gr%zm(38) =     1297.777
!              gr%zm(39) =     1335.204
!              gr%zm(40) =     1372.755
!              gr%zm(41) =     1410.431
!              gr%zm(42) =     1448.231
!              gr%zm(43) =     1486.158
!              gr%zm(44) =     1524.213
!              gr%zm(45) =     1562.395
!              gr%zm(46) =     1600.706
!              gr%zm(47) =     1639.147
!              gr%zm(48) =     1677.719
!              gr%zm(49) =     1716.423
!              gr%zm(50) =     1755.260
!              gr%zm(51) =     1794.230
!              gr%zm(52) =     1833.336
!              gr%zm(53) =     1872.576
!              gr%zm(54) =     1911.954
!              gr%zm(55) =     1951.470
!              gr%zm(56) =     1991.125
!              gr%zm(57) =     2030.919
!              gr%zm(58) =     2070.854
!              gr%zm(59) =     2110.932
!              gr%zm(60) =     2151.153
!              gr%zm(61) =     2191.518
!              gr%zm(62) =     2232.028
!              gr%zm(63) =     2272.685
!              gr%zm(64) =     2313.489
!              gr%zm(65) =     2354.442
!              gr%zm(66) =     2395.546
!              gr%zm(67) =     2436.800
!              gr%zm(68) =     2478.207
!              gr%zm(69) =     2519.779
!              gr%zm(70) =     2561.540
!              gr%zm(71) =     2603.504
!              gr%zm(72) =     2645.673
!              gr%zm(73) =     2688.050
!              gr%zm(74) =     2730.636
!              gr%zm(75) =     2773.433
!              gr%zm(76) =     2816.444
!              gr%zm(77) =     2859.670
!              gr%zm(78) =     2903.113
!              gr%zm(79) =     2946.776
!              gr%zm(80) =     2990.661
!              gr%zm(81) =     3034.770
!              gr%zm(82) =     3079.105
!              gr%zm(83) =     3123.668
!              gr%zm(84) =     3168.463
!              gr%zm(85) =     3213.491
!              gr%zm(86) =     3258.754
!              gr%zm(87) =     3304.256
!              gr%zm(88) =     3349.998
!              gr%zm(89) =     3395.982
!              gr%zm(90) =     3442.213
!              gr%zm(91) =     3488.692
!              gr%zm(92) =     3535.421
!              gr%zm(93) =     3582.403
!              gr%zm(94) =     3629.642
!              gr%zm(95) =     3677.139
!              gr%zm(96) =     3724.899
!              gr%zm(97) =     3772.922
!              gr%zm(98) =     3821.213
!              gr%zm(99) =     3869.775
!              gr%zm(100) =    3918.609
!              gr%zm(101) =    3967.720
!              gr%zm(102) =    4017.111
!              gr%zm(103) =    4066.784
!              gr%zm(104) =    4116.742
!              gr%zm(105) =    4166.990
!              gr%zm(106) =    4217.530
!              gr%zm(107) =    4268.366
!              gr%zm(108) =    4319.500
!              gr%zm(109) =    4370.938
!              gr%zm(110) =    4422.681
!
!           elseif (runtype == 2) then      ! BOMEX case
!
!              gr%zm(1) =     0.0000000
!              gr%zm(2) =      35.94881
!              gr%zm(3) =      71.99928
!              gr%zm(4) =      108.1521
!              gr%zm(5) =      144.4082
!              gr%zm(6) =      180.7679
!              gr%zm(7) =      217.2321
!              gr%zm(8) =      253.8017
!              gr%zm(9) =      290.4773
!              gr%zm(10) =     327.2596
!              gr%zm(11) =     364.1495
!              gr%zm(12) =     401.1478
!              gr%zm(13) =     438.2553
!              gr%zm(14) =     475.4726
!              gr%zm(15) =     512.8007
!              gr%zm(16) =     550.2426
!              gr%zm(17) =     587.8048
!              gr%zm(18) =     625.4881
!              gr%zm(19) =     663.2934
!              gr%zm(20) =     701.2217
!              gr%zm(21) =     739.2736
!              gr%zm(22) =     777.4504
!              gr%zm(23) =     815.7531
!              gr%zm(24) =     854.1823
!              gr%zm(25) =     892.7392
!              gr%zm(26) =     931.4247
!              gr%zm(27) =     970.2397
!              gr%zm(28) =     1009.185
!              gr%zm(29) =     1048.262
!              gr%zm(30) =     1087.471
!              gr%zm(31) =     1126.813
!              gr%zm(32) =     1166.290
!              gr%zm(33) =     1205.902
!              gr%zm(34) =     1245.650
!              gr%zm(35) =     1285.536
!              gr%zm(36) =     1325.561
!              gr%zm(37) =     1365.725
!              gr%zm(38) =     1406.029
!              gr%zm(39) =     1446.475
!              gr%zm(40) =     1487.064
!              gr%zm(41) =     1527.813
!              gr%zm(42) =     1568.729
!              gr%zm(43) =     1609.816
!              gr%zm(44) =     1651.073
!              gr%zm(45) =     1692.502
!              gr%zm(46) =     1734.104
!              gr%zm(47) =     1775.880
!              gr%zm(48) =     1817.832
!              gr%zm(49) =     1859.961
!              gr%zm(50) =     1902.268
!              gr%zm(51) =     1944.755
!              gr%zm(52) =     1987.422
!              gr%zm(53) =     2030.270
!              gr%zm(54) =     2073.292
!              gr%zm(55) =     2116.489
!              gr%zm(56) =     2159.864
!              gr%zm(57) =     2203.417
!              gr%zm(58) =     2247.150
!              gr%zm(59) =     2291.065
!              gr%zm(60) =     2335.164
!              gr%zm(61) =     2379.447
!              gr%zm(62) =     2423.917
!              gr%zm(63) =     2468.575
!              gr%zm(64) =     2513.423
!              gr%zm(65) =     2558.463
!              gr%zm(66) =     2603.696
!              gr%zm(67) =     2649.124
!              gr%zm(68) =     2694.749
!              gr%zm(69) =     2740.573
!              gr%zm(70) =     2786.596
!              gr%zm(71) =     2832.823
!              gr%zm(72) =     2879.254
!              gr%zm(73) =     2925.891
!              gr%zm(74) =     2972.735
!              gr%zm(75) =     3019.790
!
!           endif
!
!           gr%zt(1) = -0.5*(gr%zm(2)+gr%zm(1))
!           do i = 2, gr%nnzp, 1
!             gr%zt(i) = 0.5*(gr%zm(i)+gr%zm(i-1))
!           enddo
!
!!           gr%zt(1) = -0.5 * deltaz
!!           do i = 1, gr%nnzp-1, 1
!!             gr%zt(i) = 0.5*(gr%zm(i)+gr%zm(i+1))
!!           enddo
!!           gr%zt(gr%nnzp) = 0.5*(gr%zm(gr%nnzp)-gr%zm(gr%nnzp-1)) 
!!     .                    + gr%zm(gr%nnzp)


!          Original code
!          Define momentum grid points. The first momentum level is 
!          at the surface

           do i=1,gr%nnzp
             gr%zm(i) = (i-1) * deltaz
           end do

!          Define thermodynamic grid. For the time being, themodynamic 
!          grid points are located in the center point between momentum 
!          points. The first thermodynmic point is below the surface.

           gr%zt(1) = -0.5 * deltaz
           do i=2,gr%nnzp
             gr%zt(i) = 0.5 * ( gr%zm(i) + gr%zm(i-1) )
           end do

!          Define dzm. dzm is the inverse spacing between thermodynamic 
!          points

           do i=1,gr%nnzp-1
              gr%dzm(i) = 1. / ( gr%zt(i+1) - gr%zt(i) )
           end do
           gr%dzm(gr%nnzp) = gr%dzm(gr%nnzp-1)

!          Define dzt. dzt is the inverse spacing between momentum points

           do i=2,gr%nnzp
              gr%dzt(i) = 1. / ( gr%zm(i) - gr%zm(i-1) )
           end do
           gr%dzt(1) = gr%dzt(2)

           return
           end subroutine gridsetup

!          --------------------------------------------------------------
!           Function to interpolate a variable located on the 
!          thermodynamic grid (azt) to a variable on the momentum 
!          grid (azm)

           function interpolated_azm( azt )
           implicit none

!          Input

           real, dimension(1:gr%nnzp), intent(in) :: azt

!          Output

           real, dimension(1:gr%nnzp) :: interpolated_azm
        
!          Internal

           integer i

!          Do actual interpolation 

           do i=1, gr%nnzp-1
              interpolated_azm(i) = 0.5 * ( azt(i) + azt(i+1) )
           end do
!           interpolated_azm(gr%nnzp) = UNDEFINED
           interpolated_azm(gr%nnzp) = azt(gr%nnzp)

           return
           end function interpolated_azm

!          --------------------------------------------------------------
!           Function to interpolate a variable located on the
!           thermodynamic grid (azt) to a variable on the momentum
!           grid (azm). The function interpolates only for level k

           function interpolated_azmk( azt, k )
           implicit none

!           Input

           real, dimension(1:gr%nnzp), intent(in) :: azt
           integer, intent(in)                    :: k

!           Output

           real interpolated_azmk

!           Do actual interpolation

           if ( k /= gr%nnzp ) then
              interpolated_azmk = 0.5 * ( azt(k) + azt(k+1) )
           else
!              interpolated_azmk = UNDEFINED
              interpolated_azmk = azt(gr%nnzp)
           end if

           return
           end function interpolated_azmk
        
!          --------------------------------------------------------------
!           Function to interpolate a variable located on the 
!          momentum grid (azm) to a variable on the thermodynamic
!          grid (azt)

           function interpolated_azt( azm )
           implicit none

!          Input

           real, dimension(1:gr%nnzp), intent(in) :: azm

!          Output

           real, dimension(1:gr%nnzp) :: interpolated_azt
        
!          Internal

           integer i

!          Do actual interpolation 

           do i=gr%nnzp,2,-1
              interpolated_azt(i) = 0.5 * ( azm(i-1) + azm(i) )
           end do
!           interpolated_azt(1) = UNDEFINED
           interpolated_azt(1) = azm(1)

           return
           end function interpolated_azt

!          --------------------------------------------------------------
!           Function to interpolate a variable located on the 
!          momentum grid (azm) to a variable on the thermodynamic
!          grid (azt). This function interpolates only for level k.

           function interpolated_aztk( azm, k )
           implicit none

!          Input

           real, dimension(1:gr%nnzp), intent(in) :: azm
           integer, intent(in)                    :: k

!          Output

           real interpolated_aztk
        
!          Do actual interpolation 

           if ( k /= 1 ) then
             interpolated_aztk = 0.5 * ( azm(k-1) + azm(k) )
           else
!             interpolated_aztk = UNDEFINED
             interpolated_aztk = azm(1)
           end if

           return
           end function interpolated_aztk

!          --------------------------------------------------------------
!          Function to compute the vertical derivative of a variable (azm)
!          located on the momentum grid. The results is returned in an 
!          array defined on the thermodynamic grid.

           function gradzm( azm )
           implicit none

!           Input

           real, dimension(1:gr%nnzp), intent(in) :: azm

!          Output

           real, dimension(1:gr%nnzp) :: gradzm

!           Internal

           integer i

!           Compute vertical derivative

           do i=gr%nnzp, 2, -1
              gradzm(i) = ( azm(i) - azm(i-1) ) * gr%dzt(i)
           end do
!           gradzm(1) = UNDEFINED
           gradzm(1) = 0.

           return
           end function gradzm

!          --------------------------------------------------------------
!           Function to compute the vertical derivative of a variable (azt)
!          located on the thermodynamic grid. The results is returned in
!          an array defined on the momentum grid.

           function gradzt( azt )
           implicit none

!           Input

           real, dimension(1:gr%nnzp), intent(in) :: azt

!          Output

           real, dimension(1:gr%nnzp) :: gradzt

!           Internal

           integer i

!          Compute vertical derivative

           do i=1,gr%nnzp-1
              gradzt(i) = ( azt(i+1) - azt(i) ) * gr%dzm(i)
           end do
!           gradzt(gr%nnzp) = UNDEFINED
           gradzt(gr%nnzp) = 0.

           return
           end function gradzt

        end module grid_class

!------------------------------------------------------------------------
