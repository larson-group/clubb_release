!------------------------------------------------------------------------------
! $Id$
module gfdl_activation

  ! Description:
  ! Contains subroutines, functions and variables used by GFDL for their droplet
  ! activation scheme.
  !-----------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  ! Namelist parameters
  logical, public  :: nooc        = .false.   ! include organic aerosols as ccns ?
  real( kind = core_rknd ), public     :: sul_concen  = 0.1_core_rknd
  real( kind = core_rknd ), public     :: low_concen  = 0.1_core_rknd
  real( kind = core_rknd ), public     :: high_concen = 1._core_rknd
  
  !Parameters for look-up tables
  real( kind = core_rknd ), public ::  lowup  = 0.3_core_rknd !m/s
  real( kind = core_rknd ), public ::  highup = 10._core_rknd

  ! earlier values: lowup2 = 0.001_core_rknd, highmass2 = 1000._core_rknd, 
  !   highmass3 = 1000._core_rknd
  !real( kind = core_rknd) ::  lowup2=0.0001_core_rknd !m/s
  real( kind = core_rknd ), public ::  lowup2    = 0.01_core_rknd   !m/s
  real( kind = core_rknd ), public ::  highup2   = 0.3_core_rknd
  real( kind = core_rknd ), public ::  lowmass2  = 0.01_core_rknd !ug m-3
  !real( kind = core_rknd) ::  highmass2=1000._core_rknd
  real( kind = core_rknd ), public ::  highmass2 = 100._core_rknd
  real( kind = core_rknd ), public ::  lowmass3  = 0.01_core_rknd !ug m-3
  !real( kind = core_rknd) ::  highmass3=1000._core_rknd
  real( kind = core_rknd ), public ::  highmass3 = 100._core_rknd
  real( kind = core_rknd ), public ::  lowmass4  = 0.01_core_rknd !ug m-3
  real( kind = core_rknd ), public ::  highmass4 = 100._core_rknd
  real( kind = core_rknd ), public ::  lowmass5  = 0.01_core_rknd !ug m-3
  real( kind = core_rknd ), public ::  highmass5 = 100._core_rknd
  real( kind = core_rknd ), public :: lowT2      = 243.15_core_rknd !K
  real( kind = core_rknd ), public :: highT2     = 308.15_core_rknd

  real( kind = core_rknd ), public :: aeromass_value = 2.25e-12_core_rknd

  ! Subroutines
  public :: aer_act_clubb_quadrature_Gauss, Loading

  ! Functions
  private :: erff, get_unit

  private

  contains

  !-----------------------------------------------------------------------------
  SUBROUTINE aer_act_clubb_quadrature_Gauss( gr, pdf_params, p_in_Pa,            & ! Intent(in)
                                             aeromass_clubb, temp_clubb_act, & ! Intent(in)
                                             Ndrop_max )                       ! Intent(out)

        ! Description:
        ! The main subroutine used for the GFDL droplet activation.
        !=======================================================================
    use aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act_wpdf_k
    use grid_class, only: grid
    use pdf_parameter_module, only: pdf_parameter
    use clubb_precision, only: &
      core_rknd ! Variable(s)
    implicit none

    type(grid), target, intent(in) :: gr

    intrinsic :: real
    ! ---> h1g, 2010-08-24, dumping Nact
    ! removed since it was not being used  -meyern
    !type(time_type),         intent(in)      ::    Time_next
    ! <--- h1g, 2010-08-24

    type( pdf_parameter ), intent(in) :: pdf_params
    real( kind = core_rknd ),  intent(in),    dimension(:)       ::    p_in_Pa
    real( kind = core_rknd ),  intent(in),    dimension(:)       ::    temp_clubb_act
    real( kind = core_rknd ), intent(inout),  dimension(:, :)    ::    aeromass_clubb
    real( kind = core_rknd ), intent(out),    dimension(:)       ::    Ndrop_max

    real( kind = core_rknd )                ::  drop
    integer             ::  Tym, ier
    character(len=256)  ::  ermesg

    integer            :: iz_clubb

    real( kind = core_rknd ) :: P1_updraft, P2_updraft ! probability of updraft
    ! updraft probability threshold
    real( kind = core_rknd ), parameter :: P_updraft_eps = 1.e-16_core_rknd    
    real( kind = core_rknd ), parameter :: wp2_eps = 0.0001_core_rknd  ! w variance threshold

    ! temporary variables for passing into aer_ccn_act_wpdf_k, which uses reals
    real, dimension(size(aeromass_clubb(1,:))) :: aeromass_clubb_r4
    real                                       :: drop_r4
 
    !=======================================================================
    Tym = size(aeromass_clubb, 2)

    do iz_clubb = 2, gr%nz

      if( pdf_params%varnce_w_1(1,iz_clubb) > wp2_eps) then
        P1_updraft = 0.5_core_rknd + 0.5_core_rknd &
            * erff( pdf_params%w_1(1,iz_clubb) &
            / sqrt( 2.0_core_rknd*pdf_params%varnce_w_1(1,iz_clubb)) )
        P1_updraft = P1_updraft * pdf_params%mixt_frac(1,iz_clubb) &
            * pdf_params%cloud_frac_1(1,iz_clubb)
      else if( pdf_params%w_1(1,iz_clubb) > 0.0_core_rknd) then
          P1_updraft = pdf_params%mixt_frac(1,iz_clubb) * pdf_params%cloud_frac_1(1,iz_clubb)
      else
        ! Eric Raut added to remove compiler warning
          P1_updraft = 0.0_core_rknd
      end if


      if( pdf_params%varnce_w_2(1,iz_clubb) > wp2_eps) then
        P2_updraft = 0.5_core_rknd + &
                      0.5_core_rknd*erff( pdf_params%w_2(1,iz_clubb) &
                      / sqrt( 2.0_core_rknd*pdf_params%varnce_w_2(1,iz_clubb)) )
        P2_updraft = P2_updraft * ( 1.0_core_rknd-pdf_params%mixt_frac(1,iz_clubb) ) &
                      * pdf_params%cloud_frac_2(1,iz_clubb)
      else if( pdf_params%w_2(1,iz_clubb) > 0.0_core_rknd) then
           P2_updraft = ( 1.0_core_rknd-pdf_params%mixt_frac(1,iz_clubb) ) &
                          * pdf_params%cloud_frac_2(1,iz_clubb)
      else
        ! Eric Raut added to remove compiler warning
        P2_updraft = 0.0_core_rknd
      end if

      if( P1_updraft + P2_updraft   > P_updraft_eps  ) then
        P1_updraft =  P1_updraft / (  P1_updraft + P2_updraft  )
        P2_updraft =  P2_updraft / (  P1_updraft + P2_updraft  )
      else
        P1_updraft = 0.0_core_rknd
        P2_updraft = 0.0_core_rknd
      end if

      aeromass_clubb_r4 = real(aeromass_clubb(iz_clubb, :))
      ! Eric Raut initialized drop, which was previously uninitialized!!!
      drop = 0.0_core_rknd
      drop_r4 = real(drop)

      call aer_ccn_act_wpdf_k( real(temp_clubb_act(iz_clubb)), real(p_in_Pa(iz_clubb)),&!intent(in)
                              real(pdf_params%w_1(1,iz_clubb)),                      &! intent(in)
                              real(pdf_params%varnce_w_1(1,iz_clubb)),               &! intent(in)
                              aeromass_clubb_r4, Tym,             &! intent(in)
                              drop_r4,   ier,   ermesg )                        ! intent(out)
    
      Ndrop_max(iz_clubb) = drop * P1_updraft

      call aer_ccn_act_wpdf_k( real(temp_clubb_act(iz_clubb)), real(p_in_Pa(iz_clubb)),&!intent(in)
                             real(pdf_params%w_2(1,iz_clubb)),                        &! intent(in)
                             real(pdf_params%varnce_w_2(1,iz_clubb)),                 &! intent(in)
                             aeromass_clubb_r4, Tym,               &! intent(in)
                             drop_r4,   ier,   ermesg )                         ! intent(out)

      aeromass_clubb(iz_clubb, :) = real(aeromass_clubb_r4, kind = core_rknd)
      drop = real(drop, kind = core_rknd)

      ! in-cloud activated droplet concentration
      Ndrop_max(iz_clubb) = Ndrop_max(iz_clubb) + drop * P2_updraft

      ! get the layer-averaged activated droplet concentration (/cm3)
      Ndrop_max(iz_clubb) = Ndrop_max(iz_clubb) *  &
                 (  pdf_params%mixt_frac(1,iz_clubb)  * pdf_params%cloud_frac_1(1,iz_clubb) + &
                   (1._core_rknd- pdf_params%mixt_frac(1,iz_clubb)) * &
                   pdf_params%cloud_frac_2(1,iz_clubb) )

  end do
return
end subroutine aer_act_clubb_quadrature_Gauss

  !-----------------------------------------------------------------------------
  recursive function erff(x) RESULT(y)

    ! Description:
    ! Error function from Numerical Recipes.
    ! erf(x) = 1 - erfc(x)
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), intent(in) :: x
    real( kind = core_rknd ) dumerfc
    real( kind = core_rknd ) t, z, y

    z = abs(x)
    t = 1.0_core_rknd / ( 1.0_core_rknd + 0.5_core_rknd * z )

    dumerfc = t * exp(-z * z - 1.26551223_core_rknd + t *      &
              ( 1.00002368_core_rknd + t * ( 0.37409196_core_rknd + t *  &
              ( 0.09678418_core_rknd + t * (-0.18628806_core_rknd + t *  &
              ( 0.27886807_core_rknd + t * (-1.13520398_core_rknd + t *  &
              ( 1.48851587_core_rknd + t * (-0.82215223_core_rknd + t *  &
                0.17087277_core_rknd )))))))))

    if ( x < 0.0_core_rknd ) dumerfc = 2.0_core_rknd - dumerfc
 
    y = 1.0_core_rknd - dumerfc

  end function erff

  !-----------------------------------------------------------------------------
  subroutine Loading(droplets, droplets2)

    ! Description:
    ! Loads the lookup tables for droplet activation into memory from flat data files.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), dimension(:,:,:,:,:), intent(out) :: droplets, droplets2
    real( kind = core_rknd ) xx
    integer i, j, k, l, m, unit
    integer res, res2

    ! --- Begin Code ---
    res = size(droplets,1)
    res2 = size(droplets2,1)
    unit = get_unit()
    open(unit, FILE='../input/scm_activation/droplets.dat')
    do k=1,res
      do i=1,res
        do j=1, res
          do l=1, res
            do m=1, res
              read(unit,*) xx
              droplets(m,l,j,i,k)=xx
            end do
          end do
        end do
      end do
    end do
    close(unit)

    unit = get_unit()
    open(unit, FILE='../input/scm_activation/droplets2.dat')
    do k=1,res2
      do i=1,res2
        do j=1, res2
          do l=1, res2
            do m=1, res2
              read(unit,*) xx
              droplets2(m,l,j,i,k)=xx
            end do
          end do
        end do
      end do
    end do
    close(unit)

  end subroutine Loading

  !-----------------------------------------------------------------------------
  integer function get_unit()

    ! Description:
    ! Used to determine a free unit with which to open a file.
    !---------------------------------------------------------------------------

    implicit none

    integer,save :: i
    logical      :: l_open

    ! --- Begin Code ---
    !  9 is reserved for etc_unit
    do i=10,99
      inquire(unit=i,opened=l_open)
      if(.not.l_open)exit
    end do

    if(i==100)then
      ! Eric Raut added to remove compiler warning (obviously result is not used)
      get_unit = 0
      error stop "Unable to open unit."
    else
      get_unit = i
    endif

    return
  end function get_unit

end module gfdl_activation
