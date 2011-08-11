!------------------------------------------------------------------------------
! $Id$
module gfdl_activation

  ! Description:
  ! Contains subroutines, functions and variables used by GFDL for their droplet
  ! activation scheme.
  !-----------------------------------------------------------------------------

  implicit none

  ! Namelist parameters
  logical, public  :: nooc        = .false.   ! include organic aerosols as ccns ?
  real, public     :: sul_concen  = 0.1
  real, public     :: low_concen  = 0.1
  real, public     :: high_concen = 1.
  
  !Parameters for look-up tables
  real, public ::  lowup  = 0.3 !m/s
  real, public ::  highup = 10.

  ! earlier values: lowup2 = 0.001, highmass2 = 1000., highmass3 = 1000.
  !real ::  lowup2=0.0001 !m/s
  real, public ::  lowup2    = 0.01   !m/s
  real, public ::  highup2   = 0.3
  real, public ::  lowmass2  = 0.01 !ug m-3
  !real ::  highmass2=1000.
  real, public ::  highmass2 = 100.
  real, public ::  lowmass3  = 0.01 !ug m-3
  !real ::  highmass3=1000.
  real, public ::  highmass3 = 100.
  real, public ::  lowmass4  = 0.01 !ug m-3
  real, public ::  highmass4 = 100.
  real, public ::  lowmass5  = 0.01 !ug m-3
  real, public ::  highmass5 = 100.
  real, public :: lowT2      = 243.15 !K
  real, public :: highT2     = 308.15

  real, public :: aeromass_value = 2.25e-12 

  ! Subroutines
  public :: aer_act_clubb_quadrature_Gauss, Loading

  ! Functions
  private :: erff, get_unit

  private

  contains

  !-----------------------------------------------------------------------------
  SUBROUTINE aer_act_clubb_quadrature_Gauss( aeromass_clubb, temp_clubb_act, & ! Intent(in)
               Ndrop_max   )                                                   ! Intent(out)

        ! Description:
        ! The main subroutine used for the GFDL droplet activation.
        !=======================================================================
    use aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act_wpdf_k
    use variables_prognostic_module, only: pdf_params, p_in_Pa
    use grid_class, only: gr

    implicit none
    ! ---> h1g, 2010-08-24, dumping Nact
    ! removed since it was not being used  -meyern
    !type(time_type),         intent(in)      ::    Time_next
    ! <--- h1g, 2010-08-24

    real,  intent(in),    dimension(:)       ::    temp_clubb_act
    real, intent(inout),  dimension(:, :)    ::    aeromass_clubb
    real, intent(out),    dimension(:)       ::    Ndrop_max

    real                ::  drop
    integer             ::  Tym, ier
    character(len=256)  ::  ermesg

    integer            :: iz_clubb

    real               :: P1_updraft,   P2_updraft  ! probability of updraft
    real, parameter    :: P_updraft_eps = 1.e-16    ! updraft probability threshold
    real, parameter    :: wp2_eps = 0.0001          ! w variance threshold
 
    !=======================================================================
    Tym = size(aeromass_clubb, 2)

    do iz_clubb = 2, gr%nzmax

      if( pdf_params( iz_clubb)%varnce_w1 > wp2_eps) then
        P1_updraft = 0.5 + &
                      0.5*erff( pdf_params(iz_clubb)%w1/sqrt( 2.0*pdf_params(iz_clubb)%varnce_w1) )
        P1_updraft = P1_updraft * pdf_params(iz_clubb)%mixt_frac * pdf_params(iz_clubb)%cloud_frac1
      else
        if( pdf_params( iz_clubb)%w1 > 0.0) &
          P1_updraft = pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1
      end if


      if( pdf_params( iz_clubb)%varnce_w2 > wp2_eps) then
        P2_updraft = 0.5 + &
                      0.5*erff( pdf_params(iz_clubb)%w2/sqrt( 2.0*pdf_params(iz_clubb)%varnce_w2) )
        P2_updraft = P2_updraft * ( 1.0-pdf_params( iz_clubb )%mixt_frac ) &
                      * pdf_params( iz_clubb)%cloud_frac2
      else
        if( pdf_params( iz_clubb)%w2 > 0.0) &
           P2_updraft = ( 1.0-pdf_params( iz_clubb )%mixt_frac ) &
                          * pdf_params( iz_clubb)%cloud_frac2
      end if

      if( P1_updraft + P2_updraft   > P_updraft_eps  ) then
        P1_updraft =  P1_updraft / (  P1_updraft + P2_updraft  )
        P2_updraft =  P2_updraft / (  P1_updraft + P2_updraft  )
      else
        P1_updraft = 0.0
        P2_updraft = 0.0
      end if

      call aer_ccn_act_wpdf_k( temp_clubb_act(iz_clubb), p_in_Pa(iz_clubb), &! intent(in)
                              pdf_params(iz_clubb)%w1,                      &! intent(in)
                              pdf_params(iz_clubb)%varnce_w1,               &! intent(in)
                              aeromass_clubb(iz_clubb, :), Tym,             &! intent(in)
                              drop,   ier,   ermesg )                        ! intent(out)
    
      Ndrop_max(iz_clubb) = drop * P1_updraft

      call aer_ccn_act_wpdf_k( temp_clubb_act(iz_clubb),  p_in_Pa(iz_clubb), &! intent(in)
                             pdf_params(iz_clubb)%w2,                        &! intent(in)
                             pdf_params(iz_clubb)%varnce_w2,                 &! intent(in)
                             aeromass_clubb(iz_clubb, :), Tym,               &! intent(in)
                             drop,   ier,   ermesg )                         ! intent(out)
      ! in-cloud activated droplet concentration
      Ndrop_max(iz_clubb) = Ndrop_max(iz_clubb) + drop * P2_updraft

      ! get the layer-averaged activated droplet concentration (/cm3)
      Ndrop_max(iz_clubb) = Ndrop_max(iz_clubb) *  &
                 (  pdf_params(iz_clubb)%mixt_frac  * pdf_params(iz_clubb)%cloud_frac1 + &
                   (1.- pdf_params(iz_clubb)%mixt_frac) * pdf_params(iz_clubb)%cloud_frac2 )

  end do
return
end subroutine aer_act_clubb_quadrature_Gauss

  !-----------------------------------------------------------------------------
  recursive function erff(x) RESULT(y)

    ! Description:
    ! Error function from Numerical Recipes.
    ! erf(x) = 1 - erfc(x)
    !---------------------------------------------------------------------------
    implicit none

    real, intent(in) :: x
    real dumerfc
    real t, z, y

    z = abs(x)
    t = 1.0 / ( 1.0 + 0.5 * z )

    dumerfc = t * exp(-z * z - 1.26551223 + t *      &
              ( 1.00002368 + t * ( 0.37409196 + t *  &
              ( 0.09678418 + t * (-0.18628806 + t *  &
              ( 0.27886807 + t * (-1.13520398 + t *  &
              ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

    if ( x < 0.0 ) dumerfc = 2.0 - dumerfc
 
    y = 1.0 - dumerfc

  end function erff

  !-----------------------------------------------------------------------------
  subroutine Loading(droplets, droplets2)

    ! Description:
    ! Loads the lookup tables for droplet activation into memory from flat data files.
    !---------------------------------------------------------------------------

    implicit none

    real, dimension(:,:,:,:,:), intent(out) :: droplets, droplets2
    real xx
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
      stop "Unable to open unit."
    else
      get_unit = i
    endif

    return
  end function get_unit

end module gfdl_activation
