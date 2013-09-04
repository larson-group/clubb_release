!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------
module array_index

! Description:
!   Contains indices to variables in larger arrays.
!   Note that the 'ii' is necessary because 'i' is used in
!   statistics to track locations in the zt/zm/sfc derived types.

! References:
!   None
!-----------------------------------------------------------------------
  implicit none

  ! Variables
  ! Microphysics mixing ratios
  integer, public :: &
    iirrainm,    & ! Rain water mixing ratio  [kg/kg]
    iirsnowm,    & ! Snow mixing ratio        [kg/kg]
    iiricem,     & ! Ice mixing ratio         [kg/kg]
    iirgraupelm    ! Graupel mixing ratio     [kg/kg]
!$omp threadprivate(iirrainm, iirsnowm, iiricem, iirgraupelm)

  ! Microphysics concentrations
  integer, public :: &
    iiNrm,       & ! Rain drop concentration                       [num/kg]
    iiNsnowm,    & ! Snow concentration                            [num/kg]
    iiNim,       & ! Ice concentration                             [num/kg]
    iiNgraupelm, & ! Graupel concentration                         [num/kg]
    iiNcnm,      & ! Cloud nuclei concentration                    [num/kg]
    iiNcm          ! Cloud droplet concentration (not part of PDF) [num/kg]
!$omp threadprivate(iiNrm, iiNsnowm, iiNim, iiNgraupelm, iiNcnm,iiNcm)

  ! Scalar quantities
  integer, public :: & 
    iisclr_rt, iisclr_thl, iisclr_CO2, & ! [kg/kg]/[K]/[1e6 mol/mol]
    iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "
!$omp threadprivate(iisclr_rt, iisclr_thl, iisclr_CO2, &
!$omp   iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2)

  public :: hm_idx

  private ! Default Scope

  contains

  function hm_idx(iiPDF_idx) result(ii_idx)
  ! Description:
  ! Returns the position of a certain hydrometeor within the hydromet arrays
  ! according to its iiPDF index.

  ! References:
  !
  !-----------------------------------------------------------------------

    use corr_matrix_module, only: &
        iiPDF_Ncn, &
        iiPDF_rrain, &
        iiPDF_Nr, &
        iiPDF_rsnow, &
        iiPDF_Nsnow, &
        iiPDF_rice, &
        iiPDF_Ni, &
        iiPDF_rgraupel, &
        iiPDF_Ngraupel

      implicit none

    ! Input Variables
    integer, intent(in) :: iiPDF_idx

    ! Return Variable
    integer :: ii_idx

    ! Local Variables

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( iiPDF_idx == iiPDF_Ncn ) then
       ii_idx = iiNcnm
    endif

    if ( iiPDF_idx == iiPDF_rrain ) then
       ii_idx = iirrainm
    endif

    if ( iiPDF_idx == iiPDF_Nr ) then
       ii_idx = iiNrm
    endif

    if ( iiPDF_idx == iiPDF_rsnow ) then
       ii_idx = iirsnowm
    endif

    if ( iiPDF_idx == iiPDF_Nsnow ) then
       ii_idx = iiNsnowm
    endif

    if ( iiPDF_idx == iiPDF_rice ) then
       ii_idx = iiricem
    endif

    if ( iiPDF_idx == iiPDF_Ni ) then
       ii_idx = iiNim
    endif

    if ( iiPDF_idx == iiPDF_rgraupel ) then
       ii_idx = iirgraupelm
    endif

    if ( iiPDF_idx == iiPDF_Ngraupel ) then
       ii_idx = iiNgraupelm
    endif

    return

  end function hm_idx
  !-----------------------------------------------------------------------

end module array_index
!-----------------------------------------------------------------------
