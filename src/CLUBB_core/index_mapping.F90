!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module index_mapping

  ! Description:
  ! Functions to map back and forth between the PDF arrays and the hydrometeor
  ! arrays.

  ! References:
  !   None
  !-------------------------------------------------------------------------

  ! Hydrometeor array indices
  use array_index, only: &
      iirrainm,    & ! Hydrometeor array index for rain water mixing ratio, rr
      iirsnowm,    & ! Hydrometeor array index for snow mixing ratio, rs
      iiricem,     & ! Hydrometeor array index for ice mixing ratio, ri
      iirgraupelm, & ! Hydrometeor array index for graupel mixing ratio, rg
      iiNrm,       & ! Hydrometeor array index for rain drop concentration, Nr
      iiNsnowm,    & ! Hydrometeor array index for snow concentration, Ns
      iiNim,       & ! Hydrometeor array index for ice concentration, Ni
      iiNgraupelm    ! Hydrometeor array index for graupel concentration, Ng

  ! PDF array indices
  use corr_matrix_module, only: &
      iiPDF_rrain,    & ! PDF array index for rain water mixing ratio, rr
      iiPDF_rsnow,    & ! PDF array index for snow mixing ratio, rs
      iiPDF_rice,     & ! PDF array index for ice mixing ratio, ri
      iiPDF_rgraupel, & ! PDF array index for graupel mixing ratio, rg
      iiPDF_Nr,       & ! PDF array index for rain drop concentration, Nr
      iiPDF_Nsnow,    & ! PDF array index for snow concentration, Ns
      iiPDF_Ni,       & ! PDF array index for ice concentration, Ni
      iiPDF_Ngraupel    ! PDF array index for graupel concentration, Ng

  implicit none

  private ! Default Scope

  public :: pdf2hydromet_idx, &
            hydromet2pdf_idx

contains

  !=============================================================================
  function pdf2hydromet_idx( pdf_idx ) result( hydromet_idx )

    ! Description:
    ! Returns the position of a specific precipitating hydrometeor corresponding
    ! to the PDF index (pdf_idx) in the precipitating hydrometeor array
    ! (hydromet_idx).

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      pdf_idx    ! Index of a hydrometeor in the PDF array.

    ! Return Variable
    integer :: &
      hydromet_idx    ! Index of a hydrometeor in the hydromet array.


    ! Initialize hydromet_idx
    hydromet_idx = 0

    if ( pdf_idx == iiPDF_rrain ) then

       ! Index for rain water mixing ratio, rr.
       hydromet_idx = iirrainm

    elseif ( pdf_idx == iiPDF_Nr ) then

       ! Index for rain drop concentration, Nr.
       hydromet_idx = iiNrm

    elseif ( pdf_idx == iiPDF_rsnow ) then

       ! Index for snow mixing ratio, rs.
       hydromet_idx = iirsnowm

    elseif ( pdf_idx == iiPDF_Nsnow ) then

       ! Index for snow flake concentration, Ns.
       hydromet_idx = iiNsnowm

    elseif ( pdf_idx == iiPDF_rgraupel ) then

       ! Index for graupel mixing ratio, rg.
       hydromet_idx = iirgraupelm

    elseif ( pdf_idx == iiPDF_Ngraupel ) then

       ! Index for graupel concentration, Ng.
       hydromet_idx = iiNgraupelm

    elseif ( pdf_idx == iiPDF_rice ) then

       ! Index for ice mixing ratio, ri.
       hydromet_idx = iiricem

    elseif ( pdf_idx == iiPDF_Ni ) then

       ! Index for ice concentration, Ni.
       hydromet_idx = iiNim

    endif


    return

  end function pdf2hydromet_idx

  !=============================================================================
  function hydromet2pdf_idx( hydromet_idx ) result( pdf_idx )

    ! Description:
    ! Returns the position of a specific precipitating hydrometeor corresponding
    ! to the precipitating hydrometeor index (hydromet_idx) in the PDF array
    ! (pdf_idx).

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      hydromet_idx    ! Index of a hydrometeor in the hydromet array.

    ! Return Variable
    integer :: &
      pdf_idx    ! Index of a hydrometeor in the PDF array.


    ! Initialize pdf_idx.
    pdf_idx = 0

    if ( hydromet_idx == iirrainm ) then

       ! Index for rain water mixing ratio, rr.
       pdf_idx = iiPDF_rrain

    elseif ( hydromet_idx == iiNrm ) then

       ! Index for rain drop concentration, Nr.
       pdf_idx = iiPDF_Nr

    elseif ( hydromet_idx == iiricem ) then

       ! Index for ice mixing ratio, ri.
       pdf_idx = iiPDF_rice

    elseif ( hydromet_idx == iiNim ) then

       ! Index for ice concentration, Ni.
       pdf_idx = iiPDF_Ni

    elseif ( hydromet_idx == iirsnowm ) then

       ! Index for snow mixing ratio, rs.
       pdf_idx = iiPDF_rsnow

    elseif ( hydromet_idx == iiNsnowm ) then

       ! Index for snow flake concentration, Ns.
       pdf_idx = iiPDF_Nsnow

    elseif ( hydromet_idx == iirgraupelm ) then

       ! Index for graupel mixing ratio, rg.
       pdf_idx = iiPDF_rgraupel

    elseif ( hydromet_idx == iiNgraupelm ) then

       ! Index for graupel concentration, Ng.
       pdf_idx = iiPDF_Ngraupel

    endif


    return

  end function hydromet2pdf_idx

!===============================================================================

end module index_mapping
