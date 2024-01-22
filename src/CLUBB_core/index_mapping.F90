!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module index_mapping

  ! Description:
  ! Functions to map back and forth between the PDF arrays and the hydrometeor
  ! arrays.

  ! The “iiPDF” indices are used to index all PDF variates, including all
  ! hydrometeor variates.  
  ! The “ii” indices are used to index hydrometeor arrays.  
  ! The “ii” variates are a subset of the “iiPDF” variates.  
  ! Conversions between the two sets of indices are done by the 
  ! functions pdf2hydromet_idx and hydromet2pdf_idx below.
  ! 
  ! ------------------------------------------------------------------------
  ! 
  ! iiPDF indices:
  ! 
  ! Included indices:  
  ! iiPDF_chi, iiPDF_eta, iiPDF_w, iiPDF_Ncn, iiPDF_rr, & all other hydrometeors
  ! 
  ! Number of indices:  pdf_dim
  ! 
  ! Examples of arrays dimensioned by pdf_dim:
  ! mu_x_1_n, corr_array_n_cloud, . . .
  ! 
  ! Declared as module variables in module array_index
  ! 
  ! Initialized in subroutine setup_pdf_indices
  ! 
  ! ----------------------------------------------------------------------
  ! 
  ! ii indices:
  ! 
  ! Included indices:  
  ! iirr, iiNr, iiri, iiNi, iirs, iiNs, iirg, iiNg
  ! 
  ! Number of indices:  hydromet_dim
  ! 
  ! Examples of arrays dimensioned by hydromet_dim: 
  ! hydromet, wphydrometp, . . .
  ! 
  ! Declared as module variables in module array_index.
  ! 
  ! Initialized in subroutine init_microphys
  ! 
  ! -----------------------------------------------------------------------
  !
  ! References:
  !   None
  !-------------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: pdf2hydromet_idx, &
            hydromet2pdf_idx, &
            rx2Nx_hm_idx,     &
            Nx2rx_hm_idx,     &
            mvr_hm_max

contains

  !=============================================================================
  function pdf2hydromet_idx( pdf_idx, hm_metadata ) result( hydromet_idx )

    ! Description:
    ! Returns the position of a specific precipitating hydrometeor corresponding
    ! to the PDF index (pdf_idx) in the precipitating hydrometeor array
    ! (hydromet_idx).

    ! References:
    !-----------------------------------------------------------------------

    use corr_varnce_module, only: &
      hm_metadata_type

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      pdf_idx    ! Index of a hydrometeor in the PDF array.

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    ! Return Variable
    integer :: &
      hydromet_idx    ! Index of a hydrometeor in the hydromet array.


    ! Initialize hydromet_idx
    hydromet_idx = 0

    if ( pdf_idx == hm_metadata%iiPDF_rr ) then

       ! Index for rain water mixing ratio, rr.
       hydromet_idx = hm_metadata%iirr

    elseif ( pdf_idx == hm_metadata%iiPDF_Nr ) then

       ! Index for rain drop concentration, Nr.
       hydromet_idx = hm_metadata%iiNr

    elseif ( pdf_idx == hm_metadata%iiPDF_rs ) then

       ! Index for snow mixing ratio, rs.
       hydromet_idx = hm_metadata%iirs

    elseif ( pdf_idx == hm_metadata%iiPDF_Ns ) then

       ! Index for snow flake concentration, Ns.
       hydromet_idx = hm_metadata%iiNs

    elseif ( pdf_idx == hm_metadata%iiPDF_rg ) then

       ! Index for graupel mixing ratio, rg.
       hydromet_idx = hm_metadata%iirg

    elseif ( pdf_idx == hm_metadata%iiPDF_Ng ) then

       ! Index for graupel concentration, Ng.
       hydromet_idx = hm_metadata%iiNg

    elseif ( pdf_idx == hm_metadata%iiPDF_ri ) then

       ! Index for ice mixing ratio, ri.
       hydromet_idx = hm_metadata%iiri

    elseif ( pdf_idx == hm_metadata%iiPDF_Ni ) then

       ! Index for ice concentration, Ni.
       hydromet_idx = hm_metadata%iiNi

    endif

    return

  end function pdf2hydromet_idx

  !=============================================================================
  function hydromet2pdf_idx( hydromet_idx, hm_metadata ) result( pdf_idx )

    ! Description:
    ! Returns the position of a specific precipitating hydrometeor corresponding
    ! to the precipitating hydrometeor index (hydromet_idx) in the PDF array
    ! (pdf_idx).

    ! References:
    !-----------------------------------------------------------------------

    use corr_varnce_module, only: &
      hm_metadata_type

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      hydromet_idx    ! Index of a hydrometeor in the hydromet array.

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    ! Return Variable
    integer :: &
      pdf_idx    ! Index of a hydrometeor in the PDF array.


    ! Initialize pdf_idx.
    pdf_idx = 0

    if ( hydromet_idx == hm_metadata%iirr ) then

       ! Index for rain water mixing ratio, rr.
       pdf_idx = hm_metadata%iiPDF_rr

    elseif ( hydromet_idx == hm_metadata%iiNr ) then

       ! Index for rain drop concentration, Nr.
       pdf_idx = hm_metadata%iiPDF_Nr

    elseif ( hydromet_idx == hm_metadata%iiri ) then

       ! Index for ice mixing ratio, ri.
       pdf_idx = hm_metadata%iiPDF_ri

    elseif ( hydromet_idx == hm_metadata%iiNi ) then

       ! Index for ice concentration, Ni.
       pdf_idx = hm_metadata%iiPDF_Ni

    elseif ( hydromet_idx == hm_metadata%iirs ) then

       ! Index for snow mixing ratio, rs.
       pdf_idx = hm_metadata%iiPDF_rs

    elseif ( hydromet_idx == hm_metadata%iiNs ) then

       ! Index for snow flake concentration, Ns.
       pdf_idx = hm_metadata%iiPDF_Ns

    elseif ( hydromet_idx == hm_metadata%iirg ) then

       ! Index for graupel mixing ratio, rg.
       pdf_idx = hm_metadata%iiPDF_rg

    elseif ( hydromet_idx == hm_metadata%iiNg ) then

       ! Index for graupel concentration, Ng.
       pdf_idx = hm_metadata%iiPDF_Ng

    endif


    return

  end function hydromet2pdf_idx

  !=============================================================================
  function rx2Nx_hm_idx( rx_idx, hm_metadata ) result( Nx_idx )

    ! Description:
    ! Returns the position in the hydrometeor array of the specific
    ! precipitating hydrometeor concentration (Nx_idx) corresponding to the
    ! precipitating hydrometeor mixing ratio (rx_idx) of the same species of
    ! precipitating hydrometeor (rain, ice, snow, or graupel).

    ! References:
    !-----------------------------------------------------------------------

    use corr_varnce_module, only: &
        hm_metadata_type

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      rx_idx    ! Index of the mixing ratio in the hydrometeor array.

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    ! Return Variable
    integer :: &
      Nx_idx    ! Index of the concentration in the hydrometeor array.


    ! Initialize Nx_idx.
    Nx_idx = 0

    if ( rx_idx == hm_metadata%iirr ) then

       ! Index for rain drop concentration, Nr.
       Nx_idx = hm_metadata%iiNr

    elseif ( rx_idx == hm_metadata%iiri ) then

       ! Index for ice crystal concentration, Ni.
       Nx_idx = hm_metadata%iiNi

    elseif ( rx_idx == hm_metadata%iirs ) then

       ! Index for snow flake concentration, Ns.
       Nx_idx = hm_metadata%iiNs

    elseif ( rx_idx == hm_metadata%iirg ) then

       ! Index for graupel concentration, Ng.
       Nx_idx = hm_metadata%iiNg

    endif

    return

  end function rx2Nx_hm_idx

  !=============================================================================
  function Nx2rx_hm_idx( Nx_idx, hm_metadata ) result( rx_idx )

    ! Description:
    ! Returns the position in the hydrometeor array of the specific
    ! precipitating hydrometeor mixing ratio (rx_idx) corresponding to the
    ! precipitating hydrometeor concentration (Nx_idx) of the same species of
    ! precipitating hydrometeor (rain, ice, snow, or graupel).

    ! References:
    !-----------------------------------------------------------------------

    use corr_varnce_module, only: &
        hm_metadata_type

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      Nx_idx    ! Index of the concentration in the hydrometeor array.

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    ! Return Variable
    integer :: &
      rx_idx    ! Index of the mixing ratio in the hydrometeor array.


    ! Initialize rx_idx.
    rx_idx = 0

    if ( Nx_idx == hm_metadata%iiNr ) then

       ! Index for rain water mixing ratio, rr.
       rx_idx = hm_metadata%iirr

    elseif ( Nx_idx == hm_metadata%iiNi ) then

       ! Index for ice mixing ratio, ri.
       rx_idx = hm_metadata%iiri

    elseif ( Nx_idx == hm_metadata%iiNs ) then

       ! Index for snow mixing ratio, rs.
       rx_idx = hm_metadata%iirs

    elseif ( Nx_idx == hm_metadata%iiNg ) then

       ! Index for graupel mixing ratio, rg.
       rx_idx = hm_metadata%iirg

    endif


    return

  end function Nx2rx_hm_idx

  !=============================================================================
  function mvr_hm_max( hydromet_idx, hm_metadata ) result( mvr_hydromet_max )

    ! Description:
    ! Returns the maximum allowable mean volume radius of a specific
    ! precipitating hydrometeor type (rain, ice, snow, or graupel) corresponding
    ! to the precipitating hydrometeor index, whether that index is for the
    ! mixing ratio or concentration associated with that hydrometeor type.

    ! References:
    !-----------------------------------------------------------------------
    
    use corr_varnce_module, only: &
        hm_metadata_type

    use constants_clubb, only: &
        mvr_rain_max,    & ! Constant(s)
        mvr_ice_max,     &
        mvr_snow_max,    &
        mvr_graupel_max, &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      hydromet_idx    ! Index of a hydrometeor in the hydromet array.

    type (hm_metadata_type), intent(in) :: &
      hm_metadata


    ! Return Variable
    real( kind = core_rknd ) :: &
      mvr_hydromet_max    ! Maximum allowable mean volume radius    [m]


    ! Initialize mvr_hydromet_max.
    mvr_hydromet_max = zero

    if ( hydromet_idx == hm_metadata%iirr .or. hydromet_idx == hm_metadata%iiNr ) then

       ! Maximum allowable mean volume radius for rain drops.
       mvr_hydromet_max = mvr_rain_max

    elseif ( hydromet_idx == hm_metadata%iiri .or. hydromet_idx == hm_metadata%iiNi ) then

       ! Maximum allowable mean volume radius for ice crystals.
       mvr_hydromet_max = mvr_ice_max

    elseif ( hydromet_idx == hm_metadata%iirs .or. hydromet_idx == hm_metadata%iiNs ) then

       ! Maximum allowable mean volume radius for snow flakes.
       mvr_hydromet_max = mvr_snow_max

    elseif ( hydromet_idx == hm_metadata%iirg .or. hydromet_idx == hm_metadata%iiNg ) then

       ! Maximum allowable mean volume radius for graupel.
       mvr_hydromet_max = mvr_graupel_max

    endif


    return

  end function mvr_hm_max

!===============================================================================

end module index_mapping


