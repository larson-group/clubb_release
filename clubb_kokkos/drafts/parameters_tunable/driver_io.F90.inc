  subroutine f2k_proc_read_clubb_params_nl_c(input,output) bind(C)
    ! Description: Read &clubb_params_nl using initialized caller arrays and a call-local unit.
    implicit none
    type(wire_field), intent(in) :: input(*)
    type(wire_field), intent(out) :: output(*)
    integer(c_int), pointer :: ip(:)
    integer(c_int8_t), pointer :: bp(:)
    real(core_rknd), pointer :: rp(:)
    real(c_double), pointer :: dp(:)
    character(c_char), pointer :: cp(:)
    integer :: n, j, k
    integer :: unit
    ! Required iunit is retained, but NEWUNIT follows the existing call-local I/O owner.
    namelist /clubb_params_nl/ C1, &
      C1b, &
      C1c, &
      C2rt, &
      C2thl, &
      C2rtthl, &
      C4, &
      C_uu_shr, &
      C_uu_buoy, &
      C6rt, &
      C6rtb, &
      C6rtc, &
      C6thl, &
      C6thlb, &
      C6thlc, &
      C7, &
      C7b, &
      C7c, &
      C8, &
      C8b, &
      C10, &
      C11, &
      C11b, &
      C11c, &
      C12, &
      C13, &
      C14, &
      C_wp2_pr_dfsn, &
      C_wp3_pr_tp, &
      C_wp3_pr_turb, &
      C_wp3_pr_dfsn, &
      C_wp2_splat, &
      C6rt_Lscale0, &
      C6thl_Lscale0, &
      C7_Lscale0, &
      wpxp_L_thresh, &
      c_K, &
      c_K1, &
      nu1, &
      c_K2, &
      nu2, &
      c_K6, &
      nu6, &
      c_K8, &
      nu8, &
      c_K9, &
      nu9, &
      nu10, &
      c_K_hm, &
      c_K_hmb, &
      K_hm_min_coef, &
      nu_hm, &
      slope_coef_spread_DG_means_w, &
      pdf_component_stdev_factor_w, &
      coef_spread_DG_means_rt, &
      coef_spread_DG_means_thl, &
      beta, &
      gamma_coef, &
      gamma_coefb, &
      gamma_coefc, &
      lmin_coef, &
      omicron, &
      zeta_vrnce_rat, &
      upsilon_precip_frac_rat, &
      lambda0_stability_coef, &
      mult_coef, &
      taumin, &
      taumax, &
      Lscale_mu_entrain, &
      Lscale_mu_coef, &
      Lscale_pert_coef, &
      alpha_corr, &
      Skw_denom_coef, &
      c_K10, &
      c_K10h, &
      thlp2_rad_coef, &
      thlp2_rad_cloud_frac_thresh, &
      up2_sfc_coef, &
      Skw_max_mag, &
      xp3_coef_base, &
      xp3_coef_slope, &
      altitude_threshold, &
      rtp2_clip_coef, &
      C_invrs_tau_bkgnd, &
      C_invrs_tau_sfc, &
      C_invrs_tau_shear, &
      C_invrs_tau_N2, &
      C_invrs_tau_N2_wp2, &
      C_invrs_tau_N2_xp2, &
      C_invrs_tau_N2_wpxp, &
      C_invrs_tau_N2_clear_wp3, &
      C_invrs_tau_wpxp_Ri, &
      C_invrs_tau_wpxp_N2_thresh, &
      Cx_min, &
      Cx_max, &
      Richardson_num_min, &
      Richardson_num_max, &
      wpxp_Ri_exp, &
      a3_coef_min, &
      a_const, &
      bv_efold, &
      z_displace, &
      Lscale_tke_multiplier, &
      Lscale_pdf_launch_frac, &
      Lscale_up_weight
    integer(c_int64_t) :: extents(4)
    character(len=:), allocatable :: runfile
    integer(c_int) :: ngrdcol
    integer(c_int) :: iunit
    real(core_rknd), allocatable :: C1(:)
    real(core_rknd), allocatable :: C1b(:)
    real(core_rknd), allocatable :: C1c(:)
    real(core_rknd), allocatable :: C2rt(:)
    real(core_rknd), allocatable :: C2thl(:)
    real(core_rknd), allocatable :: C2rtthl(:)
    real(core_rknd), allocatable :: C4(:)
    real(core_rknd), allocatable :: C_uu_shr(:)
    real(core_rknd), allocatable :: C_uu_buoy(:)
    real(core_rknd), allocatable :: C6rt(:)
    real(core_rknd), allocatable :: C6rtb(:)
    real(core_rknd), allocatable :: C6rtc(:)
    real(core_rknd), allocatable :: C6thl(:)
    real(core_rknd), allocatable :: C6thlb(:)
    real(core_rknd), allocatable :: C6thlc(:)
    real(core_rknd), allocatable :: C7(:)
    real(core_rknd), allocatable :: C7b(:)
    real(core_rknd), allocatable :: C7c(:)
    real(core_rknd), allocatable :: C8(:)
    real(core_rknd), allocatable :: C8b(:)
    real(core_rknd), allocatable :: C10(:)
    real(core_rknd), allocatable :: C11(:)
    real(core_rknd), allocatable :: C11b(:)
    real(core_rknd), allocatable :: C11c(:)
    real(core_rknd), allocatable :: C12(:)
    real(core_rknd), allocatable :: C13(:)
    real(core_rknd), allocatable :: C14(:)
    real(core_rknd), allocatable :: C_wp2_pr_dfsn(:)
    real(core_rknd), allocatable :: C_wp3_pr_tp(:)
    real(core_rknd), allocatable :: C_wp3_pr_turb(:)
    real(core_rknd), allocatable :: C_wp3_pr_dfsn(:)
    real(core_rknd), allocatable :: C_wp2_splat(:)
    real(core_rknd), allocatable :: C6rt_Lscale0(:)
    real(core_rknd), allocatable :: C6thl_Lscale0(:)
    real(core_rknd), allocatable :: C7_Lscale0(:)
    real(core_rknd), allocatable :: wpxp_L_thresh(:)
    real(core_rknd), allocatable :: c_K(:)
    real(core_rknd), allocatable :: c_K1(:)
    real(core_rknd), allocatable :: nu1(:)
    real(core_rknd), allocatable :: c_K2(:)
    real(core_rknd), allocatable :: nu2(:)
    real(core_rknd), allocatable :: c_K6(:)
    real(core_rknd), allocatable :: nu6(:)
    real(core_rknd), allocatable :: c_K8(:)
    real(core_rknd), allocatable :: nu8(:)
    real(core_rknd), allocatable :: c_K9(:)
    real(core_rknd), allocatable :: nu9(:)
    real(core_rknd), allocatable :: nu10(:)
    real(core_rknd), allocatable :: c_K_hm(:)
    real(core_rknd), allocatable :: c_K_hmb(:)
    real(core_rknd), allocatable :: K_hm_min_coef(:)
    real(core_rknd), allocatable :: nu_hm(:)
    real(core_rknd), allocatable :: slope_coef_spread_DG_means_w(:)
    real(core_rknd), allocatable :: pdf_component_stdev_factor_w(:)
    real(core_rknd), allocatable :: coef_spread_DG_means_rt(:)
    real(core_rknd), allocatable :: coef_spread_DG_means_thl(:)
    real(core_rknd), allocatable :: beta(:)
    real(core_rknd), allocatable :: gamma_coef(:)
    real(core_rknd), allocatable :: gamma_coefb(:)
    real(core_rknd), allocatable :: gamma_coefc(:)
    real(core_rknd), allocatable :: lmin_coef(:)
    real(core_rknd), allocatable :: omicron(:)
    real(core_rknd), allocatable :: zeta_vrnce_rat(:)
    real(core_rknd), allocatable :: upsilon_precip_frac_rat(:)
    real(core_rknd), allocatable :: lambda0_stability_coef(:)
    real(core_rknd), allocatable :: mult_coef(:)
    real(core_rknd), allocatable :: taumin(:)
    real(core_rknd), allocatable :: taumax(:)
    real(core_rknd), allocatable :: Lscale_mu_entrain(:)
    real(core_rknd), allocatable :: Lscale_mu_coef(:)
    real(core_rknd), allocatable :: Lscale_pert_coef(:)
    real(core_rknd), allocatable :: alpha_corr(:)
    real(core_rknd), allocatable :: Skw_denom_coef(:)
    real(core_rknd), allocatable :: c_K10(:)
    real(core_rknd), allocatable :: c_K10h(:)
    real(core_rknd), allocatable :: thlp2_rad_coef(:)
    real(core_rknd), allocatable :: thlp2_rad_cloud_frac_thresh(:)
    real(core_rknd), allocatable :: up2_sfc_coef(:)
    real(core_rknd), allocatable :: Skw_max_mag(:)
    real(core_rknd), allocatable :: xp3_coef_base(:)
    real(core_rknd), allocatable :: xp3_coef_slope(:)
    real(core_rknd), allocatable :: altitude_threshold(:)
    real(core_rknd), allocatable :: rtp2_clip_coef(:)
    real(core_rknd), allocatable :: C_invrs_tau_bkgnd(:)
    real(core_rknd), allocatable :: C_invrs_tau_sfc(:)
    real(core_rknd), allocatable :: C_invrs_tau_shear(:)
    real(core_rknd), allocatable :: C_invrs_tau_N2(:)
    real(core_rknd), allocatable :: C_invrs_tau_N2_wp2(:)
    real(core_rknd), allocatable :: C_invrs_tau_N2_xp2(:)
    real(core_rknd), allocatable :: C_invrs_tau_N2_wpxp(:)
    real(core_rknd), allocatable :: C_invrs_tau_N2_clear_wp3(:)
    real(core_rknd), allocatable :: C_invrs_tau_wpxp_Ri(:)
    real(core_rknd), allocatable :: C_invrs_tau_wpxp_N2_thresh(:)
    real(core_rknd), allocatable :: Cx_min(:)
    real(core_rknd), allocatable :: Cx_max(:)
    real(core_rknd), allocatable :: Richardson_num_min(:)
    real(core_rknd), allocatable :: Richardson_num_max(:)
    real(core_rknd), allocatable :: wpxp_Ri_exp(:)
    real(core_rknd), allocatable :: a3_coef_min(:)
    real(core_rknd), allocatable :: a_const(:)
    real(core_rknd), allocatable :: bv_efold(:)
    real(core_rknd), allocatable :: z_displace(:)
    real(core_rknd), allocatable :: Lscale_tke_multiplier(:)
    real(core_rknd), allocatable :: Lscale_pdf_launch_frac(:)
    real(core_rknd), allocatable :: Lscale_up_weight(:)
    integer(c_int) :: read_status
    do j = 1, 109
      call empty_field(0, 0, & ! In
                       output(j)) ! Out
    end do
    allocate(character(len=max(1,input(1)%char_len)) :: runfile)
    runfile(:) = ' '
    if (input(1)%char_len > 0) then
      call c_f_pointer(input(1)%data,cp,[input(1)%char_len])
      do k = 1, input(1)%char_len
        runfile(k:k) = cp(k)
      end do
    end if
    call c_f_pointer(input(2)%data,ip,[1])
    ngrdcol = ip(1)
    call c_f_pointer(input(3)%data,ip,[1])
    iunit = ip(1)
    if (input(4)%allocated /= 0) then
      allocate(C1(int(input(4)%extent(1))))
      n = size(C1)
      if (n > 0) then
        call c_f_pointer(input(4)%data,rp,[n])
        C1 = reshape(rp,shape(C1))
      end if
    end if
    if (input(5)%allocated /= 0) then
      allocate(C1b(int(input(5)%extent(1))))
      n = size(C1b)
      if (n > 0) then
        call c_f_pointer(input(5)%data,rp,[n])
        C1b = reshape(rp,shape(C1b))
      end if
    end if
    if (input(6)%allocated /= 0) then
      allocate(C1c(int(input(6)%extent(1))))
      n = size(C1c)
      if (n > 0) then
        call c_f_pointer(input(6)%data,rp,[n])
        C1c = reshape(rp,shape(C1c))
      end if
    end if
    if (input(7)%allocated /= 0) then
      allocate(C2rt(int(input(7)%extent(1))))
      n = size(C2rt)
      if (n > 0) then
        call c_f_pointer(input(7)%data,rp,[n])
        C2rt = reshape(rp,shape(C2rt))
      end if
    end if
    if (input(8)%allocated /= 0) then
      allocate(C2thl(int(input(8)%extent(1))))
      n = size(C2thl)
      if (n > 0) then
        call c_f_pointer(input(8)%data,rp,[n])
        C2thl = reshape(rp,shape(C2thl))
      end if
    end if
    if (input(9)%allocated /= 0) then
      allocate(C2rtthl(int(input(9)%extent(1))))
      n = size(C2rtthl)
      if (n > 0) then
        call c_f_pointer(input(9)%data,rp,[n])
        C2rtthl = reshape(rp,shape(C2rtthl))
      end if
    end if
    if (input(10)%allocated /= 0) then
      allocate(C4(int(input(10)%extent(1))))
      n = size(C4)
      if (n > 0) then
        call c_f_pointer(input(10)%data,rp,[n])
        C4 = reshape(rp,shape(C4))
      end if
    end if
    if (input(11)%allocated /= 0) then
      allocate(C_uu_shr(int(input(11)%extent(1))))
      n = size(C_uu_shr)
      if (n > 0) then
        call c_f_pointer(input(11)%data,rp,[n])
        C_uu_shr = reshape(rp,shape(C_uu_shr))
      end if
    end if
    if (input(12)%allocated /= 0) then
      allocate(C_uu_buoy(int(input(12)%extent(1))))
      n = size(C_uu_buoy)
      if (n > 0) then
        call c_f_pointer(input(12)%data,rp,[n])
        C_uu_buoy = reshape(rp,shape(C_uu_buoy))
      end if
    end if
    if (input(13)%allocated /= 0) then
      allocate(C6rt(int(input(13)%extent(1))))
      n = size(C6rt)
      if (n > 0) then
        call c_f_pointer(input(13)%data,rp,[n])
        C6rt = reshape(rp,shape(C6rt))
      end if
    end if
    if (input(14)%allocated /= 0) then
      allocate(C6rtb(int(input(14)%extent(1))))
      n = size(C6rtb)
      if (n > 0) then
        call c_f_pointer(input(14)%data,rp,[n])
        C6rtb = reshape(rp,shape(C6rtb))
      end if
    end if
    if (input(15)%allocated /= 0) then
      allocate(C6rtc(int(input(15)%extent(1))))
      n = size(C6rtc)
      if (n > 0) then
        call c_f_pointer(input(15)%data,rp,[n])
        C6rtc = reshape(rp,shape(C6rtc))
      end if
    end if
    if (input(16)%allocated /= 0) then
      allocate(C6thl(int(input(16)%extent(1))))
      n = size(C6thl)
      if (n > 0) then
        call c_f_pointer(input(16)%data,rp,[n])
        C6thl = reshape(rp,shape(C6thl))
      end if
    end if
    if (input(17)%allocated /= 0) then
      allocate(C6thlb(int(input(17)%extent(1))))
      n = size(C6thlb)
      if (n > 0) then
        call c_f_pointer(input(17)%data,rp,[n])
        C6thlb = reshape(rp,shape(C6thlb))
      end if
    end if
    if (input(18)%allocated /= 0) then
      allocate(C6thlc(int(input(18)%extent(1))))
      n = size(C6thlc)
      if (n > 0) then
        call c_f_pointer(input(18)%data,rp,[n])
        C6thlc = reshape(rp,shape(C6thlc))
      end if
    end if
    if (input(19)%allocated /= 0) then
      allocate(C7(int(input(19)%extent(1))))
      n = size(C7)
      if (n > 0) then
        call c_f_pointer(input(19)%data,rp,[n])
        C7 = reshape(rp,shape(C7))
      end if
    end if
    if (input(20)%allocated /= 0) then
      allocate(C7b(int(input(20)%extent(1))))
      n = size(C7b)
      if (n > 0) then
        call c_f_pointer(input(20)%data,rp,[n])
        C7b = reshape(rp,shape(C7b))
      end if
    end if
    if (input(21)%allocated /= 0) then
      allocate(C7c(int(input(21)%extent(1))))
      n = size(C7c)
      if (n > 0) then
        call c_f_pointer(input(21)%data,rp,[n])
        C7c = reshape(rp,shape(C7c))
      end if
    end if
    if (input(22)%allocated /= 0) then
      allocate(C8(int(input(22)%extent(1))))
      n = size(C8)
      if (n > 0) then
        call c_f_pointer(input(22)%data,rp,[n])
        C8 = reshape(rp,shape(C8))
      end if
    end if
    if (input(23)%allocated /= 0) then
      allocate(C8b(int(input(23)%extent(1))))
      n = size(C8b)
      if (n > 0) then
        call c_f_pointer(input(23)%data,rp,[n])
        C8b = reshape(rp,shape(C8b))
      end if
    end if
    if (input(24)%allocated /= 0) then
      allocate(C10(int(input(24)%extent(1))))
      n = size(C10)
      if (n > 0) then
        call c_f_pointer(input(24)%data,rp,[n])
        C10 = reshape(rp,shape(C10))
      end if
    end if
    if (input(25)%allocated /= 0) then
      allocate(C11(int(input(25)%extent(1))))
      n = size(C11)
      if (n > 0) then
        call c_f_pointer(input(25)%data,rp,[n])
        C11 = reshape(rp,shape(C11))
      end if
    end if
    if (input(26)%allocated /= 0) then
      allocate(C11b(int(input(26)%extent(1))))
      n = size(C11b)
      if (n > 0) then
        call c_f_pointer(input(26)%data,rp,[n])
        C11b = reshape(rp,shape(C11b))
      end if
    end if
    if (input(27)%allocated /= 0) then
      allocate(C11c(int(input(27)%extent(1))))
      n = size(C11c)
      if (n > 0) then
        call c_f_pointer(input(27)%data,rp,[n])
        C11c = reshape(rp,shape(C11c))
      end if
    end if
    if (input(28)%allocated /= 0) then
      allocate(C12(int(input(28)%extent(1))))
      n = size(C12)
      if (n > 0) then
        call c_f_pointer(input(28)%data,rp,[n])
        C12 = reshape(rp,shape(C12))
      end if
    end if
    if (input(29)%allocated /= 0) then
      allocate(C13(int(input(29)%extent(1))))
      n = size(C13)
      if (n > 0) then
        call c_f_pointer(input(29)%data,rp,[n])
        C13 = reshape(rp,shape(C13))
      end if
    end if
    if (input(30)%allocated /= 0) then
      allocate(C14(int(input(30)%extent(1))))
      n = size(C14)
      if (n > 0) then
        call c_f_pointer(input(30)%data,rp,[n])
        C14 = reshape(rp,shape(C14))
      end if
    end if
    if (input(31)%allocated /= 0) then
      allocate(C_wp2_pr_dfsn(int(input(31)%extent(1))))
      n = size(C_wp2_pr_dfsn)
      if (n > 0) then
        call c_f_pointer(input(31)%data,rp,[n])
        C_wp2_pr_dfsn = reshape(rp,shape(C_wp2_pr_dfsn))
      end if
    end if
    if (input(32)%allocated /= 0) then
      allocate(C_wp3_pr_tp(int(input(32)%extent(1))))
      n = size(C_wp3_pr_tp)
      if (n > 0) then
        call c_f_pointer(input(32)%data,rp,[n])
        C_wp3_pr_tp = reshape(rp,shape(C_wp3_pr_tp))
      end if
    end if
    if (input(33)%allocated /= 0) then
      allocate(C_wp3_pr_turb(int(input(33)%extent(1))))
      n = size(C_wp3_pr_turb)
      if (n > 0) then
        call c_f_pointer(input(33)%data,rp,[n])
        C_wp3_pr_turb = reshape(rp,shape(C_wp3_pr_turb))
      end if
    end if
    if (input(34)%allocated /= 0) then
      allocate(C_wp3_pr_dfsn(int(input(34)%extent(1))))
      n = size(C_wp3_pr_dfsn)
      if (n > 0) then
        call c_f_pointer(input(34)%data,rp,[n])
        C_wp3_pr_dfsn = reshape(rp,shape(C_wp3_pr_dfsn))
      end if
    end if
    if (input(35)%allocated /= 0) then
      allocate(C_wp2_splat(int(input(35)%extent(1))))
      n = size(C_wp2_splat)
      if (n > 0) then
        call c_f_pointer(input(35)%data,rp,[n])
        C_wp2_splat = reshape(rp,shape(C_wp2_splat))
      end if
    end if
    if (input(36)%allocated /= 0) then
      allocate(C6rt_Lscale0(int(input(36)%extent(1))))
      n = size(C6rt_Lscale0)
      if (n > 0) then
        call c_f_pointer(input(36)%data,rp,[n])
        C6rt_Lscale0 = reshape(rp,shape(C6rt_Lscale0))
      end if
    end if
    if (input(37)%allocated /= 0) then
      allocate(C6thl_Lscale0(int(input(37)%extent(1))))
      n = size(C6thl_Lscale0)
      if (n > 0) then
        call c_f_pointer(input(37)%data,rp,[n])
        C6thl_Lscale0 = reshape(rp,shape(C6thl_Lscale0))
      end if
    end if
    if (input(38)%allocated /= 0) then
      allocate(C7_Lscale0(int(input(38)%extent(1))))
      n = size(C7_Lscale0)
      if (n > 0) then
        call c_f_pointer(input(38)%data,rp,[n])
        C7_Lscale0 = reshape(rp,shape(C7_Lscale0))
      end if
    end if
    if (input(39)%allocated /= 0) then
      allocate(wpxp_L_thresh(int(input(39)%extent(1))))
      n = size(wpxp_L_thresh)
      if (n > 0) then
        call c_f_pointer(input(39)%data,rp,[n])
        wpxp_L_thresh = reshape(rp,shape(wpxp_L_thresh))
      end if
    end if
    if (input(40)%allocated /= 0) then
      allocate(c_K(int(input(40)%extent(1))))
      n = size(c_K)
      if (n > 0) then
        call c_f_pointer(input(40)%data,rp,[n])
        c_K = reshape(rp,shape(c_K))
      end if
    end if
    if (input(41)%allocated /= 0) then
      allocate(c_K1(int(input(41)%extent(1))))
      n = size(c_K1)
      if (n > 0) then
        call c_f_pointer(input(41)%data,rp,[n])
        c_K1 = reshape(rp,shape(c_K1))
      end if
    end if
    if (input(42)%allocated /= 0) then
      allocate(nu1(int(input(42)%extent(1))))
      n = size(nu1)
      if (n > 0) then
        call c_f_pointer(input(42)%data,rp,[n])
        nu1 = reshape(rp,shape(nu1))
      end if
    end if
    if (input(43)%allocated /= 0) then
      allocate(c_K2(int(input(43)%extent(1))))
      n = size(c_K2)
      if (n > 0) then
        call c_f_pointer(input(43)%data,rp,[n])
        c_K2 = reshape(rp,shape(c_K2))
      end if
    end if
    if (input(44)%allocated /= 0) then
      allocate(nu2(int(input(44)%extent(1))))
      n = size(nu2)
      if (n > 0) then
        call c_f_pointer(input(44)%data,rp,[n])
        nu2 = reshape(rp,shape(nu2))
      end if
    end if
    if (input(45)%allocated /= 0) then
      allocate(c_K6(int(input(45)%extent(1))))
      n = size(c_K6)
      if (n > 0) then
        call c_f_pointer(input(45)%data,rp,[n])
        c_K6 = reshape(rp,shape(c_K6))
      end if
    end if
    if (input(46)%allocated /= 0) then
      allocate(nu6(int(input(46)%extent(1))))
      n = size(nu6)
      if (n > 0) then
        call c_f_pointer(input(46)%data,rp,[n])
        nu6 = reshape(rp,shape(nu6))
      end if
    end if
    if (input(47)%allocated /= 0) then
      allocate(c_K8(int(input(47)%extent(1))))
      n = size(c_K8)
      if (n > 0) then
        call c_f_pointer(input(47)%data,rp,[n])
        c_K8 = reshape(rp,shape(c_K8))
      end if
    end if
    if (input(48)%allocated /= 0) then
      allocate(nu8(int(input(48)%extent(1))))
      n = size(nu8)
      if (n > 0) then
        call c_f_pointer(input(48)%data,rp,[n])
        nu8 = reshape(rp,shape(nu8))
      end if
    end if
    if (input(49)%allocated /= 0) then
      allocate(c_K9(int(input(49)%extent(1))))
      n = size(c_K9)
      if (n > 0) then
        call c_f_pointer(input(49)%data,rp,[n])
        c_K9 = reshape(rp,shape(c_K9))
      end if
    end if
    if (input(50)%allocated /= 0) then
      allocate(nu9(int(input(50)%extent(1))))
      n = size(nu9)
      if (n > 0) then
        call c_f_pointer(input(50)%data,rp,[n])
        nu9 = reshape(rp,shape(nu9))
      end if
    end if
    if (input(51)%allocated /= 0) then
      allocate(nu10(int(input(51)%extent(1))))
      n = size(nu10)
      if (n > 0) then
        call c_f_pointer(input(51)%data,rp,[n])
        nu10 = reshape(rp,shape(nu10))
      end if
    end if
    if (input(52)%allocated /= 0) then
      allocate(c_K_hm(int(input(52)%extent(1))))
      n = size(c_K_hm)
      if (n > 0) then
        call c_f_pointer(input(52)%data,rp,[n])
        c_K_hm = reshape(rp,shape(c_K_hm))
      end if
    end if
    if (input(53)%allocated /= 0) then
      allocate(c_K_hmb(int(input(53)%extent(1))))
      n = size(c_K_hmb)
      if (n > 0) then
        call c_f_pointer(input(53)%data,rp,[n])
        c_K_hmb = reshape(rp,shape(c_K_hmb))
      end if
    end if
    if (input(54)%allocated /= 0) then
      allocate(K_hm_min_coef(int(input(54)%extent(1))))
      n = size(K_hm_min_coef)
      if (n > 0) then
        call c_f_pointer(input(54)%data,rp,[n])
        K_hm_min_coef = reshape(rp,shape(K_hm_min_coef))
      end if
    end if
    if (input(55)%allocated /= 0) then
      allocate(nu_hm(int(input(55)%extent(1))))
      n = size(nu_hm)
      if (n > 0) then
        call c_f_pointer(input(55)%data,rp,[n])
        nu_hm = reshape(rp,shape(nu_hm))
      end if
    end if
    if (input(56)%allocated /= 0) then
      allocate(slope_coef_spread_DG_means_w(int(input(56)%extent(1))))
      n = size(slope_coef_spread_DG_means_w)
      if (n > 0) then
        call c_f_pointer(input(56)%data,rp,[n])
        slope_coef_spread_DG_means_w = reshape(rp,shape(slope_coef_spread_DG_means_w))
      end if
    end if
    if (input(57)%allocated /= 0) then
      allocate(pdf_component_stdev_factor_w(int(input(57)%extent(1))))
      n = size(pdf_component_stdev_factor_w)
      if (n > 0) then
        call c_f_pointer(input(57)%data,rp,[n])
        pdf_component_stdev_factor_w = reshape(rp,shape(pdf_component_stdev_factor_w))
      end if
    end if
    if (input(58)%allocated /= 0) then
      allocate(coef_spread_DG_means_rt(int(input(58)%extent(1))))
      n = size(coef_spread_DG_means_rt)
      if (n > 0) then
        call c_f_pointer(input(58)%data,rp,[n])
        coef_spread_DG_means_rt = reshape(rp,shape(coef_spread_DG_means_rt))
      end if
    end if
    if (input(59)%allocated /= 0) then
      allocate(coef_spread_DG_means_thl(int(input(59)%extent(1))))
      n = size(coef_spread_DG_means_thl)
      if (n > 0) then
        call c_f_pointer(input(59)%data,rp,[n])
        coef_spread_DG_means_thl = reshape(rp,shape(coef_spread_DG_means_thl))
      end if
    end if
    if (input(60)%allocated /= 0) then
      allocate(beta(int(input(60)%extent(1))))
      n = size(beta)
      if (n > 0) then
        call c_f_pointer(input(60)%data,rp,[n])
        beta = reshape(rp,shape(beta))
      end if
    end if
    if (input(61)%allocated /= 0) then
      allocate(gamma_coef(int(input(61)%extent(1))))
      n = size(gamma_coef)
      if (n > 0) then
        call c_f_pointer(input(61)%data,rp,[n])
        gamma_coef = reshape(rp,shape(gamma_coef))
      end if
    end if
    if (input(62)%allocated /= 0) then
      allocate(gamma_coefb(int(input(62)%extent(1))))
      n = size(gamma_coefb)
      if (n > 0) then
        call c_f_pointer(input(62)%data,rp,[n])
        gamma_coefb = reshape(rp,shape(gamma_coefb))
      end if
    end if
    if (input(63)%allocated /= 0) then
      allocate(gamma_coefc(int(input(63)%extent(1))))
      n = size(gamma_coefc)
      if (n > 0) then
        call c_f_pointer(input(63)%data,rp,[n])
        gamma_coefc = reshape(rp,shape(gamma_coefc))
      end if
    end if
    if (input(64)%allocated /= 0) then
      allocate(lmin_coef(int(input(64)%extent(1))))
      n = size(lmin_coef)
      if (n > 0) then
        call c_f_pointer(input(64)%data,rp,[n])
        lmin_coef = reshape(rp,shape(lmin_coef))
      end if
    end if
    if (input(65)%allocated /= 0) then
      allocate(omicron(int(input(65)%extent(1))))
      n = size(omicron)
      if (n > 0) then
        call c_f_pointer(input(65)%data,rp,[n])
        omicron = reshape(rp,shape(omicron))
      end if
    end if
    if (input(66)%allocated /= 0) then
      allocate(zeta_vrnce_rat(int(input(66)%extent(1))))
      n = size(zeta_vrnce_rat)
      if (n > 0) then
        call c_f_pointer(input(66)%data,rp,[n])
        zeta_vrnce_rat = reshape(rp,shape(zeta_vrnce_rat))
      end if
    end if
    if (input(67)%allocated /= 0) then
      allocate(upsilon_precip_frac_rat(int(input(67)%extent(1))))
      n = size(upsilon_precip_frac_rat)
      if (n > 0) then
        call c_f_pointer(input(67)%data,rp,[n])
        upsilon_precip_frac_rat = reshape(rp,shape(upsilon_precip_frac_rat))
      end if
    end if
    if (input(68)%allocated /= 0) then
      allocate(lambda0_stability_coef(int(input(68)%extent(1))))
      n = size(lambda0_stability_coef)
      if (n > 0) then
        call c_f_pointer(input(68)%data,rp,[n])
        lambda0_stability_coef = reshape(rp,shape(lambda0_stability_coef))
      end if
    end if
    if (input(69)%allocated /= 0) then
      allocate(mult_coef(int(input(69)%extent(1))))
      n = size(mult_coef)
      if (n > 0) then
        call c_f_pointer(input(69)%data,rp,[n])
        mult_coef = reshape(rp,shape(mult_coef))
      end if
    end if
    if (input(70)%allocated /= 0) then
      allocate(taumin(int(input(70)%extent(1))))
      n = size(taumin)
      if (n > 0) then
        call c_f_pointer(input(70)%data,rp,[n])
        taumin = reshape(rp,shape(taumin))
      end if
    end if
    if (input(71)%allocated /= 0) then
      allocate(taumax(int(input(71)%extent(1))))
      n = size(taumax)
      if (n > 0) then
        call c_f_pointer(input(71)%data,rp,[n])
        taumax = reshape(rp,shape(taumax))
      end if
    end if
    if (input(72)%allocated /= 0) then
      allocate(Lscale_mu_entrain(int(input(72)%extent(1))))
      n = size(Lscale_mu_entrain)
      if (n > 0) then
        call c_f_pointer(input(72)%data,rp,[n])
        Lscale_mu_entrain = reshape(rp,shape(Lscale_mu_entrain))
      end if
    end if
    if (input(73)%allocated /= 0) then
      allocate(Lscale_mu_coef(int(input(73)%extent(1))))
      n = size(Lscale_mu_coef)
      if (n > 0) then
        call c_f_pointer(input(73)%data,rp,[n])
        Lscale_mu_coef = reshape(rp,shape(Lscale_mu_coef))
      end if
    end if
    if (input(74)%allocated /= 0) then
      allocate(Lscale_pert_coef(int(input(74)%extent(1))))
      n = size(Lscale_pert_coef)
      if (n > 0) then
        call c_f_pointer(input(74)%data,rp,[n])
        Lscale_pert_coef = reshape(rp,shape(Lscale_pert_coef))
      end if
    end if
    if (input(75)%allocated /= 0) then
      allocate(alpha_corr(int(input(75)%extent(1))))
      n = size(alpha_corr)
      if (n > 0) then
        call c_f_pointer(input(75)%data,rp,[n])
        alpha_corr = reshape(rp,shape(alpha_corr))
      end if
    end if
    if (input(76)%allocated /= 0) then
      allocate(Skw_denom_coef(int(input(76)%extent(1))))
      n = size(Skw_denom_coef)
      if (n > 0) then
        call c_f_pointer(input(76)%data,rp,[n])
        Skw_denom_coef = reshape(rp,shape(Skw_denom_coef))
      end if
    end if
    if (input(77)%allocated /= 0) then
      allocate(c_K10(int(input(77)%extent(1))))
      n = size(c_K10)
      if (n > 0) then
        call c_f_pointer(input(77)%data,rp,[n])
        c_K10 = reshape(rp,shape(c_K10))
      end if
    end if
    if (input(78)%allocated /= 0) then
      allocate(c_K10h(int(input(78)%extent(1))))
      n = size(c_K10h)
      if (n > 0) then
        call c_f_pointer(input(78)%data,rp,[n])
        c_K10h = reshape(rp,shape(c_K10h))
      end if
    end if
    if (input(79)%allocated /= 0) then
      allocate(thlp2_rad_coef(int(input(79)%extent(1))))
      n = size(thlp2_rad_coef)
      if (n > 0) then
        call c_f_pointer(input(79)%data,rp,[n])
        thlp2_rad_coef = reshape(rp,shape(thlp2_rad_coef))
      end if
    end if
    if (input(80)%allocated /= 0) then
      allocate(thlp2_rad_cloud_frac_thresh(int(input(80)%extent(1))))
      n = size(thlp2_rad_cloud_frac_thresh)
      if (n > 0) then
        call c_f_pointer(input(80)%data,rp,[n])
        thlp2_rad_cloud_frac_thresh = reshape(rp,shape(thlp2_rad_cloud_frac_thresh))
      end if
    end if
    if (input(81)%allocated /= 0) then
      allocate(up2_sfc_coef(int(input(81)%extent(1))))
      n = size(up2_sfc_coef)
      if (n > 0) then
        call c_f_pointer(input(81)%data,rp,[n])
        up2_sfc_coef = reshape(rp,shape(up2_sfc_coef))
      end if
    end if
    if (input(82)%allocated /= 0) then
      allocate(Skw_max_mag(int(input(82)%extent(1))))
      n = size(Skw_max_mag)
      if (n > 0) then
        call c_f_pointer(input(82)%data,rp,[n])
        Skw_max_mag = reshape(rp,shape(Skw_max_mag))
      end if
    end if
    if (input(83)%allocated /= 0) then
      allocate(xp3_coef_base(int(input(83)%extent(1))))
      n = size(xp3_coef_base)
      if (n > 0) then
        call c_f_pointer(input(83)%data,rp,[n])
        xp3_coef_base = reshape(rp,shape(xp3_coef_base))
      end if
    end if
    if (input(84)%allocated /= 0) then
      allocate(xp3_coef_slope(int(input(84)%extent(1))))
      n = size(xp3_coef_slope)
      if (n > 0) then
        call c_f_pointer(input(84)%data,rp,[n])
        xp3_coef_slope = reshape(rp,shape(xp3_coef_slope))
      end if
    end if
    if (input(85)%allocated /= 0) then
      allocate(altitude_threshold(int(input(85)%extent(1))))
      n = size(altitude_threshold)
      if (n > 0) then
        call c_f_pointer(input(85)%data,rp,[n])
        altitude_threshold = reshape(rp,shape(altitude_threshold))
      end if
    end if
    if (input(86)%allocated /= 0) then
      allocate(rtp2_clip_coef(int(input(86)%extent(1))))
      n = size(rtp2_clip_coef)
      if (n > 0) then
        call c_f_pointer(input(86)%data,rp,[n])
        rtp2_clip_coef = reshape(rp,shape(rtp2_clip_coef))
      end if
    end if
    if (input(87)%allocated /= 0) then
      allocate(C_invrs_tau_bkgnd(int(input(87)%extent(1))))
      n = size(C_invrs_tau_bkgnd)
      if (n > 0) then
        call c_f_pointer(input(87)%data,rp,[n])
        C_invrs_tau_bkgnd = reshape(rp,shape(C_invrs_tau_bkgnd))
      end if
    end if
    if (input(88)%allocated /= 0) then
      allocate(C_invrs_tau_sfc(int(input(88)%extent(1))))
      n = size(C_invrs_tau_sfc)
      if (n > 0) then
        call c_f_pointer(input(88)%data,rp,[n])
        C_invrs_tau_sfc = reshape(rp,shape(C_invrs_tau_sfc))
      end if
    end if
    if (input(89)%allocated /= 0) then
      allocate(C_invrs_tau_shear(int(input(89)%extent(1))))
      n = size(C_invrs_tau_shear)
      if (n > 0) then
        call c_f_pointer(input(89)%data,rp,[n])
        C_invrs_tau_shear = reshape(rp,shape(C_invrs_tau_shear))
      end if
    end if
    if (input(90)%allocated /= 0) then
      allocate(C_invrs_tau_N2(int(input(90)%extent(1))))
      n = size(C_invrs_tau_N2)
      if (n > 0) then
        call c_f_pointer(input(90)%data,rp,[n])
        C_invrs_tau_N2 = reshape(rp,shape(C_invrs_tau_N2))
      end if
    end if
    if (input(91)%allocated /= 0) then
      allocate(C_invrs_tau_N2_wp2(int(input(91)%extent(1))))
      n = size(C_invrs_tau_N2_wp2)
      if (n > 0) then
        call c_f_pointer(input(91)%data,rp,[n])
        C_invrs_tau_N2_wp2 = reshape(rp,shape(C_invrs_tau_N2_wp2))
      end if
    end if
    if (input(92)%allocated /= 0) then
      allocate(C_invrs_tau_N2_xp2(int(input(92)%extent(1))))
      n = size(C_invrs_tau_N2_xp2)
      if (n > 0) then
        call c_f_pointer(input(92)%data,rp,[n])
        C_invrs_tau_N2_xp2 = reshape(rp,shape(C_invrs_tau_N2_xp2))
      end if
    end if
    if (input(93)%allocated /= 0) then
      allocate(C_invrs_tau_N2_wpxp(int(input(93)%extent(1))))
      n = size(C_invrs_tau_N2_wpxp)
      if (n > 0) then
        call c_f_pointer(input(93)%data,rp,[n])
        C_invrs_tau_N2_wpxp = reshape(rp,shape(C_invrs_tau_N2_wpxp))
      end if
    end if
    if (input(94)%allocated /= 0) then
      allocate(C_invrs_tau_N2_clear_wp3(int(input(94)%extent(1))))
      n = size(C_invrs_tau_N2_clear_wp3)
      if (n > 0) then
        call c_f_pointer(input(94)%data,rp,[n])
        C_invrs_tau_N2_clear_wp3 = reshape(rp,shape(C_invrs_tau_N2_clear_wp3))
      end if
    end if
    if (input(95)%allocated /= 0) then
      allocate(C_invrs_tau_wpxp_Ri(int(input(95)%extent(1))))
      n = size(C_invrs_tau_wpxp_Ri)
      if (n > 0) then
        call c_f_pointer(input(95)%data,rp,[n])
        C_invrs_tau_wpxp_Ri = reshape(rp,shape(C_invrs_tau_wpxp_Ri))
      end if
    end if
    if (input(96)%allocated /= 0) then
      allocate(C_invrs_tau_wpxp_N2_thresh(int(input(96)%extent(1))))
      n = size(C_invrs_tau_wpxp_N2_thresh)
      if (n > 0) then
        call c_f_pointer(input(96)%data,rp,[n])
        C_invrs_tau_wpxp_N2_thresh = reshape(rp,shape(C_invrs_tau_wpxp_N2_thresh))
      end if
    end if
    if (input(97)%allocated /= 0) then
      allocate(Cx_min(int(input(97)%extent(1))))
      n = size(Cx_min)
      if (n > 0) then
        call c_f_pointer(input(97)%data,rp,[n])
        Cx_min = reshape(rp,shape(Cx_min))
      end if
    end if
    if (input(98)%allocated /= 0) then
      allocate(Cx_max(int(input(98)%extent(1))))
      n = size(Cx_max)
      if (n > 0) then
        call c_f_pointer(input(98)%data,rp,[n])
        Cx_max = reshape(rp,shape(Cx_max))
      end if
    end if
    if (input(99)%allocated /= 0) then
      allocate(Richardson_num_min(int(input(99)%extent(1))))
      n = size(Richardson_num_min)
      if (n > 0) then
        call c_f_pointer(input(99)%data,rp,[n])
        Richardson_num_min = reshape(rp,shape(Richardson_num_min))
      end if
    end if
    if (input(100)%allocated /= 0) then
      allocate(Richardson_num_max(int(input(100)%extent(1))))
      n = size(Richardson_num_max)
      if (n > 0) then
        call c_f_pointer(input(100)%data,rp,[n])
        Richardson_num_max = reshape(rp,shape(Richardson_num_max))
      end if
    end if
    if (input(101)%allocated /= 0) then
      allocate(wpxp_Ri_exp(int(input(101)%extent(1))))
      n = size(wpxp_Ri_exp)
      if (n > 0) then
        call c_f_pointer(input(101)%data,rp,[n])
        wpxp_Ri_exp = reshape(rp,shape(wpxp_Ri_exp))
      end if
    end if
    if (input(102)%allocated /= 0) then
      allocate(a3_coef_min(int(input(102)%extent(1))))
      n = size(a3_coef_min)
      if (n > 0) then
        call c_f_pointer(input(102)%data,rp,[n])
        a3_coef_min = reshape(rp,shape(a3_coef_min))
      end if
    end if
    if (input(103)%allocated /= 0) then
      allocate(a_const(int(input(103)%extent(1))))
      n = size(a_const)
      if (n > 0) then
        call c_f_pointer(input(103)%data,rp,[n])
        a_const = reshape(rp,shape(a_const))
      end if
    end if
    if (input(104)%allocated /= 0) then
      allocate(bv_efold(int(input(104)%extent(1))))
      n = size(bv_efold)
      if (n > 0) then
        call c_f_pointer(input(104)%data,rp,[n])
        bv_efold = reshape(rp,shape(bv_efold))
      end if
    end if
    if (input(105)%allocated /= 0) then
      allocate(z_displace(int(input(105)%extent(1))))
      n = size(z_displace)
      if (n > 0) then
        call c_f_pointer(input(105)%data,rp,[n])
        z_displace = reshape(rp,shape(z_displace))
      end if
    end if
    if (input(106)%allocated /= 0) then
      allocate(Lscale_tke_multiplier(int(input(106)%extent(1))))
      n = size(Lscale_tke_multiplier)
      if (n > 0) then
        call c_f_pointer(input(106)%data,rp,[n])
        Lscale_tke_multiplier = reshape(rp,shape(Lscale_tke_multiplier))
      end if
    end if
    if (input(107)%allocated /= 0) then
      allocate(Lscale_pdf_launch_frac(int(input(107)%extent(1))))
      n = size(Lscale_pdf_launch_frac)
      if (n > 0) then
        call c_f_pointer(input(107)%data,rp,[n])
        Lscale_pdf_launch_frac = reshape(rp,shape(Lscale_pdf_launch_frac))
      end if
    end if
    if (input(108)%allocated /= 0) then
      allocate(Lscale_up_weight(int(input(108)%extent(1))))
      n = size(Lscale_up_weight)
      if (n > 0) then
        call c_f_pointer(input(108)%data,rp,[n])
        Lscale_up_weight = reshape(rp,shape(Lscale_up_weight))
      end if
    end if
    open(newunit=unit,file=trim(runfile),status="old",action="read",iostat=read_status)
    if (read_status == 0) then
      read(unit,nml=clubb_params_nl,iostat=read_status)
      close(unit)
    end if
    if (allocated(C1)) then
      extents = 1_c_int64_t
      extents(1) = size(C1,1,kind=c_int64_t)
      n = size(C1)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(4)) ! Out
      if (n > 0) then
        call c_f_pointer(output(4)%data,rp,[n])
        rp = reshape(C1,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(4)) ! Out
    end if
    if (allocated(C1b)) then
      extents = 1_c_int64_t
      extents(1) = size(C1b,1,kind=c_int64_t)
      n = size(C1b)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(5)) ! Out
      if (n > 0) then
        call c_f_pointer(output(5)%data,rp,[n])
        rp = reshape(C1b,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(5)) ! Out
    end if
    if (allocated(C1c)) then
      extents = 1_c_int64_t
      extents(1) = size(C1c,1,kind=c_int64_t)
      n = size(C1c)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(6)) ! Out
      if (n > 0) then
        call c_f_pointer(output(6)%data,rp,[n])
        rp = reshape(C1c,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(6)) ! Out
    end if
    if (allocated(C2rt)) then
      extents = 1_c_int64_t
      extents(1) = size(C2rt,1,kind=c_int64_t)
      n = size(C2rt)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(7)) ! Out
      if (n > 0) then
        call c_f_pointer(output(7)%data,rp,[n])
        rp = reshape(C2rt,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(7)) ! Out
    end if
    if (allocated(C2thl)) then
      extents = 1_c_int64_t
      extents(1) = size(C2thl,1,kind=c_int64_t)
      n = size(C2thl)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(8)) ! Out
      if (n > 0) then
        call c_f_pointer(output(8)%data,rp,[n])
        rp = reshape(C2thl,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(8)) ! Out
    end if
    if (allocated(C2rtthl)) then
      extents = 1_c_int64_t
      extents(1) = size(C2rtthl,1,kind=c_int64_t)
      n = size(C2rtthl)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(9)) ! Out
      if (n > 0) then
        call c_f_pointer(output(9)%data,rp,[n])
        rp = reshape(C2rtthl,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(9)) ! Out
    end if
    if (allocated(C4)) then
      extents = 1_c_int64_t
      extents(1) = size(C4,1,kind=c_int64_t)
      n = size(C4)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(10)) ! Out
      if (n > 0) then
        call c_f_pointer(output(10)%data,rp,[n])
        rp = reshape(C4,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(10)) ! Out
    end if
    if (allocated(C_uu_shr)) then
      extents = 1_c_int64_t
      extents(1) = size(C_uu_shr,1,kind=c_int64_t)
      n = size(C_uu_shr)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(11)) ! Out
      if (n > 0) then
        call c_f_pointer(output(11)%data,rp,[n])
        rp = reshape(C_uu_shr,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(11)) ! Out
    end if
    if (allocated(C_uu_buoy)) then
      extents = 1_c_int64_t
      extents(1) = size(C_uu_buoy,1,kind=c_int64_t)
      n = size(C_uu_buoy)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(12)) ! Out
      if (n > 0) then
        call c_f_pointer(output(12)%data,rp,[n])
        rp = reshape(C_uu_buoy,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(12)) ! Out
    end if
    if (allocated(C6rt)) then
      extents = 1_c_int64_t
      extents(1) = size(C6rt,1,kind=c_int64_t)
      n = size(C6rt)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(13)) ! Out
      if (n > 0) then
        call c_f_pointer(output(13)%data,rp,[n])
        rp = reshape(C6rt,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(13)) ! Out
    end if
    if (allocated(C6rtb)) then
      extents = 1_c_int64_t
      extents(1) = size(C6rtb,1,kind=c_int64_t)
      n = size(C6rtb)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(14)) ! Out
      if (n > 0) then
        call c_f_pointer(output(14)%data,rp,[n])
        rp = reshape(C6rtb,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(14)) ! Out
    end if
    if (allocated(C6rtc)) then
      extents = 1_c_int64_t
      extents(1) = size(C6rtc,1,kind=c_int64_t)
      n = size(C6rtc)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(15)) ! Out
      if (n > 0) then
        call c_f_pointer(output(15)%data,rp,[n])
        rp = reshape(C6rtc,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(15)) ! Out
    end if
    if (allocated(C6thl)) then
      extents = 1_c_int64_t
      extents(1) = size(C6thl,1,kind=c_int64_t)
      n = size(C6thl)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(16)) ! Out
      if (n > 0) then
        call c_f_pointer(output(16)%data,rp,[n])
        rp = reshape(C6thl,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(16)) ! Out
    end if
    if (allocated(C6thlb)) then
      extents = 1_c_int64_t
      extents(1) = size(C6thlb,1,kind=c_int64_t)
      n = size(C6thlb)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(17)) ! Out
      if (n > 0) then
        call c_f_pointer(output(17)%data,rp,[n])
        rp = reshape(C6thlb,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(17)) ! Out
    end if
    if (allocated(C6thlc)) then
      extents = 1_c_int64_t
      extents(1) = size(C6thlc,1,kind=c_int64_t)
      n = size(C6thlc)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(18)) ! Out
      if (n > 0) then
        call c_f_pointer(output(18)%data,rp,[n])
        rp = reshape(C6thlc,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(18)) ! Out
    end if
    if (allocated(C7)) then
      extents = 1_c_int64_t
      extents(1) = size(C7,1,kind=c_int64_t)
      n = size(C7)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(19)) ! Out
      if (n > 0) then
        call c_f_pointer(output(19)%data,rp,[n])
        rp = reshape(C7,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(19)) ! Out
    end if
    if (allocated(C7b)) then
      extents = 1_c_int64_t
      extents(1) = size(C7b,1,kind=c_int64_t)
      n = size(C7b)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(20)) ! Out
      if (n > 0) then
        call c_f_pointer(output(20)%data,rp,[n])
        rp = reshape(C7b,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(20)) ! Out
    end if
    if (allocated(C7c)) then
      extents = 1_c_int64_t
      extents(1) = size(C7c,1,kind=c_int64_t)
      n = size(C7c)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(21)) ! Out
      if (n > 0) then
        call c_f_pointer(output(21)%data,rp,[n])
        rp = reshape(C7c,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(21)) ! Out
    end if
    if (allocated(C8)) then
      extents = 1_c_int64_t
      extents(1) = size(C8,1,kind=c_int64_t)
      n = size(C8)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(22)) ! Out
      if (n > 0) then
        call c_f_pointer(output(22)%data,rp,[n])
        rp = reshape(C8,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(22)) ! Out
    end if
    if (allocated(C8b)) then
      extents = 1_c_int64_t
      extents(1) = size(C8b,1,kind=c_int64_t)
      n = size(C8b)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(23)) ! Out
      if (n > 0) then
        call c_f_pointer(output(23)%data,rp,[n])
        rp = reshape(C8b,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(23)) ! Out
    end if
    if (allocated(C10)) then
      extents = 1_c_int64_t
      extents(1) = size(C10,1,kind=c_int64_t)
      n = size(C10)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(24)) ! Out
      if (n > 0) then
        call c_f_pointer(output(24)%data,rp,[n])
        rp = reshape(C10,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(24)) ! Out
    end if
    if (allocated(C11)) then
      extents = 1_c_int64_t
      extents(1) = size(C11,1,kind=c_int64_t)
      n = size(C11)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(25)) ! Out
      if (n > 0) then
        call c_f_pointer(output(25)%data,rp,[n])
        rp = reshape(C11,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(25)) ! Out
    end if
    if (allocated(C11b)) then
      extents = 1_c_int64_t
      extents(1) = size(C11b,1,kind=c_int64_t)
      n = size(C11b)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(26)) ! Out
      if (n > 0) then
        call c_f_pointer(output(26)%data,rp,[n])
        rp = reshape(C11b,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(26)) ! Out
    end if
    if (allocated(C11c)) then
      extents = 1_c_int64_t
      extents(1) = size(C11c,1,kind=c_int64_t)
      n = size(C11c)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(27)) ! Out
      if (n > 0) then
        call c_f_pointer(output(27)%data,rp,[n])
        rp = reshape(C11c,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(27)) ! Out
    end if
    if (allocated(C12)) then
      extents = 1_c_int64_t
      extents(1) = size(C12,1,kind=c_int64_t)
      n = size(C12)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(28)) ! Out
      if (n > 0) then
        call c_f_pointer(output(28)%data,rp,[n])
        rp = reshape(C12,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(28)) ! Out
    end if
    if (allocated(C13)) then
      extents = 1_c_int64_t
      extents(1) = size(C13,1,kind=c_int64_t)
      n = size(C13)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(29)) ! Out
      if (n > 0) then
        call c_f_pointer(output(29)%data,rp,[n])
        rp = reshape(C13,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(29)) ! Out
    end if
    if (allocated(C14)) then
      extents = 1_c_int64_t
      extents(1) = size(C14,1,kind=c_int64_t)
      n = size(C14)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(30)) ! Out
      if (n > 0) then
        call c_f_pointer(output(30)%data,rp,[n])
        rp = reshape(C14,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(30)) ! Out
    end if
    if (allocated(C_wp2_pr_dfsn)) then
      extents = 1_c_int64_t
      extents(1) = size(C_wp2_pr_dfsn,1,kind=c_int64_t)
      n = size(C_wp2_pr_dfsn)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(31)) ! Out
      if (n > 0) then
        call c_f_pointer(output(31)%data,rp,[n])
        rp = reshape(C_wp2_pr_dfsn,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(31)) ! Out
    end if
    if (allocated(C_wp3_pr_tp)) then
      extents = 1_c_int64_t
      extents(1) = size(C_wp3_pr_tp,1,kind=c_int64_t)
      n = size(C_wp3_pr_tp)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(32)) ! Out
      if (n > 0) then
        call c_f_pointer(output(32)%data,rp,[n])
        rp = reshape(C_wp3_pr_tp,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(32)) ! Out
    end if
    if (allocated(C_wp3_pr_turb)) then
      extents = 1_c_int64_t
      extents(1) = size(C_wp3_pr_turb,1,kind=c_int64_t)
      n = size(C_wp3_pr_turb)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(33)) ! Out
      if (n > 0) then
        call c_f_pointer(output(33)%data,rp,[n])
        rp = reshape(C_wp3_pr_turb,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(33)) ! Out
    end if
    if (allocated(C_wp3_pr_dfsn)) then
      extents = 1_c_int64_t
      extents(1) = size(C_wp3_pr_dfsn,1,kind=c_int64_t)
      n = size(C_wp3_pr_dfsn)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(34)) ! Out
      if (n > 0) then
        call c_f_pointer(output(34)%data,rp,[n])
        rp = reshape(C_wp3_pr_dfsn,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(34)) ! Out
    end if
    if (allocated(C_wp2_splat)) then
      extents = 1_c_int64_t
      extents(1) = size(C_wp2_splat,1,kind=c_int64_t)
      n = size(C_wp2_splat)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(35)) ! Out
      if (n > 0) then
        call c_f_pointer(output(35)%data,rp,[n])
        rp = reshape(C_wp2_splat,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(35)) ! Out
    end if
    if (allocated(C6rt_Lscale0)) then
      extents = 1_c_int64_t
      extents(1) = size(C6rt_Lscale0,1,kind=c_int64_t)
      n = size(C6rt_Lscale0)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(36)) ! Out
      if (n > 0) then
        call c_f_pointer(output(36)%data,rp,[n])
        rp = reshape(C6rt_Lscale0,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(36)) ! Out
    end if
    if (allocated(C6thl_Lscale0)) then
      extents = 1_c_int64_t
      extents(1) = size(C6thl_Lscale0,1,kind=c_int64_t)
      n = size(C6thl_Lscale0)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(37)) ! Out
      if (n > 0) then
        call c_f_pointer(output(37)%data,rp,[n])
        rp = reshape(C6thl_Lscale0,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(37)) ! Out
    end if
    if (allocated(C7_Lscale0)) then
      extents = 1_c_int64_t
      extents(1) = size(C7_Lscale0,1,kind=c_int64_t)
      n = size(C7_Lscale0)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(38)) ! Out
      if (n > 0) then
        call c_f_pointer(output(38)%data,rp,[n])
        rp = reshape(C7_Lscale0,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(38)) ! Out
    end if
    if (allocated(wpxp_L_thresh)) then
      extents = 1_c_int64_t
      extents(1) = size(wpxp_L_thresh,1,kind=c_int64_t)
      n = size(wpxp_L_thresh)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(39)) ! Out
      if (n > 0) then
        call c_f_pointer(output(39)%data,rp,[n])
        rp = reshape(wpxp_L_thresh,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(39)) ! Out
    end if
    if (allocated(c_K)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K,1,kind=c_int64_t)
      n = size(c_K)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(40)) ! Out
      if (n > 0) then
        call c_f_pointer(output(40)%data,rp,[n])
        rp = reshape(c_K,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(40)) ! Out
    end if
    if (allocated(c_K1)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K1,1,kind=c_int64_t)
      n = size(c_K1)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(41)) ! Out
      if (n > 0) then
        call c_f_pointer(output(41)%data,rp,[n])
        rp = reshape(c_K1,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(41)) ! Out
    end if
    if (allocated(nu1)) then
      extents = 1_c_int64_t
      extents(1) = size(nu1,1,kind=c_int64_t)
      n = size(nu1)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(42)) ! Out
      if (n > 0) then
        call c_f_pointer(output(42)%data,rp,[n])
        rp = reshape(nu1,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(42)) ! Out
    end if
    if (allocated(c_K2)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K2,1,kind=c_int64_t)
      n = size(c_K2)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(43)) ! Out
      if (n > 0) then
        call c_f_pointer(output(43)%data,rp,[n])
        rp = reshape(c_K2,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(43)) ! Out
    end if
    if (allocated(nu2)) then
      extents = 1_c_int64_t
      extents(1) = size(nu2,1,kind=c_int64_t)
      n = size(nu2)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(44)) ! Out
      if (n > 0) then
        call c_f_pointer(output(44)%data,rp,[n])
        rp = reshape(nu2,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(44)) ! Out
    end if
    if (allocated(c_K6)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K6,1,kind=c_int64_t)
      n = size(c_K6)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(45)) ! Out
      if (n > 0) then
        call c_f_pointer(output(45)%data,rp,[n])
        rp = reshape(c_K6,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(45)) ! Out
    end if
    if (allocated(nu6)) then
      extents = 1_c_int64_t
      extents(1) = size(nu6,1,kind=c_int64_t)
      n = size(nu6)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(46)) ! Out
      if (n > 0) then
        call c_f_pointer(output(46)%data,rp,[n])
        rp = reshape(nu6,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(46)) ! Out
    end if
    if (allocated(c_K8)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K8,1,kind=c_int64_t)
      n = size(c_K8)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(47)) ! Out
      if (n > 0) then
        call c_f_pointer(output(47)%data,rp,[n])
        rp = reshape(c_K8,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(47)) ! Out
    end if
    if (allocated(nu8)) then
      extents = 1_c_int64_t
      extents(1) = size(nu8,1,kind=c_int64_t)
      n = size(nu8)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(48)) ! Out
      if (n > 0) then
        call c_f_pointer(output(48)%data,rp,[n])
        rp = reshape(nu8,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(48)) ! Out
    end if
    if (allocated(c_K9)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K9,1,kind=c_int64_t)
      n = size(c_K9)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(49)) ! Out
      if (n > 0) then
        call c_f_pointer(output(49)%data,rp,[n])
        rp = reshape(c_K9,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(49)) ! Out
    end if
    if (allocated(nu9)) then
      extents = 1_c_int64_t
      extents(1) = size(nu9,1,kind=c_int64_t)
      n = size(nu9)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(50)) ! Out
      if (n > 0) then
        call c_f_pointer(output(50)%data,rp,[n])
        rp = reshape(nu9,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(50)) ! Out
    end if
    if (allocated(nu10)) then
      extents = 1_c_int64_t
      extents(1) = size(nu10,1,kind=c_int64_t)
      n = size(nu10)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(51)) ! Out
      if (n > 0) then
        call c_f_pointer(output(51)%data,rp,[n])
        rp = reshape(nu10,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(51)) ! Out
    end if
    if (allocated(c_K_hm)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K_hm,1,kind=c_int64_t)
      n = size(c_K_hm)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(52)) ! Out
      if (n > 0) then
        call c_f_pointer(output(52)%data,rp,[n])
        rp = reshape(c_K_hm,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(52)) ! Out
    end if
    if (allocated(c_K_hmb)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K_hmb,1,kind=c_int64_t)
      n = size(c_K_hmb)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(53)) ! Out
      if (n > 0) then
        call c_f_pointer(output(53)%data,rp,[n])
        rp = reshape(c_K_hmb,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(53)) ! Out
    end if
    if (allocated(K_hm_min_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(K_hm_min_coef,1,kind=c_int64_t)
      n = size(K_hm_min_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(54)) ! Out
      if (n > 0) then
        call c_f_pointer(output(54)%data,rp,[n])
        rp = reshape(K_hm_min_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(54)) ! Out
    end if
    if (allocated(nu_hm)) then
      extents = 1_c_int64_t
      extents(1) = size(nu_hm,1,kind=c_int64_t)
      n = size(nu_hm)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(55)) ! Out
      if (n > 0) then
        call c_f_pointer(output(55)%data,rp,[n])
        rp = reshape(nu_hm,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(55)) ! Out
    end if
    if (allocated(slope_coef_spread_DG_means_w)) then
      extents = 1_c_int64_t
      extents(1) = size(slope_coef_spread_DG_means_w,1,kind=c_int64_t)
      n = size(slope_coef_spread_DG_means_w)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(56)) ! Out
      if (n > 0) then
        call c_f_pointer(output(56)%data,rp,[n])
        rp = reshape(slope_coef_spread_DG_means_w,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(56)) ! Out
    end if
    if (allocated(pdf_component_stdev_factor_w)) then
      extents = 1_c_int64_t
      extents(1) = size(pdf_component_stdev_factor_w,1,kind=c_int64_t)
      n = size(pdf_component_stdev_factor_w)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(57)) ! Out
      if (n > 0) then
        call c_f_pointer(output(57)%data,rp,[n])
        rp = reshape(pdf_component_stdev_factor_w,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(57)) ! Out
    end if
    if (allocated(coef_spread_DG_means_rt)) then
      extents = 1_c_int64_t
      extents(1) = size(coef_spread_DG_means_rt,1,kind=c_int64_t)
      n = size(coef_spread_DG_means_rt)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(58)) ! Out
      if (n > 0) then
        call c_f_pointer(output(58)%data,rp,[n])
        rp = reshape(coef_spread_DG_means_rt,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(58)) ! Out
    end if
    if (allocated(coef_spread_DG_means_thl)) then
      extents = 1_c_int64_t
      extents(1) = size(coef_spread_DG_means_thl,1,kind=c_int64_t)
      n = size(coef_spread_DG_means_thl)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(59)) ! Out
      if (n > 0) then
        call c_f_pointer(output(59)%data,rp,[n])
        rp = reshape(coef_spread_DG_means_thl,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(59)) ! Out
    end if
    if (allocated(beta)) then
      extents = 1_c_int64_t
      extents(1) = size(beta,1,kind=c_int64_t)
      n = size(beta)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(60)) ! Out
      if (n > 0) then
        call c_f_pointer(output(60)%data,rp,[n])
        rp = reshape(beta,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(60)) ! Out
    end if
    if (allocated(gamma_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(gamma_coef,1,kind=c_int64_t)
      n = size(gamma_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(61)) ! Out
      if (n > 0) then
        call c_f_pointer(output(61)%data,rp,[n])
        rp = reshape(gamma_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(61)) ! Out
    end if
    if (allocated(gamma_coefb)) then
      extents = 1_c_int64_t
      extents(1) = size(gamma_coefb,1,kind=c_int64_t)
      n = size(gamma_coefb)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(62)) ! Out
      if (n > 0) then
        call c_f_pointer(output(62)%data,rp,[n])
        rp = reshape(gamma_coefb,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(62)) ! Out
    end if
    if (allocated(gamma_coefc)) then
      extents = 1_c_int64_t
      extents(1) = size(gamma_coefc,1,kind=c_int64_t)
      n = size(gamma_coefc)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(63)) ! Out
      if (n > 0) then
        call c_f_pointer(output(63)%data,rp,[n])
        rp = reshape(gamma_coefc,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(63)) ! Out
    end if
    if (allocated(lmin_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(lmin_coef,1,kind=c_int64_t)
      n = size(lmin_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(64)) ! Out
      if (n > 0) then
        call c_f_pointer(output(64)%data,rp,[n])
        rp = reshape(lmin_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(64)) ! Out
    end if
    if (allocated(omicron)) then
      extents = 1_c_int64_t
      extents(1) = size(omicron,1,kind=c_int64_t)
      n = size(omicron)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(65)) ! Out
      if (n > 0) then
        call c_f_pointer(output(65)%data,rp,[n])
        rp = reshape(omicron,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(65)) ! Out
    end if
    if (allocated(zeta_vrnce_rat)) then
      extents = 1_c_int64_t
      extents(1) = size(zeta_vrnce_rat,1,kind=c_int64_t)
      n = size(zeta_vrnce_rat)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(66)) ! Out
      if (n > 0) then
        call c_f_pointer(output(66)%data,rp,[n])
        rp = reshape(zeta_vrnce_rat,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(66)) ! Out
    end if
    if (allocated(upsilon_precip_frac_rat)) then
      extents = 1_c_int64_t
      extents(1) = size(upsilon_precip_frac_rat,1,kind=c_int64_t)
      n = size(upsilon_precip_frac_rat)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(67)) ! Out
      if (n > 0) then
        call c_f_pointer(output(67)%data,rp,[n])
        rp = reshape(upsilon_precip_frac_rat,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(67)) ! Out
    end if
    if (allocated(lambda0_stability_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(lambda0_stability_coef,1,kind=c_int64_t)
      n = size(lambda0_stability_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(68)) ! Out
      if (n > 0) then
        call c_f_pointer(output(68)%data,rp,[n])
        rp = reshape(lambda0_stability_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(68)) ! Out
    end if
    if (allocated(mult_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(mult_coef,1,kind=c_int64_t)
      n = size(mult_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(69)) ! Out
      if (n > 0) then
        call c_f_pointer(output(69)%data,rp,[n])
        rp = reshape(mult_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(69)) ! Out
    end if
    if (allocated(taumin)) then
      extents = 1_c_int64_t
      extents(1) = size(taumin,1,kind=c_int64_t)
      n = size(taumin)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(70)) ! Out
      if (n > 0) then
        call c_f_pointer(output(70)%data,rp,[n])
        rp = reshape(taumin,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(70)) ! Out
    end if
    if (allocated(taumax)) then
      extents = 1_c_int64_t
      extents(1) = size(taumax,1,kind=c_int64_t)
      n = size(taumax)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(71)) ! Out
      if (n > 0) then
        call c_f_pointer(output(71)%data,rp,[n])
        rp = reshape(taumax,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(71)) ! Out
    end if
    if (allocated(Lscale_mu_entrain)) then
      extents = 1_c_int64_t
      extents(1) = size(Lscale_mu_entrain,1,kind=c_int64_t)
      n = size(Lscale_mu_entrain)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(72)) ! Out
      if (n > 0) then
        call c_f_pointer(output(72)%data,rp,[n])
        rp = reshape(Lscale_mu_entrain,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(72)) ! Out
    end if
    if (allocated(Lscale_mu_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(Lscale_mu_coef,1,kind=c_int64_t)
      n = size(Lscale_mu_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(73)) ! Out
      if (n > 0) then
        call c_f_pointer(output(73)%data,rp,[n])
        rp = reshape(Lscale_mu_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(73)) ! Out
    end if
    if (allocated(Lscale_pert_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(Lscale_pert_coef,1,kind=c_int64_t)
      n = size(Lscale_pert_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(74)) ! Out
      if (n > 0) then
        call c_f_pointer(output(74)%data,rp,[n])
        rp = reshape(Lscale_pert_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(74)) ! Out
    end if
    if (allocated(alpha_corr)) then
      extents = 1_c_int64_t
      extents(1) = size(alpha_corr,1,kind=c_int64_t)
      n = size(alpha_corr)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(75)) ! Out
      if (n > 0) then
        call c_f_pointer(output(75)%data,rp,[n])
        rp = reshape(alpha_corr,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(75)) ! Out
    end if
    if (allocated(Skw_denom_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(Skw_denom_coef,1,kind=c_int64_t)
      n = size(Skw_denom_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(76)) ! Out
      if (n > 0) then
        call c_f_pointer(output(76)%data,rp,[n])
        rp = reshape(Skw_denom_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(76)) ! Out
    end if
    if (allocated(c_K10)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K10,1,kind=c_int64_t)
      n = size(c_K10)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(77)) ! Out
      if (n > 0) then
        call c_f_pointer(output(77)%data,rp,[n])
        rp = reshape(c_K10,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(77)) ! Out
    end if
    if (allocated(c_K10h)) then
      extents = 1_c_int64_t
      extents(1) = size(c_K10h,1,kind=c_int64_t)
      n = size(c_K10h)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(78)) ! Out
      if (n > 0) then
        call c_f_pointer(output(78)%data,rp,[n])
        rp = reshape(c_K10h,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(78)) ! Out
    end if
    if (allocated(thlp2_rad_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(thlp2_rad_coef,1,kind=c_int64_t)
      n = size(thlp2_rad_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(79)) ! Out
      if (n > 0) then
        call c_f_pointer(output(79)%data,rp,[n])
        rp = reshape(thlp2_rad_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(79)) ! Out
    end if
    if (allocated(thlp2_rad_cloud_frac_thresh)) then
      extents = 1_c_int64_t
      extents(1) = size(thlp2_rad_cloud_frac_thresh,1,kind=c_int64_t)
      n = size(thlp2_rad_cloud_frac_thresh)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(80)) ! Out
      if (n > 0) then
        call c_f_pointer(output(80)%data,rp,[n])
        rp = reshape(thlp2_rad_cloud_frac_thresh,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(80)) ! Out
    end if
    if (allocated(up2_sfc_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(up2_sfc_coef,1,kind=c_int64_t)
      n = size(up2_sfc_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(81)) ! Out
      if (n > 0) then
        call c_f_pointer(output(81)%data,rp,[n])
        rp = reshape(up2_sfc_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(81)) ! Out
    end if
    if (allocated(Skw_max_mag)) then
      extents = 1_c_int64_t
      extents(1) = size(Skw_max_mag,1,kind=c_int64_t)
      n = size(Skw_max_mag)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(82)) ! Out
      if (n > 0) then
        call c_f_pointer(output(82)%data,rp,[n])
        rp = reshape(Skw_max_mag,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(82)) ! Out
    end if
    if (allocated(xp3_coef_base)) then
      extents = 1_c_int64_t
      extents(1) = size(xp3_coef_base,1,kind=c_int64_t)
      n = size(xp3_coef_base)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(83)) ! Out
      if (n > 0) then
        call c_f_pointer(output(83)%data,rp,[n])
        rp = reshape(xp3_coef_base,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(83)) ! Out
    end if
    if (allocated(xp3_coef_slope)) then
      extents = 1_c_int64_t
      extents(1) = size(xp3_coef_slope,1,kind=c_int64_t)
      n = size(xp3_coef_slope)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(84)) ! Out
      if (n > 0) then
        call c_f_pointer(output(84)%data,rp,[n])
        rp = reshape(xp3_coef_slope,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(84)) ! Out
    end if
    if (allocated(altitude_threshold)) then
      extents = 1_c_int64_t
      extents(1) = size(altitude_threshold,1,kind=c_int64_t)
      n = size(altitude_threshold)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(85)) ! Out
      if (n > 0) then
        call c_f_pointer(output(85)%data,rp,[n])
        rp = reshape(altitude_threshold,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(85)) ! Out
    end if
    if (allocated(rtp2_clip_coef)) then
      extents = 1_c_int64_t
      extents(1) = size(rtp2_clip_coef,1,kind=c_int64_t)
      n = size(rtp2_clip_coef)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(86)) ! Out
      if (n > 0) then
        call c_f_pointer(output(86)%data,rp,[n])
        rp = reshape(rtp2_clip_coef,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(86)) ! Out
    end if
    if (allocated(C_invrs_tau_bkgnd)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_bkgnd,1,kind=c_int64_t)
      n = size(C_invrs_tau_bkgnd)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(87)) ! Out
      if (n > 0) then
        call c_f_pointer(output(87)%data,rp,[n])
        rp = reshape(C_invrs_tau_bkgnd,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(87)) ! Out
    end if
    if (allocated(C_invrs_tau_sfc)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_sfc,1,kind=c_int64_t)
      n = size(C_invrs_tau_sfc)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(88)) ! Out
      if (n > 0) then
        call c_f_pointer(output(88)%data,rp,[n])
        rp = reshape(C_invrs_tau_sfc,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(88)) ! Out
    end if
    if (allocated(C_invrs_tau_shear)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_shear,1,kind=c_int64_t)
      n = size(C_invrs_tau_shear)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(89)) ! Out
      if (n > 0) then
        call c_f_pointer(output(89)%data,rp,[n])
        rp = reshape(C_invrs_tau_shear,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(89)) ! Out
    end if
    if (allocated(C_invrs_tau_N2)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_N2,1,kind=c_int64_t)
      n = size(C_invrs_tau_N2)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(90)) ! Out
      if (n > 0) then
        call c_f_pointer(output(90)%data,rp,[n])
        rp = reshape(C_invrs_tau_N2,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(90)) ! Out
    end if
    if (allocated(C_invrs_tau_N2_wp2)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_N2_wp2,1,kind=c_int64_t)
      n = size(C_invrs_tau_N2_wp2)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(91)) ! Out
      if (n > 0) then
        call c_f_pointer(output(91)%data,rp,[n])
        rp = reshape(C_invrs_tau_N2_wp2,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(91)) ! Out
    end if
    if (allocated(C_invrs_tau_N2_xp2)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_N2_xp2,1,kind=c_int64_t)
      n = size(C_invrs_tau_N2_xp2)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(92)) ! Out
      if (n > 0) then
        call c_f_pointer(output(92)%data,rp,[n])
        rp = reshape(C_invrs_tau_N2_xp2,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(92)) ! Out
    end if
    if (allocated(C_invrs_tau_N2_wpxp)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_N2_wpxp,1,kind=c_int64_t)
      n = size(C_invrs_tau_N2_wpxp)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(93)) ! Out
      if (n > 0) then
        call c_f_pointer(output(93)%data,rp,[n])
        rp = reshape(C_invrs_tau_N2_wpxp,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(93)) ! Out
    end if
    if (allocated(C_invrs_tau_N2_clear_wp3)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_N2_clear_wp3,1,kind=c_int64_t)
      n = size(C_invrs_tau_N2_clear_wp3)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(94)) ! Out
      if (n > 0) then
        call c_f_pointer(output(94)%data,rp,[n])
        rp = reshape(C_invrs_tau_N2_clear_wp3,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(94)) ! Out
    end if
    if (allocated(C_invrs_tau_wpxp_Ri)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_wpxp_Ri,1,kind=c_int64_t)
      n = size(C_invrs_tau_wpxp_Ri)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(95)) ! Out
      if (n > 0) then
        call c_f_pointer(output(95)%data,rp,[n])
        rp = reshape(C_invrs_tau_wpxp_Ri,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(95)) ! Out
    end if
    if (allocated(C_invrs_tau_wpxp_N2_thresh)) then
      extents = 1_c_int64_t
      extents(1) = size(C_invrs_tau_wpxp_N2_thresh,1,kind=c_int64_t)
      n = size(C_invrs_tau_wpxp_N2_thresh)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(96)) ! Out
      if (n > 0) then
        call c_f_pointer(output(96)%data,rp,[n])
        rp = reshape(C_invrs_tau_wpxp_N2_thresh,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(96)) ! Out
    end if
    if (allocated(Cx_min)) then
      extents = 1_c_int64_t
      extents(1) = size(Cx_min,1,kind=c_int64_t)
      n = size(Cx_min)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(97)) ! Out
      if (n > 0) then
        call c_f_pointer(output(97)%data,rp,[n])
        rp = reshape(Cx_min,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(97)) ! Out
    end if
    if (allocated(Cx_max)) then
      extents = 1_c_int64_t
      extents(1) = size(Cx_max,1,kind=c_int64_t)
      n = size(Cx_max)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(98)) ! Out
      if (n > 0) then
        call c_f_pointer(output(98)%data,rp,[n])
        rp = reshape(Cx_max,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(98)) ! Out
    end if
    if (allocated(Richardson_num_min)) then
      extents = 1_c_int64_t
      extents(1) = size(Richardson_num_min,1,kind=c_int64_t)
      n = size(Richardson_num_min)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(99)) ! Out
      if (n > 0) then
        call c_f_pointer(output(99)%data,rp,[n])
        rp = reshape(Richardson_num_min,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(99)) ! Out
    end if
    if (allocated(Richardson_num_max)) then
      extents = 1_c_int64_t
      extents(1) = size(Richardson_num_max,1,kind=c_int64_t)
      n = size(Richardson_num_max)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(100)) ! Out
      if (n > 0) then
        call c_f_pointer(output(100)%data,rp,[n])
        rp = reshape(Richardson_num_max,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(100)) ! Out
    end if
    if (allocated(wpxp_Ri_exp)) then
      extents = 1_c_int64_t
      extents(1) = size(wpxp_Ri_exp,1,kind=c_int64_t)
      n = size(wpxp_Ri_exp)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(101)) ! Out
      if (n > 0) then
        call c_f_pointer(output(101)%data,rp,[n])
        rp = reshape(wpxp_Ri_exp,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(101)) ! Out
    end if
    if (allocated(a3_coef_min)) then
      extents = 1_c_int64_t
      extents(1) = size(a3_coef_min,1,kind=c_int64_t)
      n = size(a3_coef_min)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(102)) ! Out
      if (n > 0) then
        call c_f_pointer(output(102)%data,rp,[n])
        rp = reshape(a3_coef_min,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(102)) ! Out
    end if
    if (allocated(a_const)) then
      extents = 1_c_int64_t
      extents(1) = size(a_const,1,kind=c_int64_t)
      n = size(a_const)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(103)) ! Out
      if (n > 0) then
        call c_f_pointer(output(103)%data,rp,[n])
        rp = reshape(a_const,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(103)) ! Out
    end if
    if (allocated(bv_efold)) then
      extents = 1_c_int64_t
      extents(1) = size(bv_efold,1,kind=c_int64_t)
      n = size(bv_efold)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(104)) ! Out
      if (n > 0) then
        call c_f_pointer(output(104)%data,rp,[n])
        rp = reshape(bv_efold,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(104)) ! Out
    end if
    if (allocated(z_displace)) then
      extents = 1_c_int64_t
      extents(1) = size(z_displace,1,kind=c_int64_t)
      n = size(z_displace)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(105)) ! Out
      if (n > 0) then
        call c_f_pointer(output(105)%data,rp,[n])
        rp = reshape(z_displace,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(105)) ! Out
    end if
    if (allocated(Lscale_tke_multiplier)) then
      extents = 1_c_int64_t
      extents(1) = size(Lscale_tke_multiplier,1,kind=c_int64_t)
      n = size(Lscale_tke_multiplier)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(106)) ! Out
      if (n > 0) then
        call c_f_pointer(output(106)%data,rp,[n])
        rp = reshape(Lscale_tke_multiplier,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(106)) ! Out
    end if
    if (allocated(Lscale_pdf_launch_frac)) then
      extents = 1_c_int64_t
      extents(1) = size(Lscale_pdf_launch_frac,1,kind=c_int64_t)
      n = size(Lscale_pdf_launch_frac)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(107)) ! Out
      if (n > 0) then
        call c_f_pointer(output(107)%data,rp,[n])
        rp = reshape(Lscale_pdf_launch_frac,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(107)) ! Out
    end if
    if (allocated(Lscale_up_weight)) then
      extents = 1_c_int64_t
      extents(1) = size(Lscale_up_weight,1,kind=c_int64_t)
      n = size(Lscale_up_weight)
      call allocate_field(1, extents, & ! In
                          int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                          output(108)) ! Out
      if (n > 0) then
        call c_f_pointer(output(108)%data,rp,[n])
        rp = reshape(Lscale_up_weight,[n])
      end if
    else
      call empty_field(1, 0, & ! In
                       output(108)) ! Out
    end if
    extents = 1_c_int64_t
    n = 1
    call allocate_field(0, extents, int(n,c_size_t)*c_sizeof(0_c_int), 0, & ! In
                        output(109)) ! Out
    call c_f_pointer(output(109)%data,ip,[n])
    ip(1) = read_status
  end subroutine f2k_proc_read_clubb_params_nl_c
