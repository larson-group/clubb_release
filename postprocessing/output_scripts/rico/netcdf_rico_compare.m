function[] = netcdf_rico_compare();

addpath '/home/mjfalk/mexnc2/mexnc/'
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

path = '/home/matlabuser/rico/newrico_0206/';
path_old = '/home/ldgrant/RICO_submission/Apr2010/RICO_grads_files/netcdf_files/';
path_evp = '/home/ldgrant/RICO_submission/Apr2010/RICO_grads_files/netcdf_files/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file1        = netcdf([path,'rico_file1.nc'],'nowrite');
file1_old    = netcdf([path_old,'rico_file1.nc'],'nowrite');
file1_evp    = netcdf([path_evp,'rico_file1.nc'],'nowrite');

file1_z      = file1{'zf'}(:);
file1_u      = file1{'u'}(:);
file1_v      = file1{'v'}(:);
file1_thetal = file1{'thetal'}(:);
file1_qt     = file1{'qt'}(:);
file1_rho    = file1{'rho'}(:);

file1_z_old      = file1_old{'zf'}(:);
file1_u_old      = file1_old{'u'}(:);
file1_v_old      = file1_old{'v'}(:);
file1_thetal_old = file1_old{'thetal'}(:);
file1_qt_old     = file1_old{'qt'}(:);
file1_rho_old    = file1_old{'rho'}(:);

file1_z_evp      = file1_evp{'zf'}(:);
file1_u_evp      = file1_evp{'u'}(:);
file1_v_evp      = file1_evp{'v'}(:);
file1_thetal_evp = file1_evp{'thetal'}(:);
file1_qt_evp     = file1_evp{'qt'}(:);
file1_rho_evp    = file1_evp{'rho'}(:);

clf
hold on
plot (file1_u,file1_z,'-b')
plot (file1_u_old,file1_z_old,'--g')
plot (file1_u_evp,file1_z_evp,'-r')
title ('Initial u')
grid
keyboard
hold off
clf
hold on
plot (file1_v,file1_z,'-b')
plot (file1_v_old,file1_z_old,'--g')
plot (file1_v_evp,file1_z_evp,'-r')
title ('Initial v')
grid
keyboard
hold off
clf
hold on
plot (file1_thetal,file1_z,'-b')
plot (file1_thetal_old,file1_z_old,'--g')
plot (file1_thetal_evp,file1_z_evp,'-r')
title ('Initial theta-l')
grid
keyboard
hold off
clf
hold on
plot (file1_qt,file1_z,'-b')
plot (file1_qt_old,file1_z_old,'--g')
plot (file1_qt_evp,file1_z_evp,'-r')
title ('Initial qt')
grid
keyboard
hold off
clf
hold on
plot (file1_rho,file1_z,'-b')
plot (file1_rho_old,file1_z_old,'--g')
plot (file1_rho_evp,file1_z_evp,'-r')
title ('Initial rho')
grid
keyboard
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file2             = netcdf([path,'rico_file2.nc'],'nowrite');
file2_old         = netcdf([path_old,'rico_file2.nc'],'nowrite');
file2_evp         = netcdf([path_evp,'rico_file2.nc'],'nowrite');

file2_t           = file2{'t'}(:);
file2_zcb         = file2{'zcb'}(:);
file2_ztop        = file2{'ztop'}(:);
file2_zmaxcfrac   = file2{'zmaxcfrac'}(:);
file2_lwp         = file2{'lwp'}(:);
file2_cc          = file2{'cc'}(:);
file2_shf         = file2{'shf'}(:);
file2_lhf         = file2{'lhf'}(:);
file2_smf         = file2{'smf'}(:);
file2_tke         = file2{'tke'}(:);
file2_v_lowlevel  = file2{'v_lowlevel'}(:);
file2_qv_lowlevel = file2{'qv_lowlevel'}(:);
file2_th_lowlevel = file2{'th_lowlevel'}(:);
file2_prec_2000   = file2{'prec_2000'}(:);
file2_prec_1500   = file2{'prec_1500'}(:);
file2_prec_1000   = file2{'prec_1000'}(:);
file2_prec_500    = file2{'prec_500'}(:);
file2_prec_srf    = file2{'prec_srf'}(:);
file2_qv_2000     = file2{'qv_2000'}(:);
file2_qv_1500     = file2{'qv_1500'}(:);
file2_qv_1000     = file2{'qv_1000'}(:);
file2_qv_500      = file2{'qv_500'}(:);
file2_th_2000     = file2{'th_2000'}(:);
file2_th_1500     = file2{'th_1500'}(:);
file2_th_1000     = file2{'th_1000'}(:);
file2_th_500      = file2{'th_500'}(:);
file2_Kh_1250     = file2{'Kh_1250'}(:);
file2_Kh_500      = file2{'Kh_500'}(:);
file2_Kh_300      = file2{'Kh_300'}(:);

file2_t_old           = file2_old{'time'}(:);
file2_zcb_old         = file2_old{'zcb'}(:);
file2_ztop_old        = file2_old{'ztop'}(:);
file2_zmaxcfrac_old   = file2_old{'zmaxcfrac'}(:);
file2_lwp_old         = file2_old{'LWP'}(:);
file2_cc_old          = file2_old{'cc'}(:);
file2_shf_old         = file2_old{'shf'}(:);
file2_lhf_old         = file2_old{'lhf'}(:);
file2_smf_old         = file2_old{'smf'}(:);
file2_tke_old         = file2_old{'tke'}(:);
file2_v_lowlevel_old  = file2_old{'v_lowlevel'}(:);
file2_qv_lowlevel_old = file2_old{'qv_lowlevel'}(:);
file2_th_lowlevel_old = file2_old{'th_lowlevel'}(:);
file2_prec_2000_old   = file2_old{'prec_2000'}(:);
file2_prec_1500_old   = file2_old{'prec_1500'}(:);
file2_prec_1000_old   = file2_old{'prec_1000'}(:);
file2_prec_500_old    = file2_old{'prec_500'}(:);
file2_prec_srf_old    = file2_old{'prec_srf'}(:);
file2_qv_2000_old     = file2_old{'qv_2000'}(:);
file2_qv_1500_old     = file2_old{'qv_1500'}(:);
file2_qv_1000_old     = file2_old{'qv_1000'}(:);
file2_qv_500_old      = file2_old{'qv_500'}(:);
file2_th_2000_old     = file2_old{'th_2000'}(:);
file2_th_1500_old     = file2_old{'th_1500'}(:);
file2_th_1000_old     = file2_old{'th_1000'}(:);
file2_th_500_old      = file2_old{'th_500'}(:);
file2_Kh_1250_old     = file2_old{'Kh_1250'}(:);
file2_Kh_500_old      = file2_old{'Kh_500'}(:);
file2_Kh_300_old      = file2_old{'Kh_300'}(:);

file2_t_evp           = file2_evp{'time'}(:);
file2_zcb_evp         = file2_evp{'zcb'}(:);
file2_ztop_evp        = file2_evp{'ztop'}(:);
file2_zmaxcfrac_evp   = file2_evp{'zmaxcfrac'}(:);
file2_lwp_evp         = file2_evp{'LWP'}(:);
file2_cc_evp          = file2_evp{'cc'}(:);
file2_shf_evp         = file2_evp{'shf'}(:);
file2_lhf_evp         = file2_evp{'lhf'}(:);
file2_smf_evp         = file2_evp{'smf'}(:);
file2_tke_evp         = file2_evp{'tke'}(:);
file2_v_lowlevel_evp  = file2_evp{'v_lowlevel'}(:);
file2_qv_lowlevel_evp = file2_evp{'qv_lowlevel'}(:);
file2_th_lowlevel_evp = file2_evp{'th_lowlevel'}(:);
file2_prec_2000_evp   = file2_evp{'prec_2000'}(:);
file2_prec_1500_evp   = file2_evp{'prec_1500'}(:);
file2_prec_1000_evp   = file2_evp{'prec_1000'}(:);
file2_prec_500_evp    = file2_evp{'prec_500'}(:);
file2_prec_srf_evp    = file2_evp{'prec_srf'}(:);
file2_qv_2000_evp     = file2_evp{'qv_2000'}(:);
file2_qv_1500_evp     = file2_evp{'qv_1500'}(:);
file2_qv_1000_evp     = file2_evp{'qv_1000'}(:);
file2_qv_500_evp      = file2_evp{'qv_500'}(:);
file2_th_2000_evp     = file2_evp{'th_2000'}(:);
file2_th_1500_evp     = file2_evp{'th_1500'}(:);
file2_th_1000_evp     = file2_evp{'th_1000'}(:);
file2_th_500_evp      = file2_evp{'th_500'}(:);
file2_Kh_1250_evp     = file2_evp{'Kh_1250'}(:);
file2_Kh_500_evp      = file2_evp{'Kh_500'}(:);
file2_Kh_300_evp      = file2_evp{'Kh_300'}(:);

clf
hold on
plot (file2_t,file2_zcb,'-b')
plot (file2_t_old,file2_zcb_old,'--g')
plot (file2_t_evp,file2_zcb_evp,'-r')
title ('Ts zcb')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_ztop,'-b')
plot (file2_t_old,file2_ztop_old,'--g')
plot (file2_t_evp,file2_ztop_evp,'-r')
title ('Ts ztop')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_zmaxcfrac,'-b')
plot (file2_t_old,file2_zmaxcfrac_old,'--g')
plot (file2_t_evp,file2_zmaxcfrac_evp,'-r')
title ('Ts zmaxcfrac')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_lwp,'-b')
plot (file2_t_old,file2_lwp_old,'--g')
plot (file2_t_evp,file2_lwp_evp,'-r')
title ('Ts lwp')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_cc,'-b')
plot (file2_t_old,file2_cc_old,'--g')
plot (file2_t_evp,file2_cc_evp,'-r')
title ('Ts cc')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_shf,'-b')
plot (file2_t_old,file2_shf_old,'--g')
plot (file2_t_evp,file2_shf_evp,'-r')
title ('Ts shf')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_lhf,'-b')
plot (file2_t_old,file2_lhf_old,'--g')
plot (file2_t_evp,file2_lhf_evp,'-r')
title ('Ts lhf')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_smf,'-b')
plot (file2_t_old,file2_smf_old,'--g')
plot (file2_t_evp,file2_smf_evp,'-r')
title ('Ts smf')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_tke,'-b')
plot (file2_t_old,file2_tke_old,'--g')
plot (file2_t_evp,file2_tke_evp,'-r')
title ('Ts tke')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_v_lowlevel,'-b')
plot (file2_t_old,file2_v_lowlevel_old,'--g')
plot (file2_t_evp,file2_v_lowlevel_evp,'-r')
title ('Ts v_lowlevel')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_qv_lowlevel,'-b')
plot (file2_t_old,file2_qv_lowlevel_old,'--g')
plot (file2_t_evp,file2_qv_lowlevel_evp,'-r')
title ('Ts qv_lowlevel')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_th_lowlevel,'-b')
plot (file2_t_old,file2_th_lowlevel_old,'--g')
plot (file2_t_evp,file2_th_lowlevel_evp,'-r')
title ('Ts th_lowlevel')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_prec_2000,'-b')
plot (file2_t_old,file2_prec_2000_old,'--g')
plot (file2_t_evp,file2_prec_2000_evp,'-r')
title ('Ts prec_2000')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_prec_1500,'-b')
plot (file2_t_old,file2_prec_1500_old,'--g')
plot (file2_t_evp,file2_prec_1500_evp,'-r')
title ('Ts prec_1500')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_prec_1000,'-b')
plot (file2_t_old,file2_prec_1000_old,'--g')
plot (file2_t_evp,file2_prec_1000_evp,'-r')
title ('Ts prec_1000')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_prec_500,'-b')
plot (file2_t_old,file2_prec_500_old,'--g')
plot (file2_t_evp,file2_prec_500_evp,'-r')
title ('Ts prec_500')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_prec_srf,'-b')
plot (file2_t_old,file2_prec_srf_old,'--g')
plot (file2_t_evp,file2_prec_srf_evp,'-r')
title ('Ts prec_sfc')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_qv_2000,'-b')
plot (file2_t_old,file2_qv_2000_old,'--g')
plot (file2_t_evp,file2_qv_2000_evp,'-r')
title ('Ts qv_2000')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_qv_1500,'-b')
plot (file2_t_old,file2_qv_1500_old,'--g')
plot (file2_t_evp,file2_qv_1500_evp,'-r')
title ('Ts qv_1500')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_qv_1000,'-b')
plot (file2_t_old,file2_qv_1000_old,'--g')
plot (file2_t_evp,file2_qv_1000_evp,'-r')
title ('Ts qv_1000')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_qv_500,'-b')
plot (file2_t_old,file2_qv_500_old,'--g')
plot (file2_t_evp,file2_qv_500_evp,'-r')
title ('Ts qv_500')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_th_2000,'-b')
plot (file2_t_old,file2_th_2000_old,'--g')
plot (file2_t_evp,file2_th_2000_evp,'-r')
title ('Ts th_2000')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_th_1500,'-b')
plot (file2_t_old,file2_th_1500_old,'--g')
plot (file2_t_evp,file2_th_1500_evp,'-r')
title ('Ts th_1500')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_th_1000,'-b')
plot (file2_t_old,file2_th_1000_old,'--g')
plot (file2_t_evp,file2_th_1000_evp,'-r')
title ('Ts th_1000')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_th_500,'-b')
plot (file2_t_old,file2_th_500_old,'--g')
plot (file2_t_evp,file2_th_500_evp,'-r')
title ('Ts th_500')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_Kh_1250,'-b')
plot (file2_t_old,file2_Kh_1250_old,'--g')
plot (file2_t_evp,file2_Kh_1250_evp,'-r')
title ('Ts Kh_1250')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_Kh_500,'-b')
plot (file2_t_old,file2_Kh_500_old,'--g')
plot (file2_t_evp,file2_Kh_500_evp,'-r')
title ('Ts Kh_500')
grid
keyboard
hold off
clf
hold on
plot (file2_t,file2_Kh_300,'-b')
plot (file2_t_old,file2_Kh_300_old,'--g')
plot (file2_t_evp,file2_Kh_300_evp,'-r')
title ('Ts Kh_300')
grid
keyboard
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = [23 71];
t2 = [24 72];

for i=1:2
file3        = netcdf([path,'rico_file3.nc'],'nowrite');
file3_old    = netcdf([path_old,'rico_file3.nc'],'nowrite');
file3_evp    = netcdf([path_evp,'rico_file3.nc'],'nowrite');

file3_z      = file3{'zf'}(:);
file3_u      = mean(file3{'u'}(t1(i):t2(i),:,1,1));
file3_v      = mean(file3{'v'}(t1(i):t2(i),:,1,1));
file3_thetal = mean(file3{'thetal'}(t1(i):t2(i),:,1,1));
file3_qt     = mean(file3{'qt'}(t1(i):t2(i),:,1,1));
file3_qs     = mean(file3{'qs'}(t1(i):t2(i),:,1,1));
file3_ql     = mean(file3{'ql'}(t1(i):t2(i),:,1,1));
file3_qr     = mean(file3{'qr'}(t1(i):t2(i),:,1,1));
file3_cf     = mean(file3{'cf'}(t1(i):t2(i),:,1,1));
file3_rho    = mean(file3{'rho'}(t1(i):t2(i),:,1,1));
file3_wthl   = mean(file3{'wthl'}(t1(i):t2(i),:,1,1));
file3_uw     = mean(file3{'uw'}(t1(i):t2(i),:,1,1));
file3_vw     = mean(file3{'vw'}(t1(i):t2(i),:,1,1));
file3_Mf     = mean(file3{'Mf'}(t1(i):t2(i),:,1,1));
file3_prec   = mean(file3{'prec'}(t1(i):t2(i),:,1,1));

file3_z_old      = file3_old{'zf'}(:);
file3_u_old      = mean(file3_old{'u'}(t1(i):t2(i),:,1,1));
file3_v_old      = mean(file3_old{'v'}(t1(i):t2(i),:,1,1));
file3_thetal_old = mean(file3_old{'thetal'}(t1(i):t2(i),:,1,1));
file3_qt_old     = mean(file3_old{'qt'}(t1(i):t2(i),:,1,1));
file3_qs_old     = mean(file3_old{'qs'}(t1(i):t2(i),:,1,1));
file3_ql_old     = mean(file3_old{'ql'}(t1(i):t2(i),:,1,1));
file3_qr_old     = mean(file3_old{'qr'}(t1(i):t2(i),:,1,1));
file3_cf_old     = mean(file3_old{'cf'}(t1(i):t2(i),:,1,1));
file3_rho_old    = mean(file3_old{'rho'}(t1(i):t2(i),:,1,1));
file3_wthl_old   = mean(file3_old{'wthl'}(t1(i):t2(i),:,1,1));
file3_uw_old     = mean(file3_old{'uw'}(t1(i):t2(i),:,1,1));
file3_vw_old     = mean(file3_old{'vw'}(t1(i):t2(i),:,1,1));
file3_Mf_old     = mean(file3_old{'Mf'}(t1(i):t2(i),:,1,1));
file3_prec_old   = mean(file3_old{'prec'}(t1(i):t2(i),:,1,1));

file3_z_evp      = file3_evp{'zf'}(:);
file3_u_evp      = mean(file3_evp{'u'}(t1(i):t2(i),:,1,1));
file3_v_evp      = mean(file3_evp{'v'}(t1(i):t2(i),:,1,1));
file3_thetal_evp = mean(file3_evp{'thetal'}(t1(i):t2(i),:,1,1));
file3_qt_evp     = mean(file3_evp{'qt'}(t1(i):t2(i),:,1,1));
file3_qs_evp     = mean(file3_evp{'qs'}(t1(i):t2(i),:,1,1));
file3_ql_evp     = mean(file3_evp{'ql'}(t1(i):t2(i),:,1,1));
file3_qr_evp     = mean(file3_evp{'qr'}(t1(i):t2(i),:,1,1));
file3_cf_evp     = mean(file3_evp{'cf'}(t1(i):t2(i),:,1,1));
file3_rho_evp    = mean(file3_evp{'rho'}(t1(i):t2(i),:,1,1));
file3_wthl_evp   = mean(file3_evp{'wthl'}(t1(i):t2(i),:,1,1));
file3_uw_evp     = mean(file3_evp{'uw'}(t1(i):t2(i),:,1,1));
file3_vw_evp     = mean(file3_evp{'vw'}(t1(i):t2(i),:,1,1));
file3_Mf_evp     = mean(file3_evp{'Mf'}(t1(i):t2(i),:,1,1));
file3_prec_evp   = mean(file3_evp{'prec'}(t1(i):t2(i),:,1,1));

%keyboard

clf
hold on

plot (file3_u,file3_z,'b-')
plot (file3_u_old,file3_z_old,'g--')
plot (file3_u_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours u'])
grid
keyboard
hold off
clf
hold on
plot (file3_v,file3_z,'b-')
plot (file3_v_old,file3_z_old,'g--')
plot (file3_v_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours v'])
grid
keyboard
hold off
clf
hold on
plot (file3_thetal,file3_z,'b-')
plot (file3_thetal_old,file3_z_old,'g--')
plot (file3_thetal_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours thetal'])
grid
keyboard
hold off
clf
hold on
plot (file3_qt,file3_z,'b-')
plot (file3_qt_old,file3_z_old,'g--')
plot (file3_qt_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours qt'])
grid
keyboard
hold off
clf
hold on
plot (file3_qs,file3_z,'b-')
plot (file3_qs_old,file3_z_old,'g--')
plot (file3_qs_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours qs'])
grid
keyboard
hold off
clf
hold on
plot (file3_ql,file3_z,'b-')
plot (file3_ql_old,file3_z_old,'g--')
plot (file3_ql_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours ql'])
grid
keyboard
hold off
clf
hold on
plot (file3_qr,file3_z,'b-')
plot (file3_qr_old,file3_z_old,'g--')
plot (file3_qr_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours qr'])
grid
keyboard
hold off
clf
hold on
plot (file3_cf,file3_z,'b-')
plot (file3_cf_old,file3_z_old,'g--')
plot (file3_cf_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours cf'])
grid
keyboard
hold off
clf
hold on
plot (file3_rho,file3_z,'b-')
plot (file3_rho_old,file3_z_old,'g--')
plot (file3_rho_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours rho'])
grid
keyboard
hold off
clf
hold on
plot (file3_wthl,file3_z,'b-')
plot (file3_wthl_old,file3_z_old,'g--')
plot (file3_wthl_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours wthl'])
grid
keyboard
hold off
clf
hold on
plot (file3_uw,file3_z,'b-')
plot (file3_uw_old,file3_z_old,'g--')
plot (file3_uw_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours uw'])
grid
keyboard
hold off
clf
hold on
plot (file3_vw,file3_z,'b-')
plot (file3_vw_old,file3_z_old,'g--')
plot (file3_vw_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours vw'])
grid
keyboard
hold off
clf
hold on
plot (file3_Mf,file3_z,'b-')
plot (file3_Mf_old,file3_z_old,'g--')
plot (file3_Mf_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours Mf'])
grid
keyboard
hold off
clf
hold on
plot (file3_prec,file3_z,'b-')
plot (file3_prec_old,file3_z_old,'g--')
plot (file3_prec_evp,file3_z_evp,'r-')
title ([num2str(t2(i)),' hours prec'])
grid
keyboard
hold off

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:2
file4         = netcdf([path,'rico_file4.nc'],'nowrite');
file4_old     = netcdf([path_old,'rico_file4.nc'],'nowrite');
file4_evp     = netcdf([path_evp,'rico_file4.nc'],'nowrite');

file4_z       = file4{'zf'}(:);
file4_dTdt_ls = mean(file4{'dTdt_ls'}(t1(i):t2(i),:,1,1));
file4_dqdt_ls = mean(file4{'dqdt_ls'}(t1(i):t2(i),:,1,1));
file4_dudt_ls = mean(file4{'dudt_ls'}(t1(i):t2(i),:,1,1));
file4_dvdt_ls = mean(file4{'dvdt_ls'}(t1(i):t2(i),:,1,1));

file4_z_old       = file4_old{'zf'}(:);
file4_dTdt_ls_old = mean(file4_old{'dTdt_ls'}(t1(i):t2(i),:,1,1));
file4_dqdt_ls_old = mean(file4_old{'dqdt_ls'}(t1(i):t2(i),:,1,1));
file4_dudt_ls_old = mean(file4_old{'dudt_ls'}(t1(i):t2(i),:,1,1));
file4_dvdt_ls_old = mean(file4_old{'dvdt_ls'}(t1(i):t2(i),:,1,1));

file4_z_evp       = file4_evp{'zf'}(:);
file4_dTdt_ls_evp = mean(file4_evp{'dTdt_ls'}(t1(i):t2(i),:,1,1));
file4_dqdt_ls_evp = mean(file4_evp{'dqdt_ls'}(t1(i):t2(i),:,1,1));
file4_dudt_ls_evp = mean(file4_evp{'dudt_ls'}(t1(i):t2(i),:,1,1));
file4_dvdt_ls_evp = mean(file4_evp{'dvdt_ls'}(t1(i):t2(i),:,1,1));

clf
hold on
plot (file4_dTdt_ls,file4_z,'b-')
plot (file4_dTdt_ls_old,file4_z_old,'g--')
plot (file4_dTdt_ls_evp,file4_z_evp,'r-')
title ([num2str(t2(i)),' hours dTdt_ls'])
grid
keyboard
hold off
clf
hold on
plot (file4_dqdt_ls,file4_z,'b-')
plot (file4_dqdt_ls_old,file4_z_old,'g--')
plot (file4_dqdt_ls_evp,file4_z_evp,'r-')
title ([num2str(t2(i)),' hours dqdt_ls'])
grid
keyboard
hold off
clf
hold on
plot (file4_dudt_ls,file4_z,'b-')
plot (file4_dudt_ls_old,file4_z_old,'g--')
plot (file4_dudt_ls_evp,file4_z_evp,'r-')
title ([num2str(t2(i)),' hours dudt_ls'])
grid
keyboard
hold off
clf
hold on
plot (file4_dvdt_ls,file4_z,'b-')
plot (file4_dvdt_ls_old,file4_z_old,'g--')
plot (file4_dvdt_ls_evp,file4_z_evp,'r-')
title ([num2str(t2(i)),' hours dvdt_ls'])
grid
keyboard
hold off

end
