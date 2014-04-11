function[] = netcdf_rico_reader();

addpath '/home/mjfalk/mexnc2/mexnc/'
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

path = '/home/ldgrant/RICO_submission/Apr2010/RICO_grads_files/netcdf_files/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file1        = netcdf([path,'rico_file1.nc'],'nowrite');
file1_z      = file1{'zf'}(:);
file1_u      = file1{'u'}(:);
file1_v      = file1{'v'}(:);
file1_thetal= file1{'thetal'}(:);
file1_qt     = file1{'qt'}(:);
file1_rho    = file1{'rho'}(:);

keyboard

clf
plot (file1_u,file1_z)
title ('Initial u')
grid
keyboard
plot (file1_v,file1_z)
title ('Initial v')
grid
keyboard
plot (file1_thetal,file1_z)
title ('Initial theta-l')
grid
keyboard
plot (file1_qt,file1_z)
title ('Initial qt')
grid
keyboard
plot (file1_rho,file1_z)
title ('Initial rho')
grid
keyboard


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file2             = netcdf([path,'rico_file2.nc'],'nowrite');
file2_t           = file2{'time'}(:);
file2_zcb         = file2{'zcb'}(:);
file2_ztop        = file2{'ztop'}(:);
file2_zmaxcfrac   = file2{'zmaxcfrac'}(:);
file2_lwp         = file2{'LWP'}(:);
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

keyboard

plot (file2_t,file2_zcb)
title ('Ts zcb')
grid
keyboard
plot (file2_t,file2_ztop)
title ('Ts ztop')
grid
keyboard
plot (file2_t,file2_zmaxcfrac)
title ('Ts zmaxcfrac')
grid
keyboard
plot (file2_t,file2_lwp)
title ('Ts LWP')
grid
keyboard
plot (file2_t,file2_cc)
title ('Ts cc')
grid
keyboard
plot (file2_t,file2_shf)
title ('Ts shf')
grid
keyboard
plot (file2_t,file2_lhf)
title ('Ts lhf')
grid
keyboard
plot (file2_t,file2_smf)
title ('Ts smf')
grid
keyboard
plot (file2_t,file2_tke)
title ('Ts tke')
grid
keyboard
plot (file2_t,file2_v_lowlevel)
title ('Ts v_lowlevel')
grid
keyboard
plot (file2_t,file2_qv_lowlevel)
title ('Ts qv_lowlevel')
grid
keyboard
plot (file2_t,file2_th_lowlevel)
title ('Ts th_lowlevel')
grid
keyboard
plot (file2_t,file2_prec_2000)
title ('Ts prec_2000')
grid
keyboard
plot (file2_t,file2_prec_1500)
title ('Ts prec_1500')
grid
keyboard
plot (file2_t,file2_prec_1000)
title ('Ts prec_1000')
grid
keyboard
plot (file2_t,file2_prec_500)
title ('Ts prec_500')
grid
keyboard
plot (file2_t,file2_prec_srf)
title ('Ts prec_sfc')
grid
keyboard
plot (file2_t,file2_qv_2000)
title ('Ts qv_2000')
grid
keyboard
plot (file2_t,file2_qv_1500)
title ('Ts qv_1500')
grid
keyboard
plot (file2_t,file2_qv_1000)
title ('Ts qv_1000')
grid
keyboard
plot (file2_t,file2_qv_500)
title ('Ts qv_500')
grid
keyboard
plot (file2_t,file2_th_2000)
title ('Ts th_2000')
grid
keyboard
plot (file2_t,file2_th_1500)
title ('Ts th_1500')
grid
keyboard
plot (file2_t,file2_th_1000)
title ('Ts th_1000')
grid
keyboard
plot (file2_t,file2_th_500)
title ('Ts th_500')
grid
keyboard
plot (file2_t,file2_Kh_1250)
title ('Ts Kh_1250')
grid
keyboard
plot (file2_t,file2_Kh_500)
title ('Ts Kh_500')
grid
keyboard
plot (file2_t,file2_Kh_300)
title ('Ts Kh_300')
grid
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = [23 71];
t2 = [24 72];

for i=1:2
file3        = netcdf([path,'rico_file3.nc'],'nowrite');
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

keyboard

plot (file3_u,file3_z)
title ([num2str(t2(i)),' hours u'])
grid
keyboard
plot (file3_v,file3_z)
title ([num2str(t2(i)),' hours v'])
grid
keyboard
plot (file3_thetal,file3_z)
title ([num2str(t2(i)),' hours thetal'])
grid
keyboard
plot (file3_qt,file3_z)
title ([num2str(t2(i)),' hours qt'])
grid
keyboard
plot (file3_qs,file3_z)
title ([num2str(t2(i)),' hours qs'])
grid
keyboard
plot (file3_ql,file3_z)
title ([num2str(t2(i)),' hours ql'])
grid
keyboard
plot (file3_qr,file3_z)
title ([num2str(t2(i)),' hours qr'])
grid
keyboard
plot (file3_cf,file3_z)
title ([num2str(t2(i)),' hours cf'])
grid
keyboard
plot (file3_rho,file3_z)
title ([num2str(t2(i)),' hours rho'])
grid
keyboard
plot (file3_wthl,file3_z)
title ([num2str(t2(i)),' hours wthl'])
grid
keyboard
plot (file3_uw,file3_z)
title ([num2str(t2(i)),' hours uw'])
grid
keyboard
plot (file3_vw,file3_z)
title ([num2str(t2(i)),' hours vw'])
grid
keyboard
plot (file3_Mf,file3_z)
title ([num2str(t2(i)),' hours Mf'])
grid
keyboard
plot (file3_prec,file3_z)
title ([num2str(t2(i)),' hours prec'])
grid
keyboard

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:2
file4         = netcdf([path,'rico_file4.nc'],'nowrite');
file4_z       = file4{'zf'}(:);
file4_dTdt_ls = mean(file4{'dTdt_ls'}(t1(i):t2(i),:,1,1));
file4_dqdt_ls = mean(file4{'dqdt_ls'}(t1(i):t2(i),:,1,1));
file4_dudt_ls = mean(file4{'dudt_ls'}(t1(i):t2(i),:,1,1));
file4_dvdt_ls = mean(file4{'dvdt_ls'}(t1(i):t2(i),:,1,1));

plot (file4_dTdt_ls,file4_z)
title ([num2str(t2(i)),' hours dTdt_ls'])
grid
keyboard
plot (file4_dqdt_ls,file4_z)
title ([num2str(t2(i)),' hours dqdt_ls'])
grid
keyboard
plot (file4_dudt_ls,file4_z)
title ([num2str(t2(i)),' hours dudt_ls'])
grid
keyboard
plot (file4_dvdt_ls,file4_z)
title ([num2str(t2(i)),' hours dvdt_ls'])
grid
keyboard

end
