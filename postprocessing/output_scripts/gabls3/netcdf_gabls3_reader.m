function[] = netcdf_gabls3_reader();

addpath '/home/mjfalk/mexnc2/mexnc/'
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

path = '/home/matlabuser/gabls3/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file1        = netcdf([path,'gabls3_scm_UWM_CLUBB_v3.nc'],'nowrite');
file1_time      = file1{'time'}(:)
file1_ldw       = file1{'ldw'}(:)
file1_lup       = file1{'lup'}(:)
file1_qdw       = file1{'qdw'}(:)
file1_qup       = file1{'qup'}(:) 
file1_tsk       = file1{'tsk'}(:)
file1_g         = file1{'g'}(:)
file1_shf       = file1{'shf'}(:)
file1_lhf       = file1{'lhf'}(:)
file1_ustar     = file1{'ustar'}(:)
file1_hpbl      = file1{'hpbl'}(:)
file1_t2m       = file1{'t2m'}(:)
file1_q2m       = file1{'q2m'}(:) 
file1_u10m      = file1{'u10m'}(:)
file1_v10m      = file1{'v10m'}(:)
file1_cc        = file1{'cc'}(:)

keyboard
plot (file1_time,file1_ldw)
title ('ldw')
grid

keyboard
plot (file1_time,file1_lup)
title ('lup')
grid

keyboard
plot (file1_time,file1_qdw)
title ('qdw')
grid

keyboard
plot (file1_time,file1_qup)
title ('qup')
grid

keyboard
plot (file1_time,file1_tsk)
title ('tsk')
grid

keyboard
plot (file1_time,file1_g)
title ('g')
grid

keyboard
plot (file1_time,file1_shf)
 title ('shf')
 grid
 
 keyboard
plot (file1_time,file1_lhf)
 title ('lhf')
 grid
 

 keyboard
 plot (file1_time,file1_ustar)
 title ('ustar')
 grid

keyboard
plot (file1_time,file1_hpbl)
 title ('hpbl')
 grid

 keyboard
 plot (file1_time,file1_t2m)
 title ('t2m')
 grid
 
 keyboard
plot (file1_time,file1_q2m)
 title ('q2m')
 grid

 keyboard
 plot (file1_time,file1_u10m)
 title ('u10m')
 grid

 keyboard
 plot (file1_time,file1_v10m)
 title ('v10m')
 grid

  keyboard
 plot (file1_time,file1_cc)
 title ('cc')
 grid
 
% t1= [72 144] 
% for i=1:2
file1_zf        = file1{'zf'}(:)
% file1_pf        = file1{'pf'}(t1(i),:,1,1);
file1_t         = file1{'t'}(:,:,1,1);
file1_th        = file1{'th'}(:,:,1,1);
file1_q         = file1{'q'}(:,:,1,1);
file1_u         = file1{'u'}(:,:,1,1);
file1_v         = file1{'v'}(:,:,1,1);
file1_ugeo      = file1{'ugeo'}(:,:,1,1);
file1_vgeo      = file1{'vgeo'}(:,:,1,1);
file1_dudt_ls   = file1{'dudt_ls'}(:,:,1,1);
file1_dvdt_ls   = file1{'dvdt_ls'}(:,:,1,1);
file1_dtdt_ls   = file1{'dtdt_ls'}(:,:,1,1);
file1_dqdt_ls   = file1{'dqdt_ls'}(:,:,1,1);
file1_ome       = file1{'ome'}(:,:,1,1);
% file1_zh        = file1{'zh'}(:);
% file1_ph        = file1{'ph'}(t1(i),:,1,1);
% file1_wt        = file1{'wt'}(t1(i),:,1,1);
% file1_wq        = file1{'wq'}(t1(i),:,1,1);
% file1_uw        = file1{'uw'}(t1(i),:,1,1);
% file1_vw        = file1{'vw'}(t1(i),:,1,1);
% file1_Km        = file1{'Km'}(t1(i),:,1,1);
% file1_Kh        = file1{'Kh' }(t1(i),:,1,1);
% file1_TKE       = file1{'TKE'}(t1(i),:,1,1);
% file1_shear     = file1{'shear'}(t1(i),:,1,1);
% file1_buoy      = file1{'buoy'}(t1(i),:,1,1);
% file1_trans     = file1{'trans'}(t1(i),:,1,1);
% file1_dissi     = file1{'dissi'}(t1(i),:,1,1);
% file1_zs        = file1{'zs'}(:);
% file1_ts        = file1{'ts'}(t1(i),:,1,1);
% file1_ths       = file1{'ths'}(t1(i),:,1,1);
% 
% plot (file1_pf,file1_zf)
% title ([num2str(t1(i)*10),' minutes pf'])
% grid
% 
% 
% keyboard
% 
% plot (file1_t,file1_zf)
% title ([num2str(t1(i)*10),' minutes t'])
% grid
% keyboard
% 
keyboard

plot (file1_time ./ 3600.0,file1_th(:,12))
title('th at 200m')
grid
keyboard
% 
% 
plot (file1_time ./ 3600.0,file1_q(:,12))
title('q at 200m')
grid
keyboard
 
plot (file1_time ./ 3600.0,file1_u(:,12))
title('u at 200m')
grid
keyboard

 
plot (file1_time ./ 3600.0,file1_v(:,12))
title('v at 200m')
grid
keyboard
 
 
plot (file1_time ./ 3600.0,file1_ugeo(:,12))
title('ugeo at 200m')
grid
keyboard
 
plot (file1_time ./ 3600.0,file1_vgeo(:,12))
title('vgeo at 200m')
grid
keyboard
 
plot (file1_time ./ 3600.0,file1_dudt_ls(:,12))
title('dudt_ls at 200m')
grid
keyboard
 
 
plot (file1_time ./ 3600.0,file1_dvdt_ls(:,12))
title('dvdt_ls at 200m')
grid
keyboard
 
 
plot (file1_time ./ 3600.0,file1_dtdt_ls(:,12))
title('dtdt_ls at 200m')
grid
keyboard
 
 
plot (file1_time ./ 3600.0,file1_dqdt_ls(:,12))
title('dqdt_ls at 200m')
grid
keyboard
  
plot (file1_time ./ 3600.0,file1_ome(:,12))
title ('ome at 200m')
grid
keyboard


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = [72 143]
z1 = [27 152]

for i=1:2
    for j=1:2
        plot ( file1_th(t1(i),1:z1(j)) , file1_zf(1,1:z1(j)) );
        title([num2str(t1(i)*10),' minutes th'])
        grid
        keyboard
 
        plot ( file1_q(t1(i),1:z1(j)) , file1_zf(1,1:z1(j)) )
        title([num2str(t1(i)*10),' minutes q'])
        grid
        keyboard
  
        plot ( file1_u(t1(i),1:z1(j)) , file1_zf(1,1:z1(j)) )
        title([num2str(t1(i)*10),' minutes u'])
        grid
        keyboard

        plot ( file1_v(t1(i),1:z1(j)) , file1_zf(1,1:z1(j)) )
        title([num2str(t1(i)*10),' minutes v'])
        grid
        keyboard
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% 
% plot (file1_ph,file1_zh)
% title ([num2str(t1(i)*10),' minutes ph'])
% grid
% keyboard
% 
% 
% plot (file1_wt,file1_zh)
% title ([num2str(t1(i)*10),' minutes wt'])
% grid
% keyboard
% 
% plot (file1_wq,file1_zh)
% title ([num2str(t1(i)*10),' minutes wq'])
% grid
% keyboard
% 
% 
% plot (file1_uw,file1_zh)
% title ([num2str(t1(i)*10),' minutes uw'])
% grid
% keyboard
% 
% 
% plot (file1_vw,file1_zh)
% title ([num2str(t1(i)*10),' minutes vw'])
% grid
% keyboard
% 
% plot (file1_Km,file1_zh)
% title ([num2str(t1(i)*10),' minutes Km'])
% grid
% keyboard
% 
% 
% plot (file1_Kh,file1_zh)
% title ([num2str(t1(i)*10),' minutes Kh'])
% grid
% keyboard
% 
% 
% plot (file1_TKE,file1_zh)
% title ([num2str(t1(i)*10),' minutes TKE'])
% grid
% keyboard
% 
% 
% plot (file1_shear,file1_zh)
% title ([num2str(t1(i)*10),' minutes shear'])
% grid
% keyboard
% 
% 
% plot (file1_buoy,file1_zh)
% title ([num2str(t1(i)*10),' minutes buoy'])
% grid
% keyboard
% 
% 
% plot (file1_trans,file1_zh)
% title ([num2str(t1(i)*10),' minutes trans'])
% grid
% keyboard
% 
% 
% plot (file1_dissi,file1_zh)
% title ([num2str(t1(i)*10),' minutes dissi'])
% grid
% keyboard
% keyboard 
% plot (file1_ts,file1_zs)
% title ([num2str(t1(i)*10),' minutes ts'])
% grid
% keyboard
% 
% plot (file1_ths,file1_zs)
% title ([num2str(t1(i)*10),' minutes ths'])
% grid
% keyboard
% 
% 
% end

