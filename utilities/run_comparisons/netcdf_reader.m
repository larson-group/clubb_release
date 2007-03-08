function[] = netcdf_reader();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% netcdf_reader()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, February 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters: none
% Input files: LES and SCM output files as distributed with HOC
% Output parameters: none
% Output files: none
%
% Requires: MexNC, NetCDF_toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program makes plots of the GCSS and Larsongroup cases, deciphering
% the NetCDF files and plotting LES against SCM results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath '/home/mjfalk/mexnc2/mexnc/'
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

hocpath = '/home/mjfalk/bugsrad/profiles-cmp/data/hoc/config2/';
lespath = '/home/mjfalk/bugsrad/profiles-cmp/data/les_data/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

armt1 = 480; % = 8 hours
armt2 = 540; % = 9 hours

arm_scm = netcdf([hocpath,'arm_zt.nc'],'nowrite');
arm_les = netcdf([lespath,'arm_coamps_sm.nc'],'nowrite');
arm_scm_z = arm_scm{'level110'}(:);
arm_les_z = arm_les{'level110'}(:);

arm_scm_cf = mean(arm_scm{'cf'}(armt1:armt2,:,1,1));
arm_les_cf = mean(arm_les{'cf'}(armt1:armt2,:,1,1));
arm_scm_rcm = mean(arm_scm{'rcm'}(armt1:armt2,:,1,1));
arm_les_qcm = mean(arm_les{'qcm'}(armt1:armt2,:,1,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bomext1 = 240; % = 4 hours
bomext2 = 360; % = 6 hours

bomex_scm = netcdf([hocpath,'bomex_zt.nc'],'nowrite');
bomex_les = netcdf([lespath,'bomex_coamps_sm.nc'],'nowrite');
bomex_scm_z = bomex_scm{'level075'}(:);
bomex_les_z = bomex_les{'level075'}(:);

bomex_scm_cf = mean(bomex_scm{'cf'}(bomext1:bomext2,:,1,1));
bomex_les_cf = mean(bomex_les{'cf'}(bomext1:bomext2,:,1,1));
bomex_scm_rcm = mean(bomex_scm{'rcm'}(bomext1:bomext2,:,1,1));
bomex_les_qcm = mean(bomex_les{'qcm'}(bomext1:bomext2,:,1,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atext1 = 420; % = 7 hours
atext2 = 480; % = 8 hours

atex_scm = netcdf([hocpath,'atex_zt.nc'],'nowrite');
atex_les = netcdf([lespath,'atex_coamps_sm.nc'],'nowrite');
atex_scm_z = atex_scm{'level150'}(:);
atex_les_z = atex_les{'level150'}(:);

atex_scm_cf = mean(atex_scm{'cf'}(atext1:atext2,:,1,1));
atex_les_cf = mean(atex_les{'cf'}(atext1:atext2,:,1,1));
atex_scm_rcm = mean(atex_scm{'rcm'}(atext1:atext2,:,1,1));
atex_les_qcm = mean(atex_les{'qcm'}(atext1:atext2,:,1,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rf01t1 = 180; % = 3 hours
rf01t2 = 240; % = 4 hours

rf01_scm = netcdf([hocpath,'dycoms2_rf01_zt.nc'],'nowrite');
rf01_les = netcdf([lespath,'dycoms2_rf01_coamps_sm.nc'],'nowrite');
rf01_scm_z = rf01_scm{'level132'}(:);
rf01_les_z = rf01_les{'level130'}(:);

rf01_scm_cf = mean(rf01_scm{'cf'}(rf01t1:rf01t2,:,1,1));
rf01_les_cf = mean(rf01_les{'cf'}(rf01t1:rf01t2,:,1,1));
rf01_scm_rcm = mean(rf01_scm{'rcm'}(rf01t1:rf01t2,:,1,1));
rf01_les_qcm = mean(rf01_les{'qcm'}(rf01t1:rf01t2,:,1,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rf02t1 = 300; % = 5 hours
rf02t2 = 360; % = 6 hours

rf02_scm = netcdf([hocpath,'dycoms2_rf02_nd_zt.nc'],'nowrite');
rf02_les = netcdf([lespath,'dycoms2_rf02_nd_coamps_sm.nc'],'nowrite');
rf02_scm_z = rf02_scm{'level147'}(:);
rf02_les_z = rf02_les{'level147'}(:);

rf02_scm_cf = mean(rf02_scm{'cf'}(rf02t1:rf02t2,:,1,1));
rf02_les_cf = mean(rf02_les{'cf'}(rf02t1:rf02t2,:,1,1));
rf02_scm_rcm = mean(rf02_scm{'rcm'}(rf02t1:rf02t2,:,1,1));
rf02_les_qcm = mean(rf02_les{'qcm'}(rf02t1:rf02t2,:,1,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firet1 = 60; % = 1 hour
firet2 = 120; % = 2 hours

fire_scm = netcdf([hocpath,'fire_zt.nc'],'nowrite');
fire_les = netcdf([lespath,'fire_coamps_sm.nc'],'nowrite');
fire_scm_z = fire_scm{'level048'}(:);
fire_les_z = fire_les{'level048'}(:);

fire_scm_cf = mean(fire_scm{'cf'}(firet1:firet2,:,1,1));
fire_les_cf = mean(fire_les{'cf'}(firet1:firet2,:,1,1));
fire_scm_rcm = mean(fire_scm{'rcm'}(firet1:firet2,:,1,1));
fire_les_qcm = mean(fire_les{'qcm'}(firet1:firet2,:,1,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up plots, then plot results.

width = 0.38;
wborder = (1 - (2*width)) * (1./4.);
height = 0.20;
hborder = (1-(3*height)) * (1./6.);

vector3km = [0;1;2;3]
vector2km = [0;1;2]
vector1km = [0;1]

clf
p1=subplot('Position',[wborder ((2*height) + (5*hborder)) width height]) %upper left
p2=subplot('Position',[(width+(3*wborder)) ((2*height) + (5*hborder)) width height]) %upper right
p3=subplot('Position',[wborder (height + (3*hborder)) width height]) %middle left
p4=subplot('Position',[(width+(3*wborder)) (height + (3*hborder)) width height]) %middle right
p5=subplot('Position', [wborder hborder width height]) %lower left
p6=subplot('Position',[(width+(3*wborder)) hborder width height]) %lower right
correctp1=get(p1,'Position');
correctp2=get(p2,'Position');
correctp3=get(p3,'Position');
correctp4=get(p4,'Position');
correctp5=get(p5,'Position');
correctp6=get(p6,'Position');

orient tall;

p1=subplot('Position',[wborder ((2*height) + (5*hborder)) width height]) %upper left
plot(arm_scm_cf*100,arm_scm_z./1000,'b-o',arm_les_cf*100,arm_les_z./1000,'r-x')
xlabel('%')
ylabel('Height (km)')
xlim ([0 10])
ylim ([0 3.500])
set(gca,'YTick',vector3km,'YTickLabel',vector3km)
title('Cloud fraction')
p2=subplot('Position',[(width+(3*wborder)) ((2*height) + (5*hborder)) width height]) %upper right
plot(arm_scm_rcm*1000,arm_scm_z./1000,'b-o',arm_les_qcm*1000,arm_les_z./1000,'r-x')
xlabel('g kg^{-1}')
ylabel('Height (km)')
xlim ([0 0.03])
ylim ([0 3.500])
set(gca,'YTick',vector3km,'YTickLabel',vector3km)
legend('PDF Parameterization','LES','Location','SouthEast')
legend('boxoff')
title('Cloud water mixing ratio')

p3=subplot('Position',[wborder (height + (3*hborder)) width height]) %middle left
plot(bomex_scm_cf*100,bomex_scm_z./1000,'b-o',bomex_les_cf*100,bomex_les_z./1000,'r-x')
xlabel('%')
ylabel('Height (km)')
xlim ([0 10])
ylim ([0 2.200])
set(gca,'YTick',vector2km,'YTickLabel',vector2km)
title('Cloud fraction')
p4=subplot('Position',[(width+(3*wborder)) (height + (3*hborder)) width height]) %middle right
plot(bomex_scm_rcm*1000,bomex_scm_z./1000,'b-o',bomex_les_qcm*1000,bomex_les_z./1000,'r-x')
xlabel('g kg^{-1}')
ylabel('Height (km)')
xlim ([0 0.01])
ylim ([0 2.200])
set(gca,'YTick',vector2km,'YTickLabel',vector2km)
legend('PDF Parameterization','LES','Location','SouthEast')
legend('boxoff')
title('Cloud water mixing ratio')

p5=subplot('Position', [wborder hborder width height]) %lower left
plot(atex_scm_cf*100,atex_scm_z./1000,'b-o',atex_les_cf*100,atex_les_z./1000,'r-x')
xlabel('%')
ylabel('Height (km)')
xlim ([0 100])
ylim ([0 2.000])
set(gca,'YTick',vector2km,'YTickLabel',vector2km)
title('Cloud fraction')
p6=subplot('Position',[(width+(3*wborder)) hborder width height]) %lower right
plot(atex_scm_rcm*1000,atex_scm_z./1000,'b-o',atex_les_qcm*1000,atex_les_z./1000,'r-x')
xlabel('g kg^{-1}')
ylabel('Height (km)')
xlim ([0 0.2])
ylim ([0 2.000])
set(gca,'YTick',vector2km,'YTickLabel',vector2km)
legend('PDF Parameterization','LES','Location','SouthEast')
legend('boxoff')
title('Cloud water mixing ratio')

axes('Position',[0 0 1 1],'Visible','off')
toplabel = 'ARM (continental cumulus)';
h = text(.37,.965,toplabel,'FontSize',12)
midlabel = 'BOMEX (trade-wind cumulus)';
h = text(.36,.63,midlabel,'FontSize',12)
botlabel = 'ATEX (cumulus rising into stratocumulus)';
h = text(.32,.30,botlabel,'FontSize',12)

keyboard
set(p2,'Position',correctp2);
set(p4,'Position',correctp4);
set(p6,'Position',correctp6);
print('-depsc','armbomexatex')

keyboard

clf
orient tall;

p1=subplot('Position',[wborder ((2*height) + (5*hborder)) width height]) %upper left
plot(rf01_scm_cf*100,rf01_scm_z./1000,'b-o',rf01_les_cf*100,rf01_les_z./1000,'r-x')
xlabel('%')
ylabel('Height (km)')
xlim ([0 100])
ylim ([0 1.200])
set(gca,'YTick',vector1km,'YTickLabel',vector1km)
title('Cloud fraction')
p2=subplot('Position',[(width+(3*wborder)) ((2*height) + (5*hborder)) width height]) %upper right
plot(rf01_scm_rcm*1000,rf01_scm_z./1000,'b-o',rf01_les_qcm*1000,rf01_les_z./1000,'r-x')
xlabel('g kg^{-1}')
ylabel('Height (km)')
xlim ([0 0.4])
ylim ([0 1.200])
set(gca,'YTick',vector1km,'YTickLabel',vector1km)
legend('PDF Parameterization','LES','Location','SouthEast')
legend('boxoff')
title('Cloud water mixing ratio')

p3=subplot('Position',[wborder (height + (3*hborder)) width height]) %middle left
plot(rf02_scm_cf*100,rf02_scm_z./1000,'b-o',rf02_les_cf*100,rf02_les_z./1000,'r-x')
xlabel('%')
ylabel('Height (km)')
xlim ([0 100])
ylim ([0 1.200])
set(gca,'YTick',vector1km,'YTickLabel',vector1km)
title('Cloud fraction')
p4=subplot('Position',[(width+(3*wborder)) (height + (3*hborder)) width height]) %middle right
plot(rf02_scm_rcm*1000,rf02_scm_z./1000,'b-o',rf02_les_qcm*1000,rf02_les_z./1000,'r-x')
xlabel('g kg^{-1}')
ylabel('Height (km)')
xlim ([0 0.6])
ylim ([0 1.200])
set(gca,'YTick',vector1km,'YTickLabel',vector1km)
legend('PDF Parameterization','LES','Location','SouthEast')
legend('boxoff')
title('Cloud water mixing ratio')

p5=subplot('Position', [wborder hborder width height]) %lower left
plot(fire_scm_cf*100,fire_scm_z./1000,'b-o',fire_les_cf*100,fire_les_z./1000,'r-x')
xlabel('%')
ylabel('Height (km)')
xlim ([0 100])
ylim ([0 1.200])
set(gca,'YTick',vector1km,'YTickLabel',vector1km)
title('Cloud fraction')
p6=subplot('Position',[(width+(3*wborder)) hborder width height]) %lower right
plot(fire_scm_rcm*1000,fire_scm_z./1000,'b-o',fire_les_qcm*1000,fire_les_z./1000,'r-x')
xlabel('g kg^{-1}')
ylabel('Height (km)')
xlim ([0 0.4])
ylim ([0 1.200])
set(gca,'YTick',vector1km,'YTickLabel',vector1km)
legend('PDF Parameterization','LES','Location','SouthEast')
legend('boxoff')
title('Cloud water mixing ratio')

axes('Position',[0 0 1 1],'Visible','off')
toplabel = 'DYCOMS-II RF-01 (nocturnal stratocumulus)';
h = text(.31,.97,toplabel,'FontSize',12)
midlabel = 'DYCOMS-II RF-02 (nocturnal stratocumulus)';
h = text(.31,.635,midlabel,'FontSize',12)
botlabel = 'FIRE (nocturnal stratocumulus)';
h = text(.355,.295,botlabel,'FontSize',12)

keyboard
set(p2,'Position',correctp2);
set(p4,'Position',correctp4);
set(p6,'Position',correctp6);
print('-depsc','dycomsfire')

keyboard