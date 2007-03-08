function [] = make_gabls_plots();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_gabls_plots()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, January-March 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters: none
% Input files: GABLS2 GrADS files (.ctl and .dat)
% Output parameters: none
% Output files: GABLS2 images files in .eps format
%               Filenames: gabls_xxxtt.eps
%               where xxx is the field and tt is the hour.
%
% Requires: two Michael Falk utility scripts:
%           header_read.m
%           read_grads_hoc_endian.m
%
% Reference: http://people.su.se/~gsven/gabls/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program compares two GABLS2 runs (by default, the two versions
% with lKhm_aniso varying) by plotting variables of interest in GABLS
% at times when published mean fields are available for comparison.
% In addition to scrolling through using the keyboard command, you can
% view the results by using viewing the .eps files output.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scm_path = ['/home/mjfalk/hoc_results/all_lKhm_true/'];
smfile   = 'gabls2_zt.ctl';

scm_path2 = ['/home/mjfalk/hoc_results/all_lKhm_false/'];
smfile2   = 'gabls2_zt.ctl';

% set times at 23, 31, and 35 hours, so that the hour-long averages
% cover hours 24, 32, and 36.

t = [1380 1860 2100];
sizet = size(t);
sizet = max(sizet);

[filename,nz,z,ntimesteps,numvars,list_vars] = header_read([scm_path,smfile]);
kmax = max(size(z))

for i=1:numvars
    for timestep = 1:sizet
        stringtoeval = [list_vars(i,:), ' = read_grads_hoc_endian([scm_path,filename],''ieee-le'',nz,t(timestep),t(timestep)+60,i,numvars);'];
        eval(stringtoeval)
        str = list_vars(i,:);
        for j=1:nz
            arraydata(j,timestep) = eval([str,'(j)']);
        end
    end
    eval([strtrim(str),'_array_true = arraydata;']);
end

[filename,nz,z,ntimesteps,numvars,list_vars] = header_read([scm_path2,smfile2]);
kmax = max(size(z))

for i=1:numvars
    for timestep = 1:sizet
        stringtoeval = [list_vars(i,:), ' = read_grads_hoc_endian([scm_path2,filename],''ieee-le'',nz,t(timestep),t(timestep)+60,i,numvars);'];
        eval(stringtoeval)
        str = list_vars(i,:);
        for j=1:nz
            arraydata(j,timestep) = eval([str,'(j)']);
        end
    end
    eval([strtrim(str),'_array_false = arraydata;']);
end

% Set constants
g0 = 9.8;
p0 = 1e5;
R  = 287.04;
Cp = 1004.67;
Lv = 2.5e6;

% Create arrays for wind speed
U_array_true = ((um_array_true.*um_array_true) + (vm_array_true.*vm_array_true)).^0.5;
U_array_false = ((um_array_false.*um_array_false) + (vm_array_false.*vm_array_false)).^0.5;


% Make plots
hold on
plot(p_array_true(:,1),z(:),'k-');
plot(p_array_false(:,1),z(:),'b--');
title('Pressure array, hour 24');
xlabel('Pressure [Pa]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
hold off
print('-depsc','gabls_p24')
keyboard
clf
hold on
plot(wp2zt_array_true(:,1),z(:),'k-');
plot(wp2zt_array_false(:,1),z(:),'b--');
title('wp2zt array, hour 24');
xlabel('wp2zt [m^2 s^-2]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
hold off
print('-depsc','gabls_wp2zt24')
keyboard
clf
hold on
plot(um_array_true(:,1),z(:),'k-');
plot(um_array_false(:,1),z(:), 'b--');
title('um array, hour 24');
xlabel('U-velocity [m s^-1]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_um24')
hold off
keyboard
clf
hold on
plot(vm_array_true(:,1),z(:), 'k-');
plot(vm_array_false(:,1),z(:), 'b--');
title('vm array, hour 24');
xlabel('V-velocity [m s^-1]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_vm24')
hold off
keyboard
clf
hold on
plot(U_array_true(:,1),z(:), 'k-');
plot(U_array_false(:,1),z(:), 'b--');
title('windspeed array, hour 24');
xlabel('Wind speed [m s^-1]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_U24')
hold off
keyboard
clf
hold on
plot(cf_array_true(:,1),z(:), 'k-');
plot(cf_array_false(:,1),z(:), 'b--');
title('Cloud Fraction, hour 24');
xlabel('Cloud Fraction []');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_cf24')
hold off
keyboard


hold on
plot(p_array_true(:,2),z(:),'k-');
plot(p_array_false(:,2),z(:),'b--');
title('Pressure array, hour 32');
xlabel('Pressure [Pa]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_p32')
hold off
keyboard
clf
hold on
plot(wp2zt_array_true(:,2),z(:),'k-');
plot(wp2zt_array_false(:,2),z(:),'b--');
title('wp2zt array, hour 32');
xlabel('wp2zt [m^2 s^-2]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_wp2zt32')
hold off
keyboard
clf
hold on
plot(um_array_true(:,2),z(:),'k-');
plot(um_array_false(:,2),z(:), 'b--');
title('um array, hour 32');
xlabel('U-velocity [m s^-1]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_um32')
hold off
keyboard
clf
hold on
plot(vm_array_true(:,2),z(:), 'k-');
plot(vm_array_false(:,2),z(:), 'b--');
title('vm array, hour 32');
xlabel('V-velocity [m s^-1]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_vm32')
hold off
keyboard
clf
hold on
plot(U_array_true(:,2),z(:), 'k-');
plot(U_array_false(:,2),z(:), 'b--');
title('windspeed array, hour 32');
xlabel('Wind speed [m s^-1]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_U32')
hold off
keyboard
clf
hold on
plot(cf_array_true(:,2),z(:), 'k-');
plot(cf_array_false(:,2),z(:), 'b--');
title('Cloud Fraction, hour 32');
xlabel('Cloud Fraction []');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_cf32')
hold off
keyboard

hold on
plot(p_array_true(:,3),z(:),'k-');
plot(p_array_false(:,3),z(:),'b--');
title('Pressure array, hour 36');
xlabel('Pressure [Pa]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_p36')
hold off
keyboard
clf
hold on
plot(wp2zt_array_true(:,3),z(:),'k-');
plot(wp2zt_array_false(:,3),z(:),'b--');
title('wp2zt array, hour 36');
xlabel('wp2zt [m^2 s^-2]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_wp2zt36')
hold off
keyboard
clf
hold on
plot(um_array_true(:,3),z(:),'k-');
plot(um_array_false(:,3),z(:), 'b--');
title('um array, hour 36');
xlabel('U-velocity [m s^-1]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_um36')
hold off
keyboard
clf
hold on
plot(vm_array_true(:,3),z(:), 'k-');
plot(vm_array_false(:,3),z(:), 'b--');
title('vm array, hour 36');
xlabel('V-velocity [m s^-1]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_vm36')
hold off
keyboard
clf
hold on
plot(U_array_true(:,3),z(:), 'k-');
plot(U_array_false(:,3),z(:), 'b--');
title('windspeed array, hour 36');
xlabel('Wind speed [m s^-1]');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_U36')
hold off
keyboard
clf
hold on
plot(cf_array_true(:,3),z(:), 'k-');
plot(cf_array_false(:,3),z(:), 'b--');
title('Cloud Fraction, hour 36');
xlabel('Cloud Fraction []');
ylabel('Height [m]');
legend('lKhm true','lKhm false');
print('-depsc','gabls_cf36')
hold off
keyboard