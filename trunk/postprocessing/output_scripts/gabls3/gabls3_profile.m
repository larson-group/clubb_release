function [] = gables3_profile();

% Written 11 September 2008 by Joshua Fasching
% This script converts specifications for Gables3 into
% usable values for CLUBB. Used Michael Falk's 
% mpace_a_profile as a reference.
%

%%% Initial profile %%%

% Height(m)
z = [0.0 2.0 10.0 20.0 40.0 80.0 140.0 205.0 1800.0 2200.0 5000.0 12000.0 14000.0];
% Temperature(Celcius)
T = [27.0 27.0 26.4 26.2 25.9 25.5 24.8 24.3 9.0 9.0 -6.4 -61.0 -54.0];
% Specific Humidity(kg/kg)
q = [9.3e-3 9.3e-3 8.5e-3 8.4e-3 8.3e-3 8.2e-3 8.1e-3 8.0e-3 7.5e-3 2.0e-3 0.3e-3 0.01e-3 0.003e-3];


%Height(m) for Wind
zwind = [0.0 10.0 353.0 1238.0 2000.0 5000.0];

u = [0.0 -4.0 -5.5 -5.5 -2.0 -2.0];

v = [0.0 -0.4 -0.5 -0.5 2.0 2.0];



% Constants
Rd = 287.04;
Cp = 1004.67
g0 = 9.8;
p0 = 1.000e+05;
p(1) = 1024.4e02;

% Convert Specific Humidity to Water Mixing
rt = q./(1.0 - q);

% Convert Temperature to Kelvin
T_in_K = T + 273.15

% Convert to Virtual Temperature
Tv = T_in_K .* (1 + (0.61 .* rt));

%Tvs = Ts * (1 + (0.61 * 0.001 * qs));

Tvbar(1) = (T(1) + T(2))/2;

% Determine Pressure
naturalLogPressure(2) = (-g0/(Rd*Tvbar(1)))*(z(2)-z(1)) + log(p(1));

p(2) = exp( naturalLogPressure(2) );

for i=3:max(size(Tv))
	Tvbar(i-1) = (Tv(i-1) + Tv(i))/2;
	naturalLogPressure(i) = (-g0/(Rd*Tvbar(i-1)))*(z(i)-z(i-1)) + log(p(i-1));
	p(i) = exp( naturalLogPressure(i) );
end

theta = T_in_K .* (p0./p).^(Rd/Cp);
iterpu = interp1(zwind,u,z);
iterpv = interp1(zwind,v,z);

%%% Model Forcings %%%

z_F = z(1:10);

% Prescribed Geostrophic wind
ugeo_t = [-7.8 -7.8 -6.5 -5.0 -5.0 -6.5];
vgeo_t = [0.0 0.0 4.5 4.5 4.5 2.5];
t_geo = [43200.0 64800.0 82800.0 97200.0 108000.0];


% Prescribed Vertical Dynamic tendency at specified time %
omega_t = [0.12 0.12 0.00 0.00];
% Times specified for omega
t_omega = [43200.0 61200.0 68400.0 129600.0];

% Prescribed Temperature dynamic tendency
T_adv_t = [-2.5E-5 -2.5E-5 7.5E-5 7.5E-5 0.0E+0 0.0E+0];
t_adv_T = [43200.0 90000.0 90000.0 108000.0 108000.0 129600.0];

% Prescribed Specific Humidity dynamic tendency
q_adv_t = [0.0E+0 0.0E+0 8.0E-8 8.0E-8 0.0E+0 0.0E+0 -8.0E-8 -8.0E-8 0.0E+0 0.0E+0];
t_adv_q = [43200 75600 75600 86400 86400 93600 93600 97200 97200 129600];

% Prescribed winds
u_adv_t = [0.0E+0 0.0E+0 -1.5E-4 -1.5E-4 5.0E-4 5.0E-4 0.0E+0 0.0E+0];
v_adv_t = [0.0E+0 0.0E+0 1.0E-4 1.0E-4 0.0E+0 0.0E+0 0.0E+0 0.0E+0];
t_adv_wind = [43200 64800 64800 82800 82800 97200 97200 129600];


z
theta
rt
iterpu
iterpv
keyboard
