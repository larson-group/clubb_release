function [] = rico_profile();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rico_profile()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, 13 December 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters: none
% Input files: none
% Output parameters: none
% Output files: none
%
% Requires: two Michael Falk utility scripts:
%           header_read.m
%           read_grads_hoc_endian.m
%
% References: http://www.knmi.nl/samenw/rico/setup1d_composite.html
%             -A short course in cloud physics.  R. R. Rogers and M. K. Yau,
%               Butterworth/Heinemann 1989.
%             -Atmospheric science: an introductory survey.  J. M. Wallace
%               and P. V. Hobbs, Academic Press 1977.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hoc takes inputs of Theta and mixing ratio (w); RICO specifies T and
% specific humidity (qv).  This program does the necessary conversion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up constants
p_sfc = 101540;  % From RICO spec
p0 = 100000;    % From standard meteorological definitions
g0 = 9.81;      % From RICO spec
Rd = 287;       % From RICO spec
Cp = 1005;      % From RICO spec
eps = 0.622;    % From standard meteorological definitions
zeroc = 273.15; % From standard meteorological definitions

% original spec
%z = [0,10,540,700,1540, ...                              % This whole section from RICO specification
%    2100,3300,4000,10000,12000, ...                      %
%    13000,17000,20000,50000];                            %
%T = [299.8,298.907491,294,293.024,287.9, ...             %
%    285.9,279.521053,275.8,238.046154,225.461538, ...    %
%    219.169231,194,206,270];                             %
%Tc = T - zeroc;                                          %
%qv = .001 .* [0,15.974074,14.6,13.892754,10.179710, ...  %
%     7.704348,2.4,1.6,0,0, ...                           %
%     0,0,0,0];                                           %
%u = [-8.5,-8.5,-8.5,-8.5,-6.845455, ...                  %
%    -5.742424,-3.378788,-2,22,30, ...                    %
%    30,12.857143,0,0];                                   %
%v = [-3.8,-3.8,-3.8,-3.8,-3.8, ...                       %
%    -3.8,-3.8,-3.8,-3.8,-3.8, ...                        %
%    -3.8,-3.8,-3.8,-3.8];                                % End RICO specification

% revised spec, released 15 December 2006
z = [0,740,3260,4000,9000,12000, ...                     % This whole section from RICO specification
     13000,15000,17500,20000,60000];                     %
T = [299.2,292,281.177914,278,243.909091,223.454545, ... %
     216.636364,203,194,206,270];                        %
Tc = T - zeroc;                                          %
qv = .001 .* [16,13.8,2.4,1.8,0,0, ...                   %
     0,0,0,0,0];                                         %
u = [-9.9,-8.42,-3.38,-1.9,18.0375,30, ...               %
     30,21.428571,10.714286,0,0];                        %
v = [-3.8,-3.8,-3.8,-3.8,-3.8, -3.8, ...                 %
    -3.8,-3.8,-3.8,-3.8,-3.8];                           % End RICO specification

% Compute surface variables
p(1) = p_sfc;
e(1) = 611.2 * exp ((17.67 * Tc(1))/(Tc(1) + 243.5)); % Bolton Formula (Rogers & Yau 2.17) saturated at lowest level at initial time
% check that sfc uses this formulation, not 0m qv
%qv(1) = eps * ((e(1))/(p(1) - ((1-eps)*e(1)))); % Rogers & Yau 2.19
%w(1) = eps * (e(1) / (p(1)-e(1))); % Rogers & Yau 2.18
w(1) = (eps*qv(1))/(eps-qv(1)); % Rogers & Yau 2.18/2.19 -> modified by Michael Falk derivation
%Tv(1) = T(1) / (1-((e(1)/p(1))*(1-eps))); % Wallace & Hobbs 2.17
Tv(1) = T(1)*(1 + 0.61*w(1)); % Wallace & Hobbs 2.62
Theta(1) = T(1) * (p0/p_sfc)^(Rd/Cp)

% Compute upward profile
for i=2:max(size(z))
  w(i) = (eps*qv(i))/(eps-qv(i)); % Rogers & Yau 2.18/2.19 -> modified by Michael Falk derivation
%  e(i) = qv(i) / (p(i) * eps); p(i) is undefined here
%  Tv(i) = T(i) / (1-((e(i)/p(i))*(1-eps))); % Wallace & Hobbs 2.17 ; unusable because e cannot be computed before p
  Tv(i) = T(i)*(1 + 0.61*w(i)); % Wallace & Hobbs 2.62
  p(i) = p(i-1) * exp ( (g0/Rd)*(2/(Tv(i)+Tv(i-1)))*(z(i-1)-z(i)) ); % Wallace and Hobbs 2.29
  Theta(i) = T(i) * (p0/p(i))^(Rd/Cp); % Wallace & Hobbs 2.57
end

% Pause to interact and print necessary variables
keyboard
