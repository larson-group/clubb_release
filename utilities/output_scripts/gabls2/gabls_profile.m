function [] = gabls_profile();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_gabls_plots()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, 29 December 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merely converts q (specific humidity) to r (mixing ratio) from
% GABLS specification by Gunilla Svensson.
%
% Reference: http://people.su.se/~gsven/gabls/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = [.0025, .0025, .0025, .0025, .0005, .0030, .0020, .0015]; % [q] = kg/kg per GABLS2 specification

for i=1:max(size(q))
    r(i) = q(i) / (1-q(i));
end

r_g_kg = r .* 1000