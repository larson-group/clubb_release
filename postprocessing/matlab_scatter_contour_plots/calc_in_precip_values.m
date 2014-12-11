% $Id$
function [ mean_hm_ip, hmp2_ip, hmp3_ip, skewness_hm_ip ] ...
= calc_in_precip_values( mean_hm, hmp2, hmp3, precip_frac )

% Mean (in-precip) of the hydrometeor.
mean_hm_ip = mean_hm / precip_frac;

% Variance (in-precip) of the hydrometeor.
hmp2_ip = ( hmp2 + mean_hm^2 - precip_frac * mean_hm_ip^2 ) / precip_frac;

% Third-order central moment (in-precip) of the hydrometeor.
hmp3_ip = ( hmp3 + mean_hm^3 ...
            - 3.0 * precip_frac * ( mean_hm_ip - mean_hm ) * hmp2_ip ...
            - precip_frac * mean_hm_ip^3 ...
            + 3.0 * precip_frac * mean_hm_ip^2 * mean_hm ...
            - 3.0 * precip_frac * mean_hm_ip * mean_hm^2 ) / precip_frac;
        
% Skewness (in-precip) of the hydrometeor.
skewness_hm_ip = hmp3_ip / hmp2_ip^1.5;
