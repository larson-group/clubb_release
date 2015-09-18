% $Id$
function [ wp2_sam, wp3_sam, rtp2_sam, rtp3_sam, thlp2_sam, thlp3_sam, ...
           wprtp_sam, wpthlp_sam, rtpthlp_sam, ...
           Skw_sam, Skrt_sam, Skthl_sam, ...
           corr_w_rt_ov_sam, corr_w_thl_ov_sam, corr_rt_thl_ov_sam, ...
           precip_frac_sam, mean_rr_sam, rrp2_sam, rrp3_sam, Skrr_sam, ...
           mean_rr_ip_sam, rrp2_ip_sam, rrp3_ip_sam, ...
           Skrr_ip_sam, voms_ip_rr_sam, mean_lnrr_sam, ...
           lnrrp2_sam, lnrrp3_sam, Sk_lnrr_sam, ...
           mean_Nr_sam, Nrp2_sam, Nrp3_sam, SkNr_sam, ...
           mean_Nr_ip_sam, Nrp2_ip_sam, Nrp3_ip_sam, ...
           SkNr_ip_sam, voms_ip_Nr_sam, mean_lnNr_sam, ...
           lnNrp2_sam, lnNrp3_sam, Sk_lnNr_sam ] ...
= calc_stats_SAM( sam_var_lev, nx_sam, ny_sam )

% SAM LES 3D file variable indices
global idx_3D_w
global idx_3D_rt
global idx_3D_thl
global idx_3D_rr
global idx_3D_Nr

% Calculations based on SAM LES 3D data.
mean_w_sam = mean( sam_var_lev(idx_3D_w,:) );
mean_rt_sam = mean( sam_var_lev(idx_3D_rt,:) );
mean_thetal_sam = mean( sam_var_lev(idx_3D_thl,:) );
sumw2 = 0.0;
sumw3 = 0.0;
sumrt2 = 0.0;
sumrt3 = 0.0;
sumthl2 = 0.0;
sumthl3 = 0.0;
sumwrt = 0.0;
sumwthl = 0.0;
sumrtthl = 0.0;
for idx = 1:1:nx_sam*ny_sam
   sumw2   = sumw2   + ( sam_var_lev(idx_3D_w,idx)   - mean_w_sam )^2;
   sumw3   = sumw3   + ( sam_var_lev(idx_3D_w,idx)   - mean_w_sam )^3;
   sumrt2  = sumrt2  + ( sam_var_lev(idx_3D_rt,idx)  - mean_rt_sam )^2;
   sumrt3  = sumrt3  + ( sam_var_lev(idx_3D_rt,idx)  - mean_rt_sam )^3;
   sumthl2 = sumthl2 + ( sam_var_lev(idx_3D_thl,idx) - mean_thetal_sam )^2;
   sumthl3 = sumthl3 + ( sam_var_lev(idx_3D_thl,idx) - mean_thetal_sam )^3;
   sumwrt   = sumwrt ...
              + ( sam_var_lev(idx_3D_w,idx) - mean_w_sam ) ...
                * ( sam_var_lev(idx_3D_rt,idx) - mean_rt_sam );
   sumwthl  = sumwthl ...
              + ( sam_var_lev(idx_3D_w,idx) - mean_w_sam ) ...
                * ( sam_var_lev(idx_3D_thl,idx) - mean_thetal_sam );
   sumrtthl = sumrtthl ...
              + ( sam_var_lev(idx_3D_rt,idx) - mean_rt_sam ) ...
                * ( sam_var_lev(idx_3D_thl,idx) - mean_thetal_sam );
end
wp2_sam   = sumw2   / ( nx_sam * ny_sam );
wp3_sam   = sumw3   / ( nx_sam * ny_sam );
rtp2_sam  = sumrt2  / ( nx_sam * ny_sam );
rtp3_sam  = sumrt3  / ( nx_sam * ny_sam );
thlp2_sam = sumthl2 / ( nx_sam * ny_sam );
thlp3_sam = sumthl3 / ( nx_sam * ny_sam );
wprtp_sam   = sumwrt   / ( nx_sam * ny_sam );
wpthlp_sam  = sumwthl  / ( nx_sam * ny_sam );
rtpthlp_sam = sumrtthl / ( nx_sam * ny_sam );

Skw_sam   = wp3_sam   / wp2_sam^1.5;
Skrt_sam  = rtp3_sam  / rtp2_sam^1.5;
Skthl_sam = thlp3_sam / thlp2_sam^1.5;
corr_w_rt_ov_sam ...
   = wprtp_sam   / ( sqrt( wp2_sam )  * sqrt( rtp2_sam ) );
corr_w_thl_ov_sam ...
   = wpthlp_sam  / ( sqrt( wp2_sam )  * sqrt( thlp2_sam ) );
corr_rt_thl_ov_sam ...
   = rtpthlp_sam / ( sqrt( rtp2_sam ) * sqrt( thlp2_sam ) );

% Calculations for hydrometeor fields based on SAM LES 3D data.
mean_rr_sam = mean( sam_var_lev(idx_3D_rr,:) );
mean_Nr_sam = mean( sam_var_lev(idx_3D_Nr,:) );
sumrr2 = 0.0;
sumrr3 = 0.0;
sumNr2 = 0.0;
sumNr3 = 0.0;
for idx = 1:1:nx_sam*ny_sam
   sumrr2 = sumrr2 + ( sam_var_lev(idx_3D_rr,idx) - mean_rr_sam )^2;
   sumrr3 = sumrr3 + ( sam_var_lev(idx_3D_rr,idx) - mean_rr_sam )^3;
   sumNr2 = sumNr2 + ( sam_var_lev(idx_3D_Nr,idx) - mean_Nr_sam )^2;
   sumNr3 = sumNr3 + ( sam_var_lev(idx_3D_Nr,idx) - mean_Nr_sam )^3;
end
rrp2_sam = sumrr2 / ( nx_sam * ny_sam );
rrp3_sam = sumrr3 / ( nx_sam * ny_sam );
Nrp2_sam = sumNr2 / ( nx_sam * ny_sam );
Nrp3_sam = sumNr3 / ( nx_sam * ny_sam );

Skrr_sam = rrp3_sam / rrp2_sam^1.5;
SkNr_sam = Nrp3_sam / Nrp2_sam^1.5;

precipitating = zeros(nx_sam*ny_sam,1);
for idx = 1:1:nx_sam*ny_sam
   if ( sam_var_lev(idx_3D_rr,idx) > 0.0 ...
        || sam_var_lev(idx_3D_Nr,idx) > 0.0 )
      precipitating(idx) = 1.0;
   else
      precipitating(idx) = 0.0;
   end
end
precip_frac_sam = sum( precipitating ) / ( nx_sam * ny_sam );

if ( mean_rr_sam > 0.0 )
   % Rain water mixing ratio is found at this level in SAM.
   [ mean_rr_ip_sam, rrp2_ip_sam, rrp3_ip_sam, Skrr_ip_sam ] ...
   = calc_in_precip_values( mean_rr_sam, rrp2_sam, ...
                            rrp3_sam, precip_frac_sam );
   % Calculate in-precip variance-over-mean-squared for rr.
   voms_ip_rr_sam = rrp2_ip_sam / mean_rr_ip_sam^2;
else
   % Set values to 0.
   mean_rr_ip_sam = 0.0;
   rrp2_ip_sam = 0.0;
   rrp3_ip_sam = 0.0;
   Skrr_ip_sam = 0.0;
   voms_ip_rr_sam = 0.0;
end

if ( mean_Nr_sam > 0.0 )
   % Rain drop concentration is found at this level in SAM.
   [ mean_Nr_ip_sam, Nrp2_ip_sam, Nrp3_ip_sam, SkNr_ip_sam ] ...
   = calc_in_precip_values( mean_Nr_sam, Nrp2_sam, ...
                            Nrp3_sam, precip_frac_sam );
   % Calculate in-precip variance-over-mean-squared for Nr.
   voms_ip_Nr_sam = Nrp2_ip_sam / mean_Nr_ip_sam^2;
else
   % Set values to 0.
   mean_Nr_ip_sam = 0.0;
   Nrp2_ip_sam = 0.0;
   Nrp3_ip_sam = 0.0;
   SkNr_ip_sam = 0.0;
   voms_ip_Nr_sam = 0.0;
end

% Calculate the mean, variance, 3rd-order central moment, and skewness of
% ln rr for all rr in precip. (where rr > 0) for SAM LES.
count_lnrr = 0;
for idx = 1:1:nx_sam*ny_sam
   if ( sam_var_lev(idx_3D_rr,idx) > 0.0 )
      count_lnrr = count_lnrr + 1;
      sam_lnrr(count_lnrr) = log( sam_var_lev(idx_3D_rr,idx) );
   end
end
if ( count_lnrr > 0 )
   mean_lnrr_sam = mean( sam_lnrr );
   sum_lnrr2 = 0.0;
   sum_lnrr3 = 0.0;
   for idx = 1:1:count_lnrr
      sum_lnrr2 = sum_lnrr2 + ( sam_lnrr(idx) - mean_lnrr_sam )^2;
      sum_lnrr3 = sum_lnrr3 + ( sam_lnrr(idx) - mean_lnrr_sam )^3;
   end
   lnrrp2_sam = sum_lnrr2 / count_lnrr;
   lnrrp3_sam = sum_lnrr3 / count_lnrr;
   Sk_lnrr_sam = lnrrp3_sam / lnrrp2_sam^1.5;
else % count_lnrr = 0
   mean_lnrr_sam = -realmax('single');
   lnrrp2_sam = 0.0;
   lnrrp3_sam = 0.0;
   Sk_lnrr_sam = 0.0;
end % count_lnrr > 0
count_lnNr = 0;
for idx = 1:1:nx_sam*ny_sam
   if ( sam_var_lev(idx_3D_Nr,idx) > 0.0 )
      count_lnNr = count_lnNr + 1;
      sam_lnNr(count_lnNr) = log( sam_var_lev(idx_3D_Nr,idx) );
   end
end
% Calculate the mean, variance, 3rd-order central moment, and skewness of
% ln Nr for all Nr in precip. (where Nr > 0) for SAM LES.
if ( count_lnNr > 0 )
   mean_lnNr_sam = mean( sam_lnNr );
   sum_lnNr2 = 0.0;
   sum_lnNr3 = 0.0;
   for idx = 1:1:count_lnNr
      sum_lnNr2 = sum_lnNr2 + ( sam_lnNr(idx) - mean_lnNr_sam )^2;
      sum_lnNr3 = sum_lnNr3 + ( sam_lnNr(idx) - mean_lnNr_sam )^3;
   end
   lnNrp2_sam = sum_lnNr2 / count_lnNr;
   lnNrp3_sam = sum_lnNr3 / count_lnNr;
   Sk_lnNr_sam = lnNrp3_sam / lnNrp2_sam^1.5;
else % count_lnNr = 0
   mean_lnNr_sam = -realmax('single');
   lnNrp2_sam = 0.0;
   lnNrp3_sam = 0.0;
   Sk_lnNr_sam = 0.0;
end % count_lnNr > 0
