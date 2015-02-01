% $Id$
function plot_CLUBB_PDF_LES_pts_ln_L( var_x_LES, ...
                                      nx_LES_grid, ny_LES_grid, ...
                                      num_lnx_pts, ...
                                      mu_x_1_n, mu_x_2_n, ...
                                      sigma_x_1_n, sigma_x_2_n, ...
                                      precip_frac_1, precip_frac_2, ...
                                      mixt_frac, var_lnx_label, ...
                                      case_name, field_plotted, ...
                                      time_plotted, ...
                                      altitude_plotted, ...
                                      note )
                                
% Path to PDF analysis functions.
addpath ( '../../utilities/PDF_analysis', '-end' )

scr_size= get(0,'ScreenSize');
fig_height = 0.9*scr_size(4);
fig_width  = fig_height;
figure('Position',[ 0 0 fig_width fig_height ])

% Histogram and marginal for ln x, where x > 0 (in-precip), a variable
% that is distributed normally in each PDF component.
pos_count_lnx = 0;
% The variable var_lnx_LES_pos contains only ln x, where x > 0 results
% from the LES.
for i = 1:1:nx_LES_grid*ny_LES_grid
   if ( var_x_LES(i) > 0.0 )
      pos_count_lnx = pos_count_lnx + 1;
      var_lnx_LES_pos(pos_count_lnx) = log( var_x_LES(i) );
   end
end
min_lnx = min( var_lnx_LES_pos );
max_lnx = max( var_lnx_LES_pos );
num_lnx_divs = num_lnx_pts + 1;
delta_lnx = ( max_lnx - min_lnx ) / num_lnx_pts;
lnx_divs = zeros( num_lnx_divs, 1 );
lnx = zeros( num_lnx_pts, 1 );
P_lnx = zeros( num_lnx_pts, 1 );
% Calculate the ln x "division points" (edges of the bins).
for i = 1:1:num_lnx_divs
   lnx_divs(i) = min_lnx + delta_lnx * ( i - 1 );
end
% Calculate the ln x points (center of each bin).
for i = 1:1:num_lnx_pts
   lnx(i) = 0.5 * ( lnx_divs(i) + lnx_divs(i+1) );
end
% CLUBB's PDF marginal for ln x.
for i = 1:1:num_lnx_pts
   P_lnx(i) ...
   = mixt_frac * precip_frac_1 ...
     * PDF_comp_Normal( lnx(i), mu_x_1_n, sigma_x_1_n ) ...
     + ( 1.0 - mixt_frac ) * precip_frac_2 ...
       * PDF_comp_Normal( lnx(i), mu_x_2_n, sigma_x_2_n );
end
% Centerpoints and counts for each bin (from LES results) for ln x.
binranges_lnx = lnx;
[bincounts_lnx] = histc( var_lnx_LES_pos, binranges_lnx );
% Plot normalized histogram of LES results for ln x.
bar( binranges_lnx, bincounts_lnx ...
                    / ( nx_LES_grid * ny_LES_grid * delta_lnx ), 1.0, ...
     'r', 'EdgeColor', 'r' );
hold on
% Plot normalized PDF of x for CLUBB.
plot( lnx, P_lnx, '-b', 'LineWidth', 2 )
hold off
% Set the range of the plot on both the x-axis and y-axis.
xlim( [ min_lnx max_lnx ] )
ylim( [ 0 max( max(bincounts_lnx) ...
               / ( nx_LES_grid * ny_LES_grid * delta_lnx ), ...
          max(P_lnx) ) ] );
xlabel( var_lnx_label )
ylabel( [ 'P(', field_plotted, ')' ] )
legend( 'LES', 'CLUBB', 'Location', 'NorthEast' )
grid on
box on
% Figure title and other important information.
title( [ case_name, '; ', field_plotted, '; ', time_plotted, '; ', ...
         altitude_plotted, '; ', note ] );
