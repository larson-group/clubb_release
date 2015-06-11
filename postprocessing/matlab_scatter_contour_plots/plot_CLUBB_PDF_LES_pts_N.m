% $Id$
function plot_CLUBB_PDF_LES_pts_N( var_x_LES, ...
                                   nx_LES_grid, ny_LES_grid, ...
                                   num_x_pts, log_Px_plot, ...
                                   num_clubb_files, ...
                                   mu_x_1, mu_x_2, ...
                                   sigma_x_1, sigma_x_2, ...
                                   mixt_frac, var_x_label, ...
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

subplot( 'Position', [ 0.1 0.1 0.8 0.8 ] )

% Histogram and marginal for x, a variable that is distributed normally in
% each PDF component.
min_x = min( var_x_LES );
max_x = max( var_x_LES );
num_x_divs = num_x_pts + 1;
delta_x = ( max_x - min_x ) / num_x_pts;
x_divs = zeros( num_x_divs, 1 );
x = zeros( num_x_pts, 1 );
P_x = zeros( num_clubb_files, num_x_pts );
% Calculate the x "division points" (edges of the bins).
for i = 1:1:num_x_divs
   x_divs(i) = min_x + delta_x * ( i - 1 );
end
% Calculate the x points (center of each bin).
for i = 1:1:num_x_pts
   x(i) = 0.5 * ( x_divs(i) + x_divs(i+1) );
end
% CLUBB's PDF marginal for x.
for clubb_idx = 1:1:num_clubb_files
   for i = 1:1:num_x_pts
      P_x(clubb_idx,i) ...
      = mixt_frac(clubb_idx) ...
        * PDF_comp_Normal( x(i), mu_x_1(clubb_idx), ...
                           sigma_x_1(clubb_idx) ) ...
        + ( 1.0 - mixt_frac(clubb_idx) ) ...
          * PDF_comp_Normal( x(i), mu_x_2(clubb_idx), ...
                             sigma_x_2(clubb_idx) );
   end
end % clubb_idx = 1:1:num_clubb_files
% Centerpoints and counts for each bin (from LES results) for x.
binranges_x = x;
[bincounts_x] = histc( var_x_LES, binranges_x );
% Plot normalized histogram of LES results for x.
if ( log_Px_plot )
   semilogy( binranges_x, ...
             max( bincounts_x, 1.0e-300 ) ...
             / ( nx_LES_grid * ny_LES_grid * delta_x ), ...
             '-r', 'LineWidth', 3 );
else
   bar( binranges_x, ...
        bincounts_x / ( nx_LES_grid * ny_LES_grid * delta_x ), ...
        1.0, 'r', 'EdgeColor', 'r' );
end
% Plot PDF of x for CLUBB.
for clubb_idx = 1:1:num_clubb_files
   hold on
   if ( log_Px_plot )
      semilogy( x, P_x(clubb_idx,:), '-b', 'LineWidth', 2 )
   else
      plot( x, P_x(clubb_idx,:), '-b', 'LineWidth', 2 )
   end
end % clubb_idx = 1:1:num_clubb_files
hold off
% Set the range of the plot on both the x-axis and y-axis.
xlim( [ min_x max_x ] )
if ( log_Px_plot )
   min_pos_bincounts_x = max(bincounts_x);
   for i = 1:1:num_x_pts
      if ( bincounts_x(i) > 0 )
         min_pos_bincounts_x = min( min_pos_bincounts_x, bincounts_x(i) );
      end
   end
   ylim( [ 0.001 * min_pos_bincounts_x ...
                   / ( nx_LES_grid * ny_LES_grid * delta_x ) ...
           max( max(bincounts_x) ...
                / ( nx_LES_grid * ny_LES_grid * delta_x ), ...
                max(max(P_x)) ) ] );
else
   ylim( [ 0 max( max(bincounts_x) ...
                  / ( nx_LES_grid * ny_LES_grid * delta_x ), ...
                  max(max(P_x)) ) ] );
end
xlabel( var_x_label )
ylabel( [ 'P( ', field_plotted, ' )' ] )
legend( 'LES', 'CLUBB', 'Location', 'NorthEast' )
grid on
box on

% Figure title and other important information.
title( [ case_name, '; ', field_plotted, '; ', time_plotted, '; ', ...
         altitude_plotted, '; ', note ] );
