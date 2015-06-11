% $Id$
function plot_CLUBB_PDF_LES_pts_NL( var_x_LES, ...
                                    var_y_LES, ...
                                    nx_LES_grid, ny_LES_grid, ...
                                    num_x_pts, num_y_pts, num_contours, ...
                                    num_std_devs_min_contour, ... 
                                    mu_x_1, mu_x_2, mu_y_1_n, mu_y_2_n, ...
                                    sigma_x_1, sigma_x_2, sigma_y_1_n, ...
                                    sigma_y_2_n, corr_x_y_1_n, ...
                                    corr_x_y_2_n, precip_frac_1, ...
                                    precip_frac_2, mixt_frac, ...
                                    num_clubb_files, ...
                                    var_x_label, var_y_label, ...
                                    case_name, fields_plotted, ...
                                    time_plotted, ...
                                    altitude_plotted, ...
                                    note )
                                
% Path to PDF analysis functions.
addpath ( '../../utilities/PDF_analysis', '-end' )

scr_size= get(0,'ScreenSize');
fig_height = 0.9*scr_size(4);
fig_width  = fig_height;
figure('Position',[ 0 0 fig_width fig_height ])

subplot( 'Position', [ 0.1 0.75 0.55 0.2 ] )

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
bar( binranges_x, bincounts_x / ( nx_LES_grid * ny_LES_grid * delta_x ), ...
     1.0, 'r', 'EdgeColor', 'r' );
% Plot PDF of x for CLUBB.
for clubb_idx = 1:1:num_clubb_files
   hold on
   plot( x, P_x(clubb_idx,:), '-b', 'LineWidth', 2 )
end % clubb_idx = 1:1:num_clubb_files
hold off
% Set the range of the plot on both the x-axis and y-axis.
xlim( [ min_x max_x ] );
ylim( [ 0 max( max(bincounts_x) / ( nx_LES_grid * ny_LES_grid * delta_x ), ...
               max(max(P_x)) ) ] );
%xlabel( var_x_label )
legend( 'LES', 'CLUBB', 'Location', 'NorthEast' )
grid on
box on

subplot( 'Position', [ 0.75 0.1 0.2 0.55 ] )

% Histogram and marginal for y, a variable that is distributed
% delta-lognormally in each PDF component.
pos_count_y = 0;
zero_count_y = 0;
% The variable var_y_LES_pos contains only y > 0 results from the LES.
for j = 1:1:nx_LES_grid*ny_LES_grid
   if ( var_y_LES(j) > 0.0 )
      pos_count_y = pos_count_y + 1;
      var_y_LES_pos(pos_count_y) = var_y_LES(j);
   else
      zero_count_y = zero_count_y + 1;
   end
end
precip_frac_y_LES = pos_count_y / ( nx_LES_grid * ny_LES_grid );
for clubb_idx = 1:1:num_clubb_files
   precip_frac_clubb(clubb_idx) ...
   = mixt_frac(clubb_idx) * precip_frac_1(clubb_idx) ...
     + ( 1.0 - mixt_frac(clubb_idx) ) * precip_frac_2(clubb_idx);
end % clubb_idx = 1:1:num_clubb_files
min_y = min( var_y_LES_pos );
max_y = max( var_y_LES_pos );
num_y_divs = num_y_pts + 1;
delta_y = ( max_y - min_y ) / num_y_pts;
y_divs = zeros( num_y_divs, 1 );
y = zeros( num_y_pts, 1 );
P_y = zeros( num_y_pts, 1 );
% Calculate the y "division points" (edges of the bins).
for j = 1:1:num_y_divs
   y_divs(j) = min_y + delta_y * ( j - 1 );
end
% Calculate the y points (center of each bin).
for j = 1:1:num_y_pts
   y(j) = 0.5 * ( y_divs(j) + y_divs(j+1) );
end
% CLUBB's PDF marginal for y (lognormal or in-precip portion).
for clubb_idx = 1:1:num_clubb_files
   for j = 1:1:num_y_pts
      P_y(clubb_idx,j) ...
      = mixt_frac(clubb_idx) * precip_frac_1(clubb_idx) ...
        * PDF_comp_Lognormal( y(j), mu_y_1_n(clubb_idx), ...
                              sigma_y_1_n(clubb_idx) ) ...
        + ( 1.0 - mixt_frac(clubb_idx) ) * precip_frac_2(clubb_idx) ...
          * PDF_comp_Lognormal( y(j), mu_y_2_n(clubb_idx), ...
                                sigma_y_2_n(clubb_idx) );
   end
end % clubb_idx = 1:1:num_clubb_files
% Centerpoints and counts for each bin (from LES results) for y, where
% y > 0.
binranges_y = y;
[bincounts_y] = histc( var_y_LES, binranges_y );
% Plot normalized histogram of LES results for y, where y > 0.
bar( binranges_y, bincounts_y / ( nx_LES_grid * ny_LES_grid * delta_y ), ...
     1.0, 'r', 'EdgeColor', 'r' );
if ( precip_frac_y_LES < 1.0 )
   hold on
   % Plot normalized histogram of LES results for y, where y = 0.
   arrow_vec = get( gca, 'Position' );
   annotation( 'arrow', [ arrow_vec(1) arrow_vec(1)+arrow_vec(3) ], ...
               [ arrow_vec(2) arrow_vec(2) ], 'Color', 'r', ...
               'LineStyle', '-', 'LineWidth', 3, ...
               'HeadLength', 11, 'HeadWidth', 11 );
end
% Plot PDF of y for CLUBB, where y > 0.
for clubb_idx = 1:1:num_clubb_files
   hold on
   plot( y, P_y(clubb_idx,:), '-b', 'LineWidth', 2 )
end % clubb_idx = 1:1:num_clubb_files
% Plot PDF of y for CLUBB, where y = 0.
for clubb_idx = 1:1:num_clubb_files
   if ( precip_frac_clubb(clubb_idx) < 1.0 )
      hold on
      arrow_vec = get( gca, 'Position' );
      annotation( 'arrow', [ arrow_vec(1) arrow_vec(1)+arrow_vec(3) ], ...
                  [ arrow_vec(2) arrow_vec(2) ], 'Color', 'b', ...
                  'LineStyle', '-', 'LineWidth', 2, ...
                  'HeadLength', 7, 'HeadWidth', 7 );
   end
end % clubb_idx = 1:1:num_clubb_files
hold off
% Set the range of the plot on both the x-axis and y-axis.
if ( any( precip_frac_clubb < 1.0 ) || precip_frac_y_LES < 1.0 )
   xlim([0 max_y]);
else
   xlim([min_y max_y]);
end
ylim( [ 0 max( max(bincounts_y) / ( nx_LES_grid * ny_LES_grid * delta_y ), ...
               max(max(P_y)) ) ] );
%xlabel( var_y_label )
% Reverse the plot and turn the plot 90 degrees clockwise.
set(gca, 'xdir', 'reverse')
view( [90, 90] )
grid on
box on

subplot( 'Position', [ 0.1 0.1 0.55 0.55 ] )

% Scatterplot and contour plot for x and y.
P_xy = zeros( num_clubb_files, num_y_pts, num_x_pts );
% CLUBB's PDF of x and y, where y > 0.
for clubb_idx = 1:1:num_clubb_files
   for i = 1:1:num_x_pts
      for j = 1:1:num_y_pts
         P_xy(clubb_idx,j,i) ...
         = mixt_frac(clubb_idx) * precip_frac_1(clubb_idx) ...
           * PDF_comp_bivar_NL( x(i), y(j), ...
                                mu_x_1(clubb_idx), ...
                                mu_y_1_n(clubb_idx), ...
                                sigma_x_1(clubb_idx), ...
                                sigma_y_1_n(clubb_idx), ...
                                corr_x_y_1_n(clubb_idx) ) ...
           + ( 1.0 - mixt_frac(clubb_idx) ) * precip_frac_2(clubb_idx) ...
             * PDF_comp_bivar_NL( x(i), y(j), ...
                                  mu_x_2(clubb_idx), ...
                                  mu_y_2_n(clubb_idx), ...
                                  sigma_x_2(clubb_idx), ...
                                  sigma_y_2_n(clubb_idx), ...
                                  corr_x_y_2_n(clubb_idx) );
      end
   end
end % clubb_idx = 1:1:num_clubb_files
% Find the minimum contour to plot for CLUBB results.
% Sample a value that is a set number of standard deviations away from the
% bivariate mean for each variable for each component.
num_std_devs = num_std_devs_min_contour / sqrt(2.0);
for clubb_idx = 1:1:num_clubb_files
   % For the 1st PDF component, sample P_xy at sqrt(2)*num_std_devs, which
   % is num_std_devs away from the 1st PDF component bivariate mean in both
   % x and y.
   if ( sigma_x_1(clubb_idx) > 0.0 && sigma_y_1_n(clubb_idx) > 0.0 )
      if ( mu_x_1(clubb_idx) >= mu_x_2(clubb_idx) )
         x_1_test = min( mu_x_1(clubb_idx) ...
                         + num_std_devs * sigma_x_1(clubb_idx), max_x );
      else % mu_x_1(clubb_idx) < mu_x_2(clubb_idx)
         x_1_test = max( mu_x_1(clubb_idx) ...
                         - num_std_devs * sigma_x_1(clubb_idx), min_x );
      end
      mu_y_1 = exp( mu_y_1_n(clubb_idx) + 0.5 * sigma_y_1_n(clubb_idx)^2 );
      sigma_y_1 = sqrt( ( exp( sigma_y_1_n(clubb_idx)^2 ) - 1.0 ) ...
                        * exp( 2.0 * mu_y_1_n(clubb_idx) ...
                               + sigma_y_1_n(clubb_idx)^2 ) );
      y_1_test = min( mu_y_1 + num_std_devs * sigma_y_1, max_y );
      P_low_1 ...
      = mixt_frac(clubb_idx) * precip_frac_1(clubb_idx) ...
        * PDF_comp_bivar_NL( x_1_test, y_1_test, ...
                             mu_x_1(clubb_idx), mu_y_1_n(clubb_idx), ...
                             sigma_x_1(clubb_idx), ...
                             sigma_y_1_n(clubb_idx), ...
                             corr_x_y_1_n(clubb_idx) ) ...
        + ( 1.0 - mixt_frac(clubb_idx) ) * precip_frac_2(clubb_idx) ...
          * PDF_comp_bivar_NL( x_1_test, y_1_test, ...
                               mu_x_2(clubb_idx), mu_y_2_n(clubb_idx), ...
                               sigma_x_2(clubb_idx), ...
                               sigma_y_2_n(clubb_idx), ...
                               corr_x_y_2_n(clubb_idx) );
   else
      P_low_1 = -1.0;
   end
   % For the 2nd PDF component, sample P_xy at sqrt(2)*num_std_devs, which
   % is num_std_devs away from the 2nd PDF component bivariate mean in both
   % x and y.
   if ( sigma_x_2(clubb_idx) > 0.0 && sigma_y_2_n(clubb_idx) > 0.0 )
      if ( mu_x_1(clubb_idx) >= mu_x_2(clubb_idx) )
         x_2_test = max( mu_x_2(clubb_idx) ...
                         - num_std_devs * sigma_x_2(clubb_idx), min_x );
      else % mu_x_1(clubb_idx) < mu_x_2(clubb_idx)
         x_2_test = min( mu_x_2(clubb_idx) ...
                         + num_std_devs * sigma_x_2(clubb_idx), max_x );
      end
      mu_y_2 = exp( mu_y_2_n(clubb_idx) + 0.5 * sigma_y_2_n(clubb_idx)^2 );
      sigma_y_2 = sqrt( ( exp( sigma_y_2_n(clubb_idx)^2 ) - 1.0 ) ...
                        * exp( 2.0 * mu_y_2_n(clubb_idx) ...
                               + sigma_y_2_n(clubb_idx)^2 ) );
      y_2_test = min( mu_y_2 + num_std_devs * sigma_y_2, max_y );
      P_low_2 ...
      = mixt_frac(clubb_idx) * precip_frac_1(clubb_idx) ...
        * PDF_comp_bivar_NL( x_2_test, y_2_test, ...
                             mu_x_1(clubb_idx), mu_y_1_n(clubb_idx), ...
                             sigma_x_1(clubb_idx), ...
                             sigma_y_1_n(clubb_idx), ...
                             corr_x_y_1_n(clubb_idx) ) ...
        + ( 1.0 - mixt_frac(clubb_idx) ) * precip_frac_2(clubb_idx) ...
          * PDF_comp_bivar_NL( x_2_test, y_2_test, ...
                               mu_x_2(clubb_idx), mu_y_2_n(clubb_idx), ...
                               sigma_x_2(clubb_idx), ...
                               sigma_y_2_n(clubb_idx), ...
                               corr_x_y_2_n(clubb_idx) );
   else
      P_low_2 = -1.0;
   end
   % The smallest contour is the smaller of P_low_1 and P_low_2.
   if ( P_low_1 > 0.0 && P_low_2 > 0.0 )
      contour_low = min( P_low_1, P_low_2 );
   elseif ( P_low_1 > 0.0 ) % P_low_2 <= 0.0
      contour_low = P_low_1;
   elseif ( P_low_2 > 0.0 ) % P_low_1 <= 0.0
      contour_low = P_low_2;
   else % P_low_1 <= 0.0 && P_low_2 <= 0.0
      contour_low = 0.0;
   end
   % In a scenario where contour_low is so far out that it doesn't plot
   % well, use an alternative.  Sort all the values of P_xy and plot the
   % lowest contour at the 10th percentile of P_xy.
   P_xy_vals = zeros( num_x_pts*num_y_pts, 1 );
   idx = 0;
   for i = 1:1:num_x_pts
      for j = 1:1:num_y_pts
         idx = idx + 1;
         P_xy_vals(idx) = P_xy(clubb_idx,j,i);
      end
   end
   P_xy_sort = sort( P_xy_vals );
   P_xy_low = P_xy_sort( floor( 0.10*num_x_pts*num_y_pts ) );
   if ( contour_low < P_xy_low )
      contour_low = P_xy_low;
   end
   % Find the maximum contour to plot for CLUBB results.
   P_xy_high = max( max( max( P_xy ) ) );
   contour_high = P_xy_high;
   % Set the contour vector.
   delta_contour = ( contour_high - contour_low ) / ( num_contours - 1 );
   for ci = 1:1:num_contours
      contours(clubb_idx,ci) = contour_low + delta_contour * ( ci - 1 );
   end
end % clubb_idx = 1:1:num_clubb_files
% Scatterplot of LES results for x and y.
scatter( var_x_LES, var_y_LES, ...
         'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r' )
% Contour plot of the PDF of x and y for CLUBB, where y > 0.
for clubb_idx = 1:1:num_clubb_files
   hold on
   contour( x, y, reshape( P_xy(clubb_idx,:,:), num_y_pts, num_x_pts ), ...
            contours(clubb_idx,:), 'Linewidth', 1.5 )
end % clubb_idx = 1:1:num_clubb_files
% Plot of the PDF of x and y for CLUBB, where y = 0.
for clubb_idx = 1:1:num_clubb_files
   hold on
   if ( precip_frac_clubb(clubb_idx) < 1.0 )
      plot( [ min_x max_x ], [ 0 0 ], '-b', 'LineWidth', 1.5 );
   end
end % clubb_idx = 1:1:num_clubb_files
hold off
legend( 'LES', 'CLUBB', 'Location', 'NorthEast' )
% Set the range of the plot on both the x-axis and y-axis.
xlim([min_x max_x])
if ( any( precip_frac_clubb ) < 1.0 || precip_frac_y_LES < 1.0 )
   ylim([0 max_y]);
else
   ylim([min_y max_y])
end
xlabel( var_x_label )
ylabel( var_y_label )
grid on
box on

% Figure title and other important information.
figure_title_text = { case_name; fields_plotted; time_plotted; ...
                      altitude_plotted; note };

axes( 'Position', [0 0 1 1], 'Xlim', [0 1], 'Ylim', [0 1], ...
      'Box', 'off', 'Visible', 'off', 'Units', 'normalized', ...
      'clipping', 'off' );

text( 0.825, 0.95, figure_title_text, ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'top' );
