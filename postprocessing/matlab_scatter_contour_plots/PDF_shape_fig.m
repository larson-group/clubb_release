% $Id$
function PDF_shape_fig

% Path to PDF analysis functions.
addpath ( '../../utilities/PDF_analysis', '-end' )

mean_hm = 1.0e-5;
hmp2 = 5.0e-10;

mixt_frac_1 = 0.45;
mixt_frac_2 = 0.6;

precip_frac = 0.75;
precip_frac_1 = 0.9;
precip_frac_2 = ( precip_frac - mixt_frac_1 * precip_frac_1 ) / mixt_frac_2;

hm_high = 2.5e-5;
hm_low  = 1.0e-8;
num_pts = 10001;

% DDL PDF
mu_hm_1 = 1.65 * mean_hm / precip_frac
mu_hm_2 = ( mean_hm - mixt_frac_1 * precip_frac_1 * mu_hm_1 ) ...
          / ( mixt_frac_2 * precip_frac_2 )

R = ( ( hmp2 + mean_hm^2 ) ...
      / ( mixt_frac_1 * precip_frac_1 * mu_hm_1^2 ...
          + mixt_frac_2 * precip_frac_2 * mu_hm_2^2 ) ) - 1.0

sigma_hm_1 = sqrt(R) * mu_hm_1;
sigma_hm_2 = sqrt(R) * mu_hm_2;

mu_hm_1_n = log( mu_hm_1 * ( 1.0 + sigma_hm_1^2 / mu_hm_1^2 )^-0.5 );
sigma_hm_1_n = sqrt( log( 1.0 + sigma_hm_1^2 / mu_hm_1^2 ) );
mu_hm_2_n = log( mu_hm_2 * ( 1.0 + sigma_hm_2^2 / mu_hm_2^2 )^-0.5 );
sigma_hm_2_n = sqrt( log( 1.0 + sigma_hm_2^2 / mu_hm_2^2 ) );

hm_DDL  = zeros( num_pts, 1 );
P_DDL_1 = zeros( num_pts, 1 );
P_DDL_2 = zeros( num_pts, 1 );
P_DDL   = zeros( num_pts, 1 );
for iter=1:1:num_pts

   hm_DDL(iter) = hm_low + (iter-1) * ( hm_high - hm_low ) / ( num_pts - 1 );

   P_DDL_1(iter) ...
   = mixt_frac_1 * precip_frac_1 ...
     * PDF_comp_Lognormal( hm_DDL(iter), mu_hm_1_n, sigma_hm_1_n );

   P_DDL_2(iter) ...
   = mixt_frac_2 * precip_frac_2 ...
     * PDF_comp_Lognormal( hm_DDL(iter), mu_hm_2_n, sigma_hm_2_n );

   P_DDL(iter) = P_DDL_1(iter) + P_DDL_2(iter);

end % iter = 1:1:num_pts

% DL PDF
mu_hm_1 = mean_hm / precip_frac;
mu_hm_2 = mean_hm / precip_frac;

R = ( ( hmp2 + mean_hm^2 ) ...
      / ( mixt_frac_1 * precip_frac_1 * mu_hm_1^2 ...
          + mixt_frac_2 * precip_frac_2 * mu_hm_2^2 ) ) - 1.0;

sigma_hm_1 = sqrt(R) * mu_hm_1;
sigma_hm_2 = sqrt(R) * mu_hm_2;

mu_hm_1_n = log( mu_hm_1 * ( 1.0 + sigma_hm_1^2 / mu_hm_1^2 )^-0.5 );
sigma_hm_1_n = sqrt( log( 1.0 + sigma_hm_1^2 / mu_hm_1^2 ) );
mu_hm_2_n = log( mu_hm_2 * ( 1.0 + sigma_hm_2^2 / mu_hm_2^2 )^-0.5 );
sigma_hm_2_n = sqrt( log( 1.0 + sigma_hm_2^2 / mu_hm_2^2 ) );

hm_DL  = zeros( num_pts, 1 );
P_DL_1 = zeros( num_pts, 1 );
P_DL_2 = zeros( num_pts, 1 );
P_DL   = zeros( num_pts, 1 );
for iter=1:1:num_pts

   hm_DL(iter) = hm_low + (iter-1) * ( hm_high - hm_low ) / ( num_pts - 1 );

   P_DL_1(iter) ...
   = mixt_frac_1 * precip_frac_1 ...
     * PDF_comp_Lognormal( hm_DL(iter), mu_hm_1_n, sigma_hm_1_n );

   P_DL_2(iter) ...
   = mixt_frac_2 * precip_frac_2 ...
     * PDF_comp_Lognormal( hm_DL(iter), mu_hm_2_n, sigma_hm_2_n );

   P_DL(iter) = P_DL_1(iter) + P_DL_2(iter);

end % iter = 1:1:num_pts

% SL PDF
precip_frac_1 = 1.0;
precip_frac_2 = 1.0;

mu_hm_1 = mean_hm;
mu_hm_2 = mean_hm;

R = ( ( hmp2 + mean_hm^2 ) ...
      / ( mixt_frac_1 * precip_frac_1 * mu_hm_1^2 ...
          + mixt_frac_2 * precip_frac_2 * mu_hm_2^2 ) ) - 1.0;

sigma_hm_1 = sqrt(R) * mu_hm_1;
sigma_hm_2 = sqrt(R) * mu_hm_2;

mu_hm_1_n = log( mu_hm_1 * ( 1.0 + sigma_hm_1^2 / mu_hm_1^2 )^-0.5 );
sigma_hm_1_n = sqrt( log( 1.0 + sigma_hm_1^2 / mu_hm_1^2 ) );
mu_hm_2_n = log( mu_hm_2 * ( 1.0 + sigma_hm_2^2 / mu_hm_2^2 )^-0.5 );
sigma_hm_2_n = sqrt( log( 1.0 + sigma_hm_2^2 / mu_hm_2^2 ) );

hm_SL  = zeros( num_pts, 1 );
P_SL_1 = zeros( num_pts, 1 );
P_SL_2 = zeros( num_pts, 1 );
P_SL   = zeros( num_pts, 1 );
for iter=1:1:num_pts

   hm_SL(iter) = hm_low + (iter-1) * ( hm_high - hm_low ) / ( num_pts - 1 );

   P_SL_1(iter) ...
   = mixt_frac_1 * precip_frac_1 ...
     * PDF_comp_Lognormal( hm_SL(iter), mu_hm_1_n, sigma_hm_1_n );

   P_SL_2(iter) ...
   = mixt_frac_2 * precip_frac_2 ...
     * PDF_comp_Lognormal( hm_SL(iter), mu_hm_2_n, sigma_hm_2_n );

   P_SL(iter) = P_SL_1(iter) + P_SL_2(iter);

end % iter = 1:1:num_pts

% Figure Properties for screen display.
scr_size= get(0,'ScreenSize');
fig_height = 0.2*scr_size(4);
fig_width  = 4.0*fig_height;
figure('Position',[ 0 0 fig_width fig_height ])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperPosition', [ 0.0 0.0 6.0 1.5 ])

subplot(1,3,3)
plot( hm_DDL, P_DDL, '-b', 'LineWidth', 1.5 )
hold on
plot( hm_DDL, P_DDL_1, '--k' )
hold on
plot( hm_DDL, P_DDL_2, ':k' )
hold on
arrow_vec = get( gca, 'Position' );
annotation( 'arrow', [ arrow_vec(1) arrow_vec(1) ], ...
            [ arrow_vec(2) arrow_vec(2)+arrow_vec(4) ], ...
            'Color', 'b', 'LineStyle', '-', 'LineWidth', 1.5, ...
            'HeadLength', 5, 'HeadWidth', 5 );
hold on
annotation( 'arrow', [ arrow_vec(1) arrow_vec(1) ], ...
            [ arrow_vec(2) arrow_vec(2)+arrow_vec(4) ], ...
            'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5, ...
            'HeadLength', 3, 'HeadWidth', 3 );
hold on
annotation( 'arrow', [ arrow_vec(1) arrow_vec(1) ], ...
            [ arrow_vec(2) arrow_vec(2)+arrow_vec(4) ], ...
            'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5, ...
            'HeadLength', 3, 'HeadWidth', 3 );
hold off
xlim( [ 0 max( hm_DDL ) ] )
ylim( [ 0 max( P_SL ) ] )
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'XColor',[ 0.625 0.625 0.625 ])
set(gca,'YColor',[ 0.625 0.625 0.625 ])
box off

subplot(1,3,2)
plot( hm_DL, P_DL, '-g', 'LineWidth', 1.5 )
hold on
plot( hm_DL, P_DL_1, '--k' )
hold on
plot( hm_DL, P_DL_2, ':k' )
hold on
arrow_vec = get( gca, 'Position' );
annotation( 'arrow', [ arrow_vec(1) arrow_vec(1) ], ...
            [ arrow_vec(2) arrow_vec(2)+arrow_vec(4) ], ...
            'Color', 'g', 'LineStyle', '-', 'LineWidth', 1.5, ...
            'HeadLength', 5, 'HeadWidth', 5 );
hold on
annotation( 'arrow', [ arrow_vec(1) arrow_vec(1) ], ...
            [ arrow_vec(2) arrow_vec(2)+arrow_vec(4) ], ...
            'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5, ...
            'HeadLength', 3, 'HeadWidth', 3 );
hold on
annotation( 'arrow', [ arrow_vec(1) arrow_vec(1) ], ...
            [ arrow_vec(2) arrow_vec(2)+arrow_vec(4) ], ...
            'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5, ...
            'HeadLength', 3, 'HeadWidth', 3 );
hold off
xlim( [ 0 max( hm_DL ) ] )
ylim( [ 0 max( P_SL ) ] )
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'XColor',[ 0.625 0.625 0.625 ])
set(gca,'YColor',[ 0.625 0.625 0.625 ])
box off

subplot(1,3,1)
plot( hm_SL, P_SL, '-m', 'LineWidth', 1.5 )
hold on
plot( hm_SL, P_SL_1, '--k' )
hold on
plot( hm_SL, P_SL_2, ':k' )
hold off
xlim( [ 0 max( hm_SL ) ] )
ylim( [ 0 max( P_SL ) ] )
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'XColor',[ 0.625 0.625 0.625 ])
set(gca,'YColor',[ 0.625 0.625 0.625 ])
box off

print( '-depsc', 'output/PDF_shapes_SL_DL_DDL' );
