
%//////////////////////////////////////////////////////////////////////
% GraphJacobian.m
%//////////////////////////////////////////////////////////////////////
%
% Graph the first term of the Jacobian Matrix
%
% Constant terms:
% 1:C1, 2:C2, 3:C4, 4:C5, 5:C6, 6:C7, 7:C8, 8:C10, 9:C11, 10:nu1, 11:nu2, 12:nu6, 13:nu8, 14:beta,
% 15:gamme_coef, 16:c_K, & 17:mu


x = 1;   % the constant term to plot


j116    = importdata('jacobian116.txt');
j108    = importdata('jacobian108.txt');
j104    = importdata('jacobian104.txt');
j102    = importdata('jacobian102.txt');
j101    = importdata('jacobian101.txt');
j1005   = importdata('jacobian1005.txt');
j10025  = importdata('jacobian10025.txt');
j100125 = importdata('jacobian100125.txt');

jacobian = [ j100125(:,x) , j10025(:,x), j1005(:,x), j101(:,x), j102(:,x), j104(:,x), j108(:,x), j116(:,x) ];

index = [ 1.00125, 1.0025, 1.005, 1.01, 1.02, 1.04, 1.08, 1.16 ];

%index = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ]';
plot( index, j116(:,x), index, j108(:,x), index, j104(:,x), index, j102(:,x), index, j101(:,x), index, j1005(:,x), index, j10025(:,x), index, j100125(:,x) )

title('Plot of the first column of a jacobian matrix', 'FontSize', 12 )
grid on;
xlabel('\Delta C1');
ylabel('\Delta Variables');
set (gca, 'XTick', index );
set (gca, 'XTickLabel', {'cf','rcm', 'rtm','thlm','um','vm','wp3','wp2','rtp2','thlp2','rtpthlp','wprtp','wpthlp'})
legend('cf','rcm', 'rtm','\Theta_l m','u_m','v_m','w\prime^3','w\prime^2','r_t\prime^2','\Theta_l\prime^2','rt\prime\Theta_l\prime','w\prime rt\prime','w\prime\Theta_l\prime');
