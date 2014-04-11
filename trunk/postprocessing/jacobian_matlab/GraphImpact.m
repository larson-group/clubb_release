
%//////////////////////////////////////////////////////////////////////
% GraphImpact2.m
%//////////////////////////////////////////////////////////////////////
%
% Graph the first term of the Jacobian Matrix
%
% Constant terms:
% 1:C1, 2:C2, 3:C4, 4:C5, 5:C6, 6:C7, 7:C8, 8:C10, 9:C11, 10:nu1, 11:nu2, 12:nu6, 13:nu8, 14:beta,
% 15:gamme_coef, 16:c_K, & 17:mu


x = 1;   % the constant term to plot


i116    = importdata('impact116.txt');
i108    = importdata('impact108.txt');
i104    = importdata('impact104.txt');
i102    = importdata('impact102.txt');
i101    = importdata('impact101.txt');
i1005   = importdata('impact1005.txt');
i10025  = importdata('impact10025.txt');
i100125 = importdata('impact100125.txt');

impact = [ i100125(:,x) , i10025(:,x), i1005(:,x), i101(:,x), i102(:,x), i104(:,x), i108(:,x), i116(:,x) ];

index = [ 1.00125, 1.0025, 1.005, 1.01, 1.02, 1.04, 1.08, 1.16 ];


plot( index, impact, '-o' )

title('Plot of the first column of a impact matrix', 'FontSize', 12 )
grid on;

xlabel('\Delta C1');
ylabel('\Delta Variables');
legend('cf','rcm', 'rtm','\Theta_l m','u_m','v_m','w\prime^3','w\prime^2','r_t\prime^2','\Theta_l\prime^2','rt\prime\Theta_l\prime','w\prime rt\prime','w\prime\Theta_l\prime');
