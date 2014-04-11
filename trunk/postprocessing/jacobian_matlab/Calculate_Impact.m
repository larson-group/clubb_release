%//////////////////////////////////////////////////////////////////////
% CalculateImpact.m
%//////////////////////////////////////////////////////////////////////
%
% Calculate the difference between the calculation of the impact
% from the jacobian with MATLAB vs Fortran
%
% Constant terms:
% 1:C1, 2:C2, 3:C4, 4:C5, 5:C6, 6:C7, 7:C8, 8:C10, 9:C11, 10:nu1, 11:nu2, 12:nu6, 13:nu8, 14:beta,
% 15:gamme_coef, 16:c_K, & 17:mu

x = 1;  %constant term

i116    = importdata('impact116.txt');
i108    = importdata('impact108.txt');
i104    = importdata('impact104.txt');
i102    = importdata('impact102.txt');
i101    = importdata('impact101.txt');
i1005   = importdata('impact1005.txt');
i10025  = importdata('impact10025.txt');
i100125 = importdata('impact100125.txt');

j116    = importdata('jacobian116.txt');
j108    = importdata('jacobian108.txt');
j104    = importdata('jacobian104.txt');
j102    = importdata('jacobian102.txt');
j101    = importdata('jacobian101.txt');
j1005   = importdata('jacobian1005.txt');
j10025  = importdata('jacobian10025.txt');
j100125 = importdata('jacobian100125.txt');

constants = zeros(17,13);
for k =1:13
constants(:,k) = importdata('constants.txt');
end
 
j1 = abs( (j116    .* constants' .* (1.16-1)  ) - i116 );
j2 = abs( (j108    .* constants' .* (1.08-1) ) - i108 ); 
j3 = abs( (j104    .* constants' .* (1.04-1) ) - i104 ); 
j4 = abs( (j102    .* constants' .* (1.02-1) ) - i102 ); 
j5 = abs( (j101    .* constants' .* (1.01-1) ) - i101 ); 
j6 = abs( (j1005   .* constants' .* (1.005-1) ) - i1005 ); 
j7 = abs( (j10025  .* constants' .* (1.0025-1) ) - i10025 ); 
j8 = abs( (j100125 .* constants' .* (1.00125-1) ) - i100125 ); 
jall = [j1(:,x), j2(:,x), j3(:,x), j4(:,x), j5(:,x), j5(:,x), j6(:,x), j7(:,x) j8(:,x)];
index = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ]';

plot( index, jall, 'o' );
grid on;
xlabel('\Delta Variables');
ylabel('The Absolute Error between Fortran and MATLAB');
set (gca, 'XTick', index );
set (gca, 'XTickLabel', {'cf','rcm', 'rtm','thlm','um','vm','wp3','wp2','rtp2','thlp2','rtpthlp','wprtp','wpthlp'})
legend('\Delta factor = 1.16', '\Delta factor = 1.08', '\Delta factor = 1.04', '\Delta factor = 1.02', '\Delta factor = 1.01', '\Delta factor = 1.005', '\Delta factor = 1.0025','\Delta factor = 1.00125' );