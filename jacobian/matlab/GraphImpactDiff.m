
%//////////////////////////////////////////////////////////////////////
% GraphImpactDiff.m
%//////////////////////////////////////////////////////////////////////
%
% Graph difference between adjacent terms
% Assumes the use of 8 impact files and makes 7 plots.
%
% Constant terms:
% 1:C1, 2:C2, 3:C4, 4:C5, 5:C6, 6:C7, 7:C8, 8:C10, 9:C11, 10:nu1, 11:nu2, 12:nu6, 13:nu8, 14:beta,
% 15:gamme_coef, 16:c_K, & 17:mu

x = 1; % constant term

i116    = importdata('impact116.txt');
i108    = importdata('impact108.txt');
i104    = importdata('impact104.txt');
i102    = importdata('impact102.txt');
i101    = importdata('impact101.txt');
i1005   = importdata('impact1005.txt');
i10025  = importdata('impact10025.txt');
i100125 = importdata('impact100125.txt');

i1 = abs(i116(:,x)   - i108(:,x));
i2 = abs(i108(:,x)   - i104(:,x));
i3 = abs(i104(:,x)   - i102(:,x));
i4 = abs(i102(:,x)   - i101(:,x));
i5 = abs(i101(:,x)   - i1005(:,x));
i6 = abs(i1005(:,x)  - i10025(:,x));
i7 = abs(i10025(:,x) - i100125(:,x));

iall = [ i1, i2, i3, i4, i5, i6, i7];

index = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ];

%plot( index, i116(:,1), index, i108(:,1), index, i104(:,1), index, i102(:,1), index, i101(:,1), index, i1005(:,1), index, i10025(:,1), index, i100125(:,1) )
plot( index, iall, 'o' )

title('Adjacent row differences in an impact matrix', 'FontSize', 12 )
grid on;
xlabel('\Delta Variables');
ylabel('\Delta C1');
set (gca, 'XTick', index );
set (gca, 'XTickLabel', {'cf','rcm', 'rtm','thlm','um','vm','wp3','wp2','rtp2','thlp2','rtpthlp','wprtp','wpthlp'});
legend('\Delta factor = 1.16 - \Delta factor = 1.08', '\Delta factor = 1.08 - \Delta factor = 1.04', '\Delta factor = 1.04 - \Delta factor = 1.02', '\Delta factor = 1.02 - \Delta factor = 1.01', '\Delta factor = 1.01 - \Delta factor = 1.005', '\Delta factor = 1.005 - \Delta factor = 1.0025','\Delta factor = 1.0025 - \Delta factor = 1.00125' );