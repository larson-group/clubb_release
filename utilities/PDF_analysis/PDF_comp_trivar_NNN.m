function [ P_i ] ...
= PDF_comp_trivar_NNN( x1, x2, x3, mu_x1_i, mu_x2_i, mu_x3_i, ...
                       sigma_x1_i, sigma_x2_i, sigma_x3_i, ...
                       rho_x1x2_i, rho_x1x3_i, rho_x2x3_i )

% This function calculates probability density given by the Probability
% Density Function (PDF) for a trivariate normal distribution.  Since the
% overall PDF consists of mulitple components, this function gives the
% results for one component, i.

% x1:  a variable that has a normal marginal distribution in the ith
%      component of the PDF.
%
% x2:  a variable that has a normal marginal distribution in the ith
%      component of the PDF.
%
% x3:  a variable that has a normal marginal distribution in the ith
%      component of the PDF.
%
% mu_x1_i:  the mean of x1 in the ith component of the PDF.  The overall
%           mean of x1 is x1_bar.
%
% mu_x2_i:  the mean of x2 in the ith component of the PDF.  The overall
%           mean of x2 is x2_bar.
%
% mu_x3_i:  the mean of x3 in the ith component of the PDF.  The overall
%           mean of x3 is x3_bar.
%
% sigma_x1_i:  the standard deviation of x1 in the ith component of the
%              PDF.  The variance of x1 in the ith component of the PDF is
%              sigma_x1_i^2, and the overall variance of x1 is x1p2.
%
% sigma_x2_i:  the standard deviation of x2 in the ith component of the
%              PDF.  The variance of x2 in the ith component of the PDF is
%              sigma_x2_i^2, and the overall variance of x2 is x2p2.
%
% sigma_x3_i:  the standard deviation of x3 in the ith component of the
%              PDF.  The variance of x3 in the ith component of the PDF is
%              sigma_x3_i^2, and the overall variance of x3 is x3p2.
%
% rho_x1x2_i:  the correlation between variables x1 and x2 in the ith
%              component of the PDF.  The overall correlation between
%              variables x1 and x2 can be found by:
%                 x1px2p / ( sqrt ( x1p2 * x2p2 ) );
%              where x1px2p is the overall covariance of x1 and x2.
%
% rho_x1x3_i:  the correlation between variables x1 and x3 in the ith
%              component of the PDF.  The overall correlation between
%              variables x1 and x3 can be found by:
%                 x1px3p / ( sqrt ( x1p2 * x3p2 ) );
%              where x1px3p is the overall covariance of x1 and x3.
%
% rho_x2x3_i:  the correlation between variables x2 and x3 in the ith
%              component of the PDF.  The overall correlation between
%              variables x2 and x3 can be found by:
%                 x2px3p / ( sqrt ( x2p2 * x3p2 ) );
%              where x2px3p is the overall covariance of x2 and x3.

C_1 = sqrt( 1.0 - ( rho_x1x2_i^2 + rho_x1x3_i^2 + rho_x2x3_i^2 ) ...
            + 2.0 * rho_x1x2_i * rho_x1x3_i * rho_x2x3_i );

C_2 = ( 1.0 - rho_x2x3_i^2 ) / sigma_x1_i^2;

C_3 = ( 1.0 - rho_x1x3_i^2 ) / sigma_x2_i^2;

C_4 = ( 1.0 - rho_x1x2_i^2 ) / sigma_x3_i^2;

C_5 = 2.0 * ( rho_x1x3_i*rho_x2x3_i - rho_x1x2_i ) ...
      / ( sigma_x1_i * sigma_x2_i );

C_6 = 2.0 * ( rho_x1x2_i*rho_x2x3_i - rho_x1x3_i ) ...
      / ( sigma_x1_i * sigma_x3_i );

C_7 = 2.0 * ( rho_x1x2_i*rho_x1x3_i - rho_x2x3_i ) ...
      / ( sigma_x2_i * sigma_x3_i );

% A trivariate normal PDF
P_i ...
= 1.0 / ( (2.0*pi)^(3.0/2.0) ...
          * sigma_x1_i * sigma_x2_i * sigma_x3_i * C_1 ) ...
  * exp( - ( 1.0 / ( 2.0 * C_1^2 ) ) ...
           * (   C_2 * ( x1 - mu_x1_i )^2 ...
               + C_3 * ( x2 - mu_x2_i )^2 ...
               + C_4 * ( x3 - mu_x3_i )^2 ...
               + C_5 * ( x1 - mu_x1_i ) * ( x2 - mu_x2_i ) ...
               + C_6 * ( x1 - mu_x1_i ) * ( x3 - mu_x3_i ) ...
               + C_7 * ( x2 - mu_x2_i ) * ( x3 - mu_x3_i ) ...
             ) ...
       );