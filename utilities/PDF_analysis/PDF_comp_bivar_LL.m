% $Id$
function [ P_i ] ...
= PDF_comp_bivar_LL( x1, x2, mu_x1_ni, mu_x2_ni, ...
                     sigma_x1_ni, sigma_x2_ni, rho_x1x2_ni )

% This function calculates probability density given by the Probability
% Density Function (PDF) for a bivariate lognormal distribution.  Since the
% overall PDF consists of mulitple components, this function gives the
% results for one component, i.

% x1:  a variable that has a lognormal marginal distribution in the ith
%      component of the PDF.
%
% x2:  a variable that has a lognormal marginal distribution in the ith
%      component of the PDF.
%
% mu_x1_ni:  the mean of ln x1 in the ith component of the PDF; where ln x1
%            is normally distributed.  The mean of ln x1 is given by the
%            equation:
%              ln( mu_x1_i ) - (1/2)*ln( 1 + sigma_x1_i^2 / mu_x1_i^2 );
%            which can also be written as:
%              ln[ mu_x1_i * ( 1 + sigma_x1_i^2 / mu_x1_i^2 )^(-1/2) ];
%            where mu_x1_i and sigma_x1_i are the mean and standard
%            deviation of x1 in the ith component of the PDF.  The overall
%            mean of x1 is x1_bar.
%
% mu_x2_ni:  the mean of ln x2 in the ith component of the PDF; where ln x2
%            is normally distributed.  The mean of ln x2 is given by the
%            equation:
%              ln( mu_x2_i ) - (1/2)*ln( 1 + sigma_x2_i^2 / mu_x2_i^2 );
%            which can also be written as:
%              ln[ mu_x2_i * ( 1 + sigma_x2_i^2 / mu_x2_i^2 )^(-1/2) ];
%            where mu_x2_i and sigma_x2_i are the mean and standard
%            deviation of x2 in the ith component of the PDF.  The overall
%            mean of x2 is x2_bar.
%
% sigma_x1_ni:  the standard deviation of ln x1 in the ith component of the
%               PDF, where ln x1 is normally distributed.  The variance of
%               ln x1 in the ith component of the PDF is sigma_x1_ni^2, and
%               is given by the equation:
%                 ln( 1 + sigma_x1_i^2 / mu_x1_i^2 );
%               where mu_x1_i and sigma_x1_i are the mean and standard
%               deviation of x1 in the ith component of the PDF.  Of
%               course:
%                 sigma_x1_ni = sqrt( sigma_x1_ni^2 )
%                             = sqrt[ ln( 1 + sigma_x1_i^2 / mu_x1_i^2 ) ].
%               The overall variance of x1 is x1p2.
%
% sigma_x2_ni:  the standard deviation of ln x2 in the ith component of the
%               PDF, where ln x2 is normally distributed.  The variance of
%               ln x2 in the ith component of the PDF is sigma_x2_ni^2, and
%               is given by the equation:
%                 ln( 1 + sigma_x2_i^2 / mu_x2_i^2 );
%               where mu_x2_i and sigma_x2_i are the mean and standard
%               deviation of x2 in the ith component of the PDF.  Of 
%               course:
%                 sigma_x2_ni = sqrt( sigma_x2_ni^2 )
%                             = sqrt[ ln( 1 + sigma_x2_i^2 / mu_x2_i^2 ) ].
%               The overall variance of x2 is x2p2.
%
% rho_x1x2_ni:  the correlation between variables ln x1 and ln x2 in the
%               ith component of the PDF, where ln x1 and ln x2 are
%               normally distributed.  The correlation is given by the
%               equation:
%                 1.0 / ( sigma_x1_ni * sigma_x2_ni )
%                 * ln[ 1 + rho_x1x2_i 
%                           * sqrt( exp( sigma_x1_ni^2 ) - 1 )
%                           * sqrt( exp( sigma_x2_ni^2 ) - 1 ) ];
%               where rho_x1x2_i is the correlation between x1 and x2 in
%               the ith component of the PDF.  The overall correlation
%               between variables x1 and x2 can be found by:
%                 x1px2p / ( sqrt ( x1p2 * x2p2 ) );
%               where x1px2p is the overall covariance of x1 and x2.

% A bivariate lognormal PDF
P_i ...
= 1.0 / ( 2.0*pi * sigma_x1_ni * sigma_x2_ni ...
          * sqrt( 1.0 - rho_x1x2_ni^2 ) * x1 * x2 ) ...
  * exp( - ( 1.0 / ( 2.0 * ( 1.0 - rho_x1x2_ni^2 ) ) ) ...
           * (   ( 1.0 / sigma_x1_ni^2 ) ...
                  * ( log( x1 ) - mu_x1_ni )^2 ...
               - ( 2.0 * rho_x1x2_ni / ( sigma_x1_ni * sigma_x2_ni ) ) ...
                  * ( log( x1 ) - mu_x1_ni ) * ( log( x2 ) - mu_x2_ni ) ...
               + ( 1.0 / sigma_x2_ni^2 ) ...
                  * ( log( x2 ) - mu_x2_ni )^2 ...
             ) ...
       );
