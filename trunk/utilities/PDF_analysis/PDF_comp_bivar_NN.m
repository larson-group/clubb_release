% $Id$
function [ P_i ] ...
= PDF_comp_bivar_NN( x1, x2, mu_x1_i, mu_x2_i, ...
                     sigma_x1_i, sigma_x2_i, rho_x1x2_i )

% This function calculates probability density given by the Probability
% Density Function (PDF) for a bivariate normal distribution.  Since the
% overall PDF consists of mulitple components, this function gives the
% results for one component, i.

% x1:  a variable that has a normal marginal distribution in the ith
%      component of the PDF.
%
% x2:  a variable that has a normal marginal distribution in the ith
%      component of the PDF.
%
% mu_x1_i:  the mean of x1 in the ith component of the PDF.  The overall
%           mean of x1 is x1_bar.
%
% mu_x2_i:  the mean of x2 in the ith component of the PDF.  The overall
%           mean of x2 is x2_bar.
%
% sigma_x1_i:  the standard deviation of x1 in the ith component of the
%              PDF.  The variance of x1 in the ith component of the PDF is
%              sigma_x1_i^2, and the overall variance of x1 is x1p2.
%
% sigma_x2_i:  the standard deviation of x2 in the ith component of the
%              PDF.  The variance of x2 in the ith component of the PDF is
%              sigma_x2_i^2, and the overall variance of x2 is x2p2.
%
% rho_x1x2_i:  the correlation between variables x1 and x2 in the ith
%              component of the PDF.  The overall correlation between
%              variables x1 and x2 can be found by:
%                 x1px2p / ( sqrt ( x1p2 * x2p2 ) );
%              where x1px2p is the overall covariance of x1 and x2.

% A bivariate normal PDF
P_i ...
= 1.0 / ( 2.0*pi * sigma_x1_i * sigma_x2_i ...
          * sqrt( 1.0 - rho_x1x2_i^2 ) ) ...
  * exp( - ( 1.0 / ( 2.0 * ( 1.0 - rho_x1x2_i^2 ) ) ) ...
           * (   ( 1.0 / sigma_x1_i^2 ) ...
                  * ( x1 - mu_x1_i )^2 ...
               - ( 2.0 * rho_x1x2_i / ( sigma_x1_i * sigma_x2_i ) ) ...
                  * ( x1 - mu_x1_i ) * ( x2 - mu_x2_i ) ...
               + ( 1.0 / sigma_x2_i^2 ) ...
                  * ( x2 - mu_x2_i )^2 ...
             ) ...
       );
