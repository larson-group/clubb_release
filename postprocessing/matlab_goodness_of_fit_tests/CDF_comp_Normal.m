% $Id$
function [ C_i ] = CDF_comp_Normal( x, mu_x_i, sigma_x_i )

% This function calculates cumulative density given by the Cumulative
% Density Function (CDF) for a univariate normal distribution.  Since the
% overall CDF consists of mulitple components, this function gives the
% results for one component, i.

% x:  a variable that is distributed normally in the ith component of the
%     PDF.
%
% mu_x_i:  the mean of x in the ith component of the PDF.  The overall mean
%          of x is x_bar.
%
% sigma_x_i:  the standard deviation of x in the ith component of the PDF.
%             The variance of x in the ith component of the PDF is
%             sigma_x_i^2, and the overall variance of x is xp2.

% A single-variable normal CDF.
C_i = 0.5 * ( 1.0 + erf( ( x - mu_x_i ) / ( sqrt(2) * sigma_x_i ) ) );
