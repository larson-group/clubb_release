function [ P_i ] = PDF_comp_Normal( x, mu_x_i, sigma_x_i )

% This function calculates probability density given by the Probability
% Density Function (PDF) for a univariate normal distribution.  Since the
% overall PDF consists of mulitple components, this function gives the
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

% A single-variable normal PDF.
P_i ...
= 1.0 / ( sqrt( 2.0*pi ) * sigma_x_i ) ...
  * exp( - ( x - mu_x_i )^2 / ( 2.0 * sigma_x_i^2 ) );