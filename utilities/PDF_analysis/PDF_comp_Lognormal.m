function [ P_i ] = PDF_comp_Lognormal( x, mu_x_ni, sigma_x_ni )

% This function calculates probability density given by the Probability
% Density Function (PDF) for a univariate lognormal distribution.  Since
% the overall PDF consists of mulitple components, this function gives the
% results for one component, i.

% x:  a variable that is distributed lognormally in the ith component of
%     the PDF.
%
% mu_x_ni:  the mean of ln x in the ith component of the PDF; where ln x is
%           normally distributed.  The mean of ln x is given by the
%           equation:
%             ln( mu_x_i ) - (1/2)*ln( 1 + sigma_x_i^2 / mu_x_i^2 );
%           which can also be written as:
%             ln[ mu_x_i * ( 1 + sigma_x_i^2 / mu_x_i^2 )^(-1/2) ];
%           where mu_x_i and sigma_x_i are the mean and standard deviation
%           of x in the ith component of the PDF.  The overall mean of x is
%           x_bar.
%
% sigma_x_ni:  the standard deviation of ln x in the ith component of the
%              PDF, where ln x is normally distributed.  The variance of
%              ln x in the ith component of the PDF is sigma_x_ni^2, and is
%              given by the equation:
%                ln( 1 + sigma_x_i^2 / mu_x_i^2 );
%              where mu_x_i and sigma_x_i are the mean and standard
%              deviation of x in the ith component of the PDF.  Of course,
%                sigma_x_ni = sqrt( sigma_x_ni^2 )
%                           = sqrt[ ln( 1 + sigma_x_i^2 / mu_x_i^2 ) ].
%              The overall variance of x is xp2.

% A single-variable normal PDF.
P_i ...
= 1.0 / ( sqrt( 2.0*pi ) * sigma_x_ni * x ) ...
  * exp( - ( log( x ) - mu_x_ni )^2 / ( 2.0 * sigma_x_ni^2 ) );