% $Id$
function [ C_i ] = CDF_comp_Lognormal( x, mu_x_i_n, sigma_x_i_n )

% This function calculates cumulative density given by the Cumulative
% Density Function (CDF) for a univariate lognormal distribution.  Since
% the overall PDF consists of mulitple components, this function gives the
% results for one component, i.

% x:  a variable that is distributed lognormally in the ith component of
%     the PDF.
%
% mu_x_i_n:  the mean of ln x in the ith component of the PDF; where ln x
%            is normally distributed.  The mean of ln x is given by the
%            equation:
%              ln( mu_x_i ) - (1/2)*ln( 1 + sigma_x_i^2 / mu_x_i^2 );
%            which can also be written as:
%              ln[ mu_x_i * ( 1 + sigma_x_i^2 / mu_x_i^2 )^(-1/2) ];
%            where mu_x_i and sigma_x_i are the mean and standard deviation
%            of x in the ith component of the PDF.  The overall mean of x
%            is x_bar.
%
% sigma_x_i_n:  the standard deviation of ln x in the ith component of the
%               PDF, where ln x is normally distributed.  The variance of
%               ln x in the ith component of the PDF is sigma_x_i_n^2, and
%               is given by the equation:
%                 ln( 1 + sigma_x_i^2 / mu_x_i^2 );
%               where mu_x_i and sigma_x_i are the mean and standard
%               deviation of x in the ith component of the PDF.  Of course,
%                 sigma_x_i_n = sqrt( sigma_x_i_n^2 )
%                             = sqrt[ ln( 1 + sigma_x_i^2 / mu_x_i^2 ) ].
%               The overall variance of x is xp2.

% A single-variable lognormal CDF.
C_i = 0.5 * ( 1.0 + erf( ( log( x ) - mu_x_i_n ) ...
                         / ( sqrt(2) * sigma_x_i_n ) ) );
