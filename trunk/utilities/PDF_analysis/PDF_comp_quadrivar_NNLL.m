% $Id$
function [ P_i ] ...
= PDF_comp_quadrivar_NNLL( x1, x2, x3, x4, mu_x1_i, mu_x2_i, mu_x3_ni, ...
                           mu_x4_ni, sigma_x1_i, sigma_x2_i, ...
                           sigma_x3_ni, sigma_x4_ni, rho_x1x2_i, ...
                           rho_x1x3_ni, rho_x1x4_ni, rho_x2x3_ni, ...
                           rho_x2x4_ni, rho_x3x4_ni )

% This function calculates probability density given by the Probability
% Density Function (PDF) for a quadrivariate
% normal-normal-lognormal-lognormal distribution.  Since the overall PDF
% consists of mulitple components, this function gives the results for one
% component, i.

% x1:  a variable that has a normal marginal distribution in the ith
%      component of the PDF.
%
% x2:  a variable that has a normal marginal distribution in the ith
%      component of the PDF.
%
% x3:  a variable that has a lognormal marginal distribution in the ith
%      component of the PDF.
%
% x4:  a variable that has a lognormal marginal distribution in the ith
%      component of the PDF.
%
% mu_x1_i:  the mean of x1 in the ith component of the PDF.  The overall
%           mean of x1 is x1_bar.
%
% mu_x2_i:  the mean of x2 in the ith component of the PDF.  The overall
%           mean of x2 is x2_bar.
%
% mu_x3_ni:  the mean of ln x3 in the ith component of the PDF; where ln x3
%            is normally distributed.  The mean of ln x3 is given by the
%            equation:
%              ln( mu_x3_i ) - (1/2)*ln( 1 + sigma_x3_i^2 / mu_x3_i^2 );
%            which can also be written as:
%              ln[ mu_x3_i * ( 1 + sigma_x3_i^2 / mu_x3_i^2 )^(-1/2) ];
%            where mu_x3_i and sigma_x3_i are the mean and standard
%            deviation of x3 in the ith component of the PDF.  The overall
%            mean of x3 is x3_bar.
%
% mu_x4_ni:  the mean of ln x4 in the ith component of the PDF; where ln x4
%            is normally distributed.  The mean of ln x4 is given by the
%            equation:
%              ln( mu_x4_i ) - (1/2)*ln( 1 + sigma_x4_i^2 / mu_x4_i^2 );
%            which can also be written as:
%              ln[ mu_x4_i * ( 1 + sigma_x4_i^2 / mu_x4_i^2 )^(-1/2) ];
%            where mu_x4_i and sigma_x4_i are the mean and standard
%            deviation of x4 in the ith component of the PDF.  The overall
%            mean of x4 is x4_bar.
%
% sigma_x1_i:  the standard deviation of x1 in the ith component of the
%              PDF.  The variance of x1 in the ith component of the PDF is
%              sigma_x1_i^2, and the overall variance of x1 is x1p2.
%
% sigma_x2_i:  the standard deviation of x2 in the ith component of the
%              PDF.  The variance of x2 in the ith component of the PDF is
%              sigma_x2_i^2, and the overall variance of x2 is x2p2.
%
% sigma_x3_ni:  the standard deviation of ln x3 in the ith component of the
%               PDF, where ln x3 is normally distributed.  The variance of
%               ln x3 in the ith component of the PDF is sigma_x3_ni^2, and
%               is given by the equation:
%                 ln( 1 + sigma_x3_i^2 / mu_x3_i^2 );
%               where mu_x3_i and sigma_x3_i are the mean and standard
%               deviation of x3 in the ith component of the PDF.  Of 
%               course:
%                 sigma_x3_ni = sqrt( sigma_x3_ni^2 )
%                             = sqrt[ ln( 1 + sigma_x3_i^2 / mu_x3_i^2 ) ].
%               The overall variance of x3 is x3p2.
%
% sigma_x4_ni:  the standard deviation of ln x4 in the ith component of the
%               PDF, where ln x4 is normally distributed.  The variance of
%               ln x4 in the ith component of the PDF is sigma_x4_ni^2, and
%               is given by the equation:
%                 ln( 1 + sigma_x4_i^2 / mu_x4_i^2 );
%               where mu_x4_i and sigma_x4_i are the mean and standard
%               deviation of x4 in the ith component of the PDF.  Of 
%               course:
%                 sigma_x4_ni = sqrt( sigma_x4_ni^2 )
%                             = sqrt[ ln( 1 + sigma_x4_i^2 / mu_x4_i^2 ) ].
%               The overall variance of x4 is x4p2.
%
% rho_x1x2_i:  the correlation between variables x1 and x2 in the ith
%              component of the PDF.  The overall correlation between
%              variables x1 and x2 can be found by:
%                 x1px2p / ( sqrt ( x1p2 * x2p2 ) );
%              where x1px2p is the overall covariance of x1 and x2.
%
% rho_x1x3_ni:  the correlation between variables x1 and ln x3 in the ith
%               component of the PDF; where ln x3 is normally distributed
%               (and of course, so is x1).  The correlation is given by the
%               equation:
%                 rho_x1x3_i
%                 * sqrt( exp( sigma_x3_ni^2 ) - 1 ) / sigma_x3_ni;
%               where rho_x1x3_i is the correlation between x1 and x3 in
%               the ith component of the PDF.  The overall correlation
%               between variables x1 and x3 can be found by:
%                 x1px3p / ( sqrt ( x1p2 * x3p2 ) );
%               where x1px3p is the overall covariance of x1 and x3.
%
% rho_x1x4_ni:  the correlation between variables x1 and ln x4 in the ith
%               component of the PDF; where ln x4 is normally distributed
%               (and of course, so is x1).  The correlation is given by the
%               equation:
%                 rho_x1x4_i
%                 * sqrt( exp( sigma_x4_ni^2 ) - 1 ) / sigma_x4_ni;
%               where rho_x1x4_i is the correlation between x1 and x4 in
%               the ith component of the PDF.  The overall correlation
%               between variables x1 and x4 can be found by:
%                 x1px4p / ( sqrt ( x1p2 * x4p2 ) );
%               where x1px4p is the overall covariance of x1 and x4.
%
% rho_x2x3_ni:  the correlation between variables x2 and ln x3 in the ith
%               component of the PDF; where ln x3 is normally distributed
%               (and of course, so is x2).  The correlation is given by the
%               equation:
%                 rho_x2x3_i
%                 * sqrt( exp( sigma_x3_ni^2 ) - 1 ) / sigma_x3_ni;
%               where rho_x2x3_i is the correlation between x2 and x3 in
%               the ith component of the PDF.  The overall correlation
%               between variables x2 and x3 can be found by:
%                 x2px3p / ( sqrt ( x2p2 * x3p2 ) );
%               where x2px3p is the overall covariance of x2 and x3.
%
% rho_x2x4_ni:  the correlation between variables x2 and ln x4 in the ith
%               component of the PDF; where ln x4 is normally distributed
%               (and of course, so is x2).  The correlation is given by the
%               equation:
%                 rho_x2x4_i
%                 * sqrt( exp( sigma_x4_ni^2 ) - 1 ) / sigma_x4_ni;
%               where rho_x2x4_i is the correlation between x2 and x4 in
%               the ith component of the PDF.  The overall correlation
%               between variables x2 and x4 can be found by:
%                 x2px4p / ( sqrt ( x2p2 * x4p2 ) );
%               where x2px4p is the overall covariance of x2 and x4.
%
% rho_x3x4_ni:  the correlation between variables ln x3 and ln x4 in the
%               ith component of the PDF, where ln x3 and ln x4 are
%               normally distributed.  The correlation is given by the
%               equation:
%                 1.0 / ( sigma_x3_ni * sigma_x4_ni )
%                 * ln[ 1 + rho_x3x4_i 
%                           * sqrt( exp( sigma_x3_ni^2 ) - 1 )
%                           * sqrt( exp( sigma_x4_ni^2 ) - 1 ) ];
%               where rho_x3x4_i is the correlation between x3 and x4 in
%               the ith component of the PDF.  The overall correlation
%               between variables x3 and x4 can be found by:
%                 x3px4p / ( sqrt ( x3p2 * x4p2 ) );
%               where x3px4p is the overall covariance of x3 and x4.

C_1 ...
= sqrt( 1.0 ...
        - ( rho_x1x2_i^2 + rho_x1x3_ni^2 + rho_x1x4_ni^2 ...
            + rho_x2x3_ni^2 + rho_x2x4_ni^2 + rho_x3x4_ni^2 ) ...
        + 2.0 * rho_x1x2_i * rho_x1x3_ni * rho_x2x3_ni ...
        + 2.0 * rho_x1x2_i * rho_x1x4_ni * rho_x2x4_ni ...
        + 2.0 * rho_x1x3_ni * rho_x1x4_ni * rho_x3x4_ni ...
        + 2.0 * rho_x2x3_ni * rho_x2x4_ni * rho_x3x4_ni ...
        + rho_x1x2_i^2 * rho_x3x4_ni^2 ...
        + rho_x1x3_ni^2 * rho_x2x4_ni^2 ...
        + rho_x1x4_ni^2 * rho_x2x3_ni^2 ...
        - 2.0 * rho_x1x2_i * rho_x1x3_ni * rho_x2x4_ni * rho_x3x4_ni ...
        - 2.0 * rho_x1x2_i * rho_x1x4_ni * rho_x2x3_ni * rho_x3x4_ni ...
        - 2.0 * rho_x1x3_ni * rho_x1x4_ni * rho_x2x3_ni * rho_x2x4_ni ...
      );

% Coefficient for ( x1 - mu_x1_i )^2.
C_2 = ( 1.0 - ( rho_x2x3_ni^2 + rho_x2x4_ni^2 + rho_x3x4_ni^2 ) ...
        + 2.0 * rho_x2x3_ni * rho_x2x4_ni * rho_x3x4_ni ) / sigma_x1_i^2;

% Coefficient for ( x2 - mu_x2_i )^2.
C_3 = ( 1.0 - ( rho_x1x3_ni^2 + rho_x1x4_ni^2 + rho_x3x4_ni^2 ) ...
        + 2.0 * rho_x1x3_ni * rho_x1x4_ni * rho_x3x4_ni ) / sigma_x2_i^2;
    
% Coefficient for ( ln x3 - mu_x3_ni )^2.
C_4 = ( 1.0 - ( rho_x1x2_i^2 + rho_x1x4_ni^2 + rho_x2x4_ni^2 ) ...
        + 2.0 * rho_x1x2_i * rho_x1x4_ni * rho_x2x4_ni ) / sigma_x3_ni^2;
    
% Coefficient for ( ln x4 - mu_x4_ni )^2.
C_5 = ( 1.0 - ( rho_x1x2_i^2 + rho_x1x3_ni^2 + rho_x2x3_ni^2 ) ...
        + 2.0 * rho_x1x2_i * rho_x1x3_ni * rho_x2x3_ni ) / sigma_x4_ni^2;
    
% Coefficient for ( x1 - mu_x1_i )( x2 - mu_x2_i ).
C_6 = 2.0 * ( rho_x1x2_i * rho_x3x4_ni^2 ...
              - rho_x1x4_ni * rho_x2x3_ni * rho_x3x4_ni ...
              - rho_x1x3_ni * rho_x2x4_ni * rho_x3x4_ni ...
              + rho_x1x3_ni * rho_x2x3_ni + rho_x1x4_ni * rho_x2x4_ni ...
              - rho_x1x2_i ) / ( sigma_x1_i * sigma_x2_i );

% Coefficient for ( x1 - mu_x1_i )( ln x3 - mu_x3_ni ).
C_7 = 2.0 * ( rho_x1x3_ni * rho_x2x4_ni^2 ...
              - rho_x1x2_i * rho_x2x4_ni * rho_x3x4_ni ...
              - rho_x1x4_ni * rho_x2x3_ni * rho_x2x4_ni ...
              + rho_x1x2_i * rho_x2x3_ni + rho_x1x4_ni * rho_x3x4_ni ...
              - rho_x1x3_ni ) / ( sigma_x1_i * sigma_x3_ni );

% Coefficient for ( x1 - mu_x1_i )( ln x4 - mu_x4_ni ).
C_8 = 2.0 * ( rho_x1x4_ni * rho_x2x3_ni^2 ...
              - rho_x1x2_i * rho_x2x3_ni * rho_x3x4_ni ...
              - rho_x1x3_ni * rho_x2x3_ni * rho_x2x4_ni ...
              + rho_x1x2_i * rho_x2x4_ni + rho_x1x3_ni * rho_x3x4_ni ...
              - rho_x1x4_ni ) / ( sigma_x1_i * sigma_x4_ni );
    
% Coefficient for ( x2 - mu_x2_i )( ln x3 - mu_x3_ni ).
C_9 = 2.0 * ( rho_x1x4_ni^2 * rho_x2x3_ni ...
              - rho_x1x2_i * rho_x1x4_ni * rho_x3x4_ni ...
              - rho_x1x3_ni * rho_x1x4_ni * rho_x2x4_ni ...
              + rho_x2x4_ni * rho_x3x4_ni + rho_x1x2_i * rho_x1x3_ni ...
              - rho_x2x3_ni ) / ( sigma_x2_i * sigma_x3_ni );

% Coefficient for ( x2 - mu_x2_i )( ln x4 - mu_x4_ni ).
C_10 = 2.0 * ( rho_x1x3_ni^2 * rho_x2x4_ni ...
               - rho_x1x2_i * rho_x1x3_ni * rho_x3x4_ni ...
               - rho_x1x3_ni * rho_x1x4_ni * rho_x2x3_ni ...
               + rho_x2x3_ni * rho_x3x4_ni + rho_x1x2_i * rho_x1x4_ni ...
               - rho_x2x4_ni ) / ( sigma_x2_i * sigma_x4_ni );

% Coefficient for ( ln x3 - mu_x3_ni )( ln x4 - mu_x4_ni ).
C_11 = 2.0 * ( rho_x1x2_i^2 * rho_x3x4_ni ...
               - rho_x1x2_i * rho_x1x4_ni * rho_x2x3_ni ...
               - rho_x1x2_i * rho_x1x3_ni * rho_x2x4_ni ...
               + rho_x2x3_ni * rho_x2x4_ni + rho_x1x3_ni * rho_x1x4_ni ...
               - rho_x3x4_ni ) / ( sigma_x3_ni * sigma_x4_ni );

% A quadrivariate normal PDF
P_i ...
= 1.0 / ( (2.0*pi)^2 * sigma_x1_i * sigma_x2_i ...
          * sigma_x3_ni * sigma_x4_ni * C_1 * x3 * x4 ) ...
  * exp( - ( 1.0 / ( 2.0 * C_1^2 ) ) ...
           * (   C_2 * ( x1 - mu_x1_i )^2 ...
               + C_3 * ( x2 - mu_x2_i )^2 ...
               + C_4 * ( log( x3 ) - mu_x3_ni )^2 ...
               + C_5 * ( log( x4 ) - mu_x4_ni )^2 ...
               + C_6 * ( x1 - mu_x1_i ) * ( x2 - mu_x2_i ) ...
               + C_7 * ( x1 - mu_x1_i ) * ( log( x3 ) - mu_x3_ni ) ...
               + C_8 * ( x1 - mu_x1_i ) * ( log( x4 ) - mu_x4_ni ) ...
               + C_9 * ( x2 - mu_x2_i ) * ( log( x3 ) - mu_x3_ni ) ...
               + C_10 * ( x2 - mu_x2_i ) * ( log( x4 ) - mu_x4_ni ) ...
               + C_11 * ( log( x3 ) - mu_x3_ni ) ...
                      * ( log( x4 ) - mu_x4_ni ) ...
             ) ...
       );
