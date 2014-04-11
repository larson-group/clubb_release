% $Id$
function PDF_integrate_fnc_NL( x1_opt )

% Function that compares the result of an evaluated bivariate integral to
% the result obtained by numerical integration.
%
% The integral is of the bivariate form:
% INT(-inf:inf) INT(0:inf) (x1 H(x1))^alpha x2^beta P_i(x1,x2) dx2 dx1;
% where P_i(x1,x2) is a bivariate normal-lognormal PDF.  The individual
% marginal for x1 is distributed normally, while the individual marginal for
% x2 is distributed lognormally.

addpath( '../' )

format long

% The actual means, variances, and correlation of x1 and x2.
%mu_x1 = 1.0e-3;
mu_x2 = 1.0e-4;
%sigma_x1 = 1.25e-4
sigma_x2 = 2.0e-5;
rho_x1x2 = -0.05;
if ( x1_opt == 1 )
   % The individual marginal PDF of x1 is basically entirely greater than 0.
   mu_x1 = 1.0e-3;
   sigma_x1 = 1.25e-4;
elseif ( x1_opt == 2 )
   % The individual marginal PDF of x1 is mostly greater than 0,
   % but there is a significant portion that is less than 0.
   mu_x1 = 1.0e-4;
   sigma_x1 = 1.25e-4;
elseif ( x1_opt == 3 )
   % The individual marginal PDF of x1 is slightly more heavily less
   % than 0 than it is greater than 0.
   mu_x1 = -1.0e-7;
   sigma_x1 = 1.25e-4;
else % default
   % Same as x1_opt = 1.
   mu_x1 = 1.0e-3;
   sigma_x1 = 1.25e-4;
end

% normalize lognormal means, variances, and correlation.
mu_x2_n = log( mu_x2 * ( 1.0 + sigma_x2^2 / mu_x2^2 )^-0.5 );
sigma_x2_n = sqrt( log( 1.0 + sigma_x2^2 / mu_x2^2 ) );
rho_x1x2_n = rho_x1x2 * sqrt( exp( sigma_x2_n^2 ) - 1.0 ) / sigma_x2_n;

% Numerical integration.
num_std_dev_x1 = 10.0;
num_std_dev_x2 = 20.0;

% lower bound of x1
integral_lb_x1 = mu_x1 - num_std_dev_x1 * sigma_x1;
% upper bound of x1
integral_ub_x1 = mu_x1 + num_std_dev_x1 * sigma_x1;

% lower bound of x2
integral_lb_x2 = max( 0.0, mu_x2 - num_std_dev_x2 * sigma_x2 );
% upper bound of x2
integral_ub_x2 = mu_x2 + num_std_dev_x2 * sigma_x2;

% Number of partitions
num_partitions_x1 = 2000;
num_partitions_x2 = 1000;

fprintf( 'number of partitions: x1: %d; x2: %d\n', ...
         num_partitions_x1, num_partitions_x2 );

% Partition size for x1 and x2.
delta_x1 = ( integral_ub_x1 - integral_lb_x1 ) / num_partitions_x1;
delta_x2 = ( integral_ub_x2 - integral_lb_x2 ) / num_partitions_x2;

fprintf( 'partition sizes: x1: %7.5e; x2: %7.5e\n', ...
         delta_x1, delta_x2 );

% Set the exponents on the function f(x1,x2) = (x1 H(x1))^alpha * x2^beta.
alpha_exp = 1.15;
beta_exp = 1.15;

sum = 0.0;
for iter_x1 = 1:1:num_partitions_x1
   for iter_x2 = 1:1:num_partitions_x2
   
      x1 = ( 0.5 + (iter_x1-1) ) * delta_x1 + integral_lb_x1;
      x2 = ( 0.5 + (iter_x2-1) ) * delta_x2 + integral_lb_x2;
      
      sum = sum ...
            + delta_x1 * delta_x2 ...
              * x1^alpha_exp * (Heaviside(x1))^alpha_exp * x2^beta_exp ...
              * PDF_comp_bivar_NL( x1, x2, mu_x1, mu_x2_n, ...
                                   sigma_x1, sigma_x2_n, rho_x1x2_n );
      
   end
end

numerical_integral = sum;

fprintf( 'numerical integral = %17.16g\n', numerical_integral );


% Calculate the actual integral
% INT(-inf:inf) INT(0:inf) (x1 H(x1))^alpha x2^beta P_i(x1,x2) dx2 dx1,
% which also can be written as
% INT(0:inf) INT(0:inf) x1^alpha x2^beta P_i(x1,x2) dx2 dx1.
s_c = ( mu_x1 / sigma_x1 ) + rho_x1x2_n * sigma_x2_n * beta_exp;

x1_alpha_x2_beta ...
= 1.0 / sqrt( 2.0*pi ) * sigma_x1^alpha_exp ...
  * exp( mu_x2_n * beta_exp ...
         + 0.5 * sigma_x2_n^2 * beta_exp^2 - 0.25 * s_c^2 ) ...
  * gamma( alpha_exp + 1.0 ) * pu( alpha_exp + 0.5, -s_c );

result = x1_alpha_x2_beta;

fprintf( 'Integral = %17.16g\n', result );


% Percent difference between the numerical integral and the actual integral.
percent_diff ...
   = 100.0 * abs( ( numerical_integral - result ) / result );

fprintf( [ 'Percent difference between numerical integral and ', ...
           'actual integral: %7.5g%%\n' ], percent_diff );


function [ H_value ] = Heaviside( x )

if ( x > 0 )
   H_value = 1;
else
   H_value = 0;
end
