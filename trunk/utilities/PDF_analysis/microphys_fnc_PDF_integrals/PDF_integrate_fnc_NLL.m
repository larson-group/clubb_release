% $Id$
function PDF_integrate_fnc_NLL( x1_opt )

% Function that compares the result of an evaluated trivariate integral to
% the result obtained by numerical integration.
%
% The integral is of the trivariate form:
% INT(-inf:inf) INT(0:inf) INT(0:inf)
%    (x1 H(-x1))^alpha x2^beta x3^gamma P_i(x1,x2,x3) dx3 dx2 dx1;
% where P_i(x1,x2,x3) is a trivariate normal-lognormal-lognormal PDF.  The
% individual marginal for x1 is distributed normally, while the individual
% marginals for x2 and x3 are both distributed lognormally.

addpath( '../' )

format long

% The actual means, variances, and correlations of x1, x2, and x3.
%mu_x1 = -1.0e-3;
mu_x2 = 1.0e-4;
mu_x3 = 100.0;
%sigma_x1 = 1.25e-4;
sigma_x2 = 2.0e-5;
sigma_x3 = 20.0;
rho_x1x2 = 0.2;
rho_x1x3 = 0.25;
rho_x2x3 = 0.3;
if ( x1_opt == 1 )
   % The individual marginal PDF of x1 is basically entirely less than 0.
   mu_x1 = -1.0e-3;
   sigma_x1 = 1.25e-4;
elseif ( x1_opt == 2 )
   % The individual marginal PDF of x1 is mostly less than 0,
   % but there is a significant portion that is greater than 0.
   mu_x1 = -1.0e-4;
   sigma_x1 = 1.25e-4;
elseif ( x1_opt == 3 )
   % The individual marginal PDF of x1 is slightly more heavily greater
   % than 0 than it is less than 0.
   mu_x1 = 1.0e-7;
   sigma_x1 = 1.25e-4;
else % default
   % Same as x1_opt = 1.
   mu_x1 = -1.0e-3;
   sigma_x1 = 1.25e-4;
end

% normalize lognormal means, variances, and correlation.
mu_x2_n = log( mu_x2 * ( 1.0 + sigma_x2^2 / mu_x2^2 )^-0.5 );
sigma_x2_n = sqrt( log( 1.0 + sigma_x2^2 / mu_x2^2 ) );
mu_x3_n = log( mu_x3 * ( 1.0 + sigma_x3^2 / mu_x3^2 )^-0.5 );
sigma_x3_n = sqrt( log( 1.0 + sigma_x3^2 / mu_x3^2 ) );
rho_x1x2_n = rho_x1x2 * sqrt( exp( sigma_x2_n^2 ) - 1.0 ) / sigma_x2_n;
rho_x1x3_n = rho_x1x3 * sqrt( exp( sigma_x3_n^2 ) - 1.0 ) / sigma_x3_n;
rho_x2x3_n = ( 1.0 / ( sigma_x2_n * sigma_x3_n ) ) ...
             * log( 1.0 + rho_x2x3 ...
                          * sqrt( exp( sigma_x2_n^2 ) - 1.0 ) ...
                          * sqrt( exp( sigma_x3_n^2 ) - 1.0 ) );

% Numerical integration.
num_std_dev_x1 = 10.0;
num_std_dev_x2 = 15.0;
num_std_dev_x3 = 15.0;

% lower bound for x1
integral_lb_x1 = mu_x1 - num_std_dev_x1 * sigma_x1;
% upper bound for x1
integral_ub_x1 = mu_x1 + num_std_dev_x1 * sigma_x1;

% lower bound for x2
integral_lb_x2 = max( 0.0, mu_x2 - num_std_dev_x2 * sigma_x2 );
% upper bound for x2
integral_ub_x2 = mu_x2 + num_std_dev_x2 * sigma_x2;

% lower bound for x3
integral_lb_x3 = max( 0.0, mu_x3 - num_std_dev_x3 * sigma_x3 );
% upper bound for x3
integral_ub_x3 = mu_x3 + num_std_dev_x3 * sigma_x3;

% Number of partitions
num_partitions_x1 = 800;
num_partitions_x2 = 400;
num_partitions_x3 = 400;

fprintf( 'number of partitions: x1: %d; x2: %d; x3: %d\n', ...
         num_partitions_x1, num_partitions_x2, num_partitions_x3 );

% Partition size for x1, x2, and x3.
delta_x1 = ( integral_ub_x1 - integral_lb_x1 ) / num_partitions_x1;
delta_x2 = ( integral_ub_x2 - integral_lb_x2 ) / num_partitions_x2;
delta_x3 = ( integral_ub_x3 - integral_lb_x3 ) / num_partitions_x3;

fprintf( 'partition sizes: x1: %7.5e; x2: %7.5e; x3: %7.5e\n', ...
         delta_x1, delta_x2, delta_x3 );

% Set the exponents on the function
% f(x1,x2,x3) = (x1 H(-x1))^alpha * x2^beta * x3^gamma.
alpha_exp = 1.0;
beta_exp  = 1.0/3.0;
gamma_exp = 2.0/3.0;

sum = 0.0;
for iter_x1 = 1:1:num_partitions_x1
   for iter_x2 = 1:1:num_partitions_x2
      for iter_x3 = 1:1:num_partitions_x3
   
         x1 = ( 0.5 + (iter_x1-1) ) * delta_x1 + integral_lb_x1;
         x2 = ( 0.5 + (iter_x2-1) ) * delta_x2 + integral_lb_x2;
         x3 = ( 0.5 + (iter_x3-1) ) * delta_x3 + integral_lb_x3;

         sum = sum ...
               + delta_x1 * delta_x2 * delta_x3 ...
                 * x1^alpha_exp * (Heaviside(-x1))^alpha_exp ...
                   * x2^beta_exp * x3^gamma_exp ...
                 * PDF_comp_trivar_NLL ...
                      ( x1, x2, x3, mu_x1, mu_x2_n, mu_x3_n, ...
                        sigma_x1, sigma_x2_n, sigma_x3_n, ...
                        rho_x1x2_n, rho_x1x3_n, rho_x2x3_n );

      end
   end
end

numerical_integral = sum;

fprintf( 'numerical integral = %17.16g\n', numerical_integral );


% Calculate the actual integral
% INT(-inf:inf) INT(0:inf) INT(0:inf)
%    (x1 H(-x1))^alpha x2^beta x3^gamma P_i(x1,x2,x3) dx3 dx2 dx1,
% which also can be written as
% INT(-inf:0) INT(0:inf) INT(0:inf)
%    x1^alpha x2^beta x3^gamma P_i(x1,x2,x3) dx3 dx2 dx1.
s_cc = ( mu_x1 / sigma_x1 ) ...
       + rho_x1x2_n * sigma_x2_n * beta_exp ...
       + rho_x1x3_n * sigma_x3_n * gamma_exp;

x1_alpha_x2_beta_x3_gamma ...
= 1.0 / sqrt( 2.0*pi ) * ( -sigma_x1 )^alpha_exp ...
  * exp( mu_x2_n * beta_exp + mu_x3_n * gamma_exp ) ...
  * exp( 0.5 * ...
         (   ( 1.0 - rho_x1x2_n^2 ) * sigma_x2_n^2 * beta_exp^2  ...
           + ( 1.0 - rho_x1x3_n^2 ) * sigma_x3_n^2 * gamma_exp^2 ...
           + 2.0 * ( rho_x2x3_n - rho_x1x2_n * rho_x1x3_n ) ...
                 * sigma_x2_n * beta_exp * sigma_x3_n * gamma_exp ...
         ) ...
       ) ...
  * exp( 0.25 * s_cc^2 - ( mu_x1 / sigma_x1 ) * s_cc ...
         + 0.5 * ( mu_x1^2 / sigma_x1^2 ) ) ...
  * gamma( alpha_exp + 1.0 ) * pu( alpha_exp + 0.5, s_cc );

result = x1_alpha_x2_beta_x3_gamma;

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
