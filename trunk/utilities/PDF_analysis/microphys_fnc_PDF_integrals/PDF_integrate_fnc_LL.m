% $Id$
function PDF_integrate_fnc_LL

% Function that compares the result of an evaluated bivariate integral to
% the result obtained by numerical integration.
%
% The integral is of the bivariate form:
% INT(0:inf) INT(0:inf) x1^alpha x2^beta P_i(x1,x2) dx2 dx1;
% where P_i(x1,x2) is a bivariate lognormal-lognormal PDF.  The individual
% marginals for x1 and x2 are both distributed lognormally.

addpath( '../' )

format long

% The actual means, variances, and correlation of x1 and x2.
mu_x1 = 1.0e-4;
mu_x2 = 100.0;
sigma_x1 = 2.0e-5;
sigma_x2 = 20.0;
rho_x1x2 = 0.3;

% normalize lognormal means, variances, and correlation.
mu_x1_n = log( mu_x1 * ( 1.0 + sigma_x1^2 / mu_x1^2 )^-0.5 );
sigma_x1_n = sqrt( log( 1.0 + sigma_x1^2 / mu_x1^2 ) );
mu_x2_n = log( mu_x2 * ( 1.0 + sigma_x2^2 / mu_x2^2 )^-0.5 );
sigma_x2_n = sqrt( log( 1.0 + sigma_x2^2 / mu_x2^2 ) );
rho_x1x2_n = ( 1.0 / ( sigma_x1_n * sigma_x2_n ) ) ...
             * log( 1.0 + rho_x1x2 ...
                          * sqrt( exp( sigma_x1_n^2 ) - 1.0 ) ...
                          * sqrt( exp( sigma_x2_n^2 ) - 1.0 ) );

% Numerical integration.
num_std_dev_x1 = 15.0;
num_std_dev_x2 = 15.0;
                      
% lower bound of x1
integral_lb_x1 = max( 0.0, mu_x1 - num_std_dev_x1 * sigma_x1 );
% upper bound of x1
integral_ub_x1 = mu_x1 + num_std_dev_x1 * sigma_x1;

% lower bound of x2
integral_lb_x2 = max( 0.0, mu_x2 - num_std_dev_x2 * sigma_x2 );
% upper bound of x2
integral_ub_x2 = mu_x2 + num_std_dev_x2 * sigma_x2;

% Number of partitions
num_partitions_x1 = 2000;
num_partitions_x2 = 2000;

fprintf( 'number of partitions: x1: %d; x2: %d\n', ...
         num_partitions_x1, num_partitions_x2 );

% Partition size for x1 and x2.
delta_x1 = ( integral_ub_x1 - integral_lb_x1 ) / num_partitions_x1;
delta_x2 = ( integral_ub_x2 - integral_lb_x2 ) / num_partitions_x2;

fprintf( 'partition sizes: x1: %7.5e; x2: %7.5e\n', ...
         delta_x1, delta_x2 );

% Set the exponents on the function f(x1,x2) = x1^alpha * x2^beta.
alpha_exp = 1.0/3.0;
beta_exp = -1.0/3.0;

sum = 0.0;
for iter_x1 = 1:1:num_partitions_x1
   for iter_x2 = 1:1:num_partitions_x2
   
      x1 = ( 0.5 + (iter_x1-1) ) * delta_x1 + integral_lb_x1;
      x2 = ( 0.5 + (iter_x2-1) ) * delta_x2 + integral_lb_x2;
                               
      sum = sum ...
            + delta_x1 * delta_x2 ...
              * x1^alpha_exp * x2^beta_exp ...
              * PDF_comp_bivar_LL( x1, x2, mu_x1_n, mu_x2_n, ...
                                   sigma_x1_n, sigma_x2_n, rho_x1x2_n );
      
   end
end

numerical_integral = sum;

fprintf( 'numerical integral = %17.16g\n', numerical_integral );


% Calculate the actual integral
% INT(0:inf) INT(0:inf) x1^alpha x2^beta P_i(x1,x2) dx2 dx1.
x1_alpha_x2_beta ...
= exp( mu_x1_n * alpha_exp + mu_x2_n * beta_exp ...
       + 0.5 * sigma_x1_n^2 * alpha_exp^2 ...
       + 0.5 * sigma_x2_n^2 * beta_exp^2 ...
       + rho_x1x2_n * sigma_x1_n * alpha_exp * sigma_x2_n * beta_exp );

result = x1_alpha_x2_beta;

fprintf( 'Integral = %17.16g\n', result );


% Percent difference between the numerical integral and the actual integral.
percent_diff ...
   = 100.0 * abs( ( numerical_integral - result ) / result );

fprintf( [ 'Percent difference between numerical integral and ', ...
           'actual integral: %7.5g%%\n' ], percent_diff );
