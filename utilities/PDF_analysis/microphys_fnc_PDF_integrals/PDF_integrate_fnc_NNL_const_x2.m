% $Id$
function PDF_integrate_fnc_NNL_const_x2( a_exp, b_exp, x2_opt )

% Function that compares the result of an evaluated trivariate integral to
% the result obtained by numerical integration.
%
% The integral is of the trivariate form:
% INT(-inf:inf) INT(-inf:inf) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(x2))^alpha x3^beta - x2_alpha_x3_beta_bar )^b
%    P_i(x1,x2,x3) dx3 dx2 dx1;
% where P_i(x1,x2,x3) is a trivariate normal-normal-lognormal PDF.  The
% individual marginals for x1 and x2 are both distributed normally,
% while the individual marginal for x3 is distributed lognormally.
%
% In this special case, the x2 variable is held constant.

addpath( '../' )

format long

x1_bar = -0.5;
x2_alpha_x3_beta_bar = 1.0e-8; 

% The actual means, variances, and correlations of x1, x2, and x3.
mu_x1 = -1.0;
%mu_x2 = 1.0e-3;
mu_x3 = 1.0e-4;
sigma_x1 = 0.125;
sigma_x3 = 2.0e-5;
rho_x1x3 = 0.10;
if ( x2_opt == 1 )
   % Since sigma_x2 = 0, x2 is greater than 0.
   mu_x2 = 1.0e-3;
elseif ( x2_opt == 2 )
   % Since sigma_x2 = 0, x2 is greater than 0.
   mu_x2 = 1.0e-4;
elseif ( x2_opt == 3 )
   % Since sigma_x2 = 0, x2 is less than 0.
   mu_x2 = -1.0e-7;
else % default
   % Same as x2_opt = 1.
   mu_x2 = 1.0e-3;
end

% normalize lognormal means, variances, and correlation.
mu_x3_n = log( mu_x3 * ( 1.0 + sigma_x3^2 / mu_x3^2 )^-0.5 );
sigma_x3_n = sqrt( log( 1.0 + sigma_x3^2 / mu_x3^2 ) );
rho_x1x3_n = rho_x1x3 * sqrt( exp( sigma_x3_n^2 ) - 1.0 ) / sigma_x3_n;

% Numerical integration.
num_std_dev_x1 = 20.0;
num_std_dev_x3 = 20.0;

% lower bound for x1
integral_lb_x1 = mu_x1 - num_std_dev_x1 * sigma_x1;
% upper bound for x1
integral_ub_x1 = mu_x1 + num_std_dev_x1 * sigma_x1;

% lower bound for x3
integral_lb_x3 = max( 0.0, mu_x3 - num_std_dev_x3 * sigma_x3 );
% upper bound for x3
integral_ub_x3 = mu_x3 + num_std_dev_x3 * sigma_x3;

% Number of partitions
num_partitions_x1 = 100;
num_partitions_x3 = 400;

fprintf( 'number of partitions: x1: %d; x3: %d\n', ...
         num_partitions_x1, num_partitions_x3 );

% Partition size for x1 and x3.
delta_x1 = ( integral_ub_x1 - integral_lb_x1 ) / num_partitions_x1;
delta_x3 = ( integral_ub_x3 - integral_lb_x3 ) / num_partitions_x3;

fprintf( 'partition sizes: x1: %7.5e; x3: %7.5e\n', delta_x1, delta_x3 );

% Set the exponents on the function f(x1,x2) = x1^alpha * x2^beta.
alpha_exp = 1.15;
beta_exp = 1.15;

% Set the order prime exponents.
if ( a_exp < 0 )
   a_exp = 0;
else
   a_exp = floor( a_exp );
end
if ( b_exp < 0 )
   b_exp = 0;
else
   b_exp = floor( b_exp );
end

fprintf( 'a_exp = %d; b_exp = %d; x2 option: %d\n', a_exp, b_exp, x2_opt );

tot_int = 0.0;
for q = 0:1:b_exp

   prefix = factorial( b_exp ) ...
            / ( factorial( b_exp - q ) * factorial( q ) ) ...
            * ( -x2_alpha_x3_beta_bar )^(b_exp-q) ...
            * mu_x2^(alpha_exp*q) * (Heaviside(mu_x2))^(alpha_exp*q);
    
   sum = 0.0;
   for iter_x1 = 1:1:num_partitions_x1
      for iter_x3 = 1:1:num_partitions_x3
   
         x1 = ( 0.5 + (iter_x1-1) ) * delta_x1 + integral_lb_x1;
         x3 = ( 0.5 + (iter_x3-1) ) * delta_x3 + integral_lb_x3;

         sum = sum ...
               + delta_x1 * delta_x3 ...
                 * ( x1 - x1_bar )^a_exp ...
                 * x3^(beta_exp*q) ...
                 * PDF_comp_bivar_NL ...
                      ( x1, x3, mu_x1, mu_x3_n, ...
                        sigma_x1, sigma_x3_n, rho_x1x3_n );

      end
   end

   tot_int = tot_int + prefix * sum;
   
end

numerical_integral = tot_int;

fprintf( 'numerical integral = %17.16g\n', numerical_integral );


% Calculate the actual integral
% INT(-inf:inf) INT(0:inf) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(x2))^alpha x3^beta - x2_alpha_x3_beta_bar )^b
%    P_i(x1,x2,x3) dx3 dx2 dx1.
if ( mu_x2 >= 0.0 )

   sigma_sum_1 = 0.0;
   for p = 0:1:floor(a_exp/2.0)
      for q = 0:1:b_exp

         sigma_sum_1 ...
         = sigma_sum_1 ...
         +   factorial( a_exp ) ...
             / ( factorial( a_exp - 2*p ) * factorial( p ) ) ...
           * factorial( b_exp ) ...
             / ( factorial( b_exp - q ) * factorial( q ) ) ...
           * ( - x2_alpha_x3_beta_bar )^(b_exp-q) ...
           * mu_x2^(alpha_exp*q) ...
           * ( 0.5 * sigma_x1^2 )^p ...
           * ( mu_x1 - x1_bar ...
               + rho_x1x3_n * sigma_x1 ...
                 * sigma_x3_n * beta_exp * q )^(a_exp-2*p) ...
           * exp( mu_x3_n * beta_exp * q ...
                  + 0.5 * sigma_x3_n^2 * beta_exp^2 * q^2 );

      end
   end

else % mu_x2 < 0

   sigma_sum_1 = 0.0;

end

fprintf( 'Integral (x2 domain from 0 to inf) = %17.16g\n', sigma_sum_1 );

% Calculate the actual integral
% INT(-inf:inf) INT(-inf:0) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(x2))^alpha x3^beta - x2_alpha_x3_beta_bar )^b
%    P_i(x1,x2,x3) dx3 dx2 dx1.
if ( mu_x2 < 0 )

   sigma_sum_2 = 0.0;
   for p = 0:1:floor(a_exp/2.0)

      sigma_sum_2 ...
      = sigma_sum_2 ...
      + ( - x2_alpha_x3_beta_bar )^b_exp ...
        * factorial( a_exp ) ...
          / ( factorial( a_exp - 2*p ) * factorial( p ) ) ...
        * ( 0.5 * sigma_x1^2 )^p ...
        * ( mu_x1 - x1_bar )^(a_exp-2*p);

   end

else % mu_x2 >= 0

   sigma_sum_2 = 0.0;

end

fprintf( 'Integral (x2 domain from -inf to 0) = %17.16g\n', sigma_sum_2 );

% Add the two to obtain the result of the integral over the full domain:
% INT(-inf:inf) INT(-inf:inf) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(x2))^alpha x3^beta - x2_alpha_x3_beta_bar )^b
%    P_i(x1,x2,x3) dx3 dx2 dx1.
result_x1_pa_x2_alpha_x3_beta_pb = sigma_sum_1 + sigma_sum_2;

fprintf( 'Integral = %17.16g\n', result_x1_pa_x2_alpha_x3_beta_pb );


% Calculate the actual integral
% INT(-inf:inf) INT(-inf:inf) INT(0:inf)
%    ( x1 - x1_bar )
%    ( (x2 H(x2))^alpha x3^beta - x2_alpha_x3_beta_bar )
%    P_i(x1,x2,x3) dx3 dx2 dx1.
if ( ( a_exp == 1 ) && ( b_exp == 1 ) )

   if ( mu_x2 >= 0.0 )
    
      x1_p1_x2_alpha_x3_beta_p1 ...
      = mu_x2^alpha_exp ...
        * ( mu_x1 - x1_bar ...
            + rho_x1x3_n * sigma_x1 * sigma_x3_n * beta_exp ) ...
        * exp( mu_x3_n * beta_exp ...
                 + 0.5 * sigma_x3_n^2 * beta_exp^2 ) ...
        - x2_alpha_x3_beta_bar * ( mu_x1 - x1_bar );

   else % mu_x2 < 0

      x1_p1_x2_alpha_x3_beta_p1 ...
      = - x2_alpha_x3_beta_bar * ( mu_x1 - x1_bar );

   end
 
   result_x1_p1_x2_alpha_x3_beta_p1 = x1_p1_x2_alpha_x3_beta_p1;

   fprintf( 'Integral (when a = 1 and b = 1) = %17.16g\n', ...
            result_x1_p1_x2_alpha_x3_beta_p1 );
    
end


% Percent difference between the numerical integral and the actual integral.
percent_diff ...
   = 100.0 * abs( ( numerical_integral ...
                    - result_x1_pa_x2_alpha_x3_beta_pb ) ...
                  / result_x1_pa_x2_alpha_x3_beta_pb );

fprintf( [ 'Percent difference between numerical integral and ', ...
           'actual integral: %7.5g%%\n' ], percent_diff );


function [ H_value ] = Heaviside( x )

if ( x > 0 )
   H_value = 1;
else
   H_value = 0;
end
