% $Id$
function PDF_integrate_fnc_NNLL_const_x2_x3_x4( a_exp, b_exp, x2_opt )

% Function that compares the result of an evaluated quadrivariate integral to
% the result obtained by numerical integration.
%
% The integral is of the quadrivariate form:
% INT(-inf:inf) INT(-inf:inf) INT(0:inf) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(-x2))^alpha x3^beta x4^gamma - x2_alpha_x3_beta_x4_gamma_bar )^b
%    P_i(x1,x2,x3,x4) dx4 dx3 dx2 dx1;
% where P_i(x1,x2,x3,x4) is a quadrivariate normal-normal-lognormal-lognormal
% PDF.  The individual marginals for x1 and x2 are both distributed normally,
% while the individual marginals for x3 and x4 are both distributed lognormally.
%
% In this special case, the x2, x3, and x4 variables are held constant.

addpath( '../' )

format long

x1_bar = -0.5;
x2_alpha_x3_beta_x4_gamma_bar = -2.0e-4;

% The actual means, variances, and correlations of x1, x2, x3, and x4.
mu_x1 = -1.0;
%mu_x2 = -1.0e-3;
mu_x3 = 1.0e-4;
mu_x4 = 100.0;
sigma_x1 = 0.125;
if ( x2_opt == 1 )
   % Since sigma_x2 = 0, x2 is less than 0.
   mu_x2 = -1.0e-3;
elseif ( x2_opt == 2 )
   % Since sigma_x2 = 0, x2 is less than 0.
   mu_x2 = -1.0e-4;
elseif ( x2_opt == 3 )
   % Since sigma_x2 = 0, x2 is greater than 0.
   mu_x2 = 1.0e-7;
else % default
   % Same as x2_opt = 1.
   mu_x2 = -1.0e-3;
end

% Numerical integration.
num_std_dev_x1 = 10.0;

% lower bound for x1
integral_lb_x1 = mu_x1 - num_std_dev_x1 * sigma_x1;
% upper bound for x1
integral_ub_x1 = mu_x1 + num_std_dev_x1 * sigma_x1;

% Number of partitions
num_partitions_x1 = 10000;

fprintf( 'number of partitions: x1: %d\n', num_partitions_x1 );

% Partition size for x1.
delta_x1 = ( integral_ub_x1 - integral_lb_x1 ) / num_partitions_x1;

fprintf( 'partition size: x1: %7.5e\n', delta_x1 );

% Set the exponents on the function
% f(x2,x3,x4) = (x2 H(-x2))^alpha * x3^beta * x4^gamma.
alpha_exp = 1.0;
beta_exp  = 1.0/3.0;
gamma_exp = 2.0/3.0;

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

prefix = ( mu_x2^alpha_exp * (Heaviside(-mu_x2))^alpha_exp * mu_x3^beta_exp ...
           * mu_x4^gamma_exp - x2_alpha_x3_beta_x4_gamma_bar )^b_exp;

sum = 0.0;
for iter_x1 = 1:1:num_partitions_x1
   
   x1 = ( 0.5 + (iter_x1-1) ) * delta_x1 + integral_lb_x1;

   sum = sum ...
         + delta_x1 * ( x1 - x1_bar )^a_exp ...
           * PDF_comp_Normal( x1, mu_x1, sigma_x1 );

end
   
tot_int = tot_int + prefix * sum;

numerical_integral = tot_int;

fprintf( 'numerical integral = %17.16g\n', numerical_integral );


% Calculate the actual integral
% INT(-inf:inf) INT(-inf:0) INT(0:inf) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(-x2))^alpha x3^beta x4^gamma - x2_alpha_x3_beta_x4_gamma_bar )^b
%    P_i(x1,x2,x3,x4) dx4 dx3 dx2 dx1.
if ( mu_x2 <= 0.0 )

   sigma_sum_1 = 0.0;
   for p = 0:1:floor(a_exp/2.0)

      sigma_sum_1 ...
      = sigma_sum_1 ...   
      +   factorial( a_exp ) ...
          / ( factorial( a_exp - 2*p ) * factorial( p ) ) ...
        * ( 0.5 * sigma_x1^2 )^p * ( mu_x1 - x1_bar )^(a_exp-2*p) ...
        * ( mu_x2^alpha_exp * mu_x3^beta_exp * mu_x4^gamma_exp ...
            - x2_alpha_x3_beta_x4_gamma_bar )^b_exp;

   end

else % mu_x2 > 0

   sigma_sum_1 = 0.0;

end

fprintf( 'Integral (x2 domain from -inf to 0) = %17.16g\n', sigma_sum_1 );

% Calculate the actual integral
% INT(-inf:inf) INT(0:inf) INT(0:inf) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(-x2))^alpha x3^beta x4^gamma - x2_alpha_x3_beta_x4_gamma_bar )^b
%    P_i(x1,x2,x3,x4) dx4 dx3 dx2 dx1.
if ( mu_x2 > 0.0 )

   sigma_sum_2 = 0.0;
   for p = 0:1:floor(a_exp/2.0)

      sigma_sum_2 ...
      = sigma_sum_2 ...
      + ( - x2_alpha_x3_beta_x4_gamma_bar )^b_exp ...
        * factorial( a_exp ) ...
          / ( factorial( a_exp - 2*p ) * factorial( p ) ) ...
        * ( 0.5 * sigma_x1^2 )^p ...
        * ( mu_x1 - x1_bar )^(a_exp-2*p);
          
   end

else % mu_x2 <= 0

   sigma_sum_2 = 0.0;

end

fprintf( 'Integral (x2 domain from 0 to inf) = %17.16g\n', sigma_sum_2 );

% Add the two to obtain the result of the integral over the full domain:
% INT(-inf:inf) INT(-inf:inf) INT(0:inf) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(-x2))^alpha x3^beta x4^gamma - x2_alpha_x3_beta_x4_gamma_bar )^b
%    P_i(x1,x2,x3,x4) dx4 dx3 dx2 dx1.
result_x1_pa_x2_alpha_x3_beta_x4_gamma_pb = sigma_sum_1 + sigma_sum_2;

fprintf( 'Integral = %17.16g\n', result_x1_pa_x2_alpha_x3_beta_x4_gamma_pb );


% Calculate the actual integral
% INT(-inf:inf) INT(-inf:inf) INT(0:inf) INT(0:inf)
%    ( x1 - x1_bar )
%    ( (x2 H(-x2))^alpha x3^beta x4^gamma - x2_alpha_x3_beta_x4_gamma_bar )
%    P_i(x1,x2,x3,x4) dx4 dx3 dx2 dx1.
if ( ( a_exp == 1 ) && ( b_exp == 1 ) )

   if ( mu_x2 <= 0.0 )

      x1_p1_x2_alpha_x3_beta_x4_gamma_p1 ...
      = ( mu_x1 - x1_bar ) ...         
        * ( mu_x2^alpha_exp * mu_x3^beta_exp * mu_x4^gamma_exp ...
            - x2_alpha_x3_beta_x4_gamma_bar );

   else % mu_x2 > 0

      x1_p1_x2_alpha_x3_beta_x4_gamma_p1 ...
      = - x2_alpha_x3_beta_x4_gamma_bar * ( mu_x1 - x1_bar );

   end
 
   result_x1_p1_x2_alpha_x3_beta_x4_gamma_p1 ...
      = x1_p1_x2_alpha_x3_beta_x4_gamma_p1;

   fprintf( 'Integral (when a = 1 and b = 1) = %17.16g\n', ...
            result_x1_p1_x2_alpha_x3_beta_x4_gamma_p1 );
    
end


% Percent difference between the numerical integral and the actual integral.
percent_diff ...
   = 100.0 * abs( ( numerical_integral ...
                    - result_x1_pa_x2_alpha_x3_beta_x4_gamma_pb ) ...
                  / result_x1_pa_x2_alpha_x3_beta_x4_gamma_pb );

fprintf( [ 'Percent difference between numerical integral and ', ...
           'actual integral: %7.5g%%\n' ], percent_diff );


function [ H_value ] = Heaviside( x )

if ( x > 0 )
   H_value = 1;
else
   H_value = 0;
end
