% $Id$
function PDF_integrate_fnc_NNL_const_x1_x3( a_exp, b_exp, x2_opt )

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
% In this special case, the x1 and x3 variables are held constant.

addpath( '../' )

format long

x1_bar = -0.5;
x2_alpha_x3_beta_bar = 1.0e-8; 

% The actual means, variances, and correlations of x1, x2, and x3.
mu_x1 = -1.0;
%mu_x2 = 1.0e-3;
mu_x3 = 1.0e-4;
%sigma_x2 = 1.25e-4;
if ( x2_opt == 1 )
   % The individual marginal PDF of x2 is basically entirely greater than 0.
   mu_x2 = 1.0e-3;
   sigma_x2 = 1.25e-4;
elseif ( x2_opt == 2 )
   % The individual marginal PDF of x2 is mostly greater than 0,
   % but there is a significant portion that is less than 0.
   mu_x2 = 1.0e-4;
   sigma_x2 = 1.25e-4;
elseif ( x2_opt == 3 )
   % The individual marginal PDF of x2 is slightly more heavily less
   % than 0 than it is greater than 0.
   mu_x2 = -1.0e-7;
   sigma_x2 = 1.25e-4;
else % default
   % Same as x2_opt = 1.
   mu_x2 = 1.0e-3;
   sigma_x2 = 1.25e-4;
end

% Numerical integration.
num_std_dev_x2 = 10.0;

% lower bound for x2
integral_lb_x2 = mu_x2 - num_std_dev_x2 * sigma_x2;
% upper bound for x2
integral_ub_x2 = mu_x2 + num_std_dev_x2 * sigma_x2;

% Number of partitions
num_partitions_x2 = 400;

fprintf( 'number of partitions: x2: %d\n', num_partitions_x2 );

% Partition size for x2.
delta_x2 = ( integral_ub_x2 - integral_lb_x2 ) / num_partitions_x2;

fprintf( 'partition size: x2: %7.5e\n', delta_x2 );

% Set the exponents on the function
% f(x2,x3) = (x2 H(x2))^alpha * x2^beta.
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

sum = 0.0;
for iter_x2 = 1:1:num_partitions_x2
   
   x2 = ( 0.5 + (iter_x2-1) ) * delta_x2 + integral_lb_x2;
      
   sum = sum ...
         + delta_x2 ...
           * ( x2^alpha_exp * (Heaviside(x2))^alpha_exp ...
               * mu_x3^beta_exp - x2_alpha_x3_beta_bar )^b_exp ...
           * PDF_comp_Normal( x2, mu_x2, sigma_x2 );

end

numerical_integral = ( mu_x1 - x1_bar )^a_exp * sum;

fprintf( 'numerical integral = %17.16g\n', numerical_integral );


% Calculate the actual integral
% INT(-inf:inf) INT(0:inf) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(x2))^alpha x3^beta - x2_alpha_x3_beta_bar )^b
%    P_i(x1,x2,x3) dx3 dx2 dx1.
sigma_sum_1 = 0.0;
for q = 0:1:b_exp
         
   sigma_sum_1 ...
   = sigma_sum_1 ...
   + 1.0 / sqrt( 2.0*pi ) ...
     * ( mu_x1 - x1_bar )^a_exp ...
     * factorial( b_exp ) ...
       / ( factorial( b_exp - q ) * factorial( q ) ) ...
     * ( - x2_alpha_x3_beta_bar )^(b_exp-q) ...
     * sigma_x2^(alpha_exp*q) ...
     * mu_x3^(beta_exp*q) ...
     * exp( - 0.25 * ( mu_x2^2 / sigma_x2^2 ) ) ...
     * gamma( alpha_exp*q + 1 ) ...
     * pu( alpha_exp*q + 0.5, - ( mu_x2 / sigma_x2 ) );
          
end

fprintf( 'Integral (x2 domain from 0 to inf) = %17.16g\n', sigma_sum_1 );

% Calculate the actual integral
% INT(-inf:inf) INT(-inf:0) INT(0:inf)
%    ( x1 - x1_bar )^a
%    ( (x2 H(x2))^alpha x3^beta - x2_alpha_x3_beta_bar )^b
%    P_i(x1,x2,x3) dx3 dx2 dx1.
sigma_sum_2 ...
= ( mu_x1 - x1_bar )^a_exp ...
  * ( - x2_alpha_x3_beta_bar )^b_exp ...
  * 0.5 * ( 1.0 - erf( mu_x2 / ( sqrt(2) * sigma_x2 ) ) );

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
if ( ( a_exp == 1 ) & ( b_exp == 1 ) )
   
   x1_p1_x2_alpha_x3_beta_p1 ...
   = ( 1.0 / sqrt( 2.0*pi ) ) * ( mu_x1 - x1_bar ) ...
     * sigma_x2^alpha_exp ...
     * mu_x3^beta_exp ...
     * exp( - 0.25 * ( mu_x2^2 / sigma_x2^2 ) ) ...
     * gamma( alpha_exp + 1.0 ) ...
     * pu( alpha_exp + 0.5, - ( mu_x2 / sigma_x2 ) ) ...
     - x2_alpha_x3_beta_bar * ( mu_x1 - x1_bar );
   
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
