% $Id$
function PDF_integrate_fnc_NNL( a_exp, b_exp, c_exp )
% Function that compares the result of an evaluated trivariate integral to
% the result obtained by numerical integration.
%
% The integral is of the trivariate form:
% INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
%    ( x1 - x1_bar )^a ( x2 - x2_bar )^b ( x3 - x3_bar )^c
%    P_i(x1,x2,x3) dx3 dx2 dx1;
% where P_i(x1,x2,x3) is a trivariate normal PDF.  The individual marginals fo
% x1, x2, and x3 are both distributed normally.

format long

x1_bar = -0.5;
x2_bar = 1.0e-3;
x3_bar = 299.0;

% The actual means, variances, and correlations of x1, x2, and x3.
mu_x1 = 1.5;
mu_x2 = 1.5e-3;
mu_x3 = 298.0;
sigma_x1 = 0.5;
sigma_x2 = 5.0e-4;
sigma_x3 = 0.4;
rho_x1x2 = 0.25;
rho_x1x3 = 0.10;
rho_x2x3 = -0.5;

% Numerical integration.
num_std_dev_x1 = 10.0;
num_std_dev_x2 = 10.0;
num_std_dev_x3 = 10.0;

% lower bound for x1
integral_lb_x1 = mu_x1 - num_std_dev_x1 * sigma_x1;
% upper bound for x1
integral_ub_x1 = mu_x1 + num_std_dev_x1 * sigma_x1;

% lower bound for x2
integral_lb_x2 = mu_x2 - num_std_dev_x2 * sigma_x2;
% upper bound for x2
integral_ub_x2 = mu_x2 + num_std_dev_x2 * sigma_x2;

% lower bound for x3
integral_lb_x3 = mu_x3 - num_std_dev_x3 * sigma_x3;
% upper bound for x3
integral_ub_x3 = mu_x3 + num_std_dev_x3 * sigma_x3;

% Number of partitions
num_partitions_x1 = 200;
num_partitions_x2 = 200;
num_partitions_x3 = 200;

fprintf( 'number of partitions: x1: %d; x2: %d; x3: %d\n', ...
         num_partitions_x1, num_partitions_x2, num_partitions_x3 );

% Partition size for x1, x2, and x3.
delta_x1 = ( integral_ub_x1 - integral_lb_x1 ) / num_partitions_x1;
delta_x2 = ( integral_ub_x2 - integral_lb_x2 ) / num_partitions_x2;
delta_x3 = ( integral_ub_x3 - integral_lb_x3 ) / num_partitions_x3;

fprintf( 'partition sizes: x1: %7.5e; x2: %7.5e; x3: %7.5e\n', ...
         delta_x1, delta_x2, delta_x3 );

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
if ( c_exp < 0 )
   c_exp = 0;
else
   c_exp = floor( c_exp );
end

fprintf( 'a_exp = %d; b_exp = %d; c_exp = %d\n', a_exp, b_exp, c_exp );

sum = 0.0;
for iter_x1 = 1:1:num_partitions_x1
   for iter_x2 = 1:1:num_partitions_x2
      for iter_x3 = 1:1:num_partitions_x3
   
         x1 = ( 0.5 + ( iter_x1 - 1 ) ) * delta_x1 + integral_lb_x1;
         x2 = ( 0.5 + ( iter_x2 - 1 ) ) * delta_x2 + integral_lb_x2;
         x3 = ( 0.5 + ( iter_x3 - 1 ) ) * delta_x3 + integral_lb_x3;

         sum = sum ...
               + delta_x1 * delta_x2 * delta_x3 ...
                 * ( x1 - x1_bar )^a_exp * ( x2 - x2_bar )^b_exp ...
                 * ( x3 - x3_bar )^c_exp ...
                 * PDF_comp_trivar_NNN( x1, x2, x3, mu_x1, mu_x2, mu_x3, ...
                                        sigma_x1, sigma_x2, sigma_x3, ...
                                        rho_x1x2, rho_x1x3, rho_x2x3 );

      end
   end
end

numerical_integral = sum;

fprintf( 'numerical integral = %17.16g\n', numerical_integral );


% Calculate the actual integral
sigma_sum = 0.0;
for g = 0:1:floor(c_exp/2.0)
   for j = 0:1:(c_exp-2*g)
      for q = 0:1:floor((b_exp+c_exp-2*g-j)/2.0)
         for r = 0:1:(b_exp+c_exp-2*g-j-2*q)
            for h = 0:1:j
               for p = 0:1:floor((a_exp+r+j-h)/2.0)

                  sigma_sum ...
                  = sigma_sum ...
                    + factorial(c_exp) ...
                      / ( factorial(c_exp-2*g) * factorial(g) ) ...
                      * factorial(c_exp-2*g) ...
                        / ( factorial(c_exp-2*g-j) * factorial(j) ) ...
                      * factorial(b_exp+c_exp-2*g-j) ...
                        / ( factorial(b_exp+c_exp-2*g-j-2*q) ...
                            * factorial(q) ) ...
                      * factorial(b_exp+c_exp-2*g-j-2*q) ...
                        / ( factorial(b_exp+c_exp-2*g-j-2*q-r) ...
                            * factorial(r) ) ...
                      * factorial(j) / ( factorial(j-h) * factorial(h) ) ...
                      * factorial(a_exp+r+j-h) ...
                        / ( factorial(a_exp+r+j-h-2*p) * factorial(p) ) ...
                      * ( sigma_x3^2 / ( 2.0 * ( 1.0 - rho_x1x2^2 ) ) ...
                          * ( 1.0 - ( rho_x1x2^2 + rho_x1x3^2 + rho_x2x3^2 ) ...
                              + 2.0 * rho_x1x2 * rho_x1x3 * rho_x2x3 ) )^g ...
                      * ( -sigma_x3 * ( rho_x1x2 * rho_x1x3 - rho_x2x3 ) ...
                           / ( sigma_x2 ...
                               * ( 1.0 - rho_x1x2^2 ) ) )^(c_exp-2*g-j) ...
                      * ( rho_x1x2 * sigma_x2 / sigma_x1 )^r ...
                      * ( -sigma_x3 * ( rho_x1x2 * rho_x2x3 - rho_x1x3 ) ...
                           / ( sigma_x1 * ( 1.0 - rho_x1x2^2 ) ) )^(j-h) ...
                      * ( ( 1.0 - rho_x1x2^2 ) * sigma_x2^2 / 2.0 )^q ...
                      * ( mu_x2 - x2_bar ...
                          - rho_x1x2 * sigma_x2 / sigma_x1 ...
                            * ( mu_x1 - x1_bar ) )^(b_exp+c_exp-2*g-j-2*q-r) ...
                      * ( sigma_x1^2 / 2.0 )^p ...
                      * ( mu_x3 - x3_bar ...
                          + sigma_x3 * ( rho_x1x2 * rho_x2x3 - rho_x1x3 ) ...
                            / ( sigma_x1 * ( 1.0 - rho_x1x2^2 ) ) ...
                            * ( mu_x1 - x1_bar ) ...
                          + sigma_x3 * ( rho_x1x2 * rho_x1x3 - rho_x2x3 ) ...
                            / ( sigma_x2 * ( 1.0 - rho_x1x2^2 ) ) ...
                            * ( mu_x2 - x2_bar ) )^h ...
                      * ( mu_x1 - x1_bar )^(a_exp+r+j-h-2*p);

               end % p = 0:1:floor((a_exp+r+j-h)/2.0)
            end % h = 0:1:j
         end % r = 0:1:(b_exp+c_exp-2*g-j-2*q)
      end % q = 0:1:floor((b_exp+c_exp-2*g-j)/2.0)
   end % j = 0:1:(c_exp-2*g)
end % g = 0:1:floor(c_exp/2.0)

result_x1_pa_x2_pb_x3_pc = sigma_sum;

fprintf( 'Integral = %17.16g\n', result_x1_pa_x2_pb_x3_pc );


% Calculate the actual integral when a_exp = 1, b_exp = 1, and c_exp = 1.
if ( ( a_exp == 1 ) && ( b_exp == 1 ) && ( c_exp == 1 ) )   
      
   x1_p1_x2_p1_x3_p1 ...
   = ( mu_x1 - x1_bar ) * ( mu_x2 - x2_bar ) * ( mu_x3 - x3_bar ) ...
     + rho_x2x3 * sigma_x2 * sigma_x3 * ( mu_x1 - x1_bar ) ...
     + rho_x1x3 * sigma_x1 * sigma_x3 * ( mu_x2 - x2_bar ) ...
     + rho_x1x2 * sigma_x1 * sigma_x2 * ( mu_x3 - x3_bar );
   
   result_x1_p1_x2_p1_x3_p1 = x1_p1_x2_p1_x3_p1;

   fprintf( 'Integral (when a = 1, b = 1, and c = 1) = %17.16g\n', ...
            result_x1_p1_x2_p1_x3_p1 );
   
end


% Percent difference between the numerical integral and the actual integral.
percent_diff ...
   = 100.0 * abs( ( numerical_integral - result_x1_pa_x2_pb_x3_pc ) ...
                  / result_x1_pa_x2_pb_x3_pc );

fprintf( [ 'Percent difference between numerical integral and ', ...
           'actual integral: %7.5g%%\n' ], percent_diff );
