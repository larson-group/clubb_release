% $Id$
function [ var_array_out ] ...
= unit_corrector( num_vars, var_array_in, units_corrector_type, idx_rho )

var_array_out = zeros( size( var_array_in ) );

for i = 1:1:num_vars
   if ( units_corrector_type(i) == 0 )
      % Keep it as it is.
      var_array_out(i,:,:,:,:) = var_array_in(i,:,:,:,:);
   elseif ( units_corrector_type(i) == 1 )
      % Divide by rho
      % (used for concentrations -- from units of m^-3 to units of kg^-1).
      var_array_out(i,:,:,:,:) ...
      = var_array_in(i,:,:,:,:) ./ var_array_in(idx_rho,:,:,:,:);
   elseif ( units_corrector_type(i) == 2 )
      % Divide by rho^2 (used for variances of concentrations -- from units
      % of m^-6 to units of kg^-2).
      var_array_out(i,:,:,:,:) ...
      = var_array_in(i,:,:,:,:) ./ var_array_in(idx_rho,:,:,:,:).^2;
   elseif ( units_corrector_type(i) == 3 )
      % Divide by 1000 (used for mixing ratios -- from units of g/kg to
      % units of kg/kg).
      var_array_out(i,:,:,:,:) = var_array_in(i,:,:,:,:) ./ 1000.0;
   elseif ( units_corrector_type(i) == 4 )
      % Multiply by 100 (used for pressure -- from units of mb to units of
      % Pa).
      var_array_out(i,:,:,:,:) = var_array_in(i,:,:,:,:) .* 100.0;
   end
end
