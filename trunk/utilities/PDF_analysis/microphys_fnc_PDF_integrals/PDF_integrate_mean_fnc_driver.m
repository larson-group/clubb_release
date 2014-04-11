% $Id$
function PDF_integrate_mean_fnc_driver

for x1_opt = 1:1:3

   % Bivariate NL Integrals.

   fprintf( 'PDF_integrate_fnc_NL\n' )
   fprintf( '\n' )
   PDF_integrate_fnc_NL( x1_opt )
   fprintf( '\n' )
   fprintf( '\n' )
   fprintf( 'PDF_integrate_fnc_NL_const_x1\n' )
   fprintf( '\n' )
   PDF_integrate_fnc_NL_const_x1( x1_opt )
   fprintf( '\n' )
   fprintf( '\n' )

   % Trivariate NLL Integrals.

   fprintf( 'PDF_integrate_fnc_NLL\n' )
   fprintf( '\n' )
   PDF_integrate_fnc_NLL( x1_opt )
   fprintf( '\n' )
   fprintf( '\n' )
   fprintf( 'PDF_integrate_fnc_NLL_const_x1\n' )
   fprintf( '\n' )
   PDF_integrate_fnc_NLL_const_x1( x1_opt )
   fprintf( '\n' )
   fprintf( '\n' )

   fprintf( '==============================================================\n' )
   fprintf( '\n' )

end

fprintf( 'PDF_integrate_fnc_LL\n' )
fprintf( '\n' )
PDF_integrate_fnc_LL
fprintf( '\n' )
fprintf( '\n' )
