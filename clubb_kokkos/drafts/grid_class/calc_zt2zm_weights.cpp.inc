// =============================================================================
static void calc_zt2zm_weights(int nzm, int nzt, int ngrdcol, // In
                               grid& gr) {                    // InOut
  // Description:
  // Function used to help in an interpolation of a variable (var_zt) located
  // on the thermodynamic grid levels (azt) to the momentum grid levels (azm).
  // This function computes a weighting factor for both the upper thermodynamic
  // level (k) and the lower thermodynamic level (k-1) applied to the central
  // momentum level (k).  For the boundary momentum grid level k=1, a weighting
  // factor for both the thermodynamic level at k=1 and the thermodynamic level
  // at k=2 are calculated based on the use of a linear extension. Likewise,
  // for the boundary momentum grid level k=nzm, a weighting factor for both
  // the thermodynamic level at k=nzt and the thermodynamic level at k=nzt-1
  // are calculated based on the use of a linear extension.  This function
  // outputs the weighting factors at a single momentum grid level (k).  The
  // formulation used is compatible with a stretched (unevenly-spaced) grid.
  // The weights are defined as follows:
  //
  // ---var_zt(k)--------------------------------------------- t(k)
  //                       azt_weight(t_above) = factor
  // ===========var_zt(interp)================================ m(k)
  //                       azt_weight(t_below) = 1 - factor
  // ---var_zt(k-1)------------------------------------------- t(k-1)
  //
  // The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
  // zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // NOTE: For an ascending grid, the grid diagram shown above has increasing
  // indices with increasing altitude, such that zt(k) > zt(k-1). For a
  // descending grid, the grid diagram shown above has increasing indices with
  // decreasing altitude, such that zt(k) < zt(k-1).
  //
  // For all momentum levels 1 < k < nzm:
  //
  // The formula for a linear interpolation is given by:
  //
  // var_zt( interp to zm(k) )
  // = [ ( var_zt(k) - var_zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
  //     * ( zm(k) - zt(k-1) ) + var_zt(k-1);
  //
  // which can be rewritten as:
  //
  // var_zt( interp to zm(k) )
  // = [ ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
  //     * ( var_zt(k) - var_zt(k-1) ) + var_zt(k-1).
  //
  // Furthermore, the formula can be rewritten as:
  //
  // var_zt( interp to zm(k) )
  // = factor * var_zt(k) + ( 1 - factor ) * var_zt(k-1);
  //
  // where:
  //
  // factor = ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ).
  //
  // One of the important uses of this function is in situations where the
  // variable to be interpolated is being treated IMPLICITLY in an equation.
  // Usually, the variable to be interpolated is involved in a derivative (such
  // as d(var_zt)/dz in the diagram below).  For the term of the equation
  // containing the derivative, grid weights are needed for two interpolations,
  // rather than just one interpolation.  Thus, four grid weights (labeled
  // A(k), B(k), C(k), and D(k) in the diagram below) are needed.
  //
  // ---var_zt(k+1)------------------------------------------- t(k+1)
  //                                       A(k)
  // ===========var_zt(interp)================================ m(k+1)
  //                                       B(k) = 1 - A(k)
  // ---var_zt(k)-----------d(var_zt)/dz---------------------- t(k)
  //                                       C(k)
  // ===========var_zt(interp)================================ m(k)
  //                                       D(k) = 1 - C(k)
  // ---var_zt(k-1)------------------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // The grid weights, indexed around the central thermodynamic level (k), are
  // defined as follows:
  //
  // A(k) = ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) );
  //
  // which is the same as "factor" for the interpolation to momentum
  // level (k+1).  In the code, this interpolation is referenced as
  // gr%weights_zt2zm(t_above,mkp1), which can be read as "grid weight in a
  // zt2zm interpolation of the thermodynamic level above momentum level (k+1)
  // (applied to momentum level (k+1))".
  //
  // B(k) = 1 - [ ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
  //      = 1 - A(k);
  //
  // which is the same as "1 - factor" for the interpolation to momentum
  // level (k+1).  This equation can also be rewritten as:
  //
  // B(k) = ( zt(k+1) - zm(k+1) ) / ( zt(k+1) - zt(k) )
  //
  // In the code, this interpolation is referenced as
  // gr%weights_zt2zm(t_below,mkp1), which can be read as "grid weight in a
  // zt2zm interpolation of the thermodynamic level below momentum level (k+1)
  // (applied to momentum level (k+1))".
  //
  // C(k) = ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) );
  //
  // which is the same as "factor" for the interpolation to momentum
  // level (k).  In the code, this interpolation is referenced as
  // gr%weights_zt2zm(t_above,mk), which can be read as "grid weight in a
  // zt2zm interpolation of the thermodynamic level above momentum level (k)
  // (applied to momentum level (k))".
  //
  // D(k) = 1 - [ ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
  //      = 1 - C(k);
  //
  // which is the same as "1 - factor" for the interpolation to momentum
  // level (k).  This equation can also be rewritten as:
  //
  // D(k) = ( zt(k) - zm(k) ) / ( zt(k) - zt(k-1) )
  //
  // In the code, this interpolation is referenced as
  // gr%weights_zt2zm(t_below,mk), which can be read as "grid weight in a
  // zt2zm interpolation of the thermodynamic level below momentum level (k)
  // (applied to momentum level (k))".
  //
  // Additionally, as long as the central thermodynamic level (k) in the above
  // scenario is not a boundary thermodynamic level (k /= nzt and k /= 1), the
  // four weighting factors have the following relationships:  A(k) = C(k+1)
  // and B(k) = D(k+1).
  //
  //
  // Special condition for grid level, k = 1:
  //
  // The boundary momentum grid level k = 1 is outside the thermodynamic grid
  // level k = 1.  Thus, a linear extension is used at this level.
  //
  // The formula for a linear extension to momentum level k = 1 is given by:
  //
  // var_zt( extend to zm(k) )
  // = [ ( var_zt(k+1) - var_zt(k) ) / ( zt(k+1) - zt(k) ) ]
  //     * ( zm(k) - zt(k) ) + var_zt(k);
  //
  // which can be rewritten as:
  //
  // var_zt( extend to zm(k) )
  // = [ ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
  //     * ( var_zt(k+1) - var_zt(k) ) + var_zt(k).
  //
  // Furthermore, the formula can be rewritten as:
  //
  // var_zt( extend to zm(k) )
  // = factor * var_zt(k+1) + ( 1 - factor ) * var_zt(k);
  //
  // where:
  //
  // factor = ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ).
  //
  // Finally, since here k = 1:
  //
  // factor = ( zm(1) - zt(1) ) / ( zt(2) - zt(1) ).
  //
  // Due to the fact that a linear extension is being used, the value of factor
  // will be less than 0.  The weight of thermodynamic level k = 2 on momentum
  // level k = 1 equals the value of factor.  The weight of thermodynamic level
  // k = 1 on momentum level k = 1 equals 1 - factor, which is greater than 1.
  // However, the sum of the weights equals 1.
  //
  //
  // Special condition for grid level, k = nzm:
  //
  // The boundary momentum grid level k = nzm is outside the thermodynamic grid
  // level k = nzt.  Thus, a linear extension is used at this level.
  //
  // The formula for a linear extension to momentum level k = nzm
  // is given by:
  //
  // var_zt( extend to zm(k) )
  // = [ ( var_zt(k-1) - var_zt(k-2) ) / ( zt(k-1) - zt(k-2) ) ]
  //     * ( zm(k) - zt(k-2) ) + var_zt(k-2);
  //
  // which can be rewritten as:
  //
  // var_zt( extend to zm(k) )
  // = [ ( zm(k) - zt(k-2) ) / ( zt(k-1) - zt(k-2) ) ]
  //     * ( var_zt(k-1) - var_zt(k-2) ) + var_zt(k-2).
  //
  // Furthermore, the formula can be rewritten as:
  //
  // var_zt( extend to zm(k) )
  // = factor * var_zt(k-1) + ( 1 - factor ) * var_zt(k-2);
  //
  // where:
  //
  // factor = ( zm(k) - zt(k-2) ) / ( zt(k-1) - zt(k-2) ).
  //
  // Finally, since here k = nzm and nzt = nzm - 1:
  //
  // factor = ( zm(nzm) - zt(nzt-1) ) / ( zt(nzt) - zt(nzt-1) ).
  //
  // Due to the fact that a linear extension is being used, the value of factor
  // will be greater than 1.  The weight of thermodynamic level k = nzt on
  // momentum level k = nzm equals the value of factor.  The weight of
  // thermodynamic level k = nzt-1 on momentum level k = nzm equals
  // 1 - factor, which is less than 0.  However, the sum of the two weights
  // equals 1.
  //
  //
  // Brian Griffin; September 12, 2008.
  // Updated for removal of "ghost" thermodynamic level that was below the
  // surface; July 17, 2023.
  // Updated to generalize the grid for both ascending and descending grid
  // directions; September 24, 2024.
  //
  //-----------------------------------------------------------------------
  // nzm, nzt: Momentum and thermodynamic level counts [#].
  // ngrdcol: Number of columns [#].
  // gr: Heights [m] and interpolation weights [-] in the existing grid.
  const auto zm = gr.zm;
  const auto zt = gr.zt;
  const auto weights_zt2zm = gr.weights_zt2zm;
  // Interior interpolation and the two independent linear extensions share
  // one launch. Native requires at least two thermodynamic levels.
  Kokkos::parallel_for(
      "calc_zt2zm_weights", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        if (k > 0 && k < nzm - 1) {
          // Weight of upper thermodynamic level on momentum level.
          weights_zt2zm(i, k, t_above) = (zm(i, k) - zt(i, k - 1)) / (zt(i, k) - zt(i, k - 1));
          // Weight of lower thermodynamic level on momentum level.
          weights_zt2zm(i, k, t_below) = (zt(i, k) - zm(i, k)) / (zt(i, k) - zt(i, k - 1));
        } else if (k == 0) {
          // Storage level 1 is outside thermodynamic level 1: extend from
          // thermodynamic levels 1 and 2. This is the lower boundary ascending
          // and the upper boundary descending; factor is less than zero.
          weights_zt2zm(i, 0, t_above) = (zm(i, 0) - zt(i, 0)) / (zt(i, 1) - zt(i, 0));
          weights_zt2zm(i, 0, t_below) = (zt(i, 1) - zm(i, 0)) / (zt(i, 1) - zt(i, 0));
        } else {
          // Storage level nzm is outside thermodynamic level nzt: extend from
          // levels nzt and nzt-1. This is the upper boundary ascending and the
          // lower boundary descending; factor is greater than one.
          weights_zt2zm(i, nzm - 1, t_above) =
              (zm(i, nzm - 1) - zt(i, nzt - 2)) / (zt(i, nzt - 1) - zt(i, nzt - 2));
          weights_zt2zm(i, nzm - 1, t_below) =
              (zt(i, nzt - 1) - zm(i, nzm - 1)) / (zt(i, nzt - 1) - zt(i, nzt - 2));
        }
      });
}
