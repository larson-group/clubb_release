// =============================================================================
static void calc_zm2zt_weights(int nzt, int ngrdcol, // In
                               grid& gr) {           // InOut
  // Description:
  // Function used to help in an interpolation of a variable (var_zm) located
  // on the momentum grid levels (azm) to the thermodynamic grid levels (azt).
  // This function computes a weighting factor for both the upper momentum
  // level (k+1) and the lower momentum level (k) applied to the central
  // thermodynamic level (k).  This function outputs the weighting factors at
  // a single thermodynamic grid level (k).   The formulation used is
  // compatible with a stretched (unevenly-spaced) grid.  The weights are
  // defined as follows:
  //
  // ===var_zm(k+1)=========================================== m(k+1)
  //                       azm_weight(m_above) = factor
  // -----------var_zm(interp)-------------------------------- t(k)
  //                       azm_weight(m_below) = 1 - factor
  // ===var_zm(k)============================================= m(k)
  //
  // The vertical indices m(k+1), t(k), and m(k) correspond with altitudes
  // zm(k+1), zt(k), and zm(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // NOTE: For an ascending grid, the grid diagram shown above has increasing
  // indices with increasing altitude, such that zm(k+1) > zm(k). For a
  // descending grid, the grid diagram shown above has increasing indices with
  // decreasing altitude, such that zm(k+1) < zm(k).
  //
  // The formula for a linear interpolation is given by:
  //
  // var_zm( interp to zt(k) )
  // = [ ( var_zm(k+1) - var_zm(k) ) / ( zm(k+1) - zm(k) ) ]
  //     * ( zt(k) - zm(k) ) + var_zm(k);
  //
  // which can be rewritten as:
  //
  // var_zm( interp to zt(k) )
  // = [ ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) ) ]
  //     * ( var_zm(k+1) - var_zm(k) ) + var_zm(k).
  //
  // Furthermore, the formula can be rewritten as:
  //
  // var_zm( interp to zt(k) )
  // = factor * var_zm(k+1) + ( 1 - factor ) * var_zm(k);
  //
  // where:
  //
  // factor = ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) ).
  //
  // One of the important uses of this function is in situations where the
  // variable to be interpolated is being treated IMPLICITLY in an equation.
  // Usually, the variable to be interpolated is involved in a derivative (such
  // as d(var_zm)/dz in the diagram below).  For the term of the equation
  // containing the derivative, grid weights are needed for two interpolations,
  // rather than just one interpolation.  Thus, four grid weights (labeled
  // A(k), B(k), C(k), and D(k) in the diagram below) are needed.
  //
  // ===var_zm(k+1)=========================================== m(k+1)
  //                                       A(k)
  // -----------var_zm(interp)-------------------------------- t(k)
  //                                       B(k) = 1 - A(k)
  // ===var_zm(k)===========d(var_zm)/dz====================== m(k)
  //                                       C(k)
  // -----------var_zm(interp)-------------------------------- t(k-1)
  //                                       D(k) = 1 - C(k)
  // ===var_zm(k-1)=========================================== m(k-1)
  //
  // The vertical indices m(k+1), t(k), m(k), t(k-1), and m(k-1) correspond
  // with altitudes zm(k+1), zt(k), zm(k), zt(k-1), and zm(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // The grid weights, indexed around the central momentum level (k), are
  // defined as follows:
  //
  // A(k) = ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) );
  //
  // which is the same as "factor" for the interpolation to thermodynamic
  // level (k).  In the code, this interpolation is referenced as
  // gr%weights_zm2zt(m_above,tk), which can be read as "grid weight in a
  // zm2zt interpolation of the momentum level above thermodynamic
  // level (k) (applied to thermodynamic level (k))".
  //
  // B(k) = 1 - [ ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) ) ]
  //      = 1 - A(k);
  //
  // which is the same as "1 - factor" for the interpolation to thermodynamic
  // level (k).  This equation can also be rewritten as:
  //
  // B(k) = ( zm(k+1) - zt(k) ) / ( zm(k+1) - zm(k) ).
  //
  // In the code, this interpolation is referenced as
  // gr%weights_zm2zt(m_below,tk), which can be read as "grid weight in a
  // zm2zt interpolation of the momentum level below thermodynamic
  // level (k) (applied to thermodynamic level (k))".
  //
  // C(k) = ( zt(k-1) - zm(k-1) ) / ( zm(k) - zm(k-1) );
  //
  // which is the same as "factor" for the interpolation to thermodynamic
  // level (k-1).  In the code, this interpolation is referenced as
  // gr%weights_zm2zt(m_above,tkm1), which can be read as "grid weight in a
  // zm2zt interpolation of the momentum level above thermodynamic level (k-1)
  // (applied to thermodynamic level (k-1))".
  //
  // D(k) = 1 - [ ( zt(k-1) - zm(k-1) ) / ( zm(k) - zm(k-1) ) ]
  //      = 1 - C(k);
  //
  // which is the same as "1 - factor" for the interpolation to thermodynamic
  // level (k-1).  The equation can also be rewritten as:
  //
  // D(k) = ( zm(k) - zt(k-1) ) / ( zm(k) - zm(k-1) )
  //
  // In the code, this interpolation is referenced as
  // gr%weights_zm2zt(m_below,tkm1), which can be read as "grid weight in a
  // zm2zt interpolation of the momentum level below thermodynamic level (k-1)
  // (applied to thermodynamic level (k-1))".
  //
  // Additionally, as long as the central momentum level (k) in the above
  // scenario is not the lowermost momentum level or the uppermost momentum
  // level (k /= 1 and k /= gr%nz), the four weighting factors have the
  // following relationships:  A(k) = C(k+1) and B(k) = D(k+1).
  //
  //
  // Brian Griffin; September 12, 2008.
  // Updated for removal of "ghost" thermodynamic level that was below the
  // surface; July 18, 2023.
  // Updated to generalize the grid for both ascending and descending grid
  // directions; September 24, 2024.
  //
  //-----------------------------------------------------------------------
  // nzt: Number of thermodynamic levels [#].
  // ngrdcol: Number of columns [#].
  // gr: Heights [m] and interpolation weights [-] in the existing grid.
  const auto zm = gr.zm;
  const auto zt = gr.zt;
  const auto weights_zm2zt = gr.weights_zm2zt;
  // The thermodynamic levels are always between two momentum levels.
  // Linear interpolation is used.
  Kokkos::parallel_for(
      "calc_zm2zt_weights", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        // Weight of upper momentum level on thermodynamic level.
        weights_zm2zt(i, k, m_above) = (zt(i, k) - zm(i, k)) / (zm(i, k + 1) - zm(i, k));
        // Weight of lower momentum level on thermodynamic level.
        weights_zm2zt(i, k, m_below) = (zm(i, k + 1) - zt(i, k)) / (zm(i, k + 1) - zm(i, k));
      });
}
