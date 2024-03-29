function [ cost_function, err_terms, err_sums, invsigma2 ] = ...
   min_les_hoc_diff( hoc_dir, les_dir, c_total, c_name, weight_case, ...
                     t1_case, t2_case, v_total, les_varname, hoc_varname, ...
                     weight_var, linitialize_sigma, invsigma2 )

%-----------------------------------------------------------------------

for c_run = 1:1:c_total
   
   % Path and Filename for GrADS Header Files for each case.
   hoc_file_header = [ hoc_dir, '/', deblank(c_name(c_run,:)), '_zt.ctl' ];
   if ( strcmp( deblank(c_name(c_run,:)), 'jun25_altocu' ) )
      les_file_header = [ les_dir, '/', deblank(c_name(c_run,:)), ...
                          '_qc3_coamps_sm.ctl' ];
%      les_file_header = [ les_dir, '/', deblank(c_name(c_run,:)), ...
%                          '_qcbar_coamps_sm.ctl' ];
   else
      les_file_header = [ les_dir, '/', deblank(c_name(c_run,:)), ...
                          '_coamps_sm.ctl' ];
   end

   % Determine how large the GrADS input is
   [ hoc_filename, hoc_nz, hoc_z, hoc_t_time_steps, ...
     hoc_time_step_length, hoc_numvars, hoc_listofparams ] ...
   = header_read_expanded(hoc_file_header);
   [ les_filename, les_nz, les_z, les_t_time_steps, ...
     les_time_step_length, les_numvars, les_listofparams ] ...
   = header_read_expanded(les_file_header);

   % Start with first HOC & LES variables, then loop through and
   % calculate the mean squared difference for all the variables
   for i = 1:1:v_total

      % Find the varnum from the varname
      % LES
      for num = 1:1:les_numvars
         if ( strcmp( deblank(les_varname(i,:)), deblank(les_listofparams(num,:)) ) )
            les_varnum = num;
            break
         elseif ( num == les_numvars )
            [ 'There was no variable matching ', deblank(les_varname(i,:)), ...
              ' in the LES file.  Please recheck both the variable name', ...
              ' and the file.' ]
         end
      end
      % HOC
      for num = 1:1:hoc_numvars
         if ( strcmp( deblank(hoc_varname(i,:)), deblank(hoc_listofparams(num,:)) ) )
            hoc_varnum = num;
            break
         elseif ( num == hoc_numvars )
            [ 'There was no variable matching ', deblank(hoc_varname(i,:)), ...
              ' in the HOC file.  Please recheck both the variable name', ...
              ' and the file.' ]
         end
      end
      
      % Read in LES grads data for one variable, averaged
      % over specified time intervals
      les_zl = read_grads_hoc_endian([les_dir, '/', les_filename], ...
               'ieee-be', les_nz, t1_case(c_run), t2_case(c_run), ...
               les_varnum, les_numvars);

      % Read in HOC grads data for one variable, averaged
      % over specified time intervals
      hoc_zl = read_grads_hoc_endian([hoc_dir, '/', hoc_filename], ...
               'ieee-le', hoc_nz, t1_case(c_run), t2_case(c_run), ...
               hoc_varnum, hoc_numvars);

      % The same variable, but squared and then averaged
      % over specified time intervals
      hoc2_zl = read_grads_hoc_endian_2([hoc_dir, '/', hoc_filename], ...
                'ieee-le', hoc_nz, t1_case(c_run), t2_case(c_run), ...
                hoc_varnum, hoc_numvars);

%-----------------------------------------------------------------------

      % Calculate the mean squared difference between the HOC
      % and the LES variables

      % In order to deal with differences in order of magnitude
      % between the variables, the err_sum equation has been
      % modified to normalize the values with respect to the
      % the minimum and maximum in the LES. -Dave Schanen

      les_minmax = max(les_zl) - min(les_zl);

      if ( les_minmax == 0.0 ) then
         %stop "An LES var = 0 over all z-levels"
         [ 'An LES var = 0 over all z-levels' ]
      end

%            ! Old code
%!           err_sum = err_sum
%!    .              + mean_sqr_diff_zt( hoc_nz, les_nz, hoc_zl,
%!    .                                  les_zl, les_minmax )
%
%            ! Chris Golaz modification: mean_sqr_diff_2_zt was designed to try
%            ! to limit time noise in tuning simulations.
%            ! New code
%!           err_sum = err_sum
%!    .              + mean_sqr_diff_2_zt
%!    .                ( hoc_nz, les_nz, hoc_zl, hoc2_zl,
%!    .                  les_zl, les_minmax )
%            ! End new code
%            ! Modification for weighting
      err_sums(c_run,i) = mean_sqr_diff_2_zt( hoc_nz, les_nz, hoc_zl, ...
                                              hoc2_zl, les_zl, les_minmax );

      %c_terms = c_terms + 1

   end % i=1..v_total

end % end of do c_run=1, c_total

%----------------------------------------------------------------------

% Return error averaged over all cases, variables, and vertical levels
% Old Code
%!       min_les_hoc_diff = err_sum / real( c_terms )

%---------------------------------------------------------------
% Compute normalization factors to satisfy imposed weights
%---------------------------------------------------------------
if ( linitialize_sigma == 1 )
   for c_run = 1:1:c_total
      for i = 1:1:v_total
          invsigma2(c_run,i) ...
             = weight_case(c_run)*weight_var(i) / err_sums(c_run,i);
      end
   end
end

%---------------------------------------------------------------
% Compute normalized error
%---------------------------------------------------------------
% I changed the variable name for the MATLAB code in order to preserve
% the pure value of the err_sums while still allowing for a normalization
% of those values.  Brian
%err_sums = invsigma2 * err_sums
%err_sum  = sum( err_sums )
err_sums_norm = invsigma2 .* err_sums;
err_sum_norm  = 0.0;
for c_run = 1:1:c_total
   for i = 1:1:v_total
      err_sum_norm = err_sum_norm + err_sums_norm(c_run,i);
   end
end

%---------------------------------------------------------------
% Save total error and error contributions breakdown
%---------------------------------------------------------------
%err_terms = err_sums
err_terms = err_sums_norm;
%min_les_hoc_diff = err_sum
cost_function = err_sum_norm;
