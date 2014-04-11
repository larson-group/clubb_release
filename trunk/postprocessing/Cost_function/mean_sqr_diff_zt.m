function mean_sqr_diff_zt( hoc_nz, les_nz, hoc_zl, les_zl, norm_term )
%!       Description
%!       Calculate the mean squared difference between two input vectors,
%!       then normalize.
%
%!       References:
%!       None
%
%!       Notes:
%!       Configured to do interpolation on LES / HOC comparisons on the
%!       thermodynamic grid.  We use a modified version for tuning runs
%!       now, see below.
%!-----------------------------------------------------------------------
%        ! External
%        intrinsic sum
%
%        ! Input Variables
%        integer, intent(in) ::
%     .  hoc_nz, ! Vertical extent for HOC
%     .  les_nz  ! Vertical extent for the LES
%
%        real, intent(in), dimension(hoc_nz) ::
%     .  hoc_zl  ! HOC GrADS variable [units vary]
%
%        real, intent(in), dimension(les_nz) ::
%     .  les_zl  ! The LES GrADS variable [units vary]
%
%        real, intent(in) ::
%     .  norm_term ! normalization term;
%                  !typically maxval(les) - minval(les)
%
%        ! Local Variables
%        real, dimension(hoc_nz) ::
%     .  tmp_zl
%
%!----------------------------------------------------------------------

if ( ( hoc_nz - les_nz ) == 0 ) % most cases
   % Due to hoc's lower starting point, we can only use
   % (total number of z-levels) - 1 (a maximum of 74 for BOMEX).
   % The code below assumes the LES data are on an evenly spaced grid.
   % (Need to interpolate hoc to LES' levels.  Right now we just
   % compare adjacent z levels.  Vince Larson 12 Jan 2005)

   for k = 1:1:hoc_nz-1
      tmp_zl(k) = ( hoc_zl(k+1) - les_zl(k) ) / norm_term;
      tmp_zl(k) = tmp_zl(k)^2;
   end

   mean_sqr_diff_zt = 0;
   for k = 1:1:hoc_nz-1
      mean_sqr_diff_zt = mean_sqr_diff_zt + tmp_zl(k);
   end
   
elseif ( ( hoc_nz - les_nz ) == 2 ) % the DYCOMS II RF01 case
   
   for k = 1:1:les_nz
      tmp_zl(k) = ( hoc_zl(k+2) - les_zl(k) ) / norm_term;
      tmp_zl(k) = tmp_zl(k)^2;
   end

   mean_sqr_diff_zt = 0;
   for k = 1:1:les_nz
      mean_sqr_diff_zt = mean_sqr_diff_zt + tmp_zl(k);
   end
   
else

   %"HOC:", hoc_nz, "LES:", les_nz
   %stop "Not able to handle specified number of HOC z-levels"
   [ "Not able to handle specified number of HOC z-levels" ]
   
end
