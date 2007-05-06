% Compares the time-averaged profiles between any LES, HOC Golaz
% "best-ever", HOC 12/17/2005, HOC previous (prior CVS HOC), and HOC 
% current files for 18 different variables fields.  It can be easily run 
% for any case in the HOC repository.  It outputs three pages of graphs 
% per case (with six graphs per page).  The program outputs in single-page
% eps files, single-page jpeg files, or three-page ps files.
% This program was made by Brian Griffin.
function compare_plots_cases( case_name, t1_min, t2_min, graphbase, graphtop, ...
                              dir_LES, dir_cgbe, dir_1217, dir_prev, dir_curr, ...
                              cmp_les, les_type, cmp_cgbe, cmp_1217, cmp_prev, cmp_curr )

% Function call variable description:
%
% 1st entry:   Case name.
% 2nd entry:   First time of variable comparison (in minutes since the
%              beginning of the run).
% 3rd entry:   Last time of variable comparison (in minutes since the
%              beginning of the run).
% 4th entry:   The height of the bottom of the output graphs for each case.
% 5th entry:   The height of the top of the output graphs for each case.
% 6th entry:   The directory path that contains the LES files for
%              comparison. (omit the '/' at the very end.)
% 7th entry:   The directory path that contains the HOC "Chris Golaz
%              'best-ever'" files for comparison. (omit the '/' at the very
%              end.)
% 8th entry:   The directory path that contains the HOC December 17, 2005
%              files for comparison. (omit the '/' at the very end.)
% 9th entry:   The directory path that contains the HOC previous files
%              (from the prior CVS commitment) for comparison. (omit the
%              '/' at the very end.)
% 10th entry:  The directory path that contains the HOC current files (from 
%              the latest CVS commitment) for comparison. (omit the '/' at
%              the very end.)
% 11th entry:  LES file for comparison:  1 = true; 0 = false.
% 12th entry:  LES being using for comparison (supports coamps and rams).
% 13th entry:  Chris Golaz "best-ever" HOC file for comparison:  1 = true;
%                                                                0 = false.
% 14th entry:  December 17, 2005 HOC file for comparison:  1 = true;
%                                                         0 = false.
% 15th entry:  Previous CVS HOC file for comparison:  1 = true; 
%                                                    0 = false.
% 16th entry:  Current HOC file for comparison:  1 = true; 0 = false.

% A listing of cases that have LES files:
%
% ARM
% ATEX
% BOMEX
% DYCOMS2 RF01
% DYCOMS2 RF02 DO
% DYCOMS2 RF02 DS
% DYCOMS2 RF02 ND
% FIRE
% GABLS
% June 25 altocumulus QC3
% June 25 altocumulus QCBAR
% November 11 altocumulus
% RICO
% Wangara
%
% Note:  All LES cases are from COAMPS LES, with the one exception of the 
%        Wangara case, which is from RAMS LES.
%
% A listing of cases that have HOC files from Chris Golaz' "best ever" HOC
% from October 7, 2005.
%
% ARM
% ATEX
% BOMEX
% DYCOMS2 RF01
% DYCOMS2 RF02 DO
% DYCOMS2 RF02 DS
% DYCOMS2 RF02 ND
% DYCOMS2 RF02 SO
% FIRE
% November 11 altocumulus
% Wangara
%
% A listing of cases that have HOC files from the December 17, 2005
% version, which was used to submit the HOC results to the DYCOMS2 RF02 
% SCM intercomparison.
%
% ARM
% ATEX
% BOMEX
% DYCOMS2 RF01
% DYCOMS2 RF02 DO
% DYCOMS2 RF02 DS
% DYCOMS2 RF02 ND
% DYCOMS2 RF02 SO
% FIRE
% November 11 altocumulus
% Wangara
%
% Therefore, the cases that contain output from all three categories are:
%
% ARM
% ATEX
% BOMEX
% DYCOMS2 RF01
% DYCOMS2 RF02 DO
% DYCOMS2 RF02 DS
% DYCOMS2 RF02 ND
% FIRE
% November 11 altocumulus
% Wangara
%
% The cases that have LES output, but no previous benchmark HOC output are:
%
% GABLS
% June 25 altocumulus QC3
% June 25 altocumulus QCBAR
% RICO
%
% The cases that have previous benchmark HOC output, but no previous 
% LES output are:
%
% DYCOMS2 RF02 DS

% Print the case name
[ 'Case:  ', case_name ]
[ '-------------------------------' ]

% File names
if ( cmp_les == 1 )
   if ( strcmp( case_name, 'jun25_altocu' ) )
      case_name_les = [ case_name, '_qc3' ];
%      case_name_les = [ case_name, '_qcbar' ];
   else
      case_name_les = case_name;
   end
   if ( strcmp( les_type, 'coamps' ) )
      file_header_les = [dir_LES, '/', case_name_les, '_coamps_sm.ctl'];
      file_header_les2 = [dir_LES, '/', case_name_les, '_coamps_sw.ctl'];
   elseif ( strcmp( les_type, 'rams' ) )
      file_header_les = [dir_LES, '/', case_name_les, '_rams.ctl'];
   end
end
if ( cmp_cgbe == 1 )
   file_header_cgbe_zt = [dir_cgbe, '/', case_name, '_zt.ctl'];
   file_header_cgbe_zm = [dir_cgbe, '/', case_name, '_zm.ctl'];
end
if ( cmp_1217 == 1 )
   file_header_1217_zt = [dir_1217, '/', case_name, '_zt.ctl'];
   file_header_1217_zm = [dir_1217, '/', case_name, '_zm.ctl'];
end
if ( cmp_prev == 1 )
   file_header_prev_zt = [dir_prev, '/', case_name, '_zt.ctl'];
   file_header_prev_zm = [dir_prev, '/', case_name, '_zm.ctl'];
end
if ( cmp_curr == 1 )
   file_header_curr_zt = [dir_curr, '/', case_name, '_zt.ctl'];
   file_header_curr_zm = [dir_curr, '/', case_name, '_zm.ctl'];
end

% LES File
if ( cmp_les == 1 )
   % Main LES file (COAMPS sm or RAMS)
   [filename_les, nz_les, z_les, t_time_steps_les, ts_length_les, ...
   numvars_les, listofparams_les] = header_read_expanded(file_header_les);
   % Secondary LES file (COAMPS sw)
   if ( strcmp( les_type, 'coamps' ) )
      [filename_les2, nz_les2, z_les2, t_time_steps_les2, ts_length_les2, ...
      numvars_les2, listofparams_les2] = header_read_expanded(file_header_les2);
   end
end

% HOC output -- Chris Golaz "best ever", October 7, 2005
if ( cmp_cgbe == 1 )
   % Thermodynamic-level output
   [filename_cgbe_zt, nz_cgbe_zt, z_cgbe_zt, t_time_steps_cgbe_zt, ts_length_cgbe_zt, ...
   numvars_cgbe_zt, listofparams_cgbe_zt] = header_read_expanded(file_header_cgbe_zt);
   % Momentum-level output
   [filename_cgbe_zm, nz_cgbe_zm, z_cgbe_zm, t_time_steps_cgbe_zm, ts_length_cgbe_zm, ...
   numvars_cgbe_zm, listofparams_cgbe_zm] = header_read_expanded(file_header_cgbe_zm);
end

% HOC output -- December 17, 2005 version
if ( cmp_1217 == 1 )
   % Thermodynamic-level output
   [filename_1217_zt, nz_1217_zt, z_1217_zt, t_time_steps_1217_zt, ts_length_1217_zt, ... 
   numvars_1217_zt, listofparams_1217_zt] = header_read_expanded(file_header_1217_zt);
   % Momentum-level output
   [filename_1217_zm, nz_1217_zm, z_1217_zm, t_time_steps_1217_zm, ts_length_1217_zm, ...
   numvars_1217_zm, listofparams_1217_zm] = header_read_expanded(file_header_1217_zm);
end

% HOC output -- previous
if ( cmp_prev == 1 )
   % Thermodynamic-level output
   [filename_prev_zt, nz_prev_zt, z_prev_zt, t_time_steps_prev_zt, ts_length_prev_zt, ...
   numvars_prev_zt, listofparams_prev_zt] = header_read_expanded(file_header_prev_zt);
   % Momentum-level output
   [filename_prev_zm, nz_prev_zm, z_prev_zm, t_time_steps_prev_zm, ts_length_prev_zm, ...
   numvars_prev_zm, listofparams_prev_zm] = header_read_expanded(file_header_prev_zm);
end

% HOC output -- current
if ( cmp_curr == 1 )
   % Thermodynamic-level output
   [filename_curr_zt, nz_curr_zt, z_curr_zt, t_time_steps_curr_zt, ts_length_curr_zt, ...
   numvars_curr_zt, listofparams_curr_zt] = header_read_expanded(file_header_curr_zt);
   % Momentum-level output
   [filename_curr_zm, nz_curr_zm, z_curr_zm, t_time_steps_curr_zm, ts_length_curr_zm, ...
   numvars_curr_zm, listofparams_curr_zm] = header_read_expanded(file_header_curr_zm);
end

%==========================================================================

% A listing of variables that we want to compare output with.
%
% thlm
% rtm
% cf
% rcm
% wp2
% wp3
%
% wpthlp
% wprtp
% thlp2
% rtp2
% rtpthlp
% wm
%
% um
% vm
% upwp
% vpwp
% rrm
% Nrm

if ( cmp_les == 1 )

    if ( strcmp( les_type, 'coamps' ) )

      % COAMPS LES variable names (and string lengths for comparison):
      les_thlm        = 'thlm ';
      les_thlm_len    = 5;
      les_rtm         = 'qtm ';
      les_rtm_len     = 4;
      les_cf          = 'cf ';
      les_cf_len      = 3;
      les_rcm         = 'qcm ';
      les_rcm_len     = 4;
      les_wp2         = 'wp2 ';
      les_wp2_len     = 4;
      les_wp3         = 'wp3 ';
      les_wp3_len     = 4;
      les_wpthlp      = 'wpthlp ';
      les_wpthlp_len  = 7;
      les_wprtp       = 'wpqtp ';
      les_wprtp_len   = 6;
      les_thlp2       = 'thlp2 ';
      les_thlp2_len   = 6;
      les_rtp2        = 'qtp2 ';
      les_rtp2_len    = 5;
      les_rtpthlp     = 'qtpthlp ';
      les_rtpthlp_len = 8;
      les_wm          = 'wlsm ';
      les_wm_len      = 5;
      les_um          = 'um ';
      les_um_len      = 3;
      les_vm          = 'vm ';
      les_vm_len      = 3;
      les_upwp        = 'wpup ';
      les_upwp_len    = 5;
      les_vpwp        = 'wpvp ';
      les_vpwp_len    = 5;
      les_rrm         = 'qrm ';
      les_rrm_len     = 4;
      les_Nrm         = 'nrm ';
      les_Nrm_len     = 4;

   elseif ( strcmp(les_type, 'rams' ) )

      % RAMS LES variable names (and string lengths for comparison):
      les_thlm        = 'thlm ';
      les_thlm_len    = 5;
      les_rtm         = 'qtm ';
      les_rtm_len     = 4;
      les_cf          = 'cf ';
      les_cf_len      = 3;
      les_rcm         = 'qcm ';
      les_rcm_len     = 4;
      les_wp2         = 'wp2 ';
      les_wp2_len     = 4;
      les_wp3         = 'wp3 ';
      les_wp3_len     = 4;
      les_wpthlp      = 'wpthlp ';
      les_wpthlp_len  = 7;
      les_wprtp       = 'wprtp ';
      les_wprtp_len   = 6;
      les_thlp2       = 'thlp2 ';
      les_thlp2_len   = 6;
      les_rtp2        = 'qtp2 ';
      les_rtp2_len    = 5;
      les_rtpthlp     = 'rtpthlp ';
      les_rtpthlp_len = 8;
      les_wm          = 'wm ';
      les_wm_len      = 3;
      les_um          = 'um ';
      les_um_len      = 3;
      les_vm          = 'vm ';
      les_vm_len      = 3;
      les_upwp        = 'wpup ';
      les_upwp_len    = 5;
      les_vpwp        = 'wpvp ';
      les_vpwp_len    = 5;
      les_rrm         = 'qrm ';
      les_rrm_len     = 4;
      les_Nrm         = 'nrm ';
      les_Nrm_len     = 4;
   
    end

end

% Find the appropriate output timesteps in each file for the beginning and
% the ending of the averaging period.  The averaging period is sent into
% this function as a number of minutes.  The number of minutes divided by
% the output timestep for each file determines which timestep is
% appropriate.  Since HOC GrADS output measures an accumulation of
% statistics over the previous timestep period, integers for timestep
% values should always be rounded up.  This is the usefulness of the "ceil"
% command.

if ( cmp_les == 1 )

    if ( strcmp( les_type, 'coamps' ) )
        
       % COAMPS LES has the first GrADS output at the time that is one
       % GrADS output timestep after the initial runtime, so the formula
       % for elapsed time for the GrADS file is:  
       % elapsed_time = timestep * time_step_length;
       % which results in:  timestep = elapsed_time/time_step_length.

       t1_les = ceil(t1_min/ts_length_les);
       if ( t1_les < 1 )
          t1_les = 1;
       elseif ( t1_les > t_time_steps_les )
          t1_les = t_time_steps_les;
       end

       t2_les = ceil(t2_min/ts_length_les);
       if ( t2_les < 1 )
          t2_les = 1;
       elseif ( t2_les > t_time_steps_les )
          t2_les = t_time_steps_les;
       end
        
       t1_les2 = ceil(t1_min/ts_length_les2);
       if ( t1_les2 < 1 )
          t1_les2 = 1;
       elseif ( t1_les2 > t_time_steps_les2 )
          t1_les2 = t_time_steps_les2;
       end

       t2_les2 = ceil(t2_min/ts_length_les2);
       if ( t2_les2 < 1 )
          t2_les2 = 1;
       elseif ( t2_les2 > t_time_steps_les2 )
          t2_les2 = t_time_steps_les2;
       end
       
    elseif ( strcmp(les_type, 'rams' ) )
   
       % The RAMS LES files have the first GrADS output at the initial
       % runtime, so the formula for elapsed time for the GrADS file is:  
       % elapsed_time = ( timestep - 1 ) * time_step_length;
       % which results in:  timestep = (elapsed_time/time_step_length) + 1.
        
       t1_les = ceil(t1_min/ts_length_les) + 1;
       if ( t1_les < 1 )
          t1_les = 1;
       elseif ( t1_les > t_time_steps_les )
          t1_les = t_time_steps_les;
       end

       t2_les = ceil(t2_min/ts_length_les) + 1;
       if ( t2_les < 1 )
          t2_les = 1;
       elseif ( t2_les > t_time_steps_les )
          t2_les = t_time_steps_les;
       end
        
    end
    
end

if ( cmp_cgbe == 1 )

    % The HOC Chris Golaz "best-ever" files have the first GrADS output at
    % the initial runtime, so the formula for elapsed time for the GrADS
    % file is:  elapsed_time = ( timestep - 1 ) * time_step_length;
    % which results in:  timestep = (elapsed_time/time_step_length) + 1.
    
    t1_cgbe_zt = ceil(t1_min/ts_length_cgbe_zt) + 1;
    if ( t1_cgbe_zt < 1 )
       t1_cgbe_zt = 1;
    elseif ( t1_cgbe_zt > t_time_steps_cgbe_zt )
       t1_cgbe_zt = t_time_steps_cgbe_zt;
    end

    t2_cgbe_zt = ceil(t2_min/ts_length_cgbe_zt) + 1;
    if ( t2_cgbe_zt < 1 )
       t2_cgbe_zt = 1;
    elseif ( t2_cgbe_zt > t_time_steps_cgbe_zt )
       t2_cgbe_zt = t_time_steps_cgbe_zt;
    end
    
    t1_cgbe_zm = ceil(t1_min/ts_length_cgbe_zt) + 1;
    if ( t1_cgbe_zm < 1 )
       t1_cgbe_zm = 1;
    elseif ( t1_cgbe_zm > t_time_steps_cgbe_zm )
       t1_cgbe_zm = t_time_steps_cgbe_zm;
    end
    
    t2_cgbe_zm = ceil(t2_min/ts_length_cgbe_zt) + 1;
    if ( t2_cgbe_zm < 1 )
       t2_cgbe_zm = 1;
    elseif ( t2_cgbe_zm > t_time_steps_cgbe_zm )
       t2_cgbe_zm = t_time_steps_cgbe_zm;
    end

end

if ( cmp_1217 == 1 )
    
    % The HOC December 17, 2005 files have the first GrADS output at
    % the initial runtime, so the formula for elapsed time for the GrADS
    % file is:  elapsed_time = ( timestep - 1 ) * time_step_length;
    % which results in:  timestep = (elapsed_time/time_step_length) + 1.

    t1_1217_zt = ceil(t1_min/ts_length_1217_zt) + 1;
    if ( t1_1217_zt < 1 )
       t1_1217_zt = 1;
    elseif ( t1_1217_zt > t_time_steps_1217_zt )
       t1_1217_zt = t_time_steps_1217_zt;
    end

    t2_1217_zt = ceil(t2_min/ts_length_1217_zt) + 1;
    if ( t2_1217_zt < 1 )
       t2_1217_zt = 1;
    elseif ( t2_1217_zt > t_time_steps_1217_zt )
       t2_1217_zt = t_time_steps_1217_zt;
    end
    
    t1_1217_zm = ceil(t1_min/ts_length_1217_zt) + 1;
    if ( t1_1217_zm < 1 )
       t1_1217_zm = 1;
    elseif ( t1_1217_zm > t_time_steps_1217_zm )
       t1_1217_zm = t_time_steps_1217_zm;
    end
    
    t2_1217_zm = ceil(t2_min/ts_length_1217_zt) + 1;
    if ( t2_1217_zm < 1 )
       t2_1217_zm = 1;
    elseif ( t2_1217_zm > t_time_steps_1217_zm )
       t2_1217_zm = t_time_steps_1217_zm;
    end
    
end

if ( cmp_prev == 1 )
    
    % HOC now has the first GrADS output at the time that is one GrADS
    % output timestep after the initial runtime, so the formula for 
    % elapsed time for the GrADS file is:  
    % elapsed_time = timestep * time_step_length;
    % which results in:  timestep = elapsed_time/time_step_length.

    t1_prev_zt = ceil(t1_min/ts_length_prev_zt);
    if ( t1_prev_zt < 1 )
       t1_prev_zt = 1;
    elseif ( t1_prev_zt > t_time_steps_prev_zt )
       t1_prev_zt = t_time_steps_prev_zt;
    end

    t2_prev_zt = ceil(t2_min/ts_length_prev_zt);
    if ( t2_prev_zt < 1 )
       t2_prev_zt = 1;
    elseif ( t2_prev_zt > t_time_steps_prev_zt )
       t2_prev_zt = t_time_steps_prev_zt;
    end
    
    t1_prev_zm = ceil(t1_min/ts_length_prev_zt);
    if ( t1_prev_zm < 1 )
       t1_prev_zm = 1;
    elseif ( t1_prev_zm > t_time_steps_prev_zm )
       t1_prev_zm = t_time_steps_prev_zm;
    end
    
    t2_prev_zm = ceil(t2_min/ts_length_prev_zt);
    if ( t2_prev_zm < 1 )
       t2_prev_zm = 1;
    elseif ( t2_prev_zm > t_time_steps_prev_zm )
       t2_prev_zm = t_time_steps_prev_zm;
    end
    
end

if ( cmp_curr == 1 )

    % HOC now has the first GrADS output at the time that is one GrADS
    % output timestep after the initial runtime, so the formula for 
    % elapsed time for the GrADS file is:  
    % elapsed_time = timestep * time_step_length;
    % which results in:  timestep = elapsed_time/time_step_length.
    
    t1_curr_zt = ceil(t1_min/ts_length_curr_zt);
    if ( t1_curr_zt < 1 )
       t1_curr_zt = 1;
    elseif ( t1_curr_zt > t_time_steps_curr_zt )
       t1_curr_zt = t_time_steps_curr_zt;
    end

    t2_curr_zt = ceil(t2_min/ts_length_curr_zt);
    if ( t2_curr_zt < 1 )
       t2_curr_zt = 1;
    elseif ( t2_curr_zt > t_time_steps_curr_zt )
       t2_curr_zt = t_time_steps_curr_zt;
    end
    
    t1_curr_zm = ceil(t1_min/ts_length_curr_zt);
    if ( t1_curr_zm < 1 )
       t1_curr_zm = 1;
    elseif ( t1_curr_zm > t_time_steps_curr_zm )
       t1_curr_zm = t_time_steps_curr_zm;
    end
    
    t2_curr_zm = ceil(t2_min/ts_length_curr_zt);
    if ( t2_curr_zm < 1 )
       t2_curr_zm = 1;
    elseif ( t2_curr_zm > t_time_steps_curr_zm )
       t2_curr_zm = t_time_steps_curr_zm;
    end
    
end

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

% Find variable indices (from the list of params) and get results.

%--------------------------------------------------------------------------

% Liquid Water Potential Temperature (var "thlm")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_thlm_len), les_thlm ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable thlm not found in LES; value being set to 300.'
         avg_thlm_les(1:nz_les) = 300.0;
      elseif ( varfnd == 1 )
         avg_thlm_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                        'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:5), 'thlm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable thlm not found in HOC (Golaz best-ever); value being set to 300.'
         avg_thlm_cgbe(1:nz_cgbe_zt) = 300.0;
      elseif ( varfnd == 1 )
         avg_thlm_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                         'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:5), 'thlm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable thlm not found in HOC (12/17/2005); value being set to 300.'
         avg_thlm_1217(1:nz_1217_zt) = 300.0;
      elseif ( varfnd == 1 )
         avg_thlm_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                         'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:5), 'thlm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable thlm not found in HOC (previous); value being set to 300.'
         avg_thlm_prev(1:nz_prev_zt) = 300.0;
      elseif ( varfnd == 1 )
         avg_thlm_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                         'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:5), 'thlm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable thlm not found in HOC (current); value being set to 300.'
         avg_thlm_curr(1:nz_curr_zt) = 300.0;
      elseif ( varfnd == 1 )
         avg_thlm_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                         'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

%--------------------------------------------------------------------------
    
% Total Water Mixing Ratio (var "rtm")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_rtm_len), les_rtm ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable qtm not found in LES; value being set to 0.'
         avg_rtm_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtm_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                       'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:4), 'rtm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable rtm not found in HOC (Golaz best-ever); value being set to 0.'
         avg_rtm_cgbe(1:nz_cgbe_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtm_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                        'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:4), 'rtm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable rtm not found in HOC (12/17/2005); value being set to 0.'
         avg_rtm_1217(1:nz_1217_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtm_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                        'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:4), 'rtm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable rtm not found in HOC (previous); value being set to 0.'
         avg_rtm_prev(1:nz_prev_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtm_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                        'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:4), 'rtm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable rtm not found in HOC (current); value being set to 0.'
         avg_rtm_curr(1:nz_curr_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtm_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                        'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

%--------------------------------------------------------------------------

% Cloud Fraction (var "cf")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_cf_len), les_cf ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable cf not found in LES; value being set to 0.'
         avg_cf_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_cf_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                      'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:3), 'cf ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable cf not found in HOC (Golaz best-ever); value being set to 0.'
         avg_cf_cgbe(1:nz_cgbe_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_cf_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                       'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:3), 'cf ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable cf not found in HOC (12/17/2005); value being set to 0.'
         avg_cf_1217(1:nz_1217_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_cf_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                       'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:3), 'cf ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable cf not found in HOC (previous); value being set to 0.'
         avg_cf_prev(1:nz_prev_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_cf_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                       'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:3), 'cf ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable cf not found in HOC (current); value being set to 0.'
         avg_cf_curr(1:nz_curr_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_cf_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                       'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

%--------------------------------------------------------------------------

% Cloud Water Mixing Ratio (var "rcm")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_rcm_len), les_rcm ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable qcm not found in LES; value being set to 0.'
         avg_rcm_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_rcm_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                       'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:4), 'rcm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable rcm not found in HOC (Golaz best-ever); value being set to 0.'
         avg_rcm_cgbe(1:nz_cgbe_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rcm_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                        'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:4), 'rcm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable rcm not found in HOC (12.17/2005); value being set to 0.'
         avg_rcm_1217(1:nz_1217_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rcm_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                        'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:4), 'rcm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable rcm not found in HOC (previous); value being set to 0.'
         avg_rcm_prev(1:nz_prev_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rcm_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                        'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:4), 'rcm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable rcm not found in HOC (current); value being set to 0.'
         avg_rcm_curr(1:nz_curr_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rcm_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                        'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

%--------------------------------------------------------------------------

% w'^2 (var "wp2")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_wp2_len), les_wp2 ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable wp2 not found in LES; value being set to 0.'
         avg_wp2_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp2_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                       'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zm
      if ( strcmp( listofparams_cgbe_zm(i,1:4), 'wp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zm) & (varfnd == 0) )
         'variable wp2 not found in HOC (Golaz best-ever); value being set to 0.'
         avg_wp2_cgbe(1:nz_cgbe_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp2_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zm], ...
                        'ieee-le', nz_cgbe_zm, t1_cgbe_zm, t2_cgbe_zm, varnum, numvars_cgbe_zm);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zm
      if ( strcmp( listofparams_1217_zm(i,1:4), 'wp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zm) & (varfnd == 0) )
         'variable wp2 not found in HOC (12/17/2005); value being set to 0.'
         avg_wp2_1217(1:nz_1217_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp2_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zm], ...
                        'ieee-le', nz_1217_zm, t1_1217_zm, t2_1217_zm, varnum, numvars_1217_zm);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zm
      if ( strcmp( listofparams_prev_zm(i,1:4), 'wp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zm) & (varfnd == 0) )
         'variable wp2 not found in HOC (previous); value being set to 0.'
         avg_wp2_prev(1:nz_prev_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp2_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zm], ...
                        'ieee-le', nz_prev_zm, t1_prev_zm, t2_prev_zm, varnum, numvars_prev_zm);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zm
      if ( strcmp( listofparams_curr_zm(i,1:4), 'wp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zm) & (varfnd == 0) )
         'variable wp2 not found in HOC (current); value being set to 0.'
         avg_wp2_curr(1:nz_curr_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp2_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zm], ...
                        'ieee-le', nz_curr_zm, t1_curr_zm, t2_curr_zm, varnum, numvars_curr_zm);
         break
      end
   end
end

%--------------------------------------------------------------------------

% w'^3 (var "wp3")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_wp3_len), les_wp3 ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable wp3 not found in LES; value being set to 0.'
         avg_wp3_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp3_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                       'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:4), 'wp3 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable wp3 not found in HOC (Golaz best-ever); value being set to 0.'
         avg_wp3_cgbe(1:nz_cgbe_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp3_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                        'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:4), 'wp3 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable wp3 not found in HOC (12/17/2005); value being set to 0.'
         avg_wp3_1217(1:nz_1217_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp3_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                        'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:4), 'wp3 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable wp3 not found in HOC (previous); value being set to 0.'
         avg_wp3_prev(1:nz_prev_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp3_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                        'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:4), 'wp3 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable wp3 not found in HOC (current); value being set to 0.'
         avg_wp3_curr(1:nz_curr_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_wp3_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                        'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

%--------------------------------------------------------------------------

% w'thl' (var "wpthlp")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_wpthlp_len), les_wpthlp ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable wpthlp not found in LES; value being set to 0.'
         avg_wpthlp_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_wpthlp_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                          'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zm
      if ( strcmp( listofparams_cgbe_zm(i,1:7), 'wpthlp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zm) & (varfnd == 0) )
         'variable wpthlp not found in HOC (Golaz best-ever); value being set to 0.'
         avg_wpthlp_cgbe(1:nz_cgbe_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wpthlp_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zm], ...
                           'ieee-le', nz_cgbe_zm, t1_cgbe_zm, t2_cgbe_zm, varnum, numvars_cgbe_zm);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zm
      if ( strcmp( listofparams_1217_zm(i,1:7), 'wpthlp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zm) & (varfnd == 0) )
         'variable wpthlp not found in HOC (12/17/2005); value being set to 0.'
         avg_wpthlp_1217(1:nz_1217_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wpthlp_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zm], ...
                           'ieee-le', nz_1217_zm, t1_1217_zm, t2_1217_zm, varnum, numvars_1217_zm);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
  varfnd = 0;
   for i = 1:1:numvars_prev_zm
      if ( strcmp( listofparams_prev_zm(i,1:7), 'wpthlp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zm) & (varfnd == 0) )
         'variable wpthlp not found in HOC (previous); value being set to 0.'
         avg_wpthlp_prev(1:nz_prev_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wpthlp_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zm], ...
                           'ieee-le', nz_prev_zm, t1_prev_zm, t2_prev_zm, varnum, numvars_prev_zm);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
  varfnd = 0;
   for i = 1:1:numvars_curr_zm
      if ( strcmp( listofparams_curr_zm(i,1:7), 'wpthlp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zm) & (varfnd == 0) )
         'variable wpthlp not found in HOC (current); value being set to 0.'
         avg_wpthlp_curr(1:nz_curr_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wpthlp_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zm], ...
                           'ieee-le', nz_curr_zm, t1_curr_zm, t2_curr_zm, varnum, numvars_curr_zm);
         break
      end
   end
end

%--------------------------------------------------------------------------

% w'rt' (var "wprtp")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_wprtp_len), les_wprtp ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable wpqtp not found in LES; value being set to 0.'
         avg_wprtp_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_wprtp_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                         'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zm
      if ( strcmp( listofparams_cgbe_zm(i,1:6), 'wprtp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zm) & (varfnd == 0) )
         'variable wprtp not found in HOC (Golaz best-ever); value being set to 0.'
         avg_wprtp_cgbe(1:nz_cgbe_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wprtp_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zm], ...
                          'ieee-le', nz_cgbe_zm, t1_cgbe_zm, t2_cgbe_zm, varnum, numvars_cgbe_zm);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zm
      if ( strcmp( listofparams_1217_zm(i,1:6), 'wprtp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zm) & (varfnd == 0) )
         'variable wprtp not found in HOC (12/17/2005); value being set to 0.'
         avg_wprtp_1217(1:nz_1217_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wprtp_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zm], ...
                          'ieee-le', nz_1217_zm, t1_1217_zm, t2_1217_zm, varnum, numvars_1217_zm);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zm
      if ( strcmp( listofparams_prev_zm(i,1:6), 'wprtp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zm) & (varfnd == 0) )
         'variable wprtp not found in HOC (previous); value being set to 0.'
         avg_wprtp_prev(1:nz_prev_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wprtp_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zm], ...
                          'ieee-le', nz_prev_zm, t1_prev_zm, t2_prev_zm, varnum, numvars_prev_zm);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zm
      if ( strcmp( listofparams_curr_zm(i,1:6), 'wprtp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zm) & (varfnd == 0) )
         'variable wprtp not found in HOC (current); value being set to 0.'
         avg_wprtp_curr(1:nz_curr_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_wprtp_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zm], ...
                          'ieee-le', nz_curr_zm, t1_curr_zm, t2_curr_zm, varnum, numvars_curr_zm);
         break
      end
   end
end

%--------------------------------------------------------------------------

% thl'^2 (var "thlp2")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_thlp2_len), les_thlp2 ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable thlp2 not found in LES; value being set to 0.'
         avg_thlp2_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_thlp2_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                         'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zm
      if ( strcmp( listofparams_cgbe_zm(i,1:6), 'thlp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zm) & (varfnd == 0) )
         'variable thlp2 not found in HOC (Golaz best-ever); value being set to 0.'
         avg_thlp2_cgbe(1:nz_cgbe_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_thlp2_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zm], ...
                          'ieee-le', nz_cgbe_zm, t1_cgbe_zm, t2_cgbe_zm, varnum, numvars_cgbe_zm);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zm
      if ( strcmp( listofparams_1217_zm(i,1:6), 'thlp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zm) & (varfnd == 0) )
         'variable thlp2 not found in HOC (12/17/2005); value being set to 0.'
         avg_thlp2_1217(1:nz_1217_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_thlp2_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zm], ...
                          'ieee-le', nz_1217_zm, t1_1217_zm, t2_1217_zm, varnum, numvars_1217_zm);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zm
      if ( strcmp( listofparams_prev_zm(i,1:6), 'thlp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zm) & (varfnd == 0) )
         'variable thlp2 not found in HOC (previous); value being set to 0.'
         avg_thlp2_prev(1:nz_prev_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_thlp2_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zm], ...
                          'ieee-le', nz_prev_zm, t1_prev_zm, t2_prev_zm, varnum, numvars_prev_zm);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zm
      if ( strcmp( listofparams_curr_zm(i,1:6), 'thlp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zm) & (varfnd == 0) )
         'variable thlp2 not found in HOC (current); value being set to 0.'
         avg_thlp2_curr(1:nz_curr_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_thlp2_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zm], ...
                          'ieee-le', nz_curr_zm, t1_curr_zm, t2_curr_zm, varnum, numvars_curr_zm);
         break
      end
   end
end

%--------------------------------------------------------------------------

% rt'^2 (var "rtp2")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_rtp2_len), les_rtp2 ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable qtp2 not found in LES; value being set to 0.'
         avg_rtp2_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtp2_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                        'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zm
      if ( strcmp( listofparams_cgbe_zm(i,1:5), 'rtp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zm) & (varfnd == 0) )
         'variable rtp2 not found in HOC (Golaz best-ever); value being set to 0.'
         avg_rtp2_cgbe(1:nz_cgbe_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtp2_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zm], ...
                         'ieee-le', nz_cgbe_zm, t1_cgbe_zm, t2_cgbe_zm, varnum, numvars_cgbe_zm);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zm
      if ( strcmp( listofparams_1217_zm(i,1:5), 'rtp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zm) & (varfnd == 0) )
         'variable rtp2 not found in HOC (12/17/2005); value being set to 0.'
         avg_rtp2_1217(1:nz_1217_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtp2_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zm], ...
                         'ieee-le', nz_1217_zm, t1_1217_zm, t2_1217_zm, varnum, numvars_1217_zm);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zm
      if ( strcmp( listofparams_prev_zm(i,1:5), 'rtp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zm) & (varfnd == 0) )
         'variable rtp2 not found in HOC (previous); value being set to 0.'
         avg_rtp2_prev(1:nz_prev_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtp2_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zm], ...
                         'ieee-le', nz_prev_zm, t1_prev_zm, t2_prev_zm, varnum, numvars_prev_zm);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zm
      if ( strcmp( listofparams_curr_zm(i,1:5), 'rtp2 ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zm) & (varfnd == 0) )
         'variable rtp2 not found in HOC (current); value being set to 0.'
         avg_rtp2_curr(1:nz_curr_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtp2_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zm], ...
                         'ieee-le', nz_curr_zm, t1_curr_zm, t2_curr_zm, varnum, numvars_curr_zm);
         break
      end
   end
end

%--------------------------------------------------------------------------

% rt'thl' (var "rtpthlp")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_rtpthlp_len), les_rtpthlp ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable qtpthlp not found in LES; value being set to 0.'
         avg_rtpthlp_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtpthlp_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                           'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zm
      if ( strcmp( listofparams_cgbe_zm(i,1:8), 'rtpthlp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zm) & (varfnd == 0) )
         'variable rtpthlp not found in HOC (Golaz best-ever); value being set to 0.'
         avg_rtpthlp_cgbe(1:nz_cgbe_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtpthlp_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zm], ...
                            'ieee-le', nz_cgbe_zm, t1_cgbe_zm, t2_cgbe_zm, varnum, numvars_cgbe_zm);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zm
      if ( strcmp( listofparams_1217_zm(i,1:8), 'rtpthlp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zm) & (varfnd == 0) )
         'variable rtpthlp not found in HOC (12/17/2005); value being set to 0.'
         avg_rtpthlp_1217(1:nz_1217_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtpthlp_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zm], ...
                            'ieee-le', nz_1217_zm, t1_1217_zm, t2_1217_zm, varnum, numvars_1217_zm);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zm
      if ( strcmp( listofparams_prev_zm(i,1:8), 'rtpthlp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zm) & (varfnd == 0) )
         'variable rtpthlp not found in HOC (previous); value being set to 0.'
         avg_rtpthlp_prev(1:nz_prev_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtpthlp_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zm], ...
                            'ieee-le', nz_prev_zm, t1_prev_zm, t2_prev_zm, varnum, numvars_prev_zm);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zm
      if ( strcmp( listofparams_curr_zm(i,1:8), 'rtpthlp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zm) & (varfnd == 0) )
         'variable rtpthlp not found in HOC (current); value being set to 0.'
         avg_rtpthlp_curr(1:nz_curr_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_rtpthlp_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zm], ...
                            'ieee-le', nz_curr_zm, t1_curr_zm, t2_curr_zm, varnum, numvars_curr_zm);
         break
      end
   end
end

%--------------------------------------------------------------------------

% Vertical Wind Component (var "wm") (due to large-scale subsidence)
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_wm_len), les_wm ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable wm not found in LES; value being set to 0.'
         avg_wm_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_wm_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                      'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:3), 'wm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable wm not found in HOC (Golaz best-ever); value being set to 0.'
         avg_wm_cgbe(1:nz_cgbe_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_wm_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                       'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:3), 'wm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable wm not found in HOC (12/17/2005); value being set to 0.'
         avg_wm_1217(1:nz_1217_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_wm_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                       'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:3), 'wm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable wm not found in HOC (previous); value being set to 0.'
         avg_wm_prev(1:nz_prev_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_wm_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                       'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:3), 'wm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable wm not found in HOC (current); value being set to 0.'
         avg_wm_curr(1:nz_curr_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_wm_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                       'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

%--------------------------------------------------------------------------

% Zonal Wind Component (var "um")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_um_len), les_um ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable um not found in LES; value being set to 0.'
         avg_um_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_um_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                      'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:3), 'um ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable um not found in HOC (Golaz best-ever); value being set to 0.'
         avg_um_cgbe(1:nz_cgbe_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_um_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                       'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:3), 'um ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable um not found in HOC (12/17/2005); value being set to 0.'
         avg_um_1217(1:nz_1217_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_um_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                       'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:3), 'um ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable um not found in HOC (previous); value being set to 0.'
         avg_um_prev(1:nz_prev_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_um_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                       'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:3), 'um ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable um not found in HOC (current); value being set to 0.'
         avg_um_curr(1:nz_curr_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_um_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                       'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

%--------------------------------------------------------------------------

% Meridional Wind Component (var "vm")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_vm_len), les_vm ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable vm not found in LES; value being set to 0.'
         avg_vm_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_vm_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                      'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:3), 'vm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable vm not found in HOC (Golaz best-ever); value being set to 0.'
         avg_vm_cgbe(1:nz_cgbe_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_vm_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                       'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:3), 'vm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable vm not found in HOC (12/17/2005); value being set to 0.'
         avg_vm_1217(1:nz_1217_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_vm_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                       'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:3), 'vm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable vm not found in HOC (previous); value being set to 0.'
         avg_vm_prev(1:nz_prev_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_vm_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                       'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:3), 'vm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable vm not found in HOC (current); value being set to 0.'
         avg_vm_curr(1:nz_curr_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_vm_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                       'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

%--------------------------------------------------------------------------
    
% u'w' (var "upwp")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_upwp_len), les_upwp ) )
         varnum = i;
         varfnd = 1;
         avg_upwp_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                       'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         les_grph_upwp = 1;
         break
      end
   end
   if ( (varfnd == 0) & ( strcmp( les_type, 'coamps' ) ) )
      for i = 1:1:numvars_les2
         if ( strcmp( listofparams_les2(i,1:les_upwp_len), les_upwp ) )
            varnum = i;
            varfnd = 1;
         end
         if ( (i == numvars_les2) & (varfnd == 0) )
            'variable upwp not found in LES; value being set to 0.'
            avg_upwp_les(1:nz_les2) = 0.0;
            les_grph_upwp = 2;
         elseif ( varfnd == 1 )
            avg_upwp_les = read_grads_hoc_endian([dir_LES, '/', filename_les2], ...
                          'ieee-be', nz_les2, t1_les2, t2_les2, varnum, numvars_les2);
            les_grph_upwp = 2;
            break
         end
      end
   elseif ( (varfnd == 0) & ( strcmp( les_type, 'rams' ) ) )
      'variable upwp not found in LES; value being set to 0.'
      avg_upwp_les(1:nz_les) = 0.0;
      les_grph_upwp = 1;
   end
end
% HOC -- Golaz "best ever"
   if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zm
      if ( strcmp( listofparams_cgbe_zm(i,1:5), 'upwp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zm) & (varfnd == 0) )
         'variable upwp not found in HOC (Golaz best-ever); value being set to 0.'
         avg_upwp_cgbe(1:nz_cgbe_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_upwp_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zm], ...
                        'ieee-le', nz_cgbe_zm, t1_cgbe_zm, t2_cgbe_zm, varnum, numvars_cgbe_zm);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zm
      if ( strcmp( listofparams_1217_zm(i,1:5), 'upwp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zm) & (varfnd == 0) )
         'variable upwp not found in HOC (12/17/2005); value being set to 0.'
         avg_upwp_1217(1:nz_1217_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_upwp_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zm], ...
                         'ieee-le', nz_1217_zm, t1_1217_zm, t2_1217_zm, varnum, numvars_1217_zm);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zm
      if ( strcmp( listofparams_prev_zm(i,1:5), 'upwp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zm) & (varfnd == 0) )
         'variable upwp not found in HOC (previous); value being set to 0.'
         avg_upwp_prev(1:nz_prev_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_upwp_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zm], ...
                        'ieee-le', nz_prev_zm, t1_prev_zm, t2_prev_zm, varnum, numvars_prev_zm);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zm
      if ( strcmp( listofparams_curr_zm(i,1:5), 'upwp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zm) & (varfnd == 0) )
         'variable upwp not found in HOC (current); value being set to 0.'
         avg_upwp_curr(1:nz_curr_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_upwp_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zm], ...
                        'ieee-le', nz_curr_zm, t1_curr_zm, t2_curr_zm, varnum, numvars_curr_zm);
         break
      end
   end
end

%--------------------------------------------------------------------------
    
% v'w' (var "vpwp")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_vpwp_len), les_vpwp ) )
         varnum = i;
         varfnd = 1;
         avg_vpwp_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                       'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         les_grph_vpwp = 1;
         break
      end
   end
   if ( (varfnd == 0) & ( strcmp( les_type, 'coamps' ) ) )
      for i = 1:1:numvars_les2
         if ( strcmp( listofparams_les2(i,1:les_vpwp_len), les_vpwp ) )
            varnum = i;
            varfnd = 1;
         end
         if ( (i == numvars_les2) & (varfnd == 0) )
            'variable vpwp not found in LES; value being set to 0.'
            avg_vpwp_les(1:nz_les2) = 0.0;
            les_grph_vpwp = 2;
         elseif ( varfnd == 1 )
            avg_vpwp_les = read_grads_hoc_endian([dir_LES, '/', filename_les2], ...
                          'ieee-be', nz_les2, t1_les2, t2_les2, varnum, numvars_les2);
            les_grph_vpwp = 2;
            break
         end
      end
   elseif ( (varfnd == 0) & ( strcmp( les_type, 'rams' ) ) )
      'variable vpwp not found in LES; value being set to 0.'
      avg_vpwp_les(1:nz_les) = 0.0;
      les_grph_vpwp = 1;
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zm
      if ( strcmp( listofparams_cgbe_zm(i,1:5), 'vpwp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zm) & (varfnd == 0) )
         'variable vpwp not found in HOC (Golaz best-ever); value being set to 0.'
         avg_vpwp_cgbe(1:nz_cgbe_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_vpwp_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zm], ...
                        'ieee-le', nz_cgbe_zm, t1_cgbe_zm, t2_cgbe_zm, varnum, numvars_cgbe_zm);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zm
      if ( strcmp( listofparams_1217_zm(i,1:5), 'vpwp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zm) & (varfnd == 0) )
         'variable vpwp not found in HOC (12/17/2005); value being set to 0.'
         avg_vpwp_1217(1:nz_1217_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_vpwp_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zm], ...
                         'ieee-le', nz_1217_zm, t1_1217_zm, t2_1217_zm, varnum, numvars_1217_zm);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zm
      if ( strcmp( listofparams_prev_zm(i,1:5), 'vpwp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zm) & (varfnd == 0) )
         'variable vpwp not found in HOC (previous); value being set to 0.'
         avg_vpwp_prev(1:nz_prev_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_vpwp_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zm], ...
                        'ieee-le', nz_prev_zm, t1_prev_zm, t2_prev_zm, varnum, numvars_prev_zm);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zm
      if ( strcmp( listofparams_curr_zm(i,1:5), 'vpwp ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zm) & (varfnd == 0) )
         'variable vpwp not found in HOC (current); value being set to 0.'
         avg_vpwp_curr(1:nz_curr_zm) = 0.0;
      elseif ( varfnd == 1 )
         avg_vpwp_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zm], ...
                        'ieee-le', nz_curr_zm, t1_curr_zm, t2_curr_zm, varnum, numvars_curr_zm);
         break
      end
   end
end

%--------------------------------------------------------------------------

% Rain water mixing ratio (var "rrm")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_rrm_len), les_rrm ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable qrm not found in LES; value being set to 0.'
         avg_rrm_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_rrm_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                       'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:4), 'rrm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable rrm not found in HOC (Golaz best-ever); value being set to 0.'
         avg_rrm_cgbe(1:nz_cgbe_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rrm_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                        'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:4), 'rrm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable rrm not found in HOC (12/17/2005); value being set to 0.'
         avg_rrm_1217(1:nz_1217_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rrm_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                        'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:4), 'rrm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable rrm not found in HOC (previous); value being set to 0.'
         avg_rrm_prev(1:nz_prev_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rrm_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                        'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:4), 'rrm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable rrm not found in HOC (current); value being set to 0.'
         avg_rrm_curr(1:nz_curr_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_rrm_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                        'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

%--------------------------------------------------------------------------

% Rain drop concentration (var "Nrm")
% LES
if ( cmp_les == 1 )
   varfnd = 0;
   for i = 1:1:numvars_les
      if ( strcmp( listofparams_les(i,1:les_Nrm_len), les_Nrm ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_les) & (varfnd == 0) )
         'variable Nrm not found in LES; value being set to 0.'
         avg_Nrm_les(1:nz_les) = 0.0;
      elseif ( varfnd == 1 )
         avg_Nrm_les = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
                       'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
         break
      end
   end
end
% HOC -- Golaz "best ever"
if ( cmp_cgbe == 1 )
   varfnd = 0;
   for i = 1:1:numvars_cgbe_zt
      if ( strcmp( listofparams_cgbe_zt(i,1:4), 'Nrm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_cgbe_zt) & (varfnd == 0) )
         'variable Nrm not found in HOC (Golaz best-ever); value being set to 0.'
         avg_Nrm_cgbe(1:nz_cgbe_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_Nrm_cgbe = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe_zt], ...
                        'ieee-le', nz_cgbe_zt, t1_cgbe_zt, t2_cgbe_zt, varnum, numvars_cgbe_zt);
         break
      end
   end
end
% HOC -- December 17, 2005
if ( cmp_1217 == 1 )
   varfnd = 0;
   for i = 1:1:numvars_1217_zt
      if ( strcmp( listofparams_1217_zt(i,1:4), 'Nrm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_1217_zt) & (varfnd == 0) )
         'variable Nrm not found in HOC (12/17/2005); value being set to 0.'
         avg_Nrm_1217(1:nz_1217_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_Nrm_1217 = read_grads_hoc_endian([dir_1217, '/', filename_1217_zt], ...
                        'ieee-le', nz_1217_zt, t1_1217_zt, t2_1217_zt, varnum, numvars_1217_zt);
         break
      end
   end
end
% HOC -- previous
if ( cmp_prev == 1 )
   varfnd = 0;
   for i = 1:1:numvars_prev_zt
      if ( strcmp( listofparams_prev_zt(i,1:4), 'Nrm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_prev_zt) & (varfnd == 0) )
         'variable Nrm not found in HOC (previous); value being set to 0.'
         avg_Nrm_prev(1:nz_prev_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_Nrm_prev = read_grads_hoc_endian([dir_prev, '/', filename_prev_zt], ...
                        'ieee-le', nz_prev_zt, t1_prev_zt, t2_prev_zt, varnum, numvars_prev_zt);
         break
      end
   end
end
% HOC -- current
if ( cmp_curr == 1 )
   varfnd = 0;
   for i = 1:1:numvars_curr_zt
      if ( strcmp( listofparams_curr_zt(i,1:4), 'Nrm ' ) )
         varnum = i;
         varfnd = 1;
      end
      if ( (i == numvars_curr_zt) & (varfnd == 0) )
         'variable Nrm not found in HOC (current); value being set to 0.'
         avg_Nrm_curr(1:nz_curr_zt) = 0.0;
      elseif ( varfnd == 1 )
         avg_Nrm_curr = read_grads_hoc_endian([dir_curr, '/', filename_curr_zt], ...
                        'ieee-le', nz_curr_zt, t1_curr_zt, t2_curr_zt, varnum, numvars_curr_zt);
         break
      end
   end
end

% Adjustment:  COAMPS LES outputs Nrm in num/cm^3.  This factor needs to be
%              multiplied by 10^6 in order to be converted to num/m^3.
if ( cmp_les == 1 )
   if ( strcmp( les_type, 'coamps' ) )
      avg_Nrm_les = (10^6).*avg_Nrm_les;
   end
end

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

% Plot the results.

% Figure Properties for screen display.
scr_size = get(0,'ScreenSize');
fig_height = scr_size(4);
fig_width = (6.5/9.0) * fig_height;
fig_width = int16(fig_width);

% Open figure to set size.
figure('Position',[ 0 0 fig_width fig_height ])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperPosition', [ 1.0 1.0 6.5 9.0 ])

% The following two parameters are factors that increase or decrease the 
% amount of white space on the left or on the right of the inside of 
% each graph.
% The variable 'percentage' is the percentage of the x-axis that is 
% covered by the variable being graphed.  If 90% is specified, that means
% that 90% of the x-axis is covered by the variable being graphed.  The 
% remaining 10% is white space on either side.  Since the graphs are 
% always centered, 5% is white space on the left and 5% is white space on 
% the right.  Percentage must be a value that is greater than 0 and less 
% than or equal to 100.  
percentage = 95; %
if ( ( percentage <= 0 ) | ( percentage > 100 ) )
   ['Invalid percentage entered (for variable extent on the x-axis): ', ...
    percentage ]
end
% Normalize units to a factor of 1.
percentage = percentage/100;
% The variable 'equiv_space' is the pure number of units of the variable 
% being graphed to be displayed on either side of the vertical line if 
% the variable is only equal to one value in the graph.  For example, if 
% all values of thlm at every height for every data set are equal to 300 K,
% then the plot will only show vertical lines at 300 K.  The number of
% units (K, in this case) to the left or to the right (all white space)
% will be determined by 'equiv_space'.  For this case, a value of 0.01
% would make the graph range between 299.99 and 300.01.  The value of
% 'equiv_space' must be a positive number.
equiv_space = 0.01;

% thlm
subplot(3,2,1)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_thlm_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_thlm_les(1:graphtopidx));
   maxval(i) = max(avg_thlm_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_thlm_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(avg_thlm_cgbe(1:graphtopidx));
   maxval(i) = max(avg_thlm_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_thlm_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(avg_thlm_1217(1:graphtopidx));
   maxval(i) = max(avg_thlm_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_thlm_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(avg_thlm_prev(1:graphtopidx));
   maxval(i) = max(avg_thlm_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_thlm_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(avg_thlm_curr(1:graphtopidx));
   maxval(i) = max(avg_thlm_curr(1:graphtopidx));
end
hold off
% Brian's New Universal Legend (for output page 1).
legend( h, legend_text, 'Location', 'SouthEast' )
% Axis labels and graph title.
xlabel('thlm    [K]')
ylabel('Height    [m]')
title('Liquid Water Potential Temperature, \theta_l')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% rtm
subplot(3,2,2)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_rtm_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_rtm_les(1:graphtopidx));
   maxval(i) = max(avg_rtm_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_rtm_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(avg_rtm_cgbe(1:graphtopidx));
   maxval(i) = max(avg_rtm_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_rtm_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(avg_rtm_1217(1:graphtopidx));
   maxval(i) = max(avg_rtm_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_rtm_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(avg_rtm_prev(1:graphtopidx));
   maxval(i) = max(avg_rtm_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_rtm_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(avg_rtm_curr(1:graphtopidx));
   maxval(i) = max(avg_rtm_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('rtm    [kg/kg]')
ylabel('Height    [m]')
title('Total Water Mixing Ratio, r_{ t}')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% cf
subplot(3,2,3)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( 100.*avg_cf_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(100.*avg_cf_les(1:graphtopidx));
   maxval(i) = max(100.*avg_cf_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( 100.*avg_cf_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(100.*avg_cf_cgbe(1:graphtopidx));
   maxval(i) = max(100.*avg_cf_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( 100.*avg_cf_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(100.*avg_cf_1217(1:graphtopidx));
   maxval(i) = max(100.*avg_cf_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( 100.*avg_cf_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(100.*avg_cf_prev(1:graphtopidx));
   maxval(i) = max(100.*avg_cf_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( 100.*avg_cf_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(100.*avg_cf_curr(1:graphtopidx));
   maxval(i) = max(100.*avg_cf_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('cf    [%]')
ylabel('Height    [m]')
title('Cloud Fraction')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% rcm
subplot(3,2,4)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_rcm_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_rcm_les(1:graphtopidx));
   maxval(i) = max(avg_rcm_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_rcm_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(avg_rcm_cgbe(1:graphtopidx));
   maxval(i) = max(avg_rcm_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_rcm_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(avg_rcm_1217(1:graphtopidx));
   maxval(i) = max(avg_rcm_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_rcm_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(avg_rcm_prev(1:graphtopidx));
   maxval(i) = max(avg_rcm_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_rcm_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(avg_rcm_curr(1:graphtopidx));
   maxval(i) = max(avg_rcm_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('rcm    [kg/kg]')
ylabel('Height    [m]')
title('Cloud Water Mixing Ratio, r_c')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% wp2
subplot(3,2,5)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_wp2_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_wp2_les(1:graphtopidx));
   maxval(i) = max(avg_wp2_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_wp2_cgbe, z_cgbe_zm, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zm
      if ( z_cgbe_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zm )
         graphtopidx = nz_cgbe_zm;
      end
   end
   minval(i) = min(avg_wp2_cgbe(1:graphtopidx));
   maxval(i) = max(avg_wp2_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_wp2_1217, z_1217_zm, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zm
      if ( z_1217_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zm )
         graphtopidx = nz_1217_zm;
      end
   end
   minval(i) = min(avg_wp2_1217(1:graphtopidx));
   maxval(i) = max(avg_wp2_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_wp2_prev, z_prev_zm, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zm
      if ( z_prev_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zm )
         graphtopidx = nz_prev_zm;
      end
   end
   minval(i) = min(avg_wp2_prev(1:graphtopidx));
   maxval(i) = max(avg_wp2_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_wp2_curr, z_curr_zm, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zm
      if ( z_curr_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zm )
         graphtopidx = nz_curr_zm;
      end
   end
   minval(i) = min(avg_wp2_curr(1:graphtopidx));
   maxval(i) = max(avg_wp2_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('wp2    [m^2/s^2]')
ylabel('Height    [m]')
title('Variance of w')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% wp3
subplot(3,2,6)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_wp3_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_wp3_les(1:graphtopidx));
   maxval(i) = max(avg_wp3_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_wp3_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(avg_wp3_cgbe(1:graphtopidx));
   maxval(i) = max(avg_wp3_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_wp3_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(avg_wp3_1217(1:graphtopidx));
   maxval(i) = max(avg_wp3_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_wp3_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(avg_wp3_prev(1:graphtopidx));
   maxval(i) = max(avg_wp3_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_wp3_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(avg_wp3_curr(1:graphtopidx));
   maxval(i) = max(avg_wp3_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('wp3    [m^3/s^3]')
ylabel('Height    [m]')
title('Third-order Moment of w')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

% Print 1st Page for Case (.eps document)
outputfilename = [ 'output/eps/', case_name, '_page1.eps' ];
print( '-depsc2', outputfilename )
% Print 1st Page for Case (.jpg document)
outputfilename = [ 'output/jpg/', case_name, '_page1.jpg' ];
print( '-djpeg', outputfilename )
% Print 3-page Document for each Case (.ps document)
outputfilename = [ 'output/ps/', case_name, '_doc.ps' ];
print( '-dpsc', outputfilename )

% Close the current figure (so it doesn't produce a new one for each page
% and for each case).
closereq

%--------------------------------------------------------------------------

figure('Position',[ 0 0 fig_width fig_height ])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperPosition', [ 1.0 1.0 6.5 9.0 ])

% wpthlp
subplot(3,2,1)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_wpthlp_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_wpthlp_les(1:graphtopidx));
   maxval(i) = max(avg_wpthlp_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_wpthlp_cgbe, z_cgbe_zm, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zm
      if ( z_cgbe_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zm )
         graphtopidx = nz_cgbe_zm;
      end
   end
   minval(i) = min(avg_wpthlp_cgbe(1:graphtopidx));
   maxval(i) = max(avg_wpthlp_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_wpthlp_1217, z_1217_zm, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zm
      if ( z_1217_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zm )
         graphtopidx = nz_1217_zm;
      end
   end
   minval(i) = min(avg_wpthlp_1217(1:graphtopidx));
   maxval(i) = max(avg_wpthlp_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_wpthlp_prev, z_prev_zm, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zm
      if ( z_prev_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zm )
         graphtopidx = nz_prev_zm;
      end
   end
   minval(i) = min(avg_wpthlp_prev(1:graphtopidx));
   maxval(i) = max(avg_wpthlp_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_wpthlp_curr, z_curr_zm, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zm
      if ( z_curr_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zm )
         graphtopidx = nz_curr_zm;
      end
   end
   minval(i) = min(avg_wpthlp_curr(1:graphtopidx));
   maxval(i) = max(avg_wpthlp_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('wpthlp    [K m/s]')
ylabel('Height    [m]')
title('Turbulent Flux of \theta_l')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% wprtp
subplot(3,2,2)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_wprtp_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_wprtp_les(1:graphtopidx));
   maxval(i) = max(avg_wprtp_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_wprtp_cgbe, z_cgbe_zm, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zm
      if ( z_cgbe_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zm )
         graphtopidx = nz_cgbe_zm;
      end
   end
   minval(i) = min(avg_wprtp_cgbe(1:graphtopidx));
   maxval(i) = max(avg_wprtp_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_wprtp_1217, z_1217_zm, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zm
      if ( z_1217_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zm )
         graphtopidx = nz_1217_zm;
      end
   end
   minval(i) = min(avg_wprtp_1217(1:graphtopidx));
   maxval(i) = max(avg_wprtp_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_wprtp_prev, z_prev_zm, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zm
      if ( z_prev_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zm )
         graphtopidx = nz_prev_zm;
      end
   end
   minval(i) = min(avg_wprtp_prev(1:graphtopidx));
   maxval(i) = max(avg_wprtp_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_wprtp_curr, z_curr_zm, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zm
      if ( z_curr_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zm )
         graphtopidx = nz_curr_zm;
      end
   end
   minval(i) = min(avg_wprtp_curr(1:graphtopidx));
   maxval(i) = max(avg_wprtp_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('wprtp    [(kg/kg) m/s]')
ylabel('Height    [m]')
title('Turbulent Flux of r_{ t}')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% thlp2
subplot(3,2,3)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_thlp2_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_thlp2_les(1:graphtopidx));
   maxval(i) = max(avg_thlp2_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_thlp2_cgbe, z_cgbe_zm, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zm
      if ( z_cgbe_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zm )
         graphtopidx = nz_cgbe_zm;
      end
   end
   minval(i) = min(avg_thlp2_cgbe(1:graphtopidx));
   maxval(i) = max(avg_thlp2_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_thlp2_1217, z_1217_zm, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zm
      if ( z_1217_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zm )
         graphtopidx = nz_1217_zm;
      end
   end
   minval(i) = min(avg_thlp2_1217(1:graphtopidx));
   maxval(i) = max(avg_thlp2_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_thlp2_prev, z_prev_zm, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zm
      if ( z_prev_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zm )
         graphtopidx = nz_prev_zm;
      end
   end
   minval(i) = min(avg_thlp2_prev(1:graphtopidx));
   maxval(i) = max(avg_thlp2_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_thlp2_curr, z_curr_zm, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zm
      if ( z_curr_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zm )
         graphtopidx = nz_curr_zm;
      end
   end
   minval(i) = min(avg_thlp2_curr(1:graphtopidx));
   maxval(i) = max(avg_thlp2_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('thlp2    [K^2]')
ylabel('Height    [m]')
title('Variance of \theta_l')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% rtp2
subplot(3,2,4)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_rtp2_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_rtp2_les(1:graphtopidx));
   maxval(i) = max(avg_rtp2_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_rtp2_cgbe, z_cgbe_zm, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zm
      if ( z_cgbe_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zm )
         graphtopidx = nz_cgbe_zm;
      end
   end
   minval(i) = min(avg_rtp2_cgbe(1:graphtopidx));
   maxval(i) = max(avg_rtp2_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_rtp2_1217, z_1217_zm, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zm
      if ( z_1217_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zm )
         graphtopidx = nz_1217_zm;
      end
   end
   minval(i) = min(avg_rtp2_1217(1:graphtopidx));
   maxval(i) = max(avg_rtp2_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_rtp2_prev, z_prev_zm, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zm
      if ( z_prev_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zm )
         graphtopidx = nz_prev_zm;
      end
   end
   minval(i) = min(avg_rtp2_prev(1:graphtopidx));
   maxval(i) = max(avg_rtp2_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_rtp2_curr, z_curr_zm, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zm
      if ( z_curr_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zm )
         graphtopidx = nz_curr_zm;
      end
   end
   minval(i) = min(avg_rtp2_curr(1:graphtopidx));
   maxval(i) = max(avg_rtp2_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('rtp2    [(kg/kg)^2]')
ylabel('Height    [m]')
title('Variance of r_t')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% rtpthlp
subplot(3,2,5)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_rtpthlp_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_rtpthlp_les(1:graphtopidx));
   maxval(i) = max(avg_rtpthlp_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_rtpthlp_cgbe, z_cgbe_zm, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zm
      if ( z_cgbe_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zm )
         graphtopidx = nz_cgbe_zm;
      end
   end
   minval(i) = min(avg_rtpthlp_cgbe(1:graphtopidx));
   maxval(i) = max(avg_rtpthlp_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_rtpthlp_1217, z_1217_zm, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zm
      if ( z_1217_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zm )
         graphtopidx = nz_1217_zm;
      end
   end
   minval(i) = min(avg_rtpthlp_1217(1:graphtopidx));
   maxval(i) = max(avg_rtpthlp_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_rtpthlp_prev, z_prev_zm, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zm
      if ( z_prev_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zm )
         graphtopidx = nz_prev_zm;
      end
   end
   minval(i) = min(avg_rtpthlp_prev(1:graphtopidx));
   maxval(i) = max(avg_rtpthlp_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_rtpthlp_curr, z_curr_zm, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zm
      if ( z_curr_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zm )
         graphtopidx = nz_curr_zm;
      end
   end
   minval(i) = min(avg_rtpthlp_curr(1:graphtopidx));
   maxval(i) = max(avg_rtpthlp_curr(1:graphtopidx));
end
hold off
% Brian's New Universal Legend (for output page 2).
legend( h, legend_text, 'Location', 'SouthWest' )
% Axis labels and graph title.
xlabel('rtpthlp    [(kg/kg) K]')
ylabel('Height    [m]')
title('Covariance of r_t & \theta_l')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% wm
subplot(3,2,6)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_wm_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_wm_les(1:graphtopidx));
   maxval(i) = max(avg_wm_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_wm_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(avg_wm_cgbe(1:graphtopidx));
   maxval(i) = max(avg_wm_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_wm_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(avg_wm_1217(1:graphtopidx));
   maxval(i) = max(avg_wm_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_wm_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(avg_wm_prev(1:graphtopidx));
   maxval(i) = max(avg_wm_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_wm_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(avg_wm_curr(1:graphtopidx));
   maxval(i) = max(avg_wm_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('wm    [m/s]')
ylabel('Height    [m]')
title('Vertical Wind Component, w (subsidence)')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

% Print 2nd Page for Case (.eps document)
outputfilename = [ 'output/eps/', case_name, '_page2.eps' ];
print( '-depsc2', outputfilename )
% Print 2nd Page for Case (.jpg document)
outputfilename = [ 'output/jpg/', case_name, '_page2.jpg' ];
print( '-djpeg', outputfilename )
% Print 3-page Document for each Case (.ps document)
outputfilename = [ 'output/ps/', case_name, '_doc.ps' ];
print( '-dpsc', '-append', outputfilename )

% Close the current figure (so it doesn't produce a new one for each page
% and for each case).
closereq

%--------------------------------------------------------------------------

figure('Position',[ 0 0 fig_width fig_height ])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperPosition', [ 1.0 1.0 6.5 9.0 ])

% um
subplot(3,2,1)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_um_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_um_les(1:graphtopidx));
   maxval(i) = max(avg_um_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_um_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(avg_um_cgbe(1:graphtopidx));
   maxval(i) = max(avg_um_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_um_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(avg_um_1217(1:graphtopidx));
   maxval(i) = max(avg_um_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_um_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(avg_um_prev(1:graphtopidx));
   maxval(i) = max(avg_um_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_um_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(avg_um_curr(1:graphtopidx));
   maxval(i) = max(avg_um_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('um    [m/s]')
ylabel('Height    [m]')
title('Zonal Wind Component, u')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% vm
subplot(3,2,2)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_vm_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_vm_les(1:graphtopidx));
   maxval(i) = max(avg_vm_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_vm_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(avg_vm_cgbe(1:graphtopidx));
   maxval(i) = max(avg_vm_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_vm_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(avg_vm_1217(1:graphtopidx));
   maxval(i) = max(avg_vm_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_vm_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(avg_vm_prev(1:graphtopidx));
   maxval(i) = max(avg_vm_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_vm_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(avg_vm_curr(1:graphtopidx));
   maxval(i) = max(avg_vm_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('vm    [m/s]')
ylabel('Height    [m]')
title('Meridional Wind Component, v')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% upwp
subplot(3,2,3)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   if ( les_grph_upwp == 1 )
      h(i) = plot( avg_upwp_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
      nz_upwp = nz_les;
      z_upwp = z_les;
   elseif ( les_grph_upwp == 2 )
      h(i) = plot( avg_upwp_les, z_les2, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
      nz_upwp = nz_les2;
      z_upwp = z_les2;
   end
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_upwp
      if ( z_upwp(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_upwp )
         graphtopidx = nz_upwp;
      end
   end
   minval(i) = min(avg_upwp_les(1:graphtopidx));
   maxval(i) = max(avg_upwp_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_upwp_cgbe, z_cgbe_zm, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zm
      if ( z_cgbe_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zm )
         graphtopidx = nz_cgbe_zm;
      end
   end
   minval(i) = min(avg_upwp_cgbe(1:graphtopidx));
   maxval(i) = max(avg_upwp_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_upwp_1217, z_1217_zm, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zm
      if ( z_1217_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zm )
         graphtopidx = nz_1217_zm;
      end
   end
   minval(i) = min(avg_upwp_1217(1:graphtopidx));
   maxval(i) = max(avg_upwp_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_upwp_prev, z_prev_zm, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zm
      if ( z_prev_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zm )
         graphtopidx = nz_prev_zm;
      end
   end
   minval(i) = min(avg_upwp_prev(1:graphtopidx));
   maxval(i) = max(avg_upwp_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_upwp_curr, z_curr_zm, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zm
      if ( z_curr_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zm )
         graphtopidx = nz_curr_zm;
      end
   end
   minval(i) = min(avg_upwp_curr(1:graphtopidx));
   maxval(i) = max(avg_upwp_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('upwp    [m^2/s^2]')
ylabel('Height    [m]')
title('Covariance of u & w')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% vpwp
subplot(3,2,4)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   if ( les_grph_vpwp == 1 )
      h(i) = plot( avg_vpwp_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
      nz_vpwp = nz_les;
      z_vpwp = z_les;
   elseif ( les_grph_vpwp == 2 )
      h(i) = plot( avg_vpwp_les, z_les2, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
      nz_vpwp = nz_les2;
      z_vpwp = z_les2;
   end
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_vpwp
      if ( z_vpwp(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_vpwp )
         graphtopidx = nz_vpwp;
      end
   end
   minval(i) = min(avg_vpwp_les(1:graphtopidx));
   maxval(i) = max(avg_vpwp_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_vpwp_cgbe, z_cgbe_zm, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zm
      if ( z_cgbe_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zm )
         graphtopidx = nz_cgbe_zm;
      end
   end
   minval(i) = min(avg_vpwp_cgbe(1:graphtopidx));
   maxval(i) = max(avg_vpwp_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_vpwp_1217, z_1217_zm, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zm
      if ( z_1217_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zm )
         graphtopidx = nz_1217_zm;
      end
   end
   minval(i) = min(avg_vpwp_1217(1:graphtopidx));
   maxval(i) = max(avg_vpwp_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_vpwp_prev, z_prev_zm, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zm
      if ( z_prev_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zm )
         graphtopidx = nz_prev_zm;
      end
   end
   minval(i) = min(avg_vpwp_prev(1:graphtopidx));
   maxval(i) = max(avg_vpwp_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_vpwp_curr, z_curr_zm, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zm
      if ( z_curr_zm(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zm )
         graphtopidx = nz_curr_zm;
      end
   end
   minval(i) = min(avg_vpwp_curr(1:graphtopidx));
   maxval(i) = max(avg_vpwp_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('vpwp    [m^2/s^2]')
ylabel('Height    [m]')
title('Covariance of v & w')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% rrm
subplot(3,2,5)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_rrm_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_rrm_les(1:graphtopidx));
   maxval(i) = max(avg_rrm_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_rrm_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(avg_rrm_cgbe(1:graphtopidx));
   maxval(i) = max(avg_rrm_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_rrm_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(avg_rrm_1217(1:graphtopidx));
   maxval(i) = max(avg_rrm_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_rrm_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(avg_rrm_prev(1:graphtopidx));
   maxval(i) = max(avg_rrm_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_rrm_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(avg_rrm_curr(1:graphtopidx));
   maxval(i) = max(avg_rrm_curr(1:graphtopidx));
end
hold off
% Axis labels and graph title.
xlabel('rrm    [kg/kg]')
ylabel('Height    [m]')
title('Rain Water Mixing Ratio, r_r')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

%--------------------------------------------------------------------------

% Nrm
subplot(3,2,6)
i = 0;
if ( cmp_les == 1 )
   i = i + 1;
   h(i) = plot( avg_Nrm_les, z_les, '-', 'Color', [ 1, 0, 0 ], 'LineWidth', 5 );
   legend_text(i,1:15) = '\fontsize{6}LES';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_les
      if ( z_les(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_les )
         graphtopidx = nz_les;
      end
   end
   minval(i) = min(avg_Nrm_les(1:graphtopidx));
   maxval(i) = max(avg_Nrm_les(1:graphtopidx));
end
if ( cmp_cgbe == 1 )
   i = i + 1;
   h(i) = plot( avg_Nrm_cgbe, z_cgbe_zt, '-', 'Color', [ 0, 0.50, 0 ], 'LineWidth', 3.5 );
   legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_cgbe_zt
      if ( z_cgbe_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_cgbe_zt )
         graphtopidx = nz_cgbe_zt;
      end
   end
   minval(i) = min(avg_Nrm_cgbe(1:graphtopidx));
   maxval(i) = max(avg_Nrm_cgbe(1:graphtopidx));
end
if ( cmp_1217 == 1 )
   i = i + 1;
   h(i) = plot( avg_Nrm_1217, z_1217_zt, '-.', 'Color', [ 0.63, 0, 0.79 ], 'LineWidth', 3.5 );
   legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_1217_zt
      if ( z_1217_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_1217_zt )
         graphtopidx = nz_1217_zt;
      end
   end
   minval(i) = min(avg_Nrm_1217(1:graphtopidx));
   maxval(i) = max(avg_Nrm_1217(1:graphtopidx));
end
if ( cmp_prev == 1 )
   i = i + 1;
   h(i) = plot( avg_Nrm_prev, z_prev_zt, '--', 'Color', [ 0.94, 0.50, 0.16], 'LineWidth', 2 );
   legend_text(i,1:24) = '\fontsize{6}HOC previous';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_prev_zt
      if ( z_prev_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_prev_zt )
         graphtopidx = nz_prev_zt;
      end
   end
   minval(i) = min(avg_Nrm_prev(1:graphtopidx));
   maxval(i) = max(avg_Nrm_prev(1:graphtopidx));
end
if ( cmp_curr == 1 )
   i = i + 1;
   h(i) = plot( avg_Nrm_curr, z_curr_zt, '-', 'Color', [ 0, 0.63, 1 ], 'LineWidth', 2 );
   legend_text(i,1:23) = '\fontsize{6}HOC current';
   hold on
   % Find vertical level index right below top of graph.
   for k = 2:1:nz_curr_zt
      if ( z_curr_zt(k) > graphtop )
         graphtopidx = k - 1;
         break
      elseif ( k == nz_curr_zt )
         graphtopidx = nz_curr_zt;
      end
   end
   minval(i) = min(avg_Nrm_curr(1:graphtopidx));
   maxval(i) = max(avg_Nrm_curr(1:graphtopidx));
end
hold off
% Brian's New Universal Legend (for output page 3).
legend( h, legend_text, 'Location', 'NorthEast' )
% Axis labels and graph title.
xlabel('Nrm    [num/m^3]')
ylabel('Height    [m]')
title('Rain Drop Concentration, N_r')
% Extent of graph.
xmin = min(minval);
xmax = max(maxval);
if ( xmax == xmin )
   xmin = xmin - equiv_space;
   xmax = xmax + equiv_space;
else
   xdiff = xmax - xmin;
   xrange = xdiff/percentage;
   xmedian = ( xmin + xmax ) / 2;
   xmin = xmedian - xrange/2;
   xmax = xmedian + xrange/2;
end
zmin = graphbase;
%zmax = max(z_les);
zmax = graphtop;
axis([ xmin xmax zmin zmax ])

% Print 3rd Page for Case (.eps document)
outputfilename = [ 'output/eps/', case_name, '_page3.eps' ];
print( '-depsc2', outputfilename )
% Print 3rd Page for Case (.jpg document)
outputfilename = [ 'output/jpg/', case_name, '_page3.jpg' ];
print( '-djpeg', outputfilename )
% Print 3-page Document for each Case (.ps document)
outputfilename = [ 'output/ps/', case_name, '_doc.ps' ];
print( '-dpsc', '-append', outputfilename )

% Close the current figure (so it doesn't produce a new one for each page
% and for each case).
closereq

% Statement
[ 'Case ', case_name, ' has successfully been output!' ]
[ '=====================================================================' ]
