% Compares the time-averaged profiles between any LES, HOC Golaz
% "best-ever", HOC 12/17/2005, HOC previous (prior CVS HOC), and HOC 
% current files for 18 different variables fields.  It can be easily run 
% for any case in the HOC repository.  It outputs three pages of graphs 
% per case (with six graphs per page).  The program outputs in single-page
% eps files, single-page jpeg files, or three-page ps files.
% This program was made by Brian Griffin.
function compare_plots_cases( case_name, t1_min, t2_min, graphbase_in, graphtop_in, ...
                              dir_LES_in, dir_cgbe_in, dir_1217_in, dir_prev_in, dir_curr_in, ...
                              cmp_les_in, les_type, cmp_cgbe_in, cmp_1217_in, cmp_prev_in, cmp_curr_in )

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

%Setup global variables
global dir_prev
global dir_curr
global dir_LES
global dir_cgbe
global dir_1217
global filename_les
global filename_les2
global numvars_les2
global listofparams_les2
global filename_cgbe_zt
global filename_1217_zt
global filename_prev_zt
global filename_curr_zt
global cmp_les
global cmp_cgbe
global cmp_1217
global cmp_prev
global cmp_curr
global les_upwp
global les_vpwp
global les_grph_upwp
global les_grph_vpwp
global z_les
global t1_les2
global t2_les2
global nz_les
global z_les2
global nz_les2
global z_cgbe_zt
global nz_cgbe_zt
global z_1217_zt
global nz_1217_zt
global z_prev_zt
global nz_prev_zt
global z_curr_zt
global nz_curr_zt
global les_color
global les_width
global cgbe_color
global cgbe_width
global dec12_color
global dec12_width
global prev_color
global prev_width
global curr_color
global curr_width
global graphtop
global graphbase
global equiv_space
global percentage

dir_LES = dir_LES_in;
dir_cgbe = dir_cgbe_in;
dir_1217 = dir_1217_in;
dir_prev = dir_prev_in;
dir_curr = dir_curr_in;
cmp_les = cmp_les_in;
cmp_cgbe = cmp_cgbe_in;
cmp_1217 = cmp_1217_in;
cmp_prev = cmp_prev_in;
cmp_curr = cmp_curr_in;
graphtop = graphtop_in;
graphbase = graphbase_in;
%End global setup

%Initialize variables to 0
%This prevents errors and allows the create_plot fuction
%to handle missing data
filename_les = 0;
filename_cgbe_zt = 0;
filename_1217_zt = 0;
filename_prev_zt = 0;
filename_curr_zt = 0;
filename_cgbe_zm = 0;
filename_1217_zm = 0;
filename_prev_zm = 0;
filename_curr_zm = 0;
filename_cgbe_sfc = 0;
filename_1217_sfc = 0;
filename_prev_sfc = 0;
filename_curr_sfc = 0;
avg_vm_les = 0;
avg_vm_cgbe = 0;
avg_vm_1217 = 0;
avg_vm_prev = 0;
avg_vm_curr = 0;
avg_upwp_les = 0;
avg_upwp_cgbe = 0;
avg_upwp_1217 = 0;
avg_upwp_prev = 0;
avg_upwp_curr = 0;
avg_vpwp_les = 0;
avg_vpwp_cgbe = 0;
avg_vpwp_1217 = 0;
avg_vpwp_prev = 0;
avg_vpwp_curr = 0;
avg_rrainm_les = 0;
avg_rrainm_cgbe = 0;
avg_rrainm_1217 = 0;
avg_rrainm_prev = 0;
avg_rrainm_curr = 0;
avg_Nrm_les = 0;
avg_Nrm_cgbe = 0;
avg_Nrm_1217 = 0;
avg_Nrm_prev = 0;
avg_Nrm_curr = 0;
avg_thlm_les = 0;
avg_thlm_cgbe = 0;
avg_thlm_1217 = 0;
avg_thlm_prev = 0;
avg_thlm_curr = 0;
avg_rtm_les = 0;
avg_rtm_cgbe = 0;
avg_rtm_1217 = 0;
avg_rtm_prev = 0;
avg_rtm_curr = 0;
avg_cf_les = 0;
avg_cf_cgbe = 0;
avg_cf_1217 = 0;
avg_cf_prev = 0;
avg_cf_curr = 0;
avg_rcm_les = 0;
avg_rcm_cgbe = 0;
avg_rcm_1217 = 0;
avg_rcm_prev = 0;
avg_rcm_curr = 0;
avg_wp2_les = 0;
avg_wp2_cgbe = 0;
avg_wp2_1217 = 0;
avg_wp2_prev = 0;
avg_wp2_curr = 0;
avg_wp3_les = 0;
avg_wp3_cgbe = 0;
avg_wp3_1217 = 0;
avg_wp3_prev = 0;
avg_wp3_curr = 0;
avg_wpthlp_les = 0;
avg_wpthlp_cgbe = 0;
avg_wpthlp_1217 = 0;
avg_wpthlp_prev = 0;
avg_wpthlp_curr = 0;
avg_wprtp_les = 0;
avg_wprtp_cgbe = 0;
avg_wprtp_1217 = 0;
avg_wprtp_prev = 0;
avg_wprtp_curr = 0;
avg_thlp2_les = 0;
avg_thlp2_cgbe = 0;
avg_thlp2_1217 = 0;
avg_thlp2_prev = 0;
avg_thlp2_curr = 0;
avg_rtp2_les = 0;
avg_rtp2_cgbe = 0;
avg_rtp2_1217 = 0;
avg_rtp2_prev = 0;
avg_rtp2_curr = 0;
avg_rtpthlp_les = 0;
avg_rtpthlp_cgbe = 0;
avg_rtpthlp_1217 = 0;
avg_rtpthlp_prev = 0;
avg_rtpthlp_curr = 0;
avg_wm_les = 0;
avg_wm_cgbe = 0;
avg_wm_1217 = 0;
avg_wm_prev = 0;
avg_wm_curr = 0;
avg_um_les = 0;
avg_um_cgbe = 0;
avg_um_1217 = 0;
avg_um_prev = 0;
avg_um_curr = 0;
listofparams_les = 0;
listofparams_cgbe_zt = 0;
listofparams_1217_zt = 0;
listofparams_prev_zt = 0;
listofparams_curr_zt = 0;
listofparams_cgbe_zm = 0;
listofparams_1217_zm = 0;
listofparams_prev_zm = 0;
listofparams_curr_zm = 0;
listofparams_cgbe_sfc = 0;
listofparams_1217_sfc = 0;
listofparams_prev_sfc = 0;
listofparams_curr_sfc = 0;
numvars_les = 0;
numvars_cgbe_zt = 0;
numvars_1217_zt = 0;
numvars_prev_zt = 0;
numvars_curr_zt = 0;
numvars_cgbe_zm = 0;
numvars_1217_zm = 0;
numvars_prev_zm = 0;
numvars_curr_zm = 0;
numvars_cgbe_sfc = 0;
numvars_1217_sfc = 0;
numvars_prev_sfc = 0;
numvars_curr_sfc = 0;
nz_les = 0;
t1_les = 0;
t2_les = 0;
nz_cgbe_zt = 0;
t1_cgbe_zt = 0;
t2_cgbe_zt = 0;
nz_1217_zt = 0;
t1_1217_zt = 0;
t2_1217_zt = 0;		
nz_prev_zt = 0;
t1_prev_zt = 0;
t2_prev_zt = 0;
nz_curr_zt = 0;
t1_curr_zt = 0;
t2_curr_zt = 0;
nz_cgbe_zm = 0;
t1_cgbe_zm = 0;
t2_cgbe_zm = 0;
nz_1217_zm = 0;
t1_1217_zm = 0;
t2_1217_zm = 0;		
nz_prev_zm = 0;
t1_prev_zm = 0;
t2_prev_zm = 0;
nz_curr_zm = 0;
t1_curr_zm = 0;
t2_curr_zm = 0;
nz_cgbe_sfc = 0;
t1_cgbe_sfc = 0;
t2_cgbe_sfc = 0;
nz_1217_sfc = 0;
t1_1217_sfc = 0;
t2_1217_sfc = 0;		
nz_prev_sfc = 0;
t1_prev_sfc = 0;
t2_prev_sfc = 0;
nz_curr_sfc = 0;
t1_curr_sfc = 0;
t2_curr_sfc = 0;
t1_les2 = 0;
t2_les2 = 0;
les_var_len = 0;
les_var = 0;
cgbe_var_len = 0;
cgbe_var = 0;
dec17_var_len = 0;
dec17_var = 0;
prev_var_len = 0;
prev_var = 0;
curr_var_len = 0;
curr_var = 0;

%Special cases
les_grph_upwp = 1;
les_grph_vpwp = 1;
%End 0 initialization


%Define the colors used in the plots
les_color = [ 1.00, 0.00, 0.00 ]; %Red
les_width = 5;

cgbe_color = [ 0.00, 0.50, 0.00 ]; %Green
cgbe_width = 3.5;

dec12_color = [ 0.63, 0.00, 0.79 ]; %Purple
dec12_width = 3.5;

prev_color = [ 0.94, 0.50, 0.16 ]; %Orange
prev_width = 3.5;

curr_color = [ 0.00, 0.63, 1.00 ]; %Blue
curr_width = 2;

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
   file_header_cgbe_sfc = [dir_cgbe, '/', case_name, '_sfc.ctl'];
end
if ( cmp_1217 == 1 )
   file_header_1217_zt = [dir_1217, '/', case_name, '_zt.ctl'];
   file_header_1217_zm = [dir_1217, '/', case_name, '_zm.ctl'];
   file_header_1217_sfc = [dir_1217, '/', case_name, '_sfc.ctl'];
end
if ( cmp_prev == 1 )
   file_header_prev_zt = [dir_prev, '/', case_name, '_zt.ctl'];
   file_header_prev_zm = [dir_prev, '/', case_name, '_zm.ctl'];
   file_header_prev_sfc = [dir_prev, '/', case_name, '_sfc.ctl'];
end
if ( cmp_curr == 1 )
   file_header_curr_zt = [dir_curr, '/', case_name, '_zt.ctl'];
   file_header_curr_zm = [dir_curr, '/', case_name, '_zm.ctl'];
   file_header_curr_sfc = [dir_curr, '/', case_name, '_sfc.ctl'];
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
   % Surface level output
   [filename_cgbe_sfc, nz_cgbe_sfc, z_cgbe_sfc, t_time_steps_cgbe_sfc, ts_length_cgbe_sfc, ...
   numvars_cgbe_sfc, listofparams_cgbe_sfc] = header_read_expanded(file_header_cgbe_sfc);
end

% HOC output -- December 17, 2005 version
if ( cmp_1217 == 1 )
   % Thermodynamic-level output
   [filename_1217_zt, nz_1217_zt, z_1217_zt, t_time_steps_1217_zt, ts_length_1217_zt, ... 
   numvars_1217_zt, listofparams_1217_zt] = header_read_expanded(file_header_1217_zt);
   % Momentum-level output
   [filename_1217_zm, nz_1217_zm, z_1217_zm, t_time_steps_1217_zm, ts_length_1217_zm, ...
   numvars_1217_zm, listofparams_1217_zm] = header_read_expanded(file_header_1217_zm);
   % Surface level output
   [filename_1217_sfc, nz_1217_sfc, z_1217_sfc, t_time_steps_1217_sfc, ts_length_1217_sfc, ...
   numvars_1217_sfc, listofparams_1217_sfc] = header_read_expanded(file_header_1217_sfc);
end

% HOC output -- previous
if ( cmp_prev == 1 )
   % Thermodynamic-level output
   [filename_prev_zt, nz_prev_zt, z_prev_zt, t_time_steps_prev_zt, ts_length_prev_zt, ...
   numvars_prev_zt, listofparams_prev_zt] = header_read_expanded(file_header_prev_zt);
   % Momentum-level output
   [filename_prev_zm, nz_prev_zm, z_prev_zm, t_time_steps_prev_zm, ts_length_prev_zm, ...
   numvars_prev_zm, listofparams_prev_zm] = header_read_expanded(file_header_prev_zm);
   % Surface level output
   [filename_prev_sfc, nz_prev_sfc, z_prev_sfc, t_time_steps_prev_sfc, ts_length_prev_sfc, ...
   numvars_prev_sfc, listofparams_prev_sfc] = header_read_expanded(file_header_prev_sfc);
end

% HOC output -- current
if ( cmp_curr == 1 )
   % Thermodynamic-level output
   [filename_curr_zt, nz_curr_zt, z_curr_zt, t_time_steps_curr_zt, ts_length_curr_zt, ...
   numvars_curr_zt, listofparams_curr_zt] = header_read_expanded(file_header_curr_zt);
   % Momentum-level output
   [filename_curr_zm, nz_curr_zm, z_curr_zm, t_time_steps_curr_zm, ts_length_curr_zm, ...
   numvars_curr_zm, listofparams_curr_zm] = header_read_expanded(file_header_curr_zm);
   % Surface level output
   [filename_curr_sfc, nz_curr_sfc, z_curr_sfc, t_time_steps_curr_sfc, ts_length_curr_sfc, ...
   numvars_curr_sfc, listofparams_curr_sfc] = header_read_expanded(file_header_curr_sfc);
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
% rrainm
% Nrm

% COAMPS LES variable names (and string lengths for comparison):
les_thlm        = 'thlm ';
les_rtm         = 'qtm ';
les_cf          = 'cf ';
les_rcm         = 'qcm ';
les_wp2         = 'wp2 ';
les_wp3         = 'wp3 ';
les_wpthlp      = 'wpthlp ';
les_wprtp       = 'wpqtp ';
les_thlp2       = 'thlp2 ';
les_rtp2        = 'qtp2 ';
les_rtpthlp     = 'qtpthlp ';
les_wm          = 'wlsm ';
les_um          = 'um ';
les_vm          = 'vm ';
les_upwp        = 'wpup ';
les_vpwp        = 'wpvp ';
les_rrainm      = 'qrm ';
les_Nrm         = 'nrm ';
les_up2		= 'up2 ';
les_vp2		= 'vp2 ';
les_lwp		= 'lwp ';

if ( strcmp(les_type, 'rams' ) )

	% RAMS LES variable names (and string lengths for comparison):
	les_thlm        = 'thlm ';
	les_rtm         = 'qtm ';
	les_cf          = 'cf ';
	les_rcm         = 'qcm ';
	les_wp2         = 'wp2 ';
	les_wp3         = 'wp3 ';
	les_wpthlp      = 'wpthlp ';
	les_wprtp       = 'wprtp ';
	les_thlp2       = 'thlp2 ';
	les_rtp2        = 'qtp2 ';
	les_rtpthlp     = 'rtpthlp ';
	les_wm          = 'wm ';
	les_um          = 'um ';
	les_vm          = 'vm ';
	les_upwp        = 'wpup ';
	les_vpwp        = 'wpvp ';
	les_rrainm      = 'qrm ';
	les_Nrm         = 'nrm ';
	les_up2		= 'up2 ';
	les_vp2		= 'vp2 ';

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
       
	timesteps_les(1) = 0;	
	for j=2:t_time_steps_les,	
       		timesteps_les(j) = timesteps_les(j-1) + ts_length_les;
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

    t1_cgbe_sfc = ceil(t1_min/ts_length_cgbe_sfc) + 1;
    if ( t1_cgbe_sfc < 1 )
       t1_cgbe_sfc = 1;
    elseif ( t1_cgbe_sfc > t_time_steps_cgbe_sfc )
       t1_cgbe_sfc = t_time_steps_cgbe_sfc;
    end
    
    t2_cgbe_sfc = ceil(t2_min/ts_length_cgbe_sfc) + 1;
    if ( t2_cgbe_sfc < 1 )
       t2_cgbe_sfc = 1;
    elseif ( t2_cgbe_zm > t_time_steps_cgbe_sfc )
       t2_cgbe_sfc = t_time_steps_cgbe_sfc;
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
   
    t1_1217_sfc = ceil(t1_min/ts_length_1217_sfc) + 1;
    if ( t1_1217_sfc < 1 )
       t1_1217_sfc = 1;
    elseif ( t1_1217_sfc > t_time_steps_1217_sfc )
       t1_1217_sfc = t_time_steps_1217_sfc;
    end
    
    t2_1217_sfc = ceil(t2_min/ts_length_1217_sfc) + 1;
    if ( t2_1217_sfc < 1 )
       t2_1217_sfc = 1;
    elseif ( t2_1217_sfc > t_time_steps_1217_sfc )
       t2_1217_sfc = t_time_steps_1217_sfc;
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
    
	t1_prev_sfc = ceil(t1_min/ts_length_prev_sfc);
	if ( t1_prev_sfc < 1 )
		t1_prev_sfc = 1;
	elseif ( t1_prev_sfc > t_time_steps_prev_sfc )
		t1_prev_sfc = t_time_steps_prev_sfc;
	end
    
	t2_prev_sfc = ceil(t2_min/ts_length_prev_sfc);
	if ( t2_prev_sfc < 1 )
		t2_prev_sfc = 1;
	elseif ( t2_prev_sfc > t_time_steps_prev_sfc )
		t2_prev_sfc = t_time_steps_prev_sfc;
	end

	timesteps_prev(1) = 0;	
	for j=2:t_time_steps_prev_zt,	
       		timesteps_prev(j) = timesteps_prev(j-1) + ts_length_prev_zt;
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
    
	t1_curr_sfc = ceil(t1_min/ts_length_curr_sfc);
	if ( t1_curr_sfc < 1 )
		t1_curr_sfc = 1;
	elseif ( t1_curr_sfc > t_time_steps_curr_sfc )
		t1_curr_sfc = t_time_steps_curr_sfc;
	end
    
	t2_curr_sfc = ceil(t2_min/ts_length_curr_sfc);
	if ( t2_curr_sfc < 1 )
		t2_curr_sfc = 1;
	elseif ( t2_curr_sfc > t_time_steps_curr_sfc )
		t2_curr_sfc = t_time_steps_curr_sfc;
	end

	timesteps_curr(1) = 0;	
	for j=2:t_time_steps_curr_zt,	
       		timesteps_curr(j) = timesteps_curr(j-1) + ts_length_curr_zt;
	end

end

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

% Load each value and plot the results.

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

% Liquid Water Potential Temperature (var "thlm") - this plot has a legend
[avg_thlm_les, avg_thlm_cgbe, avg_thlm_1217, avg_thlm_prev, avg_thlm_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_thlm, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'thlm ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'thlm ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'thlm ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'thlm ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );

create_plot(3, 2, 1, 'Liquid Water Potential Temperature, \theta_l', 'thlm    [K]', avg_thlm_les, z_les, nz_les, avg_thlm_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_thlm_1217, z_1217_zt, nz_1217_zt, ...
		avg_thlm_prev, z_prev_zt, nz_prev_zt, avg_thlm_curr, z_curr_zt, nz_curr_zt, 1);

%--------------------------------------------------------------------------

% Total Water Mixing Ratio (var "rtm")
[avg_rtm_les, avg_rtm_cgbe, avg_rtm_1217, avg_rtm_prev, avg_rtm_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_rtm, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'rtm ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'rtm ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'rtm ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'rtm ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );


create_plot(3, 2, 2, 'Total Water Mixing Ratio, r_{ t}', 'rtm    [kg/kg]', avg_rtm_les, z_les, nz_les, avg_rtm_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_rtm_1217, z_1217_zt, nz_1217_zt, ...
		avg_rtm_prev, z_prev_zt, nz_prev_zt, avg_rtm_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% Cloud Fraction (var "cf")
[avg_cf_les, avg_cf_cgbe, avg_cf_1217, avg_cf_prev, avg_cf_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_cf, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'cf ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'cf ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'cf ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'cf ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );

create_plot(3, 2, 3, 'Cloud Fraction', 'cf    [%]', 100*avg_cf_les, z_les, nz_les, 100*avg_cf_cgbe, z_cgbe_zt, nz_cgbe_zt, 100*avg_cf_1217, z_1217_zt, nz_1217_zt, ...
		100*avg_cf_prev, z_prev_zt, nz_prev_zt, 100*avg_cf_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% Cloud Water Mixing Ratio (var "rcm")
[avg_rcm_les, avg_rcm_cgbe, avg_rcm_1217, avg_rcm_prev, avg_rcm_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_rcm, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'rcm ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'rcm ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'rcm ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'rcm ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );

create_plot(3, 2, 4, 'Cloud Water Mixing Ratio, r_c', 'rcm    [kg/kg]', avg_rcm_les, z_les, nz_les, avg_rcm_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_rcm_1217, z_1217_zt, nz_1217_zt, ...
		avg_rcm_prev, z_prev_zt, nz_prev_zt, avg_rcm_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% w'^2 (var "wp2")
[avg_wp2_les, avg_wp2_cgbe, avg_wp2_1217, avg_wp2_prev, avg_wp2_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_wp2, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'wp2 ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'wp2 ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'wp2 ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'wp2 ', nz_curr_zm, t1_curr_zm, t2_curr_zm, les_type, 0 );

create_plot(3, 2, 5, 'Variance of w', 'wp2    [m^2/s^2]', avg_wp2_les, z_les, nz_les, avg_wp2_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_wp2_1217, z_1217_zt, nz_1217_zt, ...
		avg_wp2_prev, z_prev_zt, nz_prev_zt, avg_wp2_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% w'^3 (var "wp3")
[avg_wp3_les, avg_wp3_cgbe, avg_wp3_1217, avg_wp3_prev, avg_wp3_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_wp3, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'wp3 ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'wp3 ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'wp3 ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'wp3 ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );

create_plot(3, 2, 6, 'Third-order Moment of w', 'wp3    [m^3/s^3]', avg_wp3_les, z_les, nz_les, avg_wp3_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_wp3_1217, z_1217_zt, nz_1217_zt, ...
		avg_wp3_prev, z_prev_zt, nz_prev_zt, avg_wp3_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

outputfilename = [ case_name, '_page1' ];
print_page(outputfilename, 1, 0, 1);

figure('Position',[ 0 0 fig_width fig_height ])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperPosition', [ 1.0 1.0 6.5 9.0 ])

% w'thl' (var "wpthlp")
[avg_wpthlp_les, avg_wpthlp_cgbe, avg_wpthlp_1217, avg_wpthlp_prev, avg_wpthlp_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_wpthlp, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'wpthlp ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'wpthlp ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'wpthlp ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'wpthlp ', nz_curr_zt, t1_curr_zm, t2_curr_zm, les_type, 0 );

create_plot(3, 2, 1, 'Turbulent Flux of \theta_l', 'wpthlp    [K m/s]', avg_wpthlp_les, z_les, nz_les, avg_wpthlp_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_wpthlp_1217, z_1217_zt, nz_1217_zt, ...
		avg_wpthlp_prev, z_prev_zt, nz_prev_zt, avg_wpthlp_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% w'rt' (var "wprtp")
[avg_wprtp_les, avg_wprtp_cgbe, avg_wprtp_1217, avg_wprtp_prev, avg_wprtp_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_wprtp, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'wprtp ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'wprtp ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'wprtp ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'wprtp ', nz_curr_zm, t1_curr_zm, t2_curr_zm, les_type, 0 );

create_plot(3, 2, 2, 'Turbulent Flux of r_{ t}', 'wprtp    [(kg/kg) m/s]', avg_wprtp_les, z_les, nz_les, avg_wprtp_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_wprtp_1217, z_1217_zt, nz_1217_zt, ...
		avg_wprtp_prev, z_prev_zt, nz_prev_zt, avg_wprtp_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% thl'^2 (var "thlp2")
[avg_thlp2_les, avg_thlp2_cgbe, avg_thlp2_1217, avg_thlp2_prev, avg_thlp2_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_thlp2, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'thlp2 ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'thlp2 ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'thlp2 ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'thlp2 ', nz_curr_zm, t1_curr_zm, t2_curr_zm, les_type, 0 );

create_plot(3, 2, 3, 'Variance of \theta_l', 'thlp2    [K^2]', avg_thlp2_les, z_les, nz_les, avg_thlp2_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_thlp2_1217, z_1217_zt, nz_1217_zt, ...
		avg_thlp2_prev, z_prev_zt, nz_prev_zt, avg_thlp2_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% rt'^2 (var "rtp2")
[avg_rtp2_les, avg_rtp2_cgbe, avg_rtp2_1217, avg_rtp2_prev, avg_rtp2_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_rtp2, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'rtp2 ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'rtp2 ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'rtp2 ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'rtp2 ', nz_curr_zm, t1_curr_zm, t2_curr_zm, les_type, 0 );

create_plot(3, 2, 4, 'Variance of r_t', 'rtp2    [(kg/kg)^2]', avg_rtp2_les, z_les, nz_les, avg_rtp2_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_rtp2_1217, z_1217_zt, nz_1217_zt, ...
		avg_rtp2_prev, z_prev_zt, nz_prev_zt, avg_rtp2_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% rt'thl' (var "rtpthlp") - this plot has a legend
[avg_rtpthlp_les, avg_rtpthlp_cgbe, avg_rtpthlp_1217, avg_rtpthlp_prev, avg_rtpthlp_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_rtpthlp, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'rtpthlp ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'rtpthlp ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'rtpthlp ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'rtpthlp ', nz_curr_zm, t1_curr_zm, t2_curr_zm, les_type, 0 );

create_plot(3, 2, 5, 'Covariance of r_t & \theta_l', 'rtpthlp    [(kg/kg) K]', avg_rtpthlp_les, z_les, nz_les, avg_rtpthlp_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_rtpthlp_1217, z_1217_zt, nz_1217_zt, ...
		avg_rtpthlp_prev, z_prev_zt, nz_prev_zt, avg_rtpthlp_curr, z_curr_zt, nz_curr_zt, 1);

%--------------------------------------------------------------------------

% Vertical Wind Component (var "wm") (due to large-scale subsidence)
[avg_wm_les, avg_wm_cgbe, avg_wm_1217, avg_wm_prev, avg_wm_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_wm, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'wm ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'wm ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'wm ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'wm ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );

create_plot(3, 2, 6, 'Vertical Wind Component, w (subsidence)', 'wm    [m/s]', avg_wm_les, z_les, nz_les, avg_wm_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_wm_1217, z_1217_zt, nz_1217_zt, ...
		avg_wm_prev, z_prev_zt, nz_prev_zt, avg_wm_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

outputfilename = [ case_name, '_page2' ];
print_page(outputfilename, 1, 0, 1);

figure('Position',[ 0 0 fig_width fig_height ])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperPosition', [ 1.0 1.0 6.5 9.0 ])

% Zonal Wind Component (var "um")
[avg_um_les, avg_um_cgbe, avg_um_1217, avg_um_prev, avg_um_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_um, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'um ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'um ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'um ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'um ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );

create_plot(3, 2, 1, 'Zonal Wind Component, u', 'um    [m/s]', avg_um_les, z_les, nz_les, avg_um_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_um_1217, z_1217_zt, nz_1217_zt, ...
		avg_um_prev, z_prev_zt, nz_prev_zt, avg_um_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% Meridional Wind Component (var "vm")
[avg_vm_les, avg_vm_cgbe, avg_vm_1217, avg_vm_prev, avg_vm_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_vm, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'vm ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'vm ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'vm ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'vm ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );

create_plot(3, 2, 2, 'Meridional Wind Component, v', 'vm    [m/s]', avg_vm_les, z_les, nz_les, avg_vm_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_vm_1217, z_1217_zt, nz_1217_zt, ...
		avg_vm_prev, z_prev_zt, nz_prev_zt, avg_vm_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% u'w' (var "upwp")
[avg_upwp_les, avg_upwp_cgbe, avg_upwp_1217, avg_upwp_prev, avg_upwp_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_upwp, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'upwp ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'upwp ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'upwp ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'upwp ', nz_curr_zm, t1_curr_zm, t2_curr_zm, les_type, 1 );

if ( les_grph_upwp == 1 )
	create_plot(3, 2, 3, 'Covariance of u & w', 'upwp    [m^2/s^2]', avg_upwp_les, z_les, nz_les, avg_upwp_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_upwp_1217, z_1217_zt, nz_1217_zt, ...
	avg_upwp_prev, z_prev_zt, nz_prev_zt, avg_upwp_curr, z_curr_zt, nz_curr_zt, 0);
elseif ( les_grph_upwp == 2 )
	create_plot(3, 2, 3, 'Covariance of u & w', 'upwp    [m^2/s^2]', avg_upwp_les, z_les2, nz_les2, avg_upwp_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_upwp_1217, z_1217_zt, nz_1217_zt, ...
       	avg_upwp_prev, z_prev_zt, nz_prev_zt, avg_upwp_curr, z_curr_zt, nz_curr_zt, 0);
end

%--------------------------------------------------------------------------

% v'w' (var "vpwp")
[avg_vpwp_les, avg_vpwp_cgbe, avg_vpwp_1217, avg_vpwp_prev, avg_vpwp_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_vpwp, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'vpwp ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'vpwp ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'vpwp ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'vpwp ', nz_curr_zm, t1_curr_zm, t2_curr_zm, les_type, 1 );

if ( les_grph_vpwp == 1 )
	create_plot(3, 2, 4, 'Covariance of v & w', 'vpwp    [m^2/s^2]', avg_vpwp_les, z_les, nz_les, avg_vpwp_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_vpwp_1217, z_1217_zt, nz_1217_zt, ...
       	avg_vpwp_prev, z_prev_zt, nz_prev_zt, avg_vpwp_curr, z_curr_zt, nz_curr_zt, 0);
elseif ( les_grph_vpwp == 2 )
	create_plot(3, 2, 4, 'Covariance of v & w', 'vpwp    [m^2/s^2]', avg_vpwp_les, z_les2, nz_les2, avg_vpwp_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_vpwp_1217, z_1217_zt, nz_1217_zt, ...
	avg_vpwp_prev, z_prev_zt, nz_prev_zt, avg_vpwp_curr, z_curr_zt, nz_curr_zt, 0);
end
  
%--------------------------------------------------------------------------

% u" (var "up2")
[avg_up2_les, avg_up2_cgbe, avg_up2_1217, avg_up2_prev, avg_up2_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_up2, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'up2 ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'up2 ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'up2 ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'up2 ', nz_curr_zm, t1_curr_zm, t2_curr_zm, les_type, 0 );

create_plot(3, 2, 5, 'Variance of u wind', 'up2    [m^2/s^2]', avg_up2_les, z_les, nz_les, avg_up2_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_up2_1217, z_1217_zt, nz_1217_zt, ...
		avg_up2_prev, z_prev_zt, nz_prev_zt, avg_up2_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% v" (var "vp2")
[avg_vp2_les, avg_vp2_cgbe, avg_vp2_1217, avg_vp2_prev, avg_vp2_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_vp2, nz_les, t1_les, t2_les, filename_cgbe_zm, listofparams_cgbe_zm, numvars_cgbe_zm, 'vp2 ', nz_cgbe_zm, t1_cgbe_zm, ...
		t2_cgbe_zm, filename_1217_zm, listofparams_1217_zm, numvars_1217_zm, 'vp2 ', nz_1217_zm, t1_1217_zm, t2_1217_zm, filename_prev_zm, listofparams_prev_zm, ...
		numvars_prev_zm, 'vp2 ', nz_prev_zm, t1_prev_zm, t2_prev_zm, filename_curr_zm, listofparams_curr_zm, numvars_curr_zm, ...
		'vp2 ', nz_curr_zm, t1_curr_zm, t2_curr_zm, les_type, 0 );

create_plot(3, 2, 6, 'Variance of v wind', 'vp2    [m^2/s^2]', avg_vp2_les, z_les, nz_les, avg_vp2_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_vp2_1217, z_1217_zt, nz_1217_zt, ...
		avg_vp2_prev, z_prev_zt, nz_prev_zt, avg_vp2_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% Print 3rd Page for Case (.eps document)
outputfilename = [ case_name, '_page3' ];
print_page(outputfilename, 1, 0, 1);

figure('Position',[ 0 0 fig_width fig_height ])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperPosition', [ 1.0 1.0 6.5 9.0 ])


% Rain water mixing ratio (var "rrainm")
[avg_rrainm_les, avg_rrainm_cgbe, avg_rrainm_1217, avg_rrainm_prev, avg_rrainm_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_rrainm, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'rrm ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'rrm ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'rrainm ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'rrainm ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );

create_plot(3, 2, 1, 'Rain Water Mixing Ratio, r_r', 'rrainm    [kg/kg]', avg_rrainm_les, z_les, nz_les, avg_rrainm_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_rrainm_1217, z_1217_zt, nz_1217_zt, ...
		avg_rrainm_prev, z_prev_zt, nz_prev_zt, avg_rrainm_curr, z_curr_zt, nz_curr_zt, 0);

%--------------------------------------------------------------------------

% Rain drop concentration (var "Nrm")
[avg_Nrm_les, avg_Nrm_cgbe, avg_Nrm_1217, avg_Nrm_prev, avg_Nrm_curr] = var_load( filename_les, listofparams_les, ...
		numvars_les, les_Nrm, nz_les, t1_les, t2_les, filename_cgbe_zt, listofparams_cgbe_zt, numvars_cgbe_zt, 'Nrm ', nz_cgbe_zt, t1_cgbe_zt, ...
		t2_cgbe_zt, filename_1217_zt, listofparams_1217_zt, numvars_1217_zt, 'Nrm ', nz_1217_zt, t1_1217_zt, t2_1217_zt, filename_prev_zt, listofparams_prev_zt, ...
		numvars_prev_zt, 'Nrm ', nz_prev_zt, t1_prev_zt, t2_prev_zt, filename_curr_zt, listofparams_curr_zt, numvars_curr_zt, ...
		'Nrm ', nz_curr_zt, t1_curr_zt, t2_curr_zt, les_type, 0 );

% Adjustment:  COAMPS LES outputs Nrm in num/cm^3.  This factor needs to be
%              multiplied by 10^6 in order to be converted to num/m^3.
if ( cmp_les == 1 )
   if ( strcmp( les_type, 'coamps' ) )
      avg_Nrm_les = (10^6).*avg_Nrm_les;
   end
end

create_plot(3, 2, 2, 'Rain Drop Concentration, N_r', 'Nrm    [num/m^3]', avg_Nrm_les, z_les, nz_les, avg_Nrm_cgbe, z_cgbe_zt, nz_cgbe_zt, avg_Nrm_1217, z_1217_zt, nz_1217_zt, ...
		avg_Nrm_prev, z_prev_zt, nz_prev_zt, avg_Nrm_curr, z_curr_zt, nz_curr_zt, 1);

%--------------------------------------------------------------------------

% Liquid water path (var lwp from sfc files)	
[avg_lwp_prev, avg_lwp_curr] = time_tendency_load( filename_prev_sfc, listofparams_prev_sfc, ...
		numvars_prev_sfc, 'lwp ', 1, max(size(timesteps_prev)), filename_curr_sfc, listofparams_curr_sfc, numvars_curr_sfc, ...
		'lwp ', 1, max(size(timesteps_curr)));

create_time_plot(3, 2, 3, 'Liquid Water Path', 'lwp    [kg/m^2]', ...
		avg_lwp_prev, timesteps_prev, avg_lwp_curr, timesteps_curr );

%--------------------------------------------------------------------------


% Print 4th Page for Case (.eps document)
outputfilename = [ case_name, '_page4' ];
print_page(outputfilename, 1, 0, 1);

% Statement
[ 'Case ', case_name, ' has successfully been output!' ]
[ '=====================================================================' ]

end

function [avg_les_values, avg_cgbe_values, avg_1217_values, avg_prev_values, avg_curr_values] = var_load( filename_les, listofparams_les,  ...
		numvars_les, les_var, nz_les, t1_les, t2_les, filename_cgbe, listofparams_cgbe, numvars_cgbe, cgbe_var, nz_cgbe, t1_cgbe, ...
		t2_cgbe, filename_1217, listofparams_1217, numvars_1217, dec17_var, nz_1217, t1_1217, t2_1217, filename_prev, listofparams_prev, ...
		numvars_prev, prev_var, nz_prev, t1_prev, t2_prev, filename_curr, listofparams_curr, numvars_curr, ...
		curr_var, nz_curr, t1_curr, t2_curr, les_type, var_type )		

	global dir_prev
	global dir_curr
	global dir_LES
	global dir_cgbe
	global dir_1217
	global cmp_les
	global cmp_cgbe
	global cmp_1217
	global cmp_prev
	global cmp_curr
	global les_grph_upwp 
	global les_grph_vpwp
	global nz_les2
	global t1_les2
	global t2_les2
	global numvars_les2
	global filename_les2
	global listofparams_les2
	global les_upwp
	global les_vpwp

	avg_les_values = 0;
	avg_cgbe_values = 0;
       	avg_1217_values = 0;
       	avg_prev_values = 0;
       	avg_curr_values = 0;
	
	%LES
	if ( cmp_les == 1 )

		%Calculate the LES variable length
		les_var_len = max(size(les_var));

		if ( var_type == 0)
			varfnd = 0;
			for i = 1:1:numvars_les
				if ( strcmp( listofparams_les(i,1:les_var_len), les_var ) )
					varnum = i;
					varfnd = 1;
				end
				if ( (i == numvars_les) & (varfnd == 0) )
					[ 'LES Variable: ', les_var, ' not found for case, being set to 0!' ]
					avg_les_values(1:nz_les) = 0.0;
				elseif ( varfnd == 1 )
					avg_les_values = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
					'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
					break
				end
			end
		%The var_type == 1 is exlusively for handling UPWP and VPWP
		elseif (var_type == 1)
			varfnd = 0;
   			for i = 1:1:numvars_les
				if ( strcmp( listofparams_les(i,1:les_var_len), les_var ) )
					varnum = i;
					varfnd = 1;
					avg_les_values = read_grads_hoc_endian([dir_LES, '/', filename_les], ...
					'ieee-be', nz_les, t1_les, t2_les, varnum, numvars_les);
					
					if (les_var == les_upwp)
						les_grph_upwp = 1;
					elseif (les_var == les_vpwp)
						les_grph_vpwp = 1;
					end

					break
				end
			end

   			if ( (varfnd == 0) & ( strcmp( les_type, 'coamps' ) ) )
      				for i = 1:1:numvars_les2
         				if ( strcmp( listofparams_les2(i,1:les_var_len), les_var ) )
            					varnum = i;
            					varfnd = 1;
         				end
         				if ( (i == numvars_les2) & (varfnd == 0) )
            					[ 'LES Variable: ', les_var, ' not found for case, being set to 0!' ]
            					avg_les_values(1:nz_les2) = 0.0;
            					
						if (les_var == les_upwp)
							les_grph_upwp = 2;
						elseif (les_var == les_vpwp)
							les_grph_vpwp = 2;
						end

         				elseif ( varfnd == 1 )
            					avg_les_values = read_grads_hoc_endian([dir_LES, '/', filename_les2], ...
                          			'ieee-be', nz_les2, t1_les2, t2_les2, varnum, numvars_les2);
            					
						if (les_var == les_upwp)
							les_grph_upwp = 2;
						elseif (les_var == les_vpwp)
							les_grph_vpwp = 2;
						end

            					break
         				end
      				end
   			elseif ( (varfnd == 0) & ( strcmp( les_type, 'rams' ) ) )
      				[ 'LES Variable: ', les_var, ' not found for case, being set to 0!' ]
      				avg_les_values(1:nz_les) = 0.0;
      				
				if (les_var == les_upwp)
					les_grph_upwp = 1;
				elseif (les_var == les_vpwp)
					les_grph_vpwp = 1;
				end
   			end
		end
	end	

	% HOC -- Golaz "best ever"
	if ( cmp_cgbe == 1 )
		
		cgbe_var_len = max(size(cgbe_var));

		varfnd = 0;
		for i = 1:1:numvars_cgbe
			if ( strcmp( listofparams_cgbe(i,1:cgbe_var_len), cgbe_var ) )
				varnum = i;
				varfnd = 1;
			end
			if ( (i == numvars_cgbe) & (varfnd == 0) )
				[ 'Chris Golaz Best Ever Variable: ', cgbe_var, ' not found for case, being set to 0!' ]
				avg_cgbe_values(1:nz_cgbe) = 0.0;
			elseif ( varfnd == 1 )
				avg_cgbe_values = read_grads_hoc_endian([dir_cgbe, '/', filename_cgbe], ...
				'ieee-le', nz_cgbe, t1_cgbe, t2_cgbe, varnum, numvars_cgbe);
				break
			end
		end
	end

	% HOC -- December 17, 2005
	if ( cmp_1217 == 1 )
		
		dec17_var_len = max(size(dec17_var));

		varfnd = 0;
		for i = 1:1:numvars_1217
			if ( strcmp( listofparams_1217(i,1:dec17_var_len), dec17_var ) )
				varnum = i;
				varfnd = 1;
			end
			if ( (i == numvars_1217) & (varfnd == 0) )
				[ 'Dec. 17 Variable: ', dec17_var, ' not found for case, being set to 0!' ]
				avg_1217_values(1:nz_1217) = 0.0;
			elseif ( varfnd == 1 )
				avg_1217_values = read_grads_hoc_endian([dir_1217, '/', filename_1217], ...
				'ieee-le', nz_1217, t1_1217, t2_1217, varnum, numvars_1217);
				break
			end
		end
	end

	% HOC -- previous
	if ( cmp_prev == 1 )
		
		prev_var_len = max(size(prev_var));
		
		varfnd = 0;
		for i = 1:1:numvars_prev
			if ( strcmp( listofparams_prev(i,1:prev_var_len), prev_var ) )
				varnum = i;
				varfnd = 1;
			end
			if ( (i == numvars_prev) & (varfnd == 0) )
				[ filename_prev, ' Variable: ', prev_var, ' not found for case, being set to 0!' ]
				avg_prev_values(1:nz_prev) = 0.0;
			elseif ( varfnd == 1 )
				avg_prev_values = read_grads_hoc_endian([dir_prev, '/', filename_prev], ...
				'ieee-le', nz_prev, t1_prev, t2_prev, varnum, numvars_prev);
				break
			end
		end
	end

	% HOC -- current
	if ( cmp_curr == 1 )
		
		curr_var_len = max(size(curr_var));
		
		varfnd = 0;
		for i = 1:1:numvars_curr
			if ( strcmp( listofparams_curr(i,1:curr_var_len), curr_var ) )
				varnum = i;
				varfnd = 1;
			end
			if ( (i == numvars_curr) & (varfnd == 0) )
				[ filename_curr, ' Variable: ', curr_var, ' not found for case, being set to 0!' ]
				avg_curr_values(1:nz_curr) = 0.0;
			elseif ( varfnd == 1 )
				avg_curr_values = read_grads_hoc_endian([dir_curr, '/', filename_curr], ...
				'ieee-le', nz_curr, t1_curr, t2_curr, varnum, numvars_curr);
				break
			end
		end
	end
%End var_load
end

%Use only for loading time tendencies from SFC files.
function [avg_prev_values, avg_curr_values] = time_tendency_load( filename_prev, listofparams_prev, ...
		numvars_prev, prev_var, t1_prev, t2_prev, filename_curr, listofparams_curr, numvars_curr, ...
		curr_var, t1_curr, t2_curr )		

	global dir_prev
	global dir_curr
	global cmp_prev
	global cmp_curr

       	avg_prev_values = 0;
       	avg_curr_values = 0;
	

	% HOC -- previous
	if ( cmp_prev == 1 )
		
		prev_var_len = max(size(prev_var));
		
		varfnd = 0;
		for i = 1:1:numvars_prev
			if ( strcmp( listofparams_prev(i,1:prev_var_len), prev_var ) )
				varnum = i;
				varfnd = 1;
			end
			if ( (i == numvars_prev) & (varfnd == 0) )
				[ filename_prev, ' Variable: ', prev_var, ' not found for case, being set to 0!' ]
				avg_prev_values(1:nz_prev) = 0.0;
			elseif ( varfnd == 1 )
				avg_prev_values = read_grads_hoc_sfc_endian([dir_prev, '/', filename_prev], ...
				'ieee-le', 1, t1_prev, t2_prev, varnum, numvars_prev);
				break
			end
		end
	end

	% HOC -- current
	if ( cmp_curr == 1 )
		
		curr_var_len = max(size(curr_var));
		
		varfnd = 0;
		for i = 1:1:numvars_curr
			if ( strcmp( listofparams_curr(i,1:curr_var_len), curr_var ) )
				varnum = i;
				varfnd = 1;
			end
			if ( (i == numvars_curr) & (varfnd == 0) )
				[ filename_curr, ' Variable: ', curr_var, ' not found for case, being set to 0!' ]
				avg_curr_values(1:nz_curr) = 0.0;
			elseif ( varfnd == 1 )
				avg_curr_values = read_grads_hoc_sfc_endian([dir_curr, '/', filename_curr], ...
				'ieee-le', 1, t1_curr, t2_curr, varnum, numvars_curr);
				break
			end
		end
	end
%End time_tendency_load
end

function create_plot( plot_x, plot_y, plot_z, plot_title, plot_units, les_data, z_les, nz_les, cgbe_data, z_cgbe, nz_cgbe, ...
			dec12_data, z_1217, nz_1217, prev_data, z_prev, nz_prev, curr_data, z_curr, nz_curr, add_legend )
	
	
	global dir_prev
	global dir_curr
	global cmp_les
	global cmp_cgbe
	global cmp_1217
	global cmp_prev
	global cmp_curr
	global les_color
	global les_width
	global cgbe_color
	global cgbe_width
	global dec12_color
	global dec12_width
	global prev_color
	global prev_width
	global curr_color
	global curr_width
	global graphtop
	global graphbase
	global equiv_space
	global percentage	
			
	subplot(plot_x, plot_y, plot_z)
	i = 0;
	clear h;
	
	if ( cmp_les == 1 )
   		i = i + 1;
   		
		h(i) = plot( les_data, z_les, '-', 'Color', les_color, 'LineWidth', les_width );
		nz_data = nz_les;
		z_data = z_les;
		
		legend_text(i,1:15) = '\fontsize{6}LES';
		hold on

   		% Find vertical level index right below top of graph.
   		for k = 2:1:nz_data
      			if ( z_data(k) > graphtop )
         			graphtopidx = k - 1;
        			break
      			elseif ( k == nz_data )
         			graphtopidx = nz_data;
      			end
   		end

   		minval(i) = min(les_data(1:graphtopidx));
   		maxval(i) = max(les_data(1:graphtopidx));
	end

	if ( cmp_cgbe == 1 )
   		i = i + 1;
   		h(i) = plot( cgbe_data, z_cgbe, '-', 'Color', cgbe_color, 'LineWidth', cgbe_width );
		legend_text(i,1:27) = '\fontsize{6}HOC "best-ever"';
   		hold on
   		% Find vertical level index right below top of graph.
   		for k = 2:1:nz_cgbe
      			if ( z_cgbe(k) > graphtop )
         			graphtopidx = k - 1;
         			break
      			elseif ( k == nz_cgbe )
         			graphtopidx = nz_cgbe;
      			end
   		end
		minval(i) = min(cgbe_data(1:graphtopidx));
   		maxval(i) = max(cgbe_data(1:graphtopidx));
	end

	if ( cmp_1217 == 1 )
   		i = i + 1;
   		h(i) = plot( dec12_data, z_1217, '-.', 'Color', dec12_color, 'LineWidth', dec12_width );
		legend_text(i,1:26) = '\fontsize{6}HOC 12/17/2005';
   		hold on
   		% Find vertical level index right below top of graph.
   		for k = 2:1:nz_1217
      			if ( z_1217(k) > graphtop )
         			graphtopidx = k - 1;
         			break
      			elseif ( k == nz_1217 )
         			graphtopidx = nz_1217;
      			end
   		end
   		minval(i) = min(dec12_data(1:graphtopidx));
   		maxval(i) = max(dec12_data(1:graphtopidx));
	end
	
	if ( cmp_prev == 1 )
   		i = i + 1;
   		h(i) = plot( prev_data, z_prev, '--', 'Color', prev_color, 'LineWidth', prev_width );
   		sim1_title = strcat('\fontsize{6}', dir_prev);
   		sim1_title = regexprep(sim1_title, '_', ' ');
   		legend_text(i,1:length(sim1_title)) = sim1_title;
   		hold on
   		% Find vertical level index right below top of graph.
   		for k = 2:1:nz_prev
      			if ( z_prev(k) > graphtop )
         			graphtopidx = k - 1;
         			break
      			elseif ( k == nz_prev )
         			graphtopidx = nz_prev;
      			end
   		end
   		minval(i) = min(prev_data(1:graphtopidx));
   		maxval(i) = max(prev_data(1:graphtopidx));
	end

	if ( cmp_curr == 1 )
   		i = i + 1;
   		h(i) = plot( curr_data, z_curr, '-', 'Color', curr_color, 'LineWidth', curr_width );
		sim2_title = strcat('\fontsize{6}', dir_curr);
   		sim2_title = regexprep(sim2_title, '_', ' ');
   		legend_text(i,1:length(sim2_title)) = sim2_title;
   		hold on
   		% Find vertical level index right below top of graph.
   		for k = 2:1:nz_curr
      			if ( z_curr(k) > graphtop )
         			graphtopidx = k - 1;
         			break
      			elseif ( k == nz_curr )
         			graphtopidx = nz_curr;
      			end
   		end
   		minval(i) = min(curr_data(1:graphtopidx));
   		maxval(i) = max(curr_data(1:graphtopidx));
	end

	hold off

	if ( add_legend == 1 )
		legend( h, legend_text, 'Location', 'NorthEast' )
	end	
	
	% Axis labels and graph title.
	xlabel(plot_units)
	ylabel('Height    [m]')
	title(plot_title)

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
	zmax = graphtop;
	axis([ xmin xmax zmin zmax ])
end

%Do not use for purposes other than plotting time tendencies (using the aforementioned time_tendency_load
function create_time_plot( plot_x, plot_y, plot_z, plot_title, plot_units, ...
			prev_data, t_prev, curr_data, t_curr )
	
	
	global dir_prev
	global dir_curr
	global cmp_les
	global cmp_prev
	global cmp_curr
	global les_color
	global les_width
	global prev_color
	global prev_width
	global curr_color
	global curr_width
	global graphtop
	global graphbase
	global equiv_space
	global percentage	
			
	subplot(plot_x, plot_y, plot_z)
	i = 0;
	clear h;
		
	if ( cmp_prev == 1 )
   		i = i + 1;
		h(i) = plot( t_prev, prev_data, '--', 'Color', prev_color, 'LineWidth', prev_width );
   		sim1_title = strcat('\fontsize{6}', dir_prev);
   		sim1_title = regexprep(sim1_title, '_', ' ');
   		legend_text(i,1:length(sim1_title)) = sim1_title;
   		hold on

		minval(i) = min(prev_data);
   		maxval(i) = max(prev_data);

		mintime(i) = min(t_prev);
		maxtime(i) = max(t_prev);
	end

	if ( cmp_curr == 1 )
   		i = i + 1;
   		h(i) = plot( t_curr, curr_data, '-', 'Color', curr_color, 'LineWidth', curr_width );
		sim2_title = strcat('\fontsize{6}', dir_curr);
   		sim2_title = regexprep(sim2_title, '_', ' ');
   		legend_text(i,1:length(sim2_title)) = sim2_title;
   		hold on

		minval(i) = min(curr_data);
   		maxval(i) = max(curr_data);

		mintime(i) = min(t_curr);
		maxtime(i) = max(t_curr);
	end

	hold off
	
	% Axis labels and graph title.
	xlabel('Time    [s]')
	ylabel(plot_units)
	title(plot_title)

	% Extent of graph.
	xmin = min(mintime);
	xmax = max(maxtime);
	ymin = min(minval);
	ymax = max(maxval);
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

	axis([ xmin xmax ymin ymax ])
end


function print_page( file_name, output_eps, output_jpeg, output_ps )
	if (output_eps == 1)
		eps_file_name = [ 'output/eps/', file_name, '.eps' ];
		print( '-depsc2', eps_file_name )
	end

	if (output_jpeg == 1)
		jpg_file_name = [ 'output/jpeg/', file_name, '.jpg' ];
		print( '-djpeg', jpg_file_name )
	end
	
	if (output_ps == 1)
		ps_file_name = [ 'output/ps/', file_name, '.ps' ];
		print( '-dpsc', '-append', ps_file_name )
	end

	% Close the current figure (so it doesn't produce a new one for each page
	% and for each case).
	closereq
end
