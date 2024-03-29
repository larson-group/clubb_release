function starter

% Cases
%
% Total Number of Cases
c_total = 4;
%
% Individual Case Names, Weights, Time Periods.
% Note:  The sum of all case weights should be 1.
%
% Case #1
% ARM case.
c_name(1,1:3) = 'arm';
weight_case(1) = 0.10;
t1_case(1) = 480;
t2_case(1) = 540;
%
% Case #2
% BOMEX case.
c_name(2,1:5) = 'bomex';
weight_case(2) = 0.25;
t1_case(2) = 180;
t2_case(2) = 360;
%
% Case #3
% DYCOMS2 RF01 case.
c_name(3,1:12) = 'dycoms2_rf01';
weight_case(3) = 0.25;
t1_case(3) = 180;
t2_case(3) = 240;
%
% Case #4
% DYCOMS2 RF02 DO case.
c_name(4,1:15) = 'dycoms2_rf02_do';
weight_case(4) = 0.40;
t1_case(4) = 180;
t2_case(4) = 360;


% Variables
%
% Total Number of Variables
v_total = 2;
%
% Individual Variable Names, Weights
% Note:  The sum of all variable weights should be 1.
%
% Variable #1
%les_varname(1) = 'thlm';
%hoc_varname(1) = 'thlm';
%weight_var(1) = 1.0;
% Variable #1
les_varname(1,1:3) = 'qcm';
hoc_varname(1,1:3) = 'rcm';
weight_var(1) = 0.5;
%
% Variable #2
les_varname(2,1:2) = 'cf';
hoc_varname(2,1:2) = 'cf';
weight_var(2) = 0.5;


% Initial file.
hoc_dir = [ 'HOC_init_data' ];
les_dir = [ 'LES_data' ];

% linitialize_sigma is set to 1 (for "true") for the reading of the initial
% file (the file that other files are being measured against in an at least
% three-way comparison with the LES).
linitialize_sigma = 1;

% invsigma2 is set to 0 here just to initialize the array in the proper
% shape, since it is being sent both in and out of the ensuing function.
invsigma2(1:c_total,1:v_total) = 0

[ 'Data from the comparison of the HOC initial file (control file) to the LES' ]

% Call min_les_hoc_diff for the initial file (control file).
[ cost_function, err_terms, err_sums, invsigma2 ] = ...
   min_les_hoc_diff( hoc_dir, les_dir, c_total, c_name, weight_case, ...
                     t1_case, t2_case, v_total, les_varname, hoc_varname, ...
                     weight_var, linitialize_sigma, invsigma2 )

                 
% Comparison file
hoc_dir = [ 'HOC_comp_data' ];
les_dir = [ 'LES_data' ];

% linitialize_sigma is set to 0 (for "false") for the reading of any
% comparison files.
linitialize_sigma = 0;

[ 'Data from the comparison of the HOC comparison file to the LES' ]

% Call min_les_hoc_diff for the comparision file.
[ cost_function, err_terms, err_sums, invsigma2 ] = ...
   min_les_hoc_diff( hoc_dir, les_dir, c_total, c_name, weight_case, ...
                     t1_case, t2_case, v_total, les_varname, hoc_varname, ...
                     weight_var, linitialize_sigma, invsigma2 )
