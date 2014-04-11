function starter_each_case

total = 9;

for case_num = 1:1:total

   % Cases
   %
   % Total Number of Cases
   % Note:  In this script, it is total number of cases at a time, which in
   % this case, would be one.
   c_total = 1;
   %
   % Individual Case Names, Weights, Time Periods.
   % Note:  The sum of all case weights should be 1.
   %
   if ( case_num == 1 )
      % Case #1
      % ARM case.
      c_name(1,1:3) = 'arm';
      c_name(1,4:15) = blanks(12);
      weight_case(1) = 1.0;
      t1_case(1) = 480;
      t2_case(1) = 540;
      [ '========== ARM case ==========' ]
   elseif ( case_num == 2 )
      % Case #2
      % ATEX case.
      c_name(1,1:4) = 'atex';
      c_name(1,5:15) = blanks(11);
      weight_case(1) = 1.0;
      t1_case(1) = 420;
      t2_case(1) = 480;
      [ '========== ATEX case =========' ]
   elseif ( case_num == 3 )
      % Case #3
      % BOMEX case.
      c_name(1,1:5) = 'bomex';
      c_name(1,6:15) = blanks(10);
      weight_case(1) = 1.0;
      t1_case(1) = 180;
      t2_case(1) = 360;
      [ '========= BOMEX case =========' ]
   elseif ( case_num == 4 )
      % Case #4
      % DYCOMS2 RF01 case.
      c_name(1,1:12) = 'dycoms2_rf01';
      c_name(1,13:15) = blanks(3);
      weight_cacdse(1) = 1.0;
      t1_case(1) = 180;
      t2_case(1) = 240;
      [ '====== DYCOMS2 RF01 case =====' ]
   elseif ( case_num == 5 )
      % Case #5
      % DYCOMS2 RF02 DO case.
      c_name(1,1:15) = 'dycoms2_rf02_do';
      weight_case(1) = 1.0;
      t1_case(1) = 180;
      t2_case(1) = 360;
      [ '==== DYCOMS2 RF02 DO case ====' ]
   elseif ( case_num == 6 )
      % Case #6
      % DYCOMS2 RF02 DS case.
      c_name(1,1:15) = 'dycoms2_rf02_ds';
      weight_case(1) = 1.0;
      t1_case(1) = 180;
      t2_case(1) = 360;
      [ '==== DYCOMS2 RF02 DS case ====' ]
   elseif ( case_num == 7 )
      % Case #7
      % DYCOMS2 RF02 ND case.
      c_name(1,1:15) = 'dycoms2_rf02_nd';
      weight_case(1) = 1.0;
      t1_case(1) = 180;
      t2_case(1) = 360;
      [ '==== DYCOMS2 RF02 DO case ====' ]
   elseif ( case_num == 8 )
      % Case #8
      % FIRE case.
      c_name(1,1:4) = 'fire';
      c_name(1,5:15) = blanks(11);
      weight_case(1) = 1.0;
      t1_case(1) = 60;
      t2_case(1) = 120;
      [ '========== FIRE case =========' ]
%   elseif ( case_num == 9 )
%      % Case #9
%      % June 25 Altocu case.
%      c_name(1,1:12) = 'jun25_altocu';
%      c_name(1,13:15) = blanks(3);
%      weight_case(1) = 1.0;
%      t1_case(1) = 60;
%      t2_case(1) = 120;
%      [ '===== June 25 Altocu case ====' ]
%   elseif ( case_num == 9 )
%      % Case #10
%      % Nov 11 Altocu case.
%      c_name(1,1:12) = 'nov11_altocu';
%      c_name(1,13:15) = blanks(3);
%      weight_case(1) = 1.0;
%      t1_case(1) = 60;
%      t2_case(1) = 120;
%      [ '===== Nov 11 Altocu case =====' ]
   elseif ( case_num == 9 )
      % Case #9
      % RICO case.
      c_name(1,1:4) = 'rico';
      c_name(1,5:15) = blanks(11);
      weight_case(1) = 1.0;
      t1_case(1) = 60;
      t2_case(1) = 120;
      [ '========== RICO case =========' ]
   end
   
   
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

   % linitialize_sigma is set to 1 (for "true") for the reading of the 
   % initial file (the file that other files are being measured against in 
   % an at least three-way comparison with the LES).
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

end
