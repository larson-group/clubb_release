% $Id: compare_plots_cases_driver.m,v 1.1 2007-02-28 17:54:24 dschanen Exp $
%function compare_plots_cases_driver( ...
%                                     arm_success, ...
%                                     atex_success, ...
%                                     bomex_success, ...
%                                     cobra_success, ...
%                                     rf01_success, ...
%                                     rf02_do_success, ...
%                                     rf02_ds_success, ...
%                                     rf02_nd_success, ...
%                                     rf02_so_success, ...
%                                     fire_success, ...
%                                     gabls_success, ...
%                                     jun25_success, ...
%                                     mpace_success, ...
%                                     nov11_success, ...
%                                     rico_success, ...
%                                     wangara_success ...
%                                                         )

function compare_plots_cases_driver

% Hardwire for successful runs.
arm_success = 1;
atex_success = 1;
bomex_success = 1;
cobra_success = 1;
rf01_success = 1;
rf02_do_success = 1;
rf02_ds_success = 1;
rf02_nd_success = 1;
rf02_so_success = 1;
fire_success = 1;
gabls_success = 1;
jun25_success = 1;
mpace_success = 1;
nov11_success = 1;
rico_success = 1;
wangara_success = 1;

% Function call variable description:
%
% 1st entry:   Case name.
% 2nd entry:   First time of variable comparison (in minutes since the
%              beginning of the run).
% 3rd entry:   Last time of variable comparison (in minutes since the
%              beginning of the run).
% 4th entry:   The height of the top of the output graphs for each case.
% 5th entry:   The directory path that contains the LES files for
%              comparison. (omit the '/' at the very end.)
% 6th entry:   The directory path that contains the HOC "Chris Golaz
%              'best-ever'" files for comparison. (omit the '/' at the very
%              end.)
% 7th entry:   The directory path that contains the HOC December 17, 2005
%              files for comparison. (omit the '/' at the very end.)
% 8th entry:   The directory path that contains the HOC previous files
%              (from the prior CVS commitment) for comparison. (omit the
%              '/' at the very end.)
% 9th entry:   The directory path that contains the HOC current files (from 
%              the latest CVS commitment) for comparison. (omit the '/' at
%              the very end.)
% 10th entry:  LES file for comparison:  1 = true; 0 = false.
% 11th entry:  LES being using for comparison (supports coamps and rams).
% 12th entry:  Chris Golaz "best-ever" HOC file for comparison:  1 = true;
%                                                                0 = false.
% 13th entry:  December 17, 2005 HOC file for comparison:  1 = true;
%                                                         0 = false.
% 14th entry:  Previous CVS HOC file for comparison:  1 = true; 
%                                                    0 = false.
% 15th entry:  Current HOC file for comparison:  1 = true; 0 = false.

for i = 1:1:16
   if ( i == 1 )
      if ( arm_success == 1 )
         compare_plots_cases( 'arm', 480, 540, 3500.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 1, 1, 1, 1 )
      else
         'ARM case crashed'
      end
   elseif ( i == 2 )
      if ( atex_success == 1 )
         compare_plots_cases( 'atex', 420, 480, 2500.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 1, 1, 1, 1 )
      else
         'ATEX case crashed'
      end
   elseif ( i == 3 )
      if ( bomex_success == 1 )
         compare_plots_cases( 'bomex', 180, 360, 2500.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 1, 1, 1, 1 )
      else
         'BOMEX case crashed'
      end
   elseif ( i == 4 )
      if ( cobra_success == 1 )
         compare_plots_cases( 'cobra', 240, 300, 3500.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              0, '', 0, 0, 1, 1 )
      else
         'COBRA case crashed'
      end
   elseif ( i == 5 )
      if ( rf01_success == 1 )
         compare_plots_cases( 'dycoms2_rf01', 180, 240, 1200.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 1, 1, 1, 1 )
      else
         'DYCOMS2 RF01 crashed'
      end
   elseif ( i == 6 )
      if ( rf02_do_success == 1 )
         compare_plots_cases( 'dycoms2_rf02_do', 300, 360, 1200.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 1, 1, 1, 1 )
      else
         'DYCOMS2 RF02 DO crashed'
      end
   elseif ( i == 7 )
      if ( rf02_ds_success == 1 )
         compare_plots_cases( 'dycoms2_rf02_ds', 300, 360, 1200.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 1, 1, 1, 1 )
      else
         'DYCOMS2 RF02 DS crashed'
      end
   elseif ( i == 8 )
      if ( rf02_nd_success == 1 )
         compare_plots_cases( 'dycoms2_rf02_nd', 300, 360, 1200.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 1, 1, 1, 1 )
      else
         'DYCOMS2 RF02 ND crashed'
      end
   elseif ( i == 9 )
      if ( rf02_so_success == 1 )
         compare_plots_cases( 'dycoms2_rf02_so', 300, 360, 1200.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              0, '', 1, 1, 1, 1 )
      else
         'DYCOMS2 RF02 SO crashed'
      end
   elseif ( i == 10 )
      if ( fire_success == 1 )
         compare_plots_cases( 'fire', 60, 120, 1000.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 1, 1, 1, 1 )
      else
         'FIRE crashed'
      end
   elseif ( i == 11 )
      if ( gabls_success == 1 )
         compare_plots_cases( 'gabls2', 2100, 2160, 2500.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 0, 0, 1, 1 )
      else
         'GABLS crashed'
      end
   elseif ( i == 12 )
      if ( jun25_success == 1 )
         compare_plots_cases( 'jun25_altocu', 60, 120, 2500.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 0, 0, 1, 1 )
      else
         'June 25 Altocumulus crashed'
      end
   elseif ( i == 13 )
      if ( mpace_success == 1 )
         compare_plots_cases( 'mpace', 540, 720, 2750.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              0, '', 0, 0, 1, 1 )
      else
         'MPACE crashed'
      end
   elseif ( i == 14 )
      if ( nov11_success == 1 )
         compare_plots_cases( 'nov11_altocu', 60, 120, 2000.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 1, 1, 1, 1 )
      else
         'November 11 Altocumulus crashed'
      end
   elseif ( i == 15 )
      if ( rico_success == 1 )
         compare_plots_cases( 'rico', 1260, 1440, 3500.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'coamps', 0, 0, 1, 1 )
      else
         'RICO crashed'
      end
   elseif ( i == 16 )
      if ( wangara_success == 1 )
         compare_plots_cases( 'wangara', 180, 240, 1900.0, ...
                              'LES_files', 'Chris_Golaz_best_ever', ...
                              'HOC_20051217', 'HOC_previous', 'HOC_current', ...
                              1, 'rams', 1, 1, 1, 1 )
      else
         'Wangara crashed'
      end
   end
   %keyboard
   %pause(1)
end
