%$Id: find_grads_differences.m,v 1.4 2008-06-27 21:47:48 vlarson Exp $

function [] = find_grads_differences( ctl_file, t1, t2, tol )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find_grads_differences( string, integer, integer, integer )
% Example: find_grads_differences( 'nov11_altocu.ctl', 1, 240, 4 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, December 2006-March 2007
% Modified by Joshua Fasching, June 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters: 
%       ctl_file - .ctl file for the data being compared
%       t1      - Lower time bound of data
%       t2      - Upper time bound of data
%       tol     - Digits of Precision ( Higher # means finer tolerance )
%
% Input files: Three GrADS .ctl files in separate directories
% Output parameters: none
% Output files: none
%
% Requires: three Michael Falk utility scripts:
%           header_read.m
%           read_grads_hoc_endian.m
%           get_profile.m
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes three separate model runs with identically-named
% files in separate directories (following Michael Falk's directory
% structure and naming convention) and compares them, telling you
% which fields differ among the three files.
%
% It accomplishes this by taking a time average over a given period for
% each file and comparing a vertical profile of each variable over that
% time period.  If the two profiles are not within tolerance, the script will
% print a message to the screen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The three directories containing three different ctlfile file sets.
path1  = ['/home/faschinj/hoc_v2.2_tuner/standalone/'];
path2  = ['/home/faschinj/previous/hoc_v2.2_tuner/standalone/'];
path3  = ['/home/faschinj/previous/hoc_v2.2_tuner/standalone/'];

% Read information from .ctl files
[filename1,nz1,z1,ntimesteps1,numvars1,list_vars1] = header_read([path1,ctl_file]);
[filename2,nz2,z2,ntimesteps2,numvars2,list_vars2] = header_read([path2,ctl_file]);
[filename3,nz3,z3,ntimesteps3,numvars3,list_vars3] = header_read([path3,ctl_file]);

% Loop over all variables
for var=1:numvars3
    var1 = deblank(list_vars3(var,:));
    
    % Get profile for variable averaged over t1 to t2
    sim_profile1 = get_profile(path1,filename1, ...
        numvars1,list_vars1,nz1,t1,t2,var1,1);
    sim_profile2 = get_profile(path2,filename2, ...
        numvars2,list_vars2,nz2,t1,t2,var1,1);
    sim_profile3 = get_profile(path3,filename3, ...
        numvars3,list_vars3,nz3,t1,t2,var1,1);

    % Determine the mean difference between profiles.
    % diff1, diff2, and diff3 are vectors, not scalars.
    diff1 = abs( sim_profile2 - sim_profile1 );
    diff2 = abs( sim_profile3 - sim_profile2 );
    diff3 = abs( sim_profile3 - sim_profile1 );
    
    % Determine the tolerance between profiles.
    % This tolerance will be undesirably reduced for domains with very high tops.
    tolerance1 = mean( abs( [sim_profile2 ; sim_profile1] ) * 1.0*10^(-tol) );
    tolerance2 = mean( abs( [sim_profile3 ; sim_profile2] ) * 1.0*10^(-tol) );
    tolerance3 = mean( abs( [sim_profile3 ; sim_profile1] ) * 1.0*10^(-tol) );
    
% Print results, if differences exist
    if ( any( diff1 > tolerance1 ) )
        disp ( ['File 1 and File 2 disagree in field ', var1] );
    end
    if ( any( diff2 > tolerance2 ) )
        disp ( ['File 2 and File 3 disagree in field ', var1] );
    end
    if ( any( diff3 > tolerance3 ) )
        disp ( ['File 1 and File 3 disagree in field ', var1] );
    end
end
