function [] = find_grads_differences();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find_grads_differences()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, December 2006-March 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters: none
% Input files: Three GrADS .ctl files in separate directories
% Output parameters: none
% Output files: none
%
% Requires: three Michael Falk utility scripts:
%           header_read.m
%           read_grads_hoc_endian.m
%           get_profile.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes three separate model runs with identically-named
% files in separate directories (following Michael Falk's directory
% structure and naming convention) and compares them, telling you
% which fields differ among the three files.
%
% It accomplishes this by taking a time average over a given period for
% each file and comparing a vertical profile of each variable over that
% time period.  If the two profiles are not identical, the script will
% print a message to the screen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The three directories containing three different gabls_zt file sets.
path1  = ['/home/mjfalk/hoc_working/standalone/gabls63/'];
path2  = ['/home/mjfalk/hoc_working/standalone/gabls63_test/'];
path3  = ['/home/mjfalk/hoc_working/standalone/gabls63_test2/'];
smfile = 'gabls_zt.ctl';

[filename1,nz1,z1,ntimesteps1,numvars1,list_vars1] = header_read([path1,smfile]);
[filename2,nz2,z2,ntimesteps2,numvars2,list_vars2] = header_read([path2,smfile]);
[filename3,nz3,z3,ntimesteps3,numvars3,list_vars3] = header_read([path3,smfile]);

% The times over which to average for the comparison
t1 = 1;
t2 = 140;

% Loop over all variables
for var=1:numvars3
    var1 = deblank(list_vars3(var,:));
    sim_profile1 = get_profile(path1,filename1, ...
        numvars1,list_vars1,nz1,t1,t2,var1,1);
    sim_profile2 = get_profile(path2,filename2, ...
        numvars2,list_vars2,nz2,t1,t2,var1,1);
    sim_profile3 = get_profile(path3,filename3, ...
        numvars3,list_vars3,nz3,t1,t2,var1,1);

    diff1 = (sim_profile2 - sim_profile1);
    diff2 = (sim_profile3 - sim_profile2);
    diff3 = (sim_profile3 - sim_profile1);

% Print results, if differences exist
    if (sum(diff1) ~= 0)
        disp (['File 1 and File 2 disagree in field ',var1]);
    end
    if (sum(diff2) ~= 0)
        disp (['File 2 and File 3 disagree in field ',var1]);
    end
    if (sum(diff3) ~= 0)
        disp (['File 1 and File 3 disagree in field ',var1]);
    end
end
