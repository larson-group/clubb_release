function stat_profile = get_profile(les_path,datfile,numvars,list_vars,nz,t1,t2, ...
                               variable,sim)

%GET_PROFILE   Accesses stats files, calls read_grads_hoc_endian to get profiles.
%
%   Original program: plot_les_scm.m
%   Written by Michael J. Falk (mjfalk@uwm.edu)
%   Original version written on 1-4 August 2005
%
%   Modified by Adam J. Smith (ajsmith4@uwm.edu)
%   Completed on 25 October 2005
%
%   This program takes a variable list from the parent program, and loops through it
%   until we find an entry matching the name of the variable to be observed.  The
%   program then generates and executes a call statement to read_grads_hoc_endian.m,
%   which returns the desired variable profile.  get_profile then passes the profile
%   back to the parent program for observation and manipulation.
%
%   NOTE: For the full original version of this program, please see
%         ../original_files/plot_les_scm.ORIGINAL.20050805.m.
%
%   ----------------------------------------------------------------------------------
%
%   INPUTS
%   ------
%   les_path  -> File path to the GrADS statistics files for our simulation
%   datfile   -> Name of the simulation's GrADS data file (usually stats_sm.dat)
%   numvars   -> Total number of variables in the GrADS data file
%   list_vars -> Array of all variables contained in the GrADS data file
%   nz        -> Total number of altitudes in the GrADS data file
%   t1        -> Initial time to examine in simulation
%   t2        -> Final time to examine in simulation; if t1 and t2 are different,
%                we obtain an averaged profile over the period from t1 to t2.
%   variable  -> The specific variable we wish to examine
%
%   OUTPUTS
%   -------
%   profile -> GrADS profile of the variable we wish to examine

% END OF DOCUMENTATION
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

% EXAMPLE CALL:     get_profile('c:\adam\nov11_les\nov11.20050707.control\','stats_sm.ctl',60,120,'cf')

% Combine les_path and datfile strings
filename = [les_path,datfile];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data field into MATLAB from .dat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This section loops once for each stats variable in the simulation.  If the
% variable name for an iteration matches the variable we are looking for, we
% generate and execute a call statement to read_grads_hoc_endian to get the
% corresponding profile.

for i=1:numvars
    varname=list_vars(i,:);  % Obtain the i'th variable in the varlist
    varname=deblank(varname);          % Removes excess spaces used in array
    
    % If varname is the variable we are looking for, call read_grads_hoc_endian
    % to get the profile
    if (strcmp(varname,variable) & sim == 1)
        stringtoeval = ['stat_profile = read_grads_hoc_endian(filename,''ieee-le'',nz,t1,t2,i,numvars);'];    
        eval(stringtoeval);
    elseif (strcmp(varname,variable) & sim == 0)
        stringtoeval = ['stat_profile = read_grads_hoc_endian(filename,''ieee-be'',nz,t1,t2,i,numvars);'];
        eval(stringtoeval);
    end
end

% End of get_profile program...Returning to parent program
