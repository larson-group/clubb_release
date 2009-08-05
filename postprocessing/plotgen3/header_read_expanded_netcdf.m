% This file is a slightly altered version of function header_read from file
% header_read.m.  The difference between this version and that version is
% that this version reads and outputs the GrADS output timestep length, in
% addition to everything else that the original version read and output.
% Another change is this version is for a NETCDF file. It us used for replicating a similar behavior
% to the functions used for GrADs files.
function [filename, nz, z, t_time_steps, time_step_length, numvars, listofparams] ...
   = header_read_expanded_netcdf(file_header)
%function [filename,nz,z,t_time_steps,numvars,listofparams] = header_read(file_header)
% header_read('gabls3_night.nc')
% Opens file and initializes some of the values.  
% Input: file_header       --    The header file provided by moments code
% Output: filename         --    The file containing data to be plotted
%         nz               --    The total number of z levels
%         z                --    The heights in the sounding, in vector form
%         t_time_steps     --    The total number of time steps for the run
%         time_step_length --    The length of the time steps in minutes
%         numvars          --    The total number of variables
%         listofparams     --    The names of all the variables, in vector form

[dummy, filename, dummy, dummy] = fileparts(file_header);

nc_file = netcdf.open(file_header, 'NC_NOWRITE');

[ndims,numvars,ngatts,unlimdimid] = netcdf.inq(nc_file);

z = netcdf.getVar( nc_file, 2 );

nz = size(z, 1);

time = netcdf.getVar( nc_file, 3 );

t_time_steps = size(time);

timeInfo = netcdf.getAtt( nc_file, 3, 'units' );

%Now that we have the text info about the unit for time, see if we can determine
%a proper dt
if findstr(timeInfo, 'seconds')
	dt = 1 / 60;
elseif findstr(timeInfo, 'minutes')
	dt = 1;
elseif findstr(timeInfo, 'hours')
	dt = 60;
else
	%Assume one dt is 1 minute
	dt = 1;
end

time_step_length = dt;

for i=5:numvars-1
	[varname,xtype,dimids,natts] = netcdf.inqVar( nc_file , i );
	while(size(varname,2) < 15 )
        	varname(size(varname,2) + 1) = ' ';
	end
	listofparams(i-4,1:15) = varname(1:15);
end

numvars = numvars - 5;

end
