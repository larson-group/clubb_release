% This file is a slightly altered version of function header_read from file
% header_read.m.  The difference between this version and that version is
% that this version reads and outputs the GrADS output timestep length, in
% addition to everything else that the original version read and output.
% Another change is this version is for a NETCDF file. It us used for replicating a similar behavior
% to the functions used for GrADs files.
function [filename, nzmax, z, t_time_steps, time_step_length, numvars, listofparams] ...
   = header_read_expanded_netcdf(file_header)
%function [filename,nzmax,z,t_time_steps,numvars,listofparams] = header_read(file_header)
% header_read('gabls3_night.nc')
% Opens file and initializes some of the values.  
% Input: file_header       --    The header file provided by moments code
% Output: filename         --    The file containing data to be plotted
%         nzmax               --    The total number of z levels
%         z                --    The heights in the sounding, in vector form
%         t_time_steps     --    The total number of time steps for the run
%         time_step_length --    The length of the time steps in minutes
%         numvars          --    The total number of variables
%         listofparams     --    The names of all the variables, in vector form

[dummy, filename, dummy] = fileparts(file_header);

nc_file = netcdf.open(file_header, 'NC_NOWRITE');

%Ensure the file will be closed no matter what happens
cleanupHandler = onCleanup(@()netcdf.close(nc_file));

[ndims,numvars,ngatts,unlimdimid] = netcdf.inq(nc_file);

if strfind ( filename, 'pf_gfdl' ) %Get height for pf_gfdl files
  %Get the correct variables for height
   varAlt = netcdf.inqVarID( nc_file, 'zfull' );
   alt = netcdf.getVar ( nc_file, varAlt );
   z(:,1) = alt(1,1,:,1);
elseif strfind ( filename, 'ph_gfdl' )  %Get height for ph_gfdl files
   varAlt = netcdf.inqVarID( nc_file, 'zhalf' );
   alt = netcdf.getVar ( nc_file, varAlt );
   z(:,1) = alt(1,1,:,1);
elseif strfind ( filename, 'ps_gfdl' )  %Get height for ps_gfdl files
    z = netcdf.getVar( nc_file, 3);
elseif strfind ( filename, '_cam' ) %For cam cases
   varAlt = netcdf.inqVarID( nc_file, 'Z3');
   alt = netcdf.getVar( nc_file, varAlt );
   varPhis = netcdf.inqVarID( nc_file, 'PHIS');
   PHIS = netcdf.getVar( nc_file, varPhis );
   z(:,1) = alt(1,1,:,1) - (PHIS(1,:,1)/9.81);
elseif strfind ( filename, 'cam.h0' ) %For cam cases
   varAlt = netcdf.inqVarID( nc_file, 'Z3');
   alt = netcdf.getVar( nc_file, varAlt );
   varPhis = netcdf.inqVarID( nc_file, 'PHIS');
   PHIS = netcdf.getVar( nc_file, varPhis );
   z(:,1) = alt(1,1,:,1) - (PHIS(1,:,1)/9.81);
elseif strfind ( filename, 'twp_ice' ) % For twp_ice
   zvarid = netcdf.inqVarID(nc_file, 'altitude');
   z = netcdf.getVar( nc_file, zvarid );
elseif strfind ( filename, '_zt' )
   zvarid = netcdf.inqVarID(nc_file, 'altitude');
   z = netcdf.getVar( nc_file, zvarid );
elseif strfind ( filename, '_zm' )
   zvarid = netcdf.inqVarID(nc_file, 'altitude');
   z = netcdf.getVar( nc_file, zvarid );
elseif strfind ( filename, '_sfc' )
   zvarid = netcdf.inqVarID(nc_file, 'altitude');
   z = netcdf.getVar( nc_file, zvarid );
else
   zvarid = netcdf.inqVarID(nc_file, 'z');
   z = netcdf.getVar( nc_file, zvarid );
end

nzmax = size(z, 1);


varTime = netcdf.inqVarID( nc_file, 'time' );
time = netcdf.getVar( nc_file, varTime );

t_time_steps = size(time,1);

if (t_time_steps > 1)
  % dt is defined as the difference between any two successive outputs in
  % the file. We only check one pair of such outputs here.
  dt = time(t_time_steps) - time(t_time_steps - 1);
elseif (t_time_steps == 1)
  % We cannot use two output times for comparison. Let's define dt as the
  % amount of simulation time before the first (and only) output!
  dt = time(1);
else
  % dt is undefined.
  dt = -999.0;
end

timeInfo = netcdf.getAtt( nc_file, varTime, 'units' );


%Now that we have the text info about the unit for time, see if we can determine
%a proper dt
if findstr(timeInfo, 'seconds')
	dt = (1 / 60) * dt;
elseif findstr(timeInfo, 'minutes')
	dt = 1 * dt;
elseif findstr(timeInfo, 'hours')
	dt = 60 * dt;
elseif findstr(timeInfo, 'days')
	dt = 24 * 60 * dt;
	time = 24 * 60 * time;
else
	%Assume one dt is 1 minute
	dt = 1 * dt;
end

time_step_length = dt;

for i=5:numvars-1
	[varname,xtype,dimids,natts] = netcdf.inqVar( nc_file , i );
	listofparams(i-4,1:size(varname,2)) = varname(1:size(varname,2));
end

numvars = numvars - 5;

end
