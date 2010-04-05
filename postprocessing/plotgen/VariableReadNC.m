function [varData, levData] = VariableReadNC( filePath, variableToRead, startTime, endTime, plotType )

%Read in some necessary information about the NETCDF file
[dataFileName, nz, z, numTimesteps, dt, numVars, listofparams] = header_read_expanded_netcdf(filePath);

% open NETCDF file
fid = netcdf.open(filePath,'NC_NOWRITE');

%Ensure the file will be closed no matter what happens
cleanupHandler = onCleanup(@()netcdf.close(fid));

%Calculate start and end timesteps, correct if they are invalid
t_start = ceil(startTime / dt);
if ( t_start < 1 )
	t_start = 1;
elseif ( t_start > numTimesteps )
	t_start = numTimesteps;
end

t_end = ceil(endTime / dt);
if ( t_end < 1 )
	t_end = 1;
elseif ( t_end > numTimesteps )
	t_end = numTimesteps;
end

%Set a default value if the passed in variable is not found
varData(1:nz) = 0;

varNum = netcdf.inqVarID(fid, variableToRead);
	
if strcmp(plotType, 'profile')
	varData = read_netcdf_hoc(filePath, nz, t_start, t_end, varNum, numVars);
elseif strcmp(plotType, 'timeseries')
	varData = read_netcdf_hoc_timeseries(filePath, nz, t_start, t_end, varNum, numVars);
end

%Determine if the data is in terms of ascending or descending height, adjust if necessary
%This is a kludge that will likely not work for NetCDF files defined in terms of pressure.
if z(max(size(z))) - z(1) < 0
	%Data is descending
	varData = flipud(varData);
	levData = flipud(z);
else
	%Data is asscending
	levData = z;
end
