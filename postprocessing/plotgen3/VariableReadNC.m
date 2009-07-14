function [varData] = VariableReadNC( variableToRead, filePath )

%Read in some necessary information about the GRaDS file
[dataFileName, nz, z, numTimesteps, dt, numVars, listofparams] = header_read_expanded_netcdf(filePath);

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
levData = 0;

for i = 1:numVars
	%See if the variable we found is the variable we are interested in
	if ( strcmp( strtrim(listofparams(i, :)), variableToRead ) )
		varData = read_netcdf_hoc(filePath, nz, t_start, t_end, i, numVars);
	end
end

levData = z;
