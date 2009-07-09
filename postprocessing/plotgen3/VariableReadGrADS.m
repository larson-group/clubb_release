function [varData] = VariableReadGrADS( filePath, variableToRead, startTime, endTime )

%Read in some necessary information about the GRaDS file
[dataFileName, nz, z, numTimesteps, dt, numVars, listofparams] = header_read_expanded(filePath);

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

%Append the data file name to the path of the passed in file
lastSlash =  max(regexp( filePath, '\/' ));
dataFilePath = filePath(1:lastSlash);

%Set a default value if the passed in variable is not found
varData(1:nz) = 0;

for i = 1:numVars
	%See if the variable we found is the variable we are interested in
	if ( strcmp( strtrim(listofparams(i, :)), variableToRead ) )
		varData = read_grads_hoc_endian([dataFilePath, dataFileName], 'ieee-le', nz, t_start, t_end, i, numVars);
	end
end
