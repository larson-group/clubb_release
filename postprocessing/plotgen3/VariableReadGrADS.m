function [varData, levData] = VariableReadGrADS( filePath, variableToRead, startTime, endTime, plotType )

%Read in some necessary information about the GRaDS file
[dataFileName, nz, z, numTimesteps, dt, numVars, listofparams] = header_read_expanded(filePath);

%Determine the endianness
fid = fopen(filePath, 'rt');
% While your not at the end of file, will advance line by line through the file.
% Searches keywords to extract needed values from the string using if statements.
%  Once the correct string is found, it chops it down to extract the specific value.
while feof(fid) == 0
	tline = fgetl(fid);
    
	if findstr(tline, 'OPTIONS');
		[remainder_1,tline] = strtok(tline);
		endianness = strtrim(strrep(tline,' ^',''));
	end  
end
fclose(fid);

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
lastSlash = max(regexp( filePath, '\/' ));
dataFilePath = filePath(1:lastSlash);

%Set a default value if the passed in variable is not found
varData(1:nz) = 0;
levData = 0;

for i = 1:numVars
	%See if the variable we found is the variable we are interested in
	if ( strcmp( strtrim(listofparams(i, :)), variableToRead ) )
		if strcmp(plotType, 'profile')
			if strcmp(endianness, 'BIG_ENDIAN')
				varData = read_grads_hoc_endian([dataFilePath, dataFileName], 'ieee-be', nz, t_start, t_end, i, numVars);
			elseif strcmp(endianness, 'LITTLE_ENDIAN')
				varData = read_grads_hoc_endian([dataFilePath, dataFileName], 'ieee-le', nz, t_start, t_end, i, numVars);
			end
		elseif strcmp(plotType, 'timeseries')
			if strcmp(endianness, 'BIG_ENDIAN')
				varData = read_grads_hoc_sfc_endian([dataFilePath, dataFileName], 'ieee-be', nz, t_start, t_end, i, numVars);
			elseif strcmp(endianness, 'LITTLE_ENDIAN')
				varData = read_grads_hoc_sfc_endian([dataFilePath, dataFileName], 'ieee-le', nz, t_start, t_end, i, numVars);
			end

		end

		break;
	end
end

levData = z;
