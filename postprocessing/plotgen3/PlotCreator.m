function PlotCreator( caseName, plotNum, plotType, startTime, endTime, startHeight, endHeight, varargin )

%Display the variables that were passed in (debug output)
caseName
plotNum
plotType
startTime
endTime
startHeight
endHeight

%Figure out the number of optional arguments passed in
optargin = size(varargin,2);

%Optional argument format is as follows:
%'/PATH/TO/FILE', 'timeseries or profile', 'varname', 'title', 'units', 'lineWidth', 'lineType', 'lineColor'

%This means we can easily figure out the number of lines on the plot by dividing by 8
numLines = optargin / 8;

%Create a blank plot of the proper type so we have somewhere to draw lines
if strcmp(plotType, 'profile')
	
elseif strcmp(plotType, 'timeseries')

end

%Loop through each line on the plot
for i=1:numLines
	filePath = varargin{1 * i};
	varName = varargin{2 * i};
	varExpression = varargin{3 * i};
	plotTitle = varargin{4 * i};
	varUnits = varargin{5 * i};
	lineWidth = varargin{6 * i};
	lineType = varargin{7 * i};
	lineColor = varargin{8 * i};

	%Determine the type of file being read in
	extension = DetermineExtension(filePath);

	%Determine the variables that need to be read in
	varsToRead = ParseVariablesFromExpression(varExpression);

	%Read in the necessary variables
	for j=1:size(varsToRead,2);
		%We need to convert the variable name to read from a cell array to a string
		varString = cell2mat(varsToRead(j));
		disp(['Reading variable ', varString]);

		if strcmp(extension, 'ctl')
			variableData = VariableReadGrADS(varString, filePath);
		elseif strcmp(extension, 'nc')
			variableData = VariableReadNC(varString, filePath);
		end

		%Store the read in values to the proper variable name (ex. variable rtm will be read in to the variable named rtm,
		%this allows the expression to be used as is).
		eval([varString, '= variableData;']);
	end

	%Now evaluate the expression using the read in values,
	eval(['valueToPlot =', varExpression, ';']);
	
	%At this point, the value of the expression is contained in valueToPlot


end

%Output the EPS file
output_file_name = [ 'output/', caseName, '_', plotNum, '.eps' ];
print( '-depsc2', output_file_name )
