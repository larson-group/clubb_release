function PlotCreator( plotNum, startTime, endTime, startHeight, endHeight, varargin )

%Display the variables that were passed in (debug output)
plotNum
startTime
endTime
startHeight
endHeight

%Figure out the number of optional arguments passed in
optargin = size(varargin,2);

%Optional argument format is as follows:
%'/PATH/TO/FILE', 'timeseries or profile', 'varname', 'title', 'units', 'lineWidth', 'lineType', 'lineColor'

%This means we can easily figure out the number of lines on the plot by dividing by 8
numLines = optargin / 9;

%Loop through each line
for i=1:numLines
	filePath = varargin{1 * i};
	plotType = varargin{2 * i};
	varName = varargin{3 * i};
	varExpression = varargin{4 * i};
	plotTitle = varargin{5 * i};
	varUnits = varargin{6 * i};
	lineWidth = varargin{7 * i};
	lineType = varargin{8 * i};
	lineColor = varargin{9 * i};

	%Determine the variables that need to be read in
	varsToRead = ParseVariablesFromExpression(varExpression);

	%Read in the necessary variables
	for j=1:size(varsToRead,2);
		varsToRead(j)
	end
end
