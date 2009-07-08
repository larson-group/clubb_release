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
numLines = optargin / 8;

%Loop through each line
for i=1:numLines
	filePath = varargin{1 * i};
	plotType = varargin{2 * i};
	varExpression = varargin{3 * i};
	plotTitle = varargin{4 * i};
	varUnits = varargin{5 * i};
	lineWidth = varargin{6 * i};
	lineType = varargin{7 * i};
	lineColor = varargin{8 * i};

	%Determine the variables that need to be read in
	varsToRead = ParseVariablesFromExpression(varExpression);
end
