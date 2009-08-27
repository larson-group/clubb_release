function [variables] = ParseVariablesFromExpression( expression )

varsToCheck = regexp( expression, '[A-Za-z0-9_]*', 'match' );

%We need to preallocate arrays for speed
variables(1,size(varsToCheck,2)) = cell(1);

varsFound = 0;

%Now, loop through all the variables found, if it is not a
%MATLAB keyword, assume it is a variable and add it to the list
for i=1:size(varsToCheck,2)
	%Current use contains a collection of all of the current
	%ways a variable is already in use, this keeps us from
	%interpreting something like 'sqrt' as a variable
	varToCheck = char(varsToCheck(i));
	
	currentUse = which('-all', varToCheck);

	%Make sure the parsed out variable name is both a valid variable
	%name and not already in use
	if (isvarname(varToCheck) && isempty(currentUse))
		varsFound = varsFound + 1;
		variables(varsFound) = varsToCheck(i);
	end
end

variables = variables(1:varsFound); %remove excess elements
