function [color] = ParseColorFromExpression( colorString )

if (strcmp(colorString, 'y') || strcmp(colorString, 'yellow'))
	color = [1.0, 1.0, 0.0];
elseif (strcmp(colorString, 'm') || strcmp(colorString, 'magenta'))
	color = [1.0, 0.0, 1.0];
elseif (strcmp(colorString, 'c') || strcmp(colorString, 'cyan'))
	color = [0.0, 1.0, 1.0];
elseif (strcmp(colorString, 'r') || strcmp(colorString, 'red'))
	color = [1.0, 0.0, 0.0];
elseif (strcmp(colorString, 'g') || strcmp(colorString, 'green'))
	color = [0.0, 1.0, 0.0];
elseif (strcmp(colorString, 'b') || strcmp(colorString, 'blue'))
	color = [0.0, 0.0, 1.0];
elseif (strcmp(colorString, 'w') || strcmp(colorString, 'white'))
	color = [1.0, 1.0, 1.0];
elseif (strcmp(colorString, 'k') || strcmp(colorString, 'black'))
	color = [0.0, 0.0, 0.0];
else
	%Assume the passed in string is a color array and determine
	%the proper values
	color = eval([colorString]);	
end
