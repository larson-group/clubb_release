%Given two arguments, 'default,' and 'fallback,' use the
%values contained in 'default' if available, otherwise
%use the values in 'fallback'.
%
%This works because PlotCreator.m sets variables with no
%data to a single '0' (thus having a size of 1).
function [data] = UseFirstArgIfExist( default, fallback )

%If default doesn't exist, use fallback, otherwise use
%default
if max(size(default)) == 1
	data = fallback;
else
	data = default;
end

end
