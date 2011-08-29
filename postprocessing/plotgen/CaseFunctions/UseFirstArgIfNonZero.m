%Given two arguments, 'default,' and 'fallback,' use the
%values contained in 'default' if non zero, otherwise
%use the values in 'fallback'.
%
function [data] = UseFirstArgIfNonZero( default, fallback )

%If default doesn't exist, use fallback, otherwise use
%default
if default == 0
	data = fallback;
else
	data = default;
end

end
