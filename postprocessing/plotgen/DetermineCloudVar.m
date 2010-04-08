function [cldData] = DetermineCloudVar( cldFrac, cld )

%Given two cloud variables, determines if cloud_frac exists.
%If yes, returns cloud_frac, otherwise, returns cld
if cldFrac == 0
	cldData = cld;
else
	cldData = cldFrac;
end

end
