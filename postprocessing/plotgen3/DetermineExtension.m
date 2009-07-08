function [extension] = DetermineExtension( filePath )

lastDot =  max(regexp( filePath, '\.' )) + 1;
extension = filePath(lastDot:size(filePath,2));
