function [extension] = DetermineExtension( filePath )

[dummy, dummy, extension, dummy] = fileparts(filePath);
