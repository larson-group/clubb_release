function [extension] = DetermineExtension( filePath )

[dummy, dummy, extension] = fileparts(filePath);
