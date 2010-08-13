%This appends a McClatchy profile to the top of the initial sounding.
%It is intended to be called from layer_parse.m. It is probably not
%very generalized as it has only be used with the McClatchy profile
%provided for use with the TWP_ICE case.
%
%Arguments:
%	caseName 	 - the name of the case files are being made for
%	fileName 	 - the name of the file to read
%	P_sfc	 	 - pressure at the surface in mb
%	startHeight	 - the pressure level at which to start using the McClatchy profile (in Pa)
%	u		 - value of the u wind (in m/s)
%	v		 - valaue of the v wind (in m/s)

function parse_McClatchey( caseName, fileName, P_sfc, startHeight, u, v )

%We need the unit conversions file
addpath('../../../postprocessing/matlab_include/convert_units.m');

%Open the file for reading only
fid = fopen(fileName,'r'); %Open the forcing file read only


%Read in the first two lines of header data
inputText = textscan(fid,'%s',1,'delimiter','\n');

%Parse out the variable names
varNames = regexp(inputText{1}, '\s*\([A-Za-z\/]*\)\s*', 'split');

%Read in the data
geoPot = [];
lev = [];
T = [];
q = [];
ozone = [];
for j = 1:33,
	inputText = textscan(fid,'%f',5,'delimiter','\n');

	varData = inputText{1};

	geoPot = [geoPot;varData(1)];
	lev = [lev;varData(2)]; %Convert from Pa to hPa
	T = [T;varData(3)];
	q = [q;varData(4)];
	ozone = [ozone;varData(5)];
end

%Close the input file
fclose(fid);

%Create the necessary output files

%Do some unit conversions
%Convert pressure levels to height (requires pressure at surface)
%height = convert_units.pressure_in_hPa_to_height_m(T, lev, P_sfc);
height = geoPot;

%To convert omega to velocity, we need rho, to get rho, we need exner
%Note, pressure must be converted from Pa to hPa
exner = convert_units.pressure_in_hPa_to_exner(lev / 100);
%We also need potential temperature
theta = convert_units.temperature_to_potential_temperature(T, exner);

%Now we need to create the soundings.in file 
fileName = [caseName '_sounding.in'];

%Open the output file
fout = fopen(fileName,'a');

%Now populate the columns
for j = 1:33,
	%Output the variable
	if (lev(j) < startHeight)
		fprintf(fout,'%f    	%f   %E      %f    %f    %f    %f    %f', lev(j), theta(j), q(j), u, v, 0, u, v);
		fprintf(fout,'\n');
	end
end


%Close the output file
fclose(fout);	

end
