%This script reads a ARM style setup files and creates a sounding.in
%and forcing.in files that can be used by CLUBB.
%
%Arguments:
%	caseName - the name of the case files are being made for
%	fileName - the name of the file to read
%	sndgTime - the time the initial sounding should be generated from (in seconds)
%	Psfc	 - pressure at the sfc in mb
function convert_arm_input_data( caseName, fileName, sndgTime, Psfc )

seconds_per_day = 86400;
seconds_per_hour = 3600;
g_per_kg = 1000;
Pa_per_hPa = 100;

%We need the unit conversions file
%addpath('../../../postprocessing/matlab_include/convert_units.m');

%Open the file for reading only
fid = fopen(fileName,'r'); %Open the forcing file read only

%Read in the first two lines of header data
inputText = textscan(fid,'%s',1,'delimiter','\n'); %Read a line delimited by a carriage return
inputText = textscan(fid,'%u',1,'delimiter','\n');

fieldLength = inputText{1}; % Pull the field length (the second line read)


%Read in the next two lines of header data
inputText = textscan(fid,'%s',1,'delimiter','\n'); %Read a line delimited by a carriage return
inputText = textscan(fid,'%f',1,'delimiter','\n');

numLevels = inputText{1}; % Pull the number of pressure levels


%Read in the pressure levels
inputText = textscan(fid,'%s',1,'delimiter','\n'); %This is the field descriptor

inputText = textscan(fid,'%f',numLevels,'delimiter','\n');

lev = inputText{1};


%Read in the time
inputText = textscan(fid,'%s',1,'delimiter','\n'); %This is the field descriptor

inputText = textscan(fid,'%f',fieldLength,'delimiter','\n');

time = inputText{1};

dt = time(2) - time(1); %This will be in a fraction of a day

%Convert dt from days to seconds
dt = dt * seconds_per_day;

%Figure out when the case starts
tStart = time(1) - floor(time(1));

%Convert to seconds
tStart = tStart * seconds_per_day;

%If the user requests a sounding generated from before the data starts, throw an error
if (tStart > sndgTime)
	'Sounding cannot be produced before data start.'
end

%Recreate the time matrix as seconds since midnight
tIndex = 0;
for i = 1:max(size(time)),
	time(i) = tStart + (dt * (i - 1));

	if ((time(i) >= sndgTime) && (tIndex == 0))
		tIndex = i;
	end
end

%Next fields are year, month, day, hour, and minute. These are don't cares
inputText = textscan(fid,'%s',1,'delimiter','\n'); %This is the field descriptor
inputText = textscan(fid,'%f',fieldLength,'delimiter','\n');
inputText = textscan(fid,'%s',1,'delimiter','\n'); %This is the field descriptor
inputText = textscan(fid,'%f',fieldLength,'delimiter','\n');
inputText = textscan(fid,'%s',1,'delimiter','\n'); %This is the field descriptor
inputText = textscan(fid,'%f',fieldLength,'delimiter','\n');
inputText = textscan(fid,'%s',1,'delimiter','\n'); %This is the field descriptor
inputText = textscan(fid,'%f',fieldLength,'delimiter','\n');
inputText = textscan(fid,'%s',1,'delimiter','\n'); %This is the field descriptor
inputText = textscan(fid,'%f',fieldLength,'delimiter','\n');
%End don't cares


%Read in the the number of fields
inputText = textscan(fid,'%s',1,'delimiter','\n'); %This is the field descriptor

inputText = textscan(fid,'%u',1,'delimiter','\n');

numFields = inputText{1};


%Now that we've chewed through all of the header information, we can read in the actual data
for i = 1:numFields,
	%Read in the field descriptor
	inputText = textscan(fid,'%s',1,'delimiter','\n');

	%Parse out the variable name
	varName = regexp(inputText{1}, '\[', 'split');
	varName = varName{1}{1};

	%Read in the data
	varData = [];
	for j = 1:numLevels,
		inputText = textscan(fid,'%f',fieldLength,'delimiter','\n');

		varData = [varData;inputText{1}];
	end

	%Store any variables needed for the soundings.in file
	if strcmp('T', varName)
		T = varData;
	end

	if strcmp('u', varName)
		u = varData;
	end

	if strcmp('v', varName)
		v = varData;
	end

	if strcmp('q', varName)
		q = varData;
	end

	if strcmp('omega', varName)
		omega = varData;
	end

	if strcmp('T_adv_h', varName)
		dTdt = varData;
	end

	if strcmp('q_adv_h', varName)
		dqdt = varData;
	end
	%End variable storage
end

%Close the input file
fclose(fid);

%Do some unit conversions
%Convert pressure levels to height (requires pressure at sfc)
height = convert_units.pressure_in_hPa_to_height_m(T(1:size(lev,1)), lev, Psfc);

exner = convert_units.pressure_in_hPa_to_exner(lev);
exnerRot = rot90(exner);
%We also need potential temperature
T_initial = [];
for j = 0:(numLevels - 1),
	%Output the variable
	index = (j * fieldLength) + tIndex;
	T_initial = [T_initial, T(index)];
end
theta = convert_units.temperature_to_potential_temperature(T_initial, exnerRot);


%Convert dTdt from temperature tendency [K/hr] to potential temperature tendency [K/s]
for i = 0:(numLevels - 1),
	for j = 1:fieldLength,
		index = (i * fieldLength)  + j;
		dTdt(index) = dTdt(index) * (1 / exner(i + 1)) * (1 / (seconds_per_hour));
	end
end

%Convert dqdt from tendency of humidity mixing ratio [g/kg/hr] to water vapor advective tendency [kg/kg/s]
dqdt = dqdt .* (1 / g_per_kg / seconds_per_hour);


%Now we need to create the soundings.in file 
fileName = [caseName '_sounding.in'];

%Open the output file
fout = fopen(fileName,'w');

%Print out the header row
fprintf(fout,'Press[Pa]	        thm[K]     rt[kg\\kg]      u[m\\s]       v[m\\s]	omega[Pa\\s]       ug[m\\s]       vg[m\\s]\n');

%Now populate the columns
for j = 0:(numLevels - 1),
	%Output the variable
	index = (j * fieldLength) + tIndex;
	fprintf(fout,'%f    	%f   %E      %f    %f    %f    %f    %f', lev(j + 1) * Pa_per_hPa, theta(j+1), (q(index) / g_per_kg), u(index), v(index), omega(index) * Pa_per_hPa / seconds_per_hour, u(index), v(index));
	fprintf(fout,'\n');
end

%Close the output file
fclose(fout);	

%Now, append the McClatchy profile
parse_McClatchey( caseName, 'McCProfiles.dat', Psfc, lev(numLevels) * Pa_per_hPa, u((numLevels - 1) * fieldLength + 1), v((numLevels - 1) * fieldLength + 1) );


%Now we need to create the forcings.in file 
fileName = [caseName '_forcings.in'];

%Open the output file
fout = fopen(fileName,'w');

%Print out the header row
fprintf(fout,'z[m]     thlm_f[K\\s]     rtm_f[kg\\kg\\s]     um_ref[m\\s]     vm_ref[m\\s]     um_f[m\\s^2]     vm_f[m\\s^2]     omega[mb\\hr]     ug[m\\s]     vg[m\\s]\n');

%Now populate the columns
for i = 1:fieldLength,
	%Output the variable
	timeVal = i * dt;
	%Output the header for this data block
	if timeVal > sndgTime
		
		fprintf(fout,'%f     %i\n', timeVal - sndgTime, numLevels);

		for j = 0:(numLevels - 1),
			index = i + (j * fieldLength);
			fprintf(fout,'%f     %E     %E     %f     %f     %f     %f     %f     %f     %f\n', height(j + 1), dTdt(index), dqdt(index), u(index), v(index), -999.9, -999.9, omega(index), -999.9, -999.9);
		end
	end	
end

%Close the output file
fclose(fout);	

%Now we need to create the sfc.in file 
fileName = [caseName '_sfc.in'];

%Open the output file
fout = fopen(fileName,'w');

%Print out the header row
fprintf(fout,'Time[s]     LH[W\\m^2]     SH[W\\m^2]     thlm[K]     rt[kg\\kg]     Press[Pa]\n');

%Now populate the columns
for i = 1:max(size(time)),
	%Output the variable
	timeVal = i * dt;

	fprintf(fout,'%f     %f     %f     %f     %f     %f\n', timeVal, 0, 0, 0, 0, 0);
end

%Close the output file
fclose(fout);	

quit;

end
