%This script converts an ill-formatted NC file to a NC file that 
%can be read with any new version of grads (not-gradsnc).
%Author: Brandon Nielsen,  Jun 2008.

%It should be noted that the script seems to work intermittently
%when using Matlab 2008a. This could likely be fixed if someone
%had the time to track down the exact problem.
%   Perhaps MATLAB has a memory leak: whenever a new variable
% is read, MATLAB stores a new copy of the entire file in memory.
% Perhaps a workaround is to close and re-open the file repeatedly.
% This is speculation, and we haven't tried it.


%You need both snctools and mexnc, they can be downloaded here:
%http://mexcdf.sourceforge.net/downloads/


%The important variables are:
%snctools_path->the path to where you have SNCTools installed at,
%no formal install is necessary, just unzip the archive

%mexnc_path->the path where MexNC is installed at, again, no
%formal install is necessary, just unzip the archive

%to_convert->the ill-formatted NC file to convert
%new_file->the converted nc file
%-------Program Setup Variables------
snctools_path = '/home/nielsenb/snctools';
mexnc_path = '/home/nielsenb/mexnc';

%Setup the path to include the necessary tools to read NetCDF files
path(path, snctools_path);
path(path, mexnc_path);

%to_convert = 'fire_zt.nc';
to_convert = 'NetCDFTest.nc';
new_file = 'NewFile.nc';
%-------End Program Setup-------


%Dump the original NC file so we can be sure its valid data
nc_dump(to_convert);


%------Shell Creation------
%Create a shell of a NetCDF file in the proper format
nc_create_empty(new_file);

%Add the 1 length longitude and latitude dimensions
nc_add_dimension(new_file, 'longitude', 1);
nc_add_dimension(new_file, 'latitude', 1);

%Determine altitude information from the original file
altitude = nc_varget(to_convert, 'z');
alt_size = max(size(altitude));

%Add the newly calculated altitude dimension
nc_add_dimension(new_file, 'altitude', alt_size);

%Determine time information from the original file
time = nc_varget(to_convert, 'time');
time_size = max(size(time));

%Add the newly calculated time dimension
nc_add_dimension(new_file, 'time', time_size);
%------End Shell Creation------


%------Variable creation------
%Create the longitude variable
lon_varstruct.Name = 'longitude';
lon_varstruct.Nctype = nc_float;
lon_varstruct.Dimension = {'longitude'};
nc_addvar(new_file, lon_varstruct);
nc_varput(new_file, 'longitude', 0);

%Create the latitude variable
lat_varstruct.Name = 'latitude';
lat_varstruct.Nctype = nc_float;
lat_varstruct.Dimension = {'latitude'};
nc_addvar(new_file, lat_varstruct);
nc_varput(new_file, 'latitude', 0);

%Create the altitude variable
alt_varstruct.Name = 'altitude';
alt_varstruct.Nctype = nc_float;
alt_varstruct.Dimension = {'altitude'};
%Retrieve the altitude data so we can populate the new file
alt_data = nc_varget(to_convert, 'z');
nc_addvar(new_file, alt_varstruct);
%Populate the new variable
nc_varput(new_file, 'altitude', alt_data);

%Create the time variable
time_varstruct.Name = 'time';
time_varstruct.Nctype = nc_double;
time_varstruct.Dimension = {'time'};
%Retrieve the time data so we can populate the new file
time_data = nc_varget(to_convert, 'time');
nc_addvar(new_file, time_varstruct);
%Populate the new variable
%Iterate through the array of time values and convert to minutes
i = 1;
while(i <= max(size(time_data)))
    %time(i) = time(i) * 864000000
    time_data(i) = i;
    i = i + 1;
end
nc_varput(new_file, 'time', time_data);
%------End Variable Creation------


%------Attribute Configuration------
%Set up the X variable
nc_attput(new_file, 'longitude', 'cartesian_axis', 'X');
nc_attput(new_file, 'longitude', 'units', 'degrees_E');
nc_attput(new_file, 'longitude', 'ipositive', '1');

%Set up the Y variable
nc_attput(new_file, 'latitude', 'cartesian_axis', 'Y');
nc_attput(new_file, 'latitude', 'units', 'degrees_N');
nc_attput(new_file, 'latitude', 'ipositive', '1');

%Set up the Z variable
nc_attput(new_file, 'altitude', 'cartesian_axis', 'Z');
nc_attput(new_file, 'altitude', 'units', 'meters');
nc_attput(new_file, 'altitude', 'positive', 'up');
nc_attput(new_file, 'altitude', 'ipositive', '1');

%Set up the time variable
nc_attput(new_file, 'time', 'cartesian_axis', 'T');
nc_attput(new_file, 'time', 'units', 'minutes since 1987-07-07 00:01:00.0');
nc_attput(new_file, 'time', 'ipositive', '1');
nc_attput(new_file, 'time', 'calendar_type', 'Gregorian');
%------End Attribute Configuration------


%------Copy variables from old file to new file------
%Open file "to_convert" as NC file ncid_1
[ncid_1, status] = mexnc('open', to_convert, nc_nowrite_mode);
[ncid_2, status] = mexnc('open', new_file, nc_share_mode);

%Determine where we need to start copying variables
[start_varid, status] = mexnc('inq_varid', ncid_1, 'time');
start_varid = start_varid + 1;

%Loop through all variables
clear i;
i = start_varid;
while(status == 0)
    %Clear used variables
    clear new_data;
    
    mexnc('redef', ncid_2);
    
    %Extract the variable name and type from the original file
    [var_name, status] = mexnc('inq_varname', ncid_1, i);
    [var_type, status] = mexnc('inq_vartype', ncid_1, i);
   
    
    %Echo the variable name so we know the script hasn't hung
    var_name
    new_varstruct.Name = var_name;
    %Recreate the variable in the new file
    new_varstruct.Nctype = var_type;
    %Get the number of dimensions of the variable so we know how many
    %dimensions the new variable should have
    [num_dims, status] = mexnc('inq_varndims', ncid_1, i);
    if num_dims == 1
        new_varstruct.Dimension = {'time'};
    else
        new_varstruct.Dimension = {'time', 'altitude'};
    end
    nc_addvar(new_file, new_varstruct);
    
    %Copy the attributes from the old file to the new one
    nc_attput(new_file, var_name, 'valid_range', '-3.402823e+38f, 3.402823e+38f');
    
    %Copy the name attribute
    [att_name, status] = mexnc('inq_attname', ncid_1, i, 0);
    [att_value, status] = mexnc('get_att_text', ncid_1, i, att_name);
    nc_attput(new_file, var_name, 'title', att_value);
    
    %Copy the units attribute
    [att_name, status] = mexnc('inq_attname', ncid_1, i, 1);
    [att_value, status] = mexnc('get_att_text', ncid_1, i, att_name);
    nc_attput(new_file, var_name, att_name, att_value);
    
    mexnc('enddef', ncid_2);
    
    %Copy the data from the old file to the new one
    %new_data = nc_varget(to_convert, var_name);
    %Get the proper varid
    %[var_id, status] = mexnc('inq_varid', ncid_2, var_name);
    %Make sure to use the proper variable write method
    %if var_type == NC_BYTE
    %    mexnc('put_var_byte', ncid_2, var_id, new_data);
    %elseif var_type == NC_CHAR
    %    mexnc('put_var_char', ncid_2, var_id, new_data);
    %elseif var_type == NC_SHORT
    %    mexnc('put_var_short', ncid_2, var_id, new_data);
    %elseif var_type == NC_INT
    %    mexnc('put_var_int', ncid_2, var_id, new_data);
    %elseif var_type == NC_FLOAT
    %    mexnc('put_var_float', ncid_2, var_id, new_data);
    %elseif var_type == NC_DOUBLE
    %    mexnc('put_var_double', ncid_2, var_id, new_data);
    %end
    %nc_varput(new_file, var_name, new_data);
    
    %Check to see if there is another variable to move    
    i = i + 1;
    [var_name, status] = mexnc('inq_varname', ncid_1, i);
end
%------End Variable Copy------

%Close any open files
mexnc('close', ncid_1);
mexnc('close', ncid_2);

%------Copy variables from old file to new file------
%Open file "to_convert" as NC file ncid_1
[ncid_1, status] = mexnc('open', to_convert, nc_nowrite_mode);

%Loop through all variables
clear i;
i = start_varid;
while(status == 0)
    %Clear used variables
    clear new_data;
    
    %Extract the variable name and type from the original file
    [var_name, status] = mexnc('inq_varname', ncid_1, i);
    
    %Echo the variable name so we know the script hasn't hung
    var_name
    
    %Copy the data from the old file to the new one
    new_data = nc_varget(to_convert, var_name);
    nc_varput(new_file, var_name, new_data);
    
    %Check to see if there is another variable to move    
    i = i + 1;
    [var_name, status] = mexnc('inq_varname', ncid_1, i);
end
%------End Variable Copy------

mexnc('close', ncid_1);

%Dump the modified file to ensure it is correct
nc_dump(new_file);
