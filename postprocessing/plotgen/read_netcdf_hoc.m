
function avg_field = read_netcdf_hoc(filename,nzmax,t1,t2,varnum,numvars)

% Reads and time-averages profiles from 4-D netCDF *.nc files 
%   with 2 singleton dimensions and averages of the time dimension
% e.g.: thlm = read_netcdf_hoc('tune/arm_zt.nc',110,1,1,1,28)
% nzmax = number of z levels in profile
% t1 = beginning timestep to look at
% t2 = ending timestep to look at
% varnum = which variable to read
% numvars = total number of variables in the netCDF file (can be determined with
%   the commands ncdump -h <filename> from the commend line.

% open NETCDF file
fid = netcdf.open(filename,'NC_NOWRITE');

%Ensure the file will be closed no matter what happens
cleanupHandler = onCleanup(@()netcdf.close(fid));

num_timesteps = (t2-t1) + 1;

%Determine number of dimensions the variable has
[varname, xtype, dimids, numatts] = netcdf.inqVar(fid,varnum);

% Read in and average profiles over all timesteps
avg_field = zeros(nzmax,1);
field = netcdf.getVar(fid,varnum);
for t=t1:t2
   %If there are only 2 dimensions, there is no need to squeeze (used for cloud_feedback LES)
   if max(size(dimids)) == 2
     new_field(1:nzmax,1) = field(1:nzmax, t);
   else
     %Sqeeze the array if there are more than 2 dimensions
     new_field(1:nzmax,1) = squeeze(field(1,1,1:nzmax,t));
   end
   avg_field = avg_field + new_field;
end
avg_field = avg_field/num_timesteps;

end
