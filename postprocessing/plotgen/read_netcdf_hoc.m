
function avg_field = read_netcdf_hoc(filename,nz,t1,t2,varnum,numvars)

% Reads and time-averages profiles from 3D NETCDF *.nc files.
% thlm = read_netcdf_hoc('tune/arm_zt.dat',110,1,1,1,28)
% nz = number of z levels in profile
% t1 = beginning timestep to look at
% t2 = ending timestep to look at
% varnum = which variable to read (see .ctl file)
% numvars = total number of variables in grads file (see .ctl file)

% open NETCDF file
fid = netcdf.open(filename,'NC_NOWRITE');

varnum = varnum+4;

num_timesteps = (t2-t1) + 1;

% Read in and average profiles over all timesteps
avg_field = zeros(nz,1);
for t=t1:t2
   field = netcdf.getVar(fid,varnum);
   new_field(1:nz,1) = squeeze(field(1,1,1:nz,t));
   avg_field = avg_field + new_field;
end
avg_field = avg_field/num_timesteps;


% close NETCDF file
netcdf.close(fid);
end
