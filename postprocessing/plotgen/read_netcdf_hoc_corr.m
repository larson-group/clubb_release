
function avg_field = read_netcdf_hoc_corr(filename,nzmax,t1,t2,varnum,numvars)

% Reads and time-averages profiles of correlations from 4-D netCDF *.nc files
%   with 2 singleton dimensions and averages of the time dimension
% e.g.: corr_rrNr = read_netcdf_hoc('tune/arm_zt.nc',110,1,1,27,28)
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

%Determine number of dimensions the variable has
[varname, xtype, dimids, numatts] = netcdf.inqVar(fid,varnum);

% For correlations, initialize the number of timesteps with valid correlation
% values (not undefined) to 0.
num_timesteps = zeros(nzmax,1);

% Read in and average profiles of correlations.
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
   % Loop over every vertical level.  When the value of a correlation is valid
   % (not undefined), add its value to the running sum and add 1 to the time
   % step counter.  When a correlation is undefined (indicated by a value of
   % -999.0), ignore it in the averaging.
   for k = 1:1:nzmax
      if ( new_field(k) ~= -999.0 )
         avg_field(k) = avg_field(k) + new_field(k);
         num_timesteps(k) = num_timesteps(k) + 1;
      end % new_field(k) ~= -999.0
   end % k = 1:1:nzmax
end
% Produce the time-averaged result for the correlation.
% When all correlation values at a vertical level are undefined throughout the
% entire averaging period, simply output a NaN so that a line isn't plotted at
% those levels (alternatively, output a 0 value so that the plots aren't
% difficult to read).
for k = 1:1:nzmax
   if ( num_timesteps(k) > 0 )
      avg_field(k) = avg_field(k)/num_timesteps(k);
   else % num_timesteps(k) = 0
      avg_field(k) = NaN;
   end % num_timesteps(k) > 0
end % k = 1:1:nzmax
