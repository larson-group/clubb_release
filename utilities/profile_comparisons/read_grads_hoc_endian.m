% $Id: read_grads_hoc_endian.m,v 1.1 2007-02-28 17:54:24 dschanen Exp $
function avg_field = read_grads_hoc_endian(filename,MachineFormat,nz,t1,t2,varnum,numvars)

% Reads and time-averages profiles from 1D GrADS *.dat files.
% thlm = read_grads_hoc('tune/arm_zt.dat','ieee-le',110,1,1,1,28)
%
% MachineFormat = string that specifies big-endian ('ieee-be', Linux COAMPS output) 
%                    or little-endian ('ieee-le', WINDOWS PC)
% nz = number of z levels in profile
% t1 = beginning timestep to look at
% t2 = ending timestep to look at
% varnum = which variable to read (see .ctl file)
% numvars = total number of variables in grads file (see .ctl file)

% open GrADS file
fid = fopen(filename,'r',MachineFormat);

num_timesteps = (t2-t1) + 1;

% Read in and average profiles over all timesteps
avg_field = zeros(nz,1);
for t=t1:t2
   % 4 bytes per 32 bit float
   byte_position = 4*( (varnum-1)*nz+numvars*nz*(t-1) );
   status = fseek(fid,byte_position,'bof');
   field = fread(fid,nz,'float32');
   avg_field = avg_field + field;
end
avg_field = avg_field/num_timesteps;


% close GrADS file
status = fclose(fid);
