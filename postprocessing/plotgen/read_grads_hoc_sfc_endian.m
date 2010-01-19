
function avg_field = read_grads_hoc_sfc_endian(filename,MachineFormat,nz,t1,t2,varnum,numvars)

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

%Ensure the file will be closed no matter what happens
cleanupHandler = onCleanup(@()fclose(fid));

%Preallocate array for speed
avg_field(t1:t2) = 0.0;

% Read in and average profiles over all timesteps
for t=t1:t2
   % 4 bytes per 32 bit float
   byte_position = 4*( (varnum-1)*nz+numvars*nz*(t-1) );
   status = fseek(fid,byte_position,'bof');
   field = fread(fid,nz,'float32');
   avg_field(t) = field;
end
