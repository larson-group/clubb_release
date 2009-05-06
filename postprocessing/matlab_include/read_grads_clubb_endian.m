
function avg_field = read_grads_clubb_endian(filename,MachineFormat,nz,t1,t2,varnum,numvars)

% READ_GRADS_CLUBB_ENDIAN   Reads and time-averages profiles from 1D GrADS *.dat files.
%
%   Written by Vince Larson (vlarson@uwm.edu)
%
%   This program access GrADS *.dat files and obtains a single-time or time-averaged profile.
%   The resulting profile is passed to the parent program for manipulation or observation.
%
% ------------------------------------------------------------------------------
%
%   INPUTS
%   ------
%   filename      -> Full filepath to the GRADS data file
%   MachineFormat -> String that specifies big-endian ('ieee-be', Linux COAMPS output) 
%                    or little-endian ('ieee-le', WINDOWS PC)
%   nz            -> Number of z levels in profile
%   t1            -> Beginning timestep to look at
%   t2            -> Ending timestep to look at
%   varnum        -> Which variable to read (see .ctl file)
%   numvars       -> Total number of variables in grads file (see .ctl file)
%
%   OUTPUTS
%   -------
%   avg_field -> GrADS profile of the variable we wish to examine

% END OF DOCUMENTATION
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

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

% End of read_grads_clubb_endian program...Returning to parent program
