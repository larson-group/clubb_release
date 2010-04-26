function[] = netcdf_rico_writer_1();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% netcdf_rico_writer_1()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, February-March 2007
% Modified by Leah Grant, April 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters: none
% Input files: RICO GrADS files (.ctl and .dat)
% Output parameters: none
% Output files: RICO NetCDF files
%
% Requires: mexnc and netcdf_toolbox (installed on condella)
%
% Reference: http://www.knmi.nl/samenw/rico/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program creates the required NetCDF output for the RICO case.
% Specifically these are the profiles at the initial time (file 1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath '/home/mjfalk/mexnc2/mexnc/'
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end


% File setup
% Original submission
%scm_path = ['/home/mjfalk/hoc_working/standalone/rico0104/'];
% Revised submission 7 March 2007
%scm_path = ['/home/mjfalk/hoc_results/rico/rico0120/'];
% Revised submission 29 August 2007
%scm_path = ['/home/mjfalk/hoc_results/rico/rico0206/'];
% Revised submission 23 April 2010 --ldgrant
scm_path = ['/home/ldgrant/RICO_submission/Apr2010/RICO_grads_files/'];
smfile   = 'rico_zt.ctl';
swfile   = 'rico_zm.ctl';


% Read RICO data from GrADS files
t = 1:864; % RICO data is now output to GrADS every 5 mins
sizet = size(t);
sizet = max(sizet);

[filename,nz,z,ntimesteps,numvars,list_vars] = header_read([scm_path,smfile]);

% instantaneous "average", for initial time
for i=1:numvars
    timestep = 1;
    stringtoeval = [list_vars(i,:), ' = read_grads_hoc_endian([scm_path,filename],''ieee-le'',nz,t(timestep),t(timestep),i,numvars);'];
    eval(stringtoeval)
    str = list_vars(i,:);
    for j=1:nz
        arraydata(j,timestep) = eval([str,'(j)']);
    end
    eval([strtrim(str),'_array = arraydata;']);
end

% Set up constants

g0 = 9.8;
p0 = 1e5;
R  = 287.04;
Cp = 1004.67;
Lv = 2.5e6;

% Set up derived variables
qtm_array = rtm_array ./ (1 + rtm_array);
theta_array = thlm_array + (Lv/Cp).*rcm_array;
T_array = theta_array .* exner_array;

% Create NetCDF files

[ncid,status] = mexnc('create', 'rico_file1.nc', nc_clobber_mode)
status = mexnc('close', ncid)

[ncid,status] = mexnc('open', 'rico_file1.nc', nc_write_mode)
status = mexnc('redef',ncid)

% Write global attributes

% syntax: status = mexnc('ATTPUT', cdfid, varid, 'name', datatype, len, value)
status = mexnc('ATTPUT',ncid,'NC_GLOBAL','Contact_Person','CHAR',-1,'Vince Larson (vlarson@uwm.edu), Michael Falk (mjf@e-falk.com), and Leah Grant (ldgrant@uwm.edu)')
status = mexnc('ATTPUT',ncid,'NC_GLOBAL','Vertical_Resolution','CHAR',-1,'128-level, 27.5 km stretched grid')
status = mexnc('ATTPUT',ncid,'NC_GLOBAL','Model_Timestep','CHAR',-1,'5 minutes')
status = mexnc('ATTPUT',ncid,'NC_GLOBAL','Run_type','CHAR',-1,'composite')
% Define dimension

[zdimid,status] = mexnc('def_dim',ncid,'zf',nz)

% Define variables

[zfvarid,status] = mexnc('def_var',ncid,'zf','float',1,[zdimid])
[uvarid,status] = mexnc('def_var',ncid,'u','float',1,[zdimid])
[vvarid,status] = mexnc('def_var',ncid,'v','float',1,[zdimid])
[thetalvarid,status] = mexnc('def_var',ncid,'thetal','float',1,[zdimid])
[qtvarid,status] = mexnc('def_var',ncid,'qt','float',1,[zdimid])
[rhovarid,status] = mexnc('def_var',ncid,'rho','float',1,[zdimid])

status = mexnc('end_def',ncid)

% Write data

for k=1:nz
    % k-1 comes from NetCDF starting variables at 0 and MATLAB starting
    % them at 1.
    % start, count, value, autoscale
    status = mexnc('VARPUT',ncid,zfvarid,k-1,1,z(k),0)
    status = mexnc('VARPUT',ncid,uvarid,k-1,1,um_array(k,1),0)
    status = mexnc('VARPUT',ncid,vvarid,k-1,1,vm_array(k,1),0)
    status = mexnc('VARPUT',ncid,thetalvarid,k-1,1,thlm_array(k,1),0)
    status = mexnc('VARPUT',ncid,qtvarid,k-1,1,qtm_array(k,1)*1000,0)
    status = mexnc('VARPUT',ncid,rhovarid,k-1,1,rho_array(k,1),0)
end


%% File inquiry, if you want to verify that everything got written properly.
[ndims,nvars,ngatts,unlimdim,status] = mexnc('inq',ncid)

% Close file
status = mexnc('close',ncid)
