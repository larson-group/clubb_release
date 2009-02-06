function[] = netcdf_rico_writer_4();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% netcdf_rico_writer_4()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, February-March 2007
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
% Specifically these are the time series of large-scale forcings (file 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath '/home/mjfalk/mexnc2/mexnc/'
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end


clear % because this requires LOTS of memory.
% File setup
% Original submission
%scm_path = ['/home/mjfalk/hoc_working/standalone/rico0104/'];
% Revised submission 7 March 2007
scm_path = ['/home/mjfalk/hoc_results/rico/rico0119/'];
smfile   = 'rico_zt.ctl';
swfile   = 'rico_zm.ctl';
sfcfile  = 'rico_sfc.ctl';

% Read RICO data from GrADS files
tmax = 4320; % should be a multiple of 60
t = 0:60:tmax; % this is why
t(1) = 1;
sizet = size(t);
sizet = max(sizet);

%mass file
[filename,nz,z,ntimesteps,numvars,list_vars] = header_read([scm_path,smfile]);

for i=1:numvars
    for timestep = 1:sizet-1
        stringtoeval = [list_vars(i,:), ' = read_grads_hoc_endian([scm_path,filename],''ieee-le'',nz,t(timestep),t(timestep+1),i,numvars);'];
        eval(stringtoeval)
        str = list_vars(i,:);
        for j=1:nz
            arraydata(j,timestep) = eval([str,'(j)']);
        end
        eval([strtrim(str),'_array = arraydata;']);
    end
end

% Set up constants
g0 = 9.8;
p0 = 1e5;
R  = 287.04;
Cp = 1004.67;
Lv = 2.5e6;
f = 4.506704e-5; % (2 omega sin 18 degrees)

% Set up derived variables
theta_array = thlm_array + (Lv/Cp).*rcm_array;
T_array = theta_array .* exner_array;
rho_array = p_array ./ (R .* T_array);

% centered difference of thlm
for k=2:nz-1
    dthlm_dz(k,:) = (thlm_array(k+1,:) - thlm_array(k-1,:)) ./ (z(k+1) - z(k-1));
end
dthlm_dz(1,:) = (thlm_array(2,:) - thlm_array(1,:)) ./ (z(2) - z(1));
dthlm_dz(nz,:) = (thlm_array(nz,:) - thlm_array(nz-1,:)) ./ (z(nz) - z(nz-1));

t_subsidence = -wm_array .* dthlm_dz .* exner_array;
t_forcing = (-2.51 / 86400) + ((-2.18 + 2.51 ) / (86400 * 4000) .* z); ... % due to radiation and advection

q_forcing = rtm_f_array ./ (1 + rtm_f_array);


% Create NetCDF files

[ncid,status] = mexnc('create', 'rico_file4.nc', nc_clobber_mode)
status = mexnc('close', ncid)

[ncid,status] = mexnc('open', 'rico_file4.nc', nc_write_mode)

% Define dimensions

status = mexnc('redef',ncid)
[tdimid,status] = mexnc('def_dim',ncid,'t',(tmax/60))
[zfdimid,status] = mexnc('def_dim',ncid,'zf',nz)

% Define variables

[zfvarid,status] = mexnc('def_var',ncid,'zf','float',2,[tdimid zfdimid])
[dTdt_varid,status] = mexnc('def_var',ncid,'dTdt_ls','float',2,[tdimid zfdimid])
[dqdt_varid,status] = mexnc('def_var',ncid,'dqdt_ls','float',2,[tdimid zfdimid])
[dudt_varid,status] = mexnc('def_var',ncid,'dudt_ls','float',2,[tdimid zfdimid])
[dvdt_varid,status] = mexnc('def_var',ncid,'dvdt_ls','float',2,[tdimid zfdimid])

status = mexnc('end_def',ncid)

% Write data

for k=1:nz
    for i=1:sizet-1
    % k-1 comes from NetCDF starting variables at 0 and MATLAB starting
    % them at 1.
    % start, count, value, autoscale
        status = mexnc('VARPUT',ncid,zfvarid,[i-1 k-1],[1 1],z(k),0)
        status = mexnc('VARPUT',ncid,dTdt_varid,[i-1 k-1],[1 1],t_forcing(k) + t_subsidence(k,i),0)
        status = mexnc('VARPUT',ncid,dqdt_varid,[i-1 k-1],[1 1],q_forcing(k,i)*1000,0)
        status = mexnc('VARPUT',ncid,dudt_varid,[i-1 k-1],[1 1],-f*(ug_array(k,i) - um_array(k,i)),0)
        status = mexnc('VARPUT',ncid,dvdt_varid,[i-1 k-1],[1 1],f*(vg_array(k,i) - vm_array(k,i)),0)
    end
end

%% File inquiry, if you want to verify that everything got written
%% properly.  I usually leave this commented.
%[ndims,nvars,ngatts,unlimdim,status] = mexnc('inq',ncid)

% Close file
status = mexnc('close',ncid)