function[] = netcdf_rico_writer_3();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% netcdf_rico_writer_3()
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
% Specifically these are the time series of mean profiles (file 3).
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
t = 0:60:tmax; % this is why it must be a multiple of 60
t(1) = 1;
sizet = size(t);
sizet = max(sizet);

% mass file
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

% momentum file
[w_filename,w_nz,w_z,w_ntimesteps,w_numvars,w_list_vars] = header_read([scm_path,swfile]);

for i=1:w_numvars
    for timestep = 1:sizet-1
        stringtoeval = [w_list_vars(i,:), ' = read_grads_hoc_endian([scm_path,w_filename],''ieee-le'',w_nz,t(timestep),t(timestep+1),i,w_numvars);'];
        eval(stringtoeval)
        str = w_list_vars(i,:);
        for j=1:w_nz
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

% Set up derived variables
theta_array = thlm_array + (Lv/Cp).*rcm_array;
T_array = theta_array .* exner_array;
rho_array = p_array ./ (R .* T_array);

qtm_array = rtm_array ./ (1 + rtm_array);
qsm_array = rsm_array ./ (1 + rsm_array);
qcm_array = rcm_array ./ (1 + rcm_array);
qrm_array = rrm_array ./ (1 + rrm_array);

wpqtp_array = wprtp_array ./ (1 + rtm_array);

% Create NetCDF files

[ncid,status] = mexnc('create', 'rico_file3.nc', nc_clobber_mode)
status = mexnc('close', ncid)

[ncid,status] = mexnc('open', 'rico_file3.nc', nc_write_mode)

% Define dimensions

status = mexnc('redef',ncid)
[tdimid,status] = mexnc('def_dim',ncid,'t',(tmax/60))
[zfdimid,status] = mexnc('def_dim',ncid,'zf',nz)
[zhdimid,status] = mexnc('def_dim',ncid,'zh',w_nz)

% Define variables

[zfvarid,status] = mexnc('def_var',ncid,'zf','float',2,[tdimid zfdimid])
[zhvarid,status] = mexnc('def_var',ncid,'zh','float',2,[tdimid zhdimid])
[uvarid,status] = mexnc('def_var',ncid,'u','float',2,[tdimid zfdimid])
[vvarid,status] = mexnc('def_var',ncid,'v','float',2,[tdimid zfdimid])
[thetalvarid,status] = mexnc('def_var',ncid,'thetal','float',2,[tdimid zfdimid])
[qtvarid,status] = mexnc('def_var',ncid,'qt','float',2,[tdimid zfdimid])
[qsvarid,status] = mexnc('def_var',ncid,'qs','float',2,[tdimid zfdimid])
[qlvarid,status] = mexnc('def_var',ncid,'ql','float',2,[tdimid zfdimid])
[qrvarid,status] = mexnc('def_var',ncid,'qr','float',2,[tdimid zfdimid])
[cfvarid,status] = mexnc('def_var',ncid,'cf','float',2,[tdimid zfdimid])
[rhovarid,status] = mexnc('def_var',ncid,'rho','float',2,[tdimid zfdimid])
[wthlvarid,status] = mexnc('def_var',ncid,'wthl','float',2,[tdimid zhdimid])
[wqtvarid,status] = mexnc('def_var',ncid,'wqt','float',2,[tdimid zhdimid])
[uwvarid,status] = mexnc('def_var',ncid,'uw','float',2,[tdimid zhdimid])
[vwvarid,status] = mexnc('def_var',ncid,'vw','float',2,[tdimid zhdimid])
[Mfvarid,status] = mexnc('def_var',ncid,'Mf','float',2,[tdimid zhdimid])
[precvarid,status] = mexnc('def_var',ncid,'prec','float',2,[tdimid zhdimid])

[w_upvarid,status] = mexnc('def_var',ncid,'w_up','float',2,[tdimid zfdimid])
[thl_upvarid,status] = mexnc('def_var',ncid,'thl_up','float',2,[tdimid zfdimid])
[qt_upvarid,status] = mexnc('def_var',ncid,'qt_up','float',2,[tdimid zfdimid])
[ql_upvarid,status] = mexnc('def_var',ncid,'ql_up','float',2,[tdimid zfdimid])
[thv_upvarid,status] = mexnc('def_var',ncid,'thv_up','float',2,[tdimid zfdimid])


status = mexnc('end_def',ncid)

% Write data

for k=1:nz
    for i=1:sizet-1
    % k-1 comes from NetCDF starting variables at 0 and MATLAB starting
    % them at 1.
    % start, count, value, autoscale
        status = mexnc('VARPUT',ncid,zfvarid,[i-1 k-1],[1 1],z(k),0);
        status = mexnc('VARPUT',ncid,zhvarid,[i-1 k-1],[1 1],w_z(k),0);
        status = mexnc('VARPUT',ncid,uvarid,[i-1 k-1],[1 1],um_array(k,i),0);
        status = mexnc('VARPUT',ncid,vvarid,[i-1 k-1],[1 1],vm_array(k,i),0);
        status = mexnc('VARPUT',ncid,thetalvarid,[i-1 k-1],[1 1],thlm_array(k,i),0);
        status = mexnc('VARPUT',ncid,qtvarid,[i-1 k-1],[1 1],qtm_array(k,i)*1000,0);
        status = mexnc('VARPUT',ncid,qsvarid,[i-1 k-1],[1 1],qsm_array(k,i)*1000,0);
        status = mexnc('VARPUT',ncid,qlvarid,[i-1 k-1],[1 1],qcm_array(k,i)*1000,0);
        status = mexnc('VARPUT',ncid,qrvarid,[i-1 k-1],[1 1],qrm_array(k,i)*1000,0);

        status = mexnc('VARPUT',ncid,cfvarid,[i-1 k-1],[1 1],cf_array(k,i),0);
        status = mexnc('VARPUT',ncid,rhovarid,[i-1 k-1],[1 1],rho_array(k,i),0);
        status = mexnc('VARPUT',ncid,wthlvarid,[i-1 k-1],[1 1],wpthlp_array(k,i),0);
        status = mexnc('VARPUT',ncid,wqtvarid,[i-1 k-1],[1 1],wpqtp_array(k,i),0);
        status = mexnc('VARPUT',ncid,uwvarid,[i-1 k-1],[1 1],upwp_array(k,i)*rhom(k),0);
        status = mexnc('VARPUT',ncid,vwvarid,[i-1 k-1],[1 1],vpwp_array(k,i)*rhom(k),0);
        status = mexnc('VARPUT',ncid,Mfvarid,[i-1 k-1],[1 1],0,0);
        status = mexnc('VARPUT',ncid,precvarid,[i-1 k-1],[1 1],28.94*rain_rate_array(k,i),0);

        status = mexnc('VARPUT',ncid,w_upvarid,[i-1 k-1],[1 1],0,0); % These are all zeros.
        status = mexnc('VARPUT',ncid,thl_upvarid,[i-1 k-1],[1 1],0,0);
        status = mexnc('VARPUT',ncid,qt_upvarid,[i-1 k-1],[1 1],0,0);
        status = mexnc('VARPUT',ncid,ql_upvarid,[i-1 k-1],[1 1],0,0);
        status = mexnc('VARPUT',ncid,thv_upvarid,[i-1 k-1],[1 1],0,0);
    end
end

%% File inquiry, if you want to verify that everything got written
%% properly.  I usually leave this commented.
%[ndims,nvars,ngatts,unlimdim,status] = mexnc('inq',ncid)

% Close file
status = mexnc('close',ncid)
