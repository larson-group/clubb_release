function[] = netcdf_rico_writer_2();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% netcdf_rico_writer_2()
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
% Specifically these are the time series of scalars (file 2).
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
tmax = 4320; % must be a multiple of 5, if you change it
t = 0:5:tmax; % this is why it must be a multiple of 5
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

% surface file
[sfc_filename,sfc_nz,sfc_z,sfc_ntimesteps,sfc_numvars,sfc_list_vars] = header_read([scm_path,sfcfile]);

for i=1:sfc_numvars
    for timestep = 1:sizet-1
        stringtoeval = [sfc_list_vars(i,:), ' = read_grads_hoc_endian([scm_path,sfc_filename],''ieee-le'',sfc_nz,t(timestep),t(timestep+1),i,sfc_numvars);'];
        eval(stringtoeval)
        str = sfc_list_vars(i,:);
        eval([strtrim(str),'_array(timestep) = ',str,';']);
%        for j=1:sfc_nz
%            arraydata(j,timestep) = eval([str,'(j)']);
%        end
%        eval([strtrim(str),'_array = arraydata;']);
    end
end

% Set up constants
g0 = 9.8
p0 = 1e5;
R  = 287.04;
Cp = 1004.67;
Lv = 2.5e6;
Cm = 0.001229; % per RICO 3D spec

% Set up derived variables
theta_array = thlm_array + (Lv/Cp).*rcm_array;
T_array = theta_array .* exner_array;
rho_array = p_array ./ (R .* T_array);
smf_array = (um_array(1,:) .* vm_array(1,:)) .* Cm;

tke_array = sum(rhom_array .* em_array * (w_z(2) - w_z(1)));

qtm_array = rtm_array ./ (1 + rtm_array);
qvm_array = (rtm_array - rcm_array) ./ (1 + (rtm_array - rcm_array));

qvm_2000_array = 0.5 * (qvm_array(51,:) + qvm_array(52,:));
qvm_1500_array = qvm_array(39,:);
qvm_1000_array = 0.5 * (qvm_array(26,:) + qvm_array(27,:));
qvm_500_array  = qvm_array(14,:);

thm_2000_array = 0.5 * (theta_array(51,:) + theta_array(52,:));
thm_1500_array = theta_array(39,:);
thm_1000_array = 0.5 * (theta_array(26,:) + theta_array(27,:));
thm_500_array  = theta_array(14,:);

Kht_1250_array = 0.25 * (Kht_array(32,:) + 3*Kht_array(33,:));
Kht_500_array  = Kht_array(14,:);
Kht_300_array  = Kht_array(9,:);

% solve for cloud base
for timestep = 1:(sizet-1)
    zcb_unset = 1;
    for k = 1:nz
        if ((zcb_unset) & (cf_array(k,timestep) > 0.))
            zcb(timestep) = z(k);
            zcb_unset = 0;
        end
    end
end

% solve for cloud top
for timestep = 1:(sizet-1)
    ztop_unset = 1;
    for k = nz:-1:1
        if ((ztop_unset) & (cf_array(k,timestep) > 0.))
            ztop(timestep) = z(k);
            ztop_unset = 0;
        end
    end
end

% solve for cf max
for timestep = 1:(sizet-1)
    cf_max = 0.0
    cf_z   = 0.0
    for k = 1:nz
        if ((cf_array(k,timestep) > cf_max))
            cf_max = cf_array(k,timestep);
            cf_z   = z(k);
        end
    end
    cf_max_z(timestep) = cf_z;
    cf_max_cf(timestep) = cf_max;
end

% Create NetCDF files

[ncid,status] = mexnc('create', 'rico_file2.nc', nc_clobber_mode)
status = mexnc('close', ncid)

[ncid,status] = mexnc('open', 'rico_file2.nc', nc_write_mode)
status = mexnc('redef',ncid)

% Define dimension

[tdimid,status] = mexnc('def_dim',ncid,'t',tmax/5)

% Define variables

[tvarid,status] = mexnc('def_var',ncid,'t','float',1,[tdimid])
[zcbvarid,status] = mexnc('def_var',ncid,'zcb','float',1,[tdimid])
[ztopvarid,status] = mexnc('def_var',ncid,'ztop','float',1,[tdimid])
[zmaxcfracvarid,status] = mexnc('def_var',ncid,'zmaxcfrac','float',1,[tdimid])

[lwpvarid,status] = mexnc('def_var',ncid,'lwp','float',1,[tdimid])
[ccvarid,status] = mexnc('def_var',ncid,'cc','float',1,[tdimid])

[shfvarid,status] = mexnc('def_var',ncid,'shf','float',1,[tdimid])
[lhfvarid,status] = mexnc('def_var',ncid,'lhf','float',1,[tdimid])
[smfvarid,status] = mexnc('def_var',ncid,'smf','float',1,[tdimid])
[tkevarid,status] = mexnc('def_var',ncid,'tke','float',1,[tdimid])

[v_lowlevelvarid,status] = mexnc('def_var',ncid,'v_lowlevel','float',1,[tdimid])
[qv_lowlevelvarid,status] = mexnc('def_var',ncid,'qv_lowlevel','float',1,[tdimid])
[th_lowlevelvarid,status] = mexnc('def_var',ncid,'th_lowlevel','float',1,[tdimid])

[prec_2000varid,status] = mexnc('def_var',ncid,'prec_2000','float',1,[tdimid])
[prec_1500varid,status] = mexnc('def_var',ncid,'prec_1500','float',1,[tdimid])
[prec_1000varid,status] = mexnc('def_var',ncid,'prec_1000','float',1,[tdimid])
[prec_500varid,status]  = mexnc('def_var',ncid,'prec_500','float',1,[tdimid])
[prec_srfvarid,status]  = mexnc('def_var',ncid,'prec_srf','float',1,[tdimid])

[qv_2000varid,status] = mexnc('def_var',ncid,'qv_2000','float',1,[tdimid])
[qv_1500varid,status] = mexnc('def_var',ncid,'qv_1500','float',1,[tdimid])
[qv_1000varid,status] = mexnc('def_var',ncid,'qv_1000','float',1,[tdimid])
[qv_500varid,status]  = mexnc('def_var',ncid,'qv_500','float',1,[tdimid])

[th_2000varid,status] = mexnc('def_var',ncid,'th_2000','float',1,[tdimid])
[th_1500varid,status] = mexnc('def_var',ncid,'th_1500','float',1,[tdimid])
[th_1000varid,status] = mexnc('def_var',ncid,'th_1000','float',1,[tdimid])
[th_500varid,status]  = mexnc('def_var',ncid,'th_500','float',1,[tdimid])

[Kh_1250varid,status] = mexnc('def_var',ncid,'Kh_1250','float',1,[tdimid])
[Kh_500varid,status]  = mexnc('def_var',ncid,'Kh_500','float',1,[tdimid])
[Kh_300varid,status]  = mexnc('def_var',ncid,'Kh_300','float',1,[tdimid])

status = mexnc('end_def',ncid)

% Write data

for i=1:sizet-1
    % k-1 comes from NetCDF starting variables at 0 and MATLAB starting
    % them at 1.
    % start, count, value, autoscale
    status = mexnc('VARPUT',ncid,tvarid,i-1,1,t(i+1)*60,0)
    status = mexnc('VARPUT',ncid,zcbvarid,i-1,1,zcb(i),0)
    status = mexnc('VARPUT',ncid,ztopvarid,i-1,1,ztop(i),0)
    status = mexnc('VARPUT',ncid,zmaxcfracvarid,i-1,1,cf_max_z(i),0)
    status = mexnc('VARPUT',ncid,lwpvarid,i-1,1,lwp_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,ccvarid,i-1,1,cc_array(i),0)
    status = mexnc('VARPUT',ncid,shfvarid,i-1,1,sh_array(i),0)
    status = mexnc('VARPUT',ncid,lhfvarid,i-1,1,lh_array(i),0)
    status = mexnc('VARPUT',ncid,smfvarid,i-1,1,smf_array(i),0)
    status = mexnc('VARPUT',ncid,tkevarid,i-1,1,tke_array(i),0)
    status = mexnc('VARPUT',ncid,v_lowlevelvarid,i-1,1,vm_array(2,i),0)
    status = mexnc('VARPUT',ncid,qv_lowlevelvarid,i-1,1,qvm_array(2,i)*1000,0)
    status = mexnc('VARPUT',ncid,th_lowlevelvarid,i-1,1,theta_array(2,i),0)
    status = mexnc('VARPUT',ncid,prec_2000varid,i-1,1,28.94*0.5*(rain_rate_array(51,i) + rain_rate_array(52,i)))
    status = mexnc('VARPUT',ncid,prec_1500varid,i-1,1,28.94*rain_rate_array(39,i),0)
    status = mexnc('VARPUT',ncid,prec_1000varid,i-1,1,28.94*0.5*(rain_rate_array(26,i) + rain_rate_array(27,i)),0)
    status = mexnc('VARPUT',ncid,prec_500varid,i-1,1,28.94*rain_rate_array(14,i),0)
    status = mexnc('VARPUT',ncid,prec_srfvarid,i-1,1,28.94*rain_rate_array(2,i),0) % 28.94 W m^-2 = 1 mm day^-1
    status = mexnc('VARPUT',ncid,qv_2000varid,i-1,1,qvm_2000_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,qv_1500varid,i-1,1,qvm_1500_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,qv_1000varid,i-1,1,qvm_1000_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,qv_500varid,i-1,1,qvm_500_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,th_2000varid,i-1,1,thm_2000_array(i),0)
    status = mexnc('VARPUT',ncid,th_1500varid,i-1,1,thm_1500_array(i),0)
    status = mexnc('VARPUT',ncid,th_1000varid,i-1,1,thm_1000_array(i),0)
    status = mexnc('VARPUT',ncid,th_500varid,i-1,1,thm_500_array(i),0)
    status = mexnc('VARPUT',ncid,Kh_1250varid,i-1,1,Kht_1250_array(i),0)
    status = mexnc('VARPUT',ncid,Kh_500varid,i-1,1,Kht_500_array(i),0)
    status = mexnc('VARPUT',ncid,Kh_300varid,i-1,1,Kht_300_array(i),0)
end

%% File inquiry, if you want to verify that everything got written
%% properly.  I usually leave this commented.
%[ndims,nvars,ngatts,unlimdim,status] = mexnc('inq',ncid)

% Close file
status = mexnc('close',ncid)