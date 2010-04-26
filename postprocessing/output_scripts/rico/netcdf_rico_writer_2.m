function[] = netcdf_rico_writer_2();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% netcdf_rico_writer_2()
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
%scm_path = ['/home/mjfalk/hoc_results/rico/rico0120/'];
% Revised submission 29 August 2007
%scm_path = ['/home/mjfalk/hoc_results/rico/rico0206/']
% Revised submission 23 April 2010 --ldgrant
scm_path = ['/home/ldgrant/RICO_submission/Apr2010/RICO_grads_files/'];
smfile   = 'rico_zt.ctl';
swfile   = 'rico_zm.ctl';
sfcfile  = 'rico_sfc.ctl';

% Read RICO data from GrADS files
t = 1:864; % Note: code changed because RICO data is now output to GrADS every 5 minutes
sizet = size(t);
sizet = max(sizet);

% mass (zt) file
[filename,nz,z,ntimesteps,numvars,list_vars] = header_read([scm_path,smfile]);
for i=1:numvars
    i
    % Note: RICO data is now output to GrADS every 5 minutes.
    for timestep = 1:sizet
        stringtoeval = [list_vars(i,:), ' = read_grads_hoc_endian([scm_path,filename],''ieee-le'',nz,t(timestep),t(timestep),i,numvars);'];
        eval(stringtoeval)
        str = list_vars(i,:);
        for j=1:nz
            arraydata(j,timestep) = eval([str,'(j)']);
        end
        eval([strtrim(str),'_array = arraydata;']);
    end
end

% momentum (zm) file
[w_filename,w_nz,w_z,w_ntimesteps,w_numvars,w_list_vars] = header_read([scm_path,swfile]);

for i=1:w_numvars
    i
    % Note: RICO data is now output to GrADS every 5 minutes.
    for timestep = 1:sizet
        stringtoeval = [w_list_vars(i,:), ' = read_grads_hoc_endian([scm_path,w_filename],''ieee-le'',w_nz,t(timestep),t(timestep),i,w_numvars);'];
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
    i
    % Note: RICO data is now output to GrADS every 5 minutes.
    for timestep = 1:sizet
        stringtoeval = [sfc_list_vars(i,:), ' = read_grads_hoc_endian([scm_path,sfc_filename],''ieee-le'',sfc_nz,t(timestep),t(timestep),i,sfc_numvars);'];
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
g0 = 9.8;
p0 = 1e5;
R  = 287.04;
Cp = 1004.67;
Lv = 2.5e6;
Cm = 0.001229; % per RICO 3D spec

% Set up derived variables
theta_array = thlm_array + (Lv/Cp).*rcm_array;
T_array = theta_array .* exner_array;
rho_array = p_in_Pa_array ./ (R .* T_array);
smf_array = (um_array(1,:) .* vm_array(1,:)) .* Cm;

% tke for RICO spec. = Vertically integrated rho*TKE [kg/s^2]
for k = 2:w_nz
    tke_2d_array(k-1,:) = rho_zm_array(k,:) .* em_array(k,:) * ( w_z(k) - w_z(k-1) );
end
tke_array = sum(tke_2d_array);

qtm_array = rtm_array ./ (1 + rtm_array);
qvm_array = (rtm_array - rcm_array) ./ (1 + (rtm_array - rcm_array));

% RICO spec. wants output at specific levels.  E.g. qvm_2000_array is qvm at 2000m.
% Because RICO now uses the 128 level stretched grid, there are 33 total levels below
% 3960m, the RICO max altitude.  The index represents the level closest to the specified
% level; e.g. qvm_array(22,:) accesses qvm at level 22, which corresponds to 1952.63 m,
% the closest level to 2000m on the 128 level stretched grid.
qvm_2000_array = qvm_array(22,:);
qvm_1500_array = qvm_array(19,:);
qvm_1000_array = qvm_array(15,:);
qvm_500_array  = qvm_array(10,:);

thm_2000_array = theta_array(22,:);
thm_1500_array = theta_array(19,:);
thm_1000_array = theta_array(15,:);
thm_500_array  = theta_array(10,:);

Kh_1250_array = Kh_zt_array(17,:);
Kh_500_array  = Kh_zt_array(10,:);
Kh_300_array  = Kh_zt_array(7,:);

% solve for cloud base
% Note: RICO data is now output to GrADS every 5 minutes.
for timestep = 1:(sizet)
    zcb_unset = 1;
    for k = 1:nz
        if ((zcb_unset) & (cloud_frac_array(k,timestep) > 0.))
            zcb(timestep) = z(k);
            zcb_unset = 0;
        end
    end
end

% solve for cloud top
% Note: RICO data is now output to GrADS every 5 minutes.
for timestep = 1:(sizet)
    ztop_unset = 1;
    for k = nz:-1:1
        if ((ztop_unset) & (cloud_frac_array(k,timestep) > 0.))
            ztop(timestep) = z(k);
            ztop_unset = 0;
        end
    end
end

% solve for cloud_frac max
% Note: RICO data is now output to GrADS every 5 minutes.
for timestep = 1:(sizet)
    cloud_frac_max = 0.0
    cloud_frac_z   = 0.0
    for k = 1:nz
        if ((cloud_frac_array(k,timestep) > cloud_frac_max))
            cloud_frac_max = cloud_frac_array(k,timestep);
            cloud_frac_z   = z(k);
        end
    end
    cloud_frac_max_z(timestep) = cloud_frac_z;
    cloud_frac_max_cloud_frac(timestep) = cloud_frac_max;
end

% Create NetCDF files

[ncid,status] = mexnc('create', 'rico_file2.nc', nc_clobber_mode)
status = mexnc('close', ncid)

[ncid,status] = mexnc('open', 'rico_file2.nc', nc_write_mode)
status = mexnc('redef',ncid)

% Write global attributes

% syntax: status = mexnc('ATTPUT', cdfid, varid, 'name', datatype, len, value)
status = mexnc('ATTPUT',ncid,'NC_GLOBAL','Contact_Person','CHAR',-1,'Vince Larson (vlarson@uwm.edu), Michael Falk (mjf@e-falk.com), and Leah Grant (ldgrant@uwm.edu)')
status = mexnc('ATTPUT',ncid,'NC_GLOBAL','Vertical_Resolution','CHAR',-1,'128-level, 27.5 km stretched grid')
status = mexnc('ATTPUT',ncid,'NC_GLOBAL','Model_Timestep','CHAR',-1,'5 minutes')
status = mexnc('ATTPUT',ncid,'NC_GLOBAL','Run_type','CHAR',-1,'composite')

% Define dimension

[tdimid,status] = mexnc('def_dim',ncid,'time',sizet)

% Define variables

[tvarid,status] = mexnc('def_var',ncid,'time','float',1,[tdimid])
[zcbvarid,status] = mexnc('def_var',ncid,'zcb','float',1,[tdimid])
[ztopvarid,status] = mexnc('def_var',ncid,'ztop','float',1,[tdimid])
[zmaxcfracvarid,status] = mexnc('def_var',ncid,'zmaxcfrac','float',1,[tdimid])
[m_basevarid] = mexnc('def_var',ncid,'M_base','float',1,[tdimid])

[lwpvarid,status] = mexnc('def_var',ncid,'LWP','float',1,[tdimid])
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

[Kh_300varid,status]  = mexnc('def_var',ncid,'Kh_300','float',1,[tdimid])
[Kh_500varid,status]  = mexnc('def_var',ncid,'Kh_500','float',1,[tdimid])
[Kh_1250varid,status] = mexnc('def_var',ncid,'Kh_1250','float',1,[tdimid])

status = mexnc('end_def',ncid)

% Write data

% Note: RICO data is now output to GrADS every 5 minutes.
for i=1:sizet
    % k-1 comes from NetCDF starting variables at 0 and MATLAB starting
    % them at 1.
    % start, count, value, autoscale
    status = mexnc('VARPUT',ncid,tvarid,i-1,1,t(i)*300,0)
    status = mexnc('VARPUT',ncid,zcbvarid,i-1,1,zcb(i),0)
    status = mexnc('VARPUT',ncid,ztopvarid,i-1,1,ztop(i),0)
    status = mexnc('VARPUT',ncid,zmaxcfracvarid,i-1,1,cloud_frac_max_z(i),0)
% Output zeroes for M_base (Mass flux at cloud base) because CLUBB does not use a mass flux scheme like other models.
    status = mexnc('VARPUT',ncid,m_basevarid,i-1,1,0,0)
    status = mexnc('VARPUT',ncid,lwpvarid,i-1,1,lwp_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,ccvarid,i-1,1,cc_array(i),0)
    status = mexnc('VARPUT',ncid,shfvarid,i-1,1,sh_array(i),0)
    status = mexnc('VARPUT',ncid,lhfvarid,i-1,1,lh_array(i),0)
    status = mexnc('VARPUT',ncid,smfvarid,i-1,1,smf_array(i),0)
    status = mexnc('VARPUT',ncid,tkevarid,i-1,1,tke_array(i),0)
    status = mexnc('VARPUT',ncid,v_lowlevelvarid,i-1,1,vm_array(2,i),0)
    status = mexnc('VARPUT',ncid,qv_lowlevelvarid,i-1,1,qvm_array(2,i)*1000,0)
    status = mexnc('VARPUT',ncid,th_lowlevelvarid,i-1,1,theta_array(2,i),0)
% rain_rate * 28.94 is virtually the same as Fprec but on the mass grid
    status = mexnc('VARPUT',ncid,prec_2000varid,i-1,1,28.94*rain_rate_array(22,i),0)
    status = mexnc('VARPUT',ncid,prec_1500varid,i-1,1,28.94*rain_rate_array(19,i),0)
    status = mexnc('VARPUT',ncid,prec_1000varid,i-1,1,28.94*rain_rate_array(15,i),0)
    status = mexnc('VARPUT',ncid,prec_500varid,i-1,1,28.94*rain_rate_array(10,i),0)
    status = mexnc('VARPUT',ncid,prec_srfvarid,i-1,1,28.94*rain_rate_array(2,i),0) % 28.94 W m^-2 = 1 mm day^-1
    status = mexnc('VARPUT',ncid,qv_2000varid,i-1,1,qvm_2000_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,qv_1500varid,i-1,1,qvm_1500_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,qv_1000varid,i-1,1,qvm_1000_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,qv_500varid,i-1,1,qvm_500_array(i)*1000,0)
    status = mexnc('VARPUT',ncid,th_2000varid,i-1,1,thm_2000_array(i),0)
    status = mexnc('VARPUT',ncid,th_1500varid,i-1,1,thm_1500_array(i),0)
    status = mexnc('VARPUT',ncid,th_1000varid,i-1,1,thm_1000_array(i),0)
    status = mexnc('VARPUT',ncid,th_500varid,i-1,1,thm_500_array(i),0)
% Output zeroes for eddy diffusivity because RICO spec. wants eddy diffusivity for heat, but CLUBB
% uses a different parameterization for w'theta_l' and only has background eddy diffusivity for heat.
    status = mexnc('VARPUT',ncid,Kh_300varid,i-1,1,0,0)
    status = mexnc('VARPUT',ncid,Kh_500varid,i-1,1,0,0)
    status = mexnc('VARPUT',ncid,Kh_1250varid,i-1,1,0,0)
end

%% File inquiry, if you want to verify that everything got written properly.
[ndims,nvars,ngatts,unlimdim,status] = mexnc('inq',ncid)

% Close file
status = mexnc('close',ncid)
