function[] = astex_a209_output_creator( scm_path, output_path )
% GABLS3_OUTPUT_CREATOR This function creates netCDF files required by the GABLS 3 intercomparison. It uses CLUBB output files as source information.
%
%   This file is also meant to be an example for future MATLAB scripts to
%   convert CLUBB output for data submission in netCDF format.
%   Essentially a script for such a conversion will break down into these
%   sections:
%
%       File Input Section -
%           This is where the input files that are to be converted are
%           specified and read into MATLAB.
%
%       Definition Section -
%           This is where netCDF definitions will be written to the output
%           file. This includes information such as variable names,
%           variable attributes, and global attributes.
%
%       Conversion Section -
%           Since the input information produced by CLUBB may not match one
%           for one with the specifications of the output file, this
%           section is needed. Here all conversions of information will
%           occur such as converting temperature into potential
%           temperature.
%
%       Output File Section -
%           This section is respondsible for writing variable data to the
%           output file.
%


% Necessary include
addpath '../../matlab_include/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File Input Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define some constants for the case
model_timestep = 300;
average_time_interval = 3600;

curr_case = 'astex_a209';

% zt Grid
smfile   = strcat(curr_case, '_zt.ctl');

% zm Grid
swfile   = strcat(curr_case, '_zm.ctl');

% sfc Grid
sfcfile  = strcat(curr_case, '_sfc.ctl');

% rad zt Grid
radztfile = strcat(curr_case, '_rad_zt.ctl');

% rad zm Grid
radzmfile = strcat(curr_case, '_rad_zm.ctl');

% Reading the header from zt file
[filename,nzmax,z,ntimesteps,numvars,list_vars] = header_read([scm_path,smfile]);

% Used to navigate through time in all of the files. At the time this
% script was written, all of the GrADS output files generated by CLUBB used
% the same timestep.
t = 1:ntimesteps;
sizet = ntimesteps;

% Read in zt file's variables into MATLAB.
% Variables will be usable in the form <GrADS Variable Name>_array.
for i=1:numvars
    for timestep = 1:sizet
    	stringtoeval = [list_vars(i,:), ' = read_grads_clubb_endian([scm_path,filename],''ieee-le'',nzmax,t(timestep),t(timestep),i,numvars);'];
    	eval(stringtoeval);
    	str = list_vars(i,:);
        arraydata(1:nzmax,timestep) = eval([str,'(1:nzmax)']);
    	eval([strtrim(str),'_array = arraydata;']);
    end
    disp(i);
end

% Reading the header from zm file
[w_filename,w_nz,w_z,w_ntimesteps,w_numvars,w_list_vars] = header_read([scm_path,swfile]);

% Read in zm file's variables into MATLAB.
% Variables will be usable in the form <GrADS Variable Name>_array
for i=1:w_numvars
     for timestep = 1:sizet
         stringtoeval = [w_list_vars(i,:), ' = read_grads_clubb_endian([scm_path,w_filename],''ieee-le'',w_nz,t(timestep),t(timestep),i,w_numvars);'];
         eval(stringtoeval)
         str = w_list_vars(i,:);
         arraydata(1:w_nz,timestep) = eval([str,'(1:w_nz)']);
         eval([strtrim(str),'_array = arraydata;']);
     end
     disp(i);
end

% Reading the header from rad zm file
[rm_filename,rm_nz,rm_z,rm_ntimesteps,rm_numvars,rm_list_vars] = header_read([scm_path,radzmfile]);
 
% Read in rad file's variables into MATLAB.
% Variables will be usable in the form <GrADS Variable Name>_array
for i=1:rm_numvars
     for timestep = 1:sizet
         stringtoeval = [rm_list_vars(i,:), ' = read_grads_clubb_endian([scm_path,rm_filename],''ieee-le'',rm_nz,t(timestep),t(timestep),i,rm_numvars);'];
         eval(stringtoeval)
         str = rm_list_vars(i,:);
         arraydata(1:rm_nz,timestep) = eval([str,'(1:rm_nz)']);
         eval([strtrim(str),'_array = arraydata;']);
     end
     disp(i);
end

% Reading the header from rad file
[rt_filename,rt_nz,rt_z,rt_ntimesteps,rt_numvars,rt_list_vars] = header_read([scm_path,radztfile]);

% Read in rad file's variables into MATLAB.
% Variables will be usable in the form <GrADS Variable Name>_array
for i=1:rt_numvars
     for timestep = 1:sizet
         stringtoeval = [rt_list_vars(i,:), ' = read_grads_clubb_endian([scm_path,rt_filename],''ieee-le'',rt_nz,t(timestep),t(timestep),i,rt_numvars);'];
         eval(stringtoeval)
         str = rt_list_vars(i,:);
         arraydata(1:rt_nz,timestep) = eval([str,'(1:rt_nz)']);
         eval([strtrim(str),'_array = arraydata;']);
     end
     disp(i);
end

% Reading the header from the sfc file
[sfc_filename,sfc_nz,sfc_z,sfc_ntimesteps,sfc_numvars,sfc_list_vars] = header_read([scm_path,sfcfile]);
 
% Read in sfc file's variables into MATLAB.
% Variables will be usable in the form <GrADS Variable Name>_array
for i=1:sfc_numvars
    for timestep = 1:sizet
        stringtoeval = [sfc_list_vars(i,:), ' = read_grads_clubb_endian([scm_path,sfc_filename],''ieee-le'',sfc_nz,t(timestep),t(timestep),i,sfc_numvars);'];
        eval(stringtoeval)
        str = sfc_list_vars(i,:);
        arraydata(1:sfc_nz,timestep) = eval([str,'(1:sfc_nz)']);
        eval([strtrim(str),'_array = arraydata(1:sfc_nz,:);']);
    end
    disp(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Conversion Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform Necessary conversions
p_array = p_in_Pa_array * 0.01;

time_out = 1:sizet;
for i=1:sizet
    time_out(i) =  i;
end

full_z  = convert_units.create_time_height_series( z, sizet );		% zt grid
full_w_z = convert_units.create_time_height_series( w_z, sizet );	% zm grid
full_sfc_z = convert_units.create_time_height_series( sfc_z, sizet );
time_array = 300 .* time_out;
time_3600_array = 1:40;
for i = 1:40
	time_3600_array(i) = 3600*i;
end
time_3600_array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Definition Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the new files. By default they are in definition mode.
ncid_initial = netcdf.create( strcat( output_path, 'larson_profiles_ini.nc' ),'NC_WRITE' );	% File 1
ncid_scalars = netcdf.create( strcat( output_path, 'larson_scalars.nc' ),'NC_WRITE');		% File 2
ncid_avg = netcdf.create( strcat( output_path, 'larson_profiles_avg.nc'),'NC_WRITE');		% File 3
%ncid_forcing = netcdf.create( strcat( output_path, 'larson_profiles_forc.nc'),'NC_WRITE');	% File 4

% Define Global Attributes

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % General
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% % Contact person
netcdf.putAtt(ncid_initial,netcdf.getConstant('NC_GLOBAL'),'contact_person','Vince Larson (vlarson@uwm.edu) and Nate Meyer (meyern@uwm.edu)');
netcdf.putAtt(ncid_scalars,netcdf.getConstant('NC_GLOBAL'),'contact_person','Vince Larson (vlarson@uwm.edu) and Nate Meyer (meyern@uwm.edu)');
netcdf.putAtt(ncid_avg,netcdf.getConstant('NC_GLOBAL'),'contact_person','Vince Larson (vlarson@uwm.edu) and Nate Meyer (meyern@uwm.edu)');
%netcdf.putAtt(ncid_forcing,netcdf.getConstant('NC_GLOBAL'),'contact_person','Vince Larson (vlarson@uwm.edu) and Keith White (kcwhite@uwm.edu)');

% % Type of model where the SCM is derived from (climate model, mesoscale
% weather prediction model, regional model) ?
%netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'type_of_model_where_the_SCM_is_derived_from','Standalone SCM');

v_i_tke = zeros( 1, sizet );
for i=1:sizet
	s = 0;
	for j=1:(nzmax - 1)	% leave off the top point since we only have 41 dz's
		s = s + ( rho_ds_zm_array(j,i) .* em_array(j,i) .* ( z(1,j+1) - z(1,j) ));
	end
	v_i_tke(i) = s;
end

% Create an array of albedo values. In the astex case, albedo is defined to be .07 in the model.in file.
SSA_array = zeros( 1, sizet );
SSA_array(1,:) = .07


%%%%%%%%%%%%%%
% FILE 1
%%%%%%%%%%%%%%
% Define dimensions
% Full Levels(zt)
levfdimid = netcdf.defdim( ncid_initial, 'zf', nzmax );
levhdimid = netcdf.defdim( ncid_initial, 'zh', nzmax );

% Define variables
zfvarid = define_variable( 'zf', 'Altitude of layer mid-points (full level)', 'm', levfdimid, ncid_initial );
uvarid = define_variable( 'u', 'Zonal wind', 'm/s', levfdimid, ncid_initial );
vvarid = define_variable( 'v', 'Meridonal wind', 'm/s', levfdimid, ncid_initial );
thetalvarid = define_variable( 'thetal', 'Liquid water potential temperature', 'K', levfdimid, ncid_initial );
qtvarid = define_variable( 'qt', 'Total water', 'g/kg', levfdimid, ncid_initial );
rhovarid = define_variable( 'rho', 'Density', 'kg/m^3', levfdimid, ncid_initial );

% Set the file to add values
netcdf.setFill(ncid_initial,'NC_FILL');
netcdf.endDef(ncid_initial);

% Output to the file
netcdf.putVar( ncid_initial, zfvarid, full_z(:,1) );
netcdf.putVar( ncid_initial, uvarid, um_array(:,1) );
netcdf.putVar( ncid_initial, vvarid, vm_array(:,1) );
netcdf.putVar( ncid_initial, thetalvarid, thlm_array(:,1) );
netcdf.putVar( ncid_initial, qtvarid, 1000 * rtm_array(:,1) );
netcdf.putVar( ncid_initial, rhovarid, rho_array(:,1) );

%%%%%%%%%%%%%%
% FILE 2
%%%%%%%%%%%%%%
% Define dimensions
% Time
tdimid = netcdf.defdim(ncid_scalars,'time', sizet);

% Define variables
timevarid = define_variable( 'time', 'Time', 's', tdimid, ncid_scalars );
zcbvarid = define_variable( 'zcb', 'Cloud base height', 'm', tdimid, ncid_scalars );
ztopvarid = define_variable( 'ztop', 'Cloud top height', 'm', tdimid, ncid_scalars );
zmaxcfracvarid = define_variable( 'zmaxcfrac', 'Height of model with largest cloud fraction', 'm', tdimid, ncid_scalars );
% skipping M_base
LWPvarid = define_variable( 'LWP', 'Liquid water path', 'g/m^2', tdimid, ncid_scalars );
precwvarid = define_variable( 'precw', 'Precipitable water', 'g/m^2', tdimid, ncid_scalars );
ccvarid = define_variable( 'cc', 'Total cloud cover', '0-1', tdimid, ncid_scalars );
shfvarid = define_variable( 'shf', 'Surface sensible heat flux', 'W/m^2', tdimid, ncid_scalars );
lhfvarid = define_variable( 'lhf', 'Surface latent heat flux', 'W/m^2', tdimid, ncid_scalars );
smfvarid = define_variable( 'smf', 'Surface momentum flux', 'm^2/s^2', tdimid, ncid_scalars );
tkevarid = define_variable( 'tke', 'Vertically integrated rho*TKE', 'kg/s^2', tdimid, ncid_scalars );
v_lowlevelvarid = define_variable( 'v_lowlevel', 'Absolute velocity at the lowest model level', 'm/s', tdimid, ncid_scalars );
qv_lowlevelvarid = define_variable( 'qv_lowlevel', 'Specific humidity at the lowest model level', 'g/kg', tdimid, ncid_scalars );
th_lowlevelvarid = define_variable( 'th_lowlevel', 'Potential temperature at the lowest model level', 'm/s', tdimid, ncid_scalars );
prec_300varid = define_variable( 'prec_300', 'Precipitation flux at 300 m', 'W/m^2', tdimid, ncid_scalars );
prec_150varid = define_variable( 'prec_150', 'Precipitation flux at 150 m', 'W/m^2', tdimid, ncid_scalars );
prec_srfvarid = define_variable( 'prec_srf', 'Precipitation flux at the surface', 'W/m^2', tdimid, ncid_scalars );
Kh_150varid = define_variable( 'Kh_150', 'Eddy diffusivity for heat at 150 m', 'm^2/s^2', tdimid, ncid_scalars );
Kh_300varid = define_variable( 'Kh_300', 'Eddy diffusivity for heat at 300 m', 'm^2/s^2', tdimid, ncid_scalars );
Kh_500varid = define_variable( 'Kh_500', 'Eddy diffusivity for heat at 500 m', 'm^2/s^2', tdimid, ncid_scalars );
Kh_1250varid = define_variable( 'Kh_1250', 'Eddy diffusivity for heat at 1250 m', 'm^2/s^2', tdimid, ncid_scalars );
SSTvarid = define_variable( 'SST', 'Sea surface temperature', 'K', tdimid, ncid_scalars );
SSAvarid = define_variable( 'SSA', 'Sea surface albedo', '0-1', tdimid, ncid_scalars );
tsairvarid = define_variable( 'tsair', 'Surface air temperature', 'K', tdimid, ncid_scalars );
psvarid = define_variable( 'ps', 'Surface pressure', 'hPa', tdimid, ncid_scalars );
% skipping pblh
fsntcvarid = define_variable( 'fsntc', 'TOA SW net downward clear-sky radiation', 'W/m^2', tdimid, ncid_scalars );
fsntvarid = define_variable( 'fsnt', 'TOA SW net downward total-sky radiation', 'W/m^2', tdimid, ncid_scalars );
flntcvarid = define_variable( 'flntc', 'TOA LW clear-sky upward radiation', 'W/m^2', tdimid, ncid_scalars );
flntvarid = define_variable( 'flnt', 'TOA LW total-sky upward radiation', 'W/m^2', tdimid, ncid_scalars );
fsnscvarid = define_variable( 'fsnsc', 'Surface SW net downward clear-sky radiation', 'W/m^2', tdimid, ncid_scalars );
fsnsvarid = define_variable( 'fsns', 'Surface SW net downward total-sky radiation', 'W/m^2', tdimid, ncid_scalars );
flnscvarid = define_variable( 'flnsc', 'Surface LW net upward clear-sky radiation', 'W/m^2', tdimid, ncid_scalars );
flnsvarid = define_variable( 'flns', 'Surface SW net upward total-sky radiation', 'W/m^2', tdimid, ncid_scalars );

% Set the file to add values
netcdf.setFill(ncid_scalars,'NC_FILL');
netcdf.endDef(ncid_scalars);

% Calculations
rvm_array = rtm_array(1,:) - rcm_array(1,:);

% Output to the file
netcdf.putVar( ncid_scalars, timevarid, time_array );
netcdf.putVar( ncid_scalars, zcbvarid, z_cloud_base_array );
netcdf.putVar( ncid_scalars, ztopvarid, non_zero_cloud( cloud_frac_array, full_z, nzmax, sizet ) );
netcdf.putVar( ncid_scalars, zmaxcfracvarid, max_in_column( cloud_frac_array, full_z, nzmax, sizet ) );
% skipping M_base
netcdf.putVar( ncid_scalars, LWPvarid, lwp_array .* 1000 );
netcdf.putVar( ncid_scalars, precwvarid, vwp_array .* 1000 );
netcdf.putVar( ncid_scalars, ccvarid, cc_array );
netcdf.putVar( ncid_scalars, shfvarid, sh_array );
netcdf.putVar( ncid_scalars, lhfvarid, lh_array );
netcdf.putVar( ncid_scalars, smfvarid, sqrt( ( upwp_sfc_array .* upwp_sfc_array ) + ( vpwp_sfc_array .* vpwp_sfc_array ) ) );
netcdf.putVar( ncid_scalars, tkevarid, v_i_tke );
netcdf.putVar( ncid_scalars, v_lowlevelvarid, sqrt( ( um_array(1,:) .* um_array(1,:) ) + ( vm_array(1,:) .* vm_array(1,:) ) ) );
netcdf.putVar( ncid_scalars, qv_lowlevelvarid, 1000 * ( rvm_array ./ ( 1 + rvm_array ) ) );
netcdf.putVar( ncid_scalars, th_lowlevelvarid, thlm_array(1,:));
netcdf.putVar( ncid_scalars, prec_300varid, Fprec_array(7,:));	% 330.26 m
netcdf.putVar( ncid_scalars, prec_150varid, Fprec_array(4,:));	% 132.89 m
netcdf.putVar( ncid_scalars, prec_srfvarid, Fprec_array(1,:));  % 0 m
netcdf.putVar( ncid_scalars, Kh_150varid, Kh_zt_array(5,:));	% 160.53 m
netcdf.putVar( ncid_scalars, Kh_300varid, Kh_zt_array(7,:));	% 292.11 m
netcdf.putVar( ncid_scalars, Kh_500varid, Kh_zt_array(10,:));	% 531.58 m
netcdf.putVar( ncid_scalars, Kh_1250varid,Kh_zt_array(17,:));	% 1265.79 m
netcdf.putVar( ncid_scalars, tsairvarid, T_in_K_array(1,:));
netcdf.putVar( ncid_scalars, SSTvarid, T_sfc_array );
netcdf.putVar( ncid_scalars, SSAvarid, SSA_array );
netcdf.putVar( ncid_scalars, psvarid, p_array(1,:));
netcdf.putVar( ncid_scalars, fsntvarid, Frad_SW_down_ra_array(rm_nz,:) - Frad_SW_up_rad_array(rm_nz,:));
netcdf.putVar( ncid_scalars, flntvarid, Frad_LW_up_rad_array(rm_nz,:)  - Frad_LW_down_ra_array(rm_nz,:));
netcdf.putVar( ncid_scalars, fsnsvarid, Frad_SW_down_ra_array(1,:) - Frad_SW_up_rad_array(1,:));
netcdf.putVar( ncid_scalars, flnsvarid, Frad_LW_up_rad_array(1,:) - Frad_LW_down_ra_array(1,:));
netcdf.putVar( ncid_scalars, flnscvarid, fulwcl_array(1,:) - fdlwcl_array(1,:));
netcdf.putVar( ncid_scalars, fsnscvarid, fdswcl_array(1,:) - fuswcl_array(1,:));
netcdf.putVar( ncid_scalars, fsntcvarid, fdswcl_array(end,:) - fuswcl_array(end,:));
netcdf.putVar( ncid_scalars, flntcvarid, fulwcl_array(end,:) - fdlwcl_array(end,:));

%%%%%%%%%%%%%%
% FILE 3
%%%%%%%%%%%%%%
% Define dimensions
% Full Levels(zt)
levfdimid = netcdf.defdim(ncid_avg,'zf', nzmax );
% Half Levels (zm)
levhdimid = netcdf.defdim(ncid_avg,'zh', w_nz);
% Time
tdimid = netcdf.defdim(ncid_avg,'time', sizet / ( average_time_interval / model_timestep ));

timevarid = define_variable( 'time', 'Time', 's', tdimid, ncid_avg );
zfvarid = define_variable( 'zf', 'Altitude of layer mid-points (full level)', 'm', [levfdimid tdimid], ncid_avg );
zhvarid = define_variable( 'zh', 'Altitude of layer bottom points (half level)', 'm', [levhdimid tdimid], ncid_avg);
presvarid = define_variable( 'pres', 'Pressure', 'Pa', [levfdimid tdimid], ncid_avg )
uvarid = define_variable( 'u', 'Zonal wind', 'm/s', [levfdimid tdimid], ncid_avg );
vvarid = define_variable( 'v', 'Meridonal wind', 'm/s', [levfdimid tdimid], ncid_avg );
thetalvarid = define_variable( 'thetal', 'Liquid water potential temperature', 'K', [levfdimid tdimid], ncid_avg );
qtvarid = define_variable( 'qt', 'Total water', 'g/kg', [levfdimid tdimid], ncid_avg );
qsvarid = define_variable( 'qs', 'Saturation specific humidity', 'g/kg', [levfdimid tdimid], ncid_avg);
qlvarid = define_variable( 'ql', 'Condensed water (liquid plus ice)', 'g/kg', [levfdimid tdimid], ncid_avg);
qrvarid = define_variable( 'qr', 'Rain water', 'g/kg', [levfdimid tdimid], ncid_avg);
cfvarid = define_variable( 'cf', 'Cloud fraction', '0-1', [levfdimid tdimid], ncid_avg);
rhovarid = define_variable( 'rho', 'Density', 'kg/m^3', [levfdimid tdimid], ncid_avg );
wthlvarid = define_variable( 'wthl', 'thetal flux', 'W/m^2', [levhdimid tdimid], ncid_avg);
wqtvarid = define_variable( 'wqt', 'qt flux', 'W/m^2', [levhdimid tdimid], ncid_avg);
uwvarid = define_variable( 'uw', 'Zonal momentum flux', 'kg/(m s^2)', [levhdimid tdimid], ncid_avg);
vwvarid = define_variable( 'vw', 'Meridonal momentum flux', 'kg/(m s^2)', [levhdimid tdimid], ncid_avg);
% skipping Mf
precvarid = define_variable( 'prec', 'Precipitation flux (positive downward)', 'W/m^2', [levhdimid tdimid], ncid_avg);
% skipping w_up
% skipping thl_up
% skipping qt_up
% skipping ql_up
% skipping thv_up
TKEvarid = define_variable( 'TKE', 'Turbulent kinetic energy', 'm^2/s^2', [levhdimid tdimid], ncid_avg);
SW_upvarid = define_variable( 'SW_up', 'Upward shortwave radiation', 'W/m^2', [levhdimid tdimid], ncid_avg );
SW_dnvarid = define_variable( 'SW_dn', 'Downward shortwave radiation', 'W/m^2', [levhdimid tdimid], ncid_avg );
LW_upvarid = define_variable( 'LW_up', 'Upward longwave radiation', 'W/m^2', [levhdimid tdimid], ncid_avg );
LW_dnvarid = define_variable( 'LW_dn', 'Downward longwave radiation', 'W/m^2', [levhdimid tdimid], ncid_avg );

% Set the file to add values
netcdf.setFill(ncid_avg,'NC_FILL');
netcdf.endDef(ncid_avg);

% Output to the file
netcdf.putVar( ncid_avg, timevarid, time_3600_array );
netcdf.putVar( ncid_avg, zfvarid, time_average( full_z, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, zhvarid, time_average( full_w_z, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, presvarid, time_average( p_in_Pa_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, uvarid, time_average( um_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, vvarid, time_average( vm_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, thetalvarid, time_average( thlm_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, qtvarid, 1000 * time_average( rtm_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, qsvarid, 1000 * time_average( rsat_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, qlvarid, 1000 * time_average( rcm_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, qrvarid, 1000 * time_average( rrm_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, cfvarid, time_average( cloud_frac_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, rhovarid, time_average( rho_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, wthlvarid, time_average( ( rho_zm_array .* convert_units.Cp .* wpthlp_array), model_timestep, average_time_interval, sizet, nzmax ) ); % check
netcdf.putVar( ncid_avg, wqtvarid, time_average( ( ( rho_zm_array .* convert_units.Lv .* wprtp_array ) ./ ( 1+rtm_array ) ), model_timestep, average_time_interval, sizet, nzmax ) ); % check
netcdf.putVar( ncid_avg, uwvarid, time_average( ( rho_zm_array .* upwp_array ), model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, vwvarid, time_average( ( rho_zm_array .* vpwp_array ), model_timestep, average_time_interval, sizet, nzmax ) );
% skipping Mf
netcdf.putVar( ncid_avg, precvarid, time_average( Fprec_array, model_timestep, average_time_interval, sizet, nzmax ) );
% skipping w_up
% skipping thl_up
% skipping qt_up
% skipping ql_up
% skipping thv_up
netcdf.putVar( ncid_avg, TKEvarid, time_average( em_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, SW_upvarid, time_average( Frad_SW_up_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, SW_dnvarid, time_average( Frad_SW_down_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, LW_upvarid, time_average( Frad_LW_up_array, model_timestep, average_time_interval, sizet, nzmax ) );
netcdf.putVar( ncid_avg, LW_dnvarid, time_average( Frad_LW_down_array, model_timestep, average_time_interval, sizet, nzmax ) );

%%%%%%%%%%%%%%
% FILE 4
%%%%%%%%%%%%%%
% zfvarid is defined in file 1 section
% dTdt_ls
% dqdt_ls
% dudt_ls
% dvdt_ls

% Output to the file

% Close files
netcdf.close(ncid_initial);
netcdf.close(ncid_scalars);
netcdf.close(ncid_avg);
%netcdf.close(ncid_forcing);
end

function varid = define_variable( shrt_name, long_name, units, dim_ids, file_id )

varid = netcdf.defVar(file_id, shrt_name, 'NC_FLOAT',dim_ids);
netcdf.putAtt(file_id, varid,'unit',units);
netcdf.putAtt(file_id, varid,'long_name',long_name);

end

% Returns an array of time-averaged values, depending on the desired interval and the model timestep.  The number of values
% averaged for a data point will be the desired interval (average_time_interval) divided by the model timestep (model_timestep).
% The dimensions of the returned array will be the height of the profile (nzmax), the rows, by the number of timesteps (num_timesteps)
% divided by the number of columns being averaged (step), the columns.
function avg_array = time_average( array_to_average, model_timestep, average_time_interval, num_timesteps, nzmax )
	step = average_time_interval / model_timestep;	% number of columns per interval for averaging
	num_columns = num_timesteps / step;	% number of columns in return array
	avg = zeros(nzmax, 1);			% temp array to hold totals for averaging
	avg_array = zeros(nzmax, num_columns);	% array to return
	index = 1;
	for i=1:step:num_timesteps
		for j=i:i+step-1
			avg = avg + array_to_average(:,j);	% add all the columns over the step interval
		end
		avg_array(:,index) = avg ./ step;
		index = index + 1;
		avg = zeros(nzmax,1);		% reset the temp array
	end
		

end

% This function finds the max value in each column and translates it to the associated height.
% Note that this will find the highest, or last, maximum it encounters.
function max_array = max_in_column( input_array, height_array, nzmax, sizet )
	max_array = zeros( 1, sizet );
	for i=1:sizet		% columns
		compare = 0;
		row = 0;
		column =0;
		for j=1:nzmax	% rows
			if input_array(j,i) >= compare
				compare = input_array(j,i);
				row = j;
				column = i;
			end
		end
		max_array(i) = height_array(row, column);
	end
end

% This function finds the highest non-zero value in each coulumn and translates it to the associated height.
function non_zero_cf = non_zero_cloud( input_array, height_array, nzmax, sizet )
	non_zero_cf = zeros( 1, sizet );
	for i=1:sizet		% columns
		row = 0;
		column =0;
		for j=1:nzmax	% rows
			if input_array(j,i) > 0
				row = j;
				column = i;
			end
		end
		if row >0
			height_array(row, column)
			non_zero_cf(i) = height_array(row, column);
		else
			non_zero_cf(i) = 0;
		end
	end
end