function[] = twp_ice_timeseries_creator(ensembleNumber)
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

%Define any needed constants
sec_per_day = 86400;
mm_per_m = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File Input Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Source files of the GABLS3 case

case_name = 'twp_ice';

clubb_path = ['/home/nielsenb/clubb/'];

% Path of the GrADS input files
scm_path = [clubb_path, 'run_scripts/ensemble_run/output/ensemble', sprintf('%02d',ensembleNumber), '/'];
output_path = [clubb_path, 'postprocessing/output_scripts/twp_ice/timeseries.CLUBB_p', sprintf('%02d',ensembleNumber), '.nc'];

% zt Grid
smfile   = [case_name, '_zt.ctl'];

% zm Grid
swfile   = [case_name, '_zm.ctl'];

% sfc Grid
sfcfile  = [case_name, '_sfc.ctl'];

% Reading the header from zt file
[filename,nzmax,z,ntimesteps,numvars,list_vars] = header_read([scm_path,smfile]);

% Used to navigate through time in all of the files. At the time this
% script was written, all of the GrADS output files generated by CLUBB used
% the same timestep.
t_start = 1; %This defines the starting timestep
sizet = ntimesteps; %This defines the end timestep
t = 1:ntimesteps;

% Read in zt file's variables into MATLAB.
% Variables will be usable in the form <GrADS Variable Name>_array.
for i=1:numvars
    t_index = 1;
    for timestep = t_start:sizet
    	stringtoeval = [list_vars(i,:), ' = read_grads_clubb_endian([scm_path,filename],''ieee-le'',nzmax,t(timestep),t(timestep),i,numvars);'];
    	eval(stringtoeval);
    	str = list_vars(i,:);
        arraydata(1:nzmax,t_index) = eval([str,'(1:nzmax)']);
    	eval([strtrim(str),'_array = arraydata;']);
	t_index = t_index + 1;
    end
    disp(i);
end

% Reading the header from zm file
[w_filename,w_nz,w_z,w_ntimesteps,w_numvars,w_list_vars] = header_read([scm_path,swfile]);
 
% Read in zm file's variables into MATLAB.
% Variables will be usable in the form <GrADS Variable Name>_array
for i=1:w_numvars
     t_index = 1;
     for timestep = t_start:sizet
         stringtoeval = [w_list_vars(i,:), ' = read_grads_clubb_endian([scm_path,w_filename],''ieee-le'',w_nz,t(timestep),t(timestep),i,w_numvars);'];
         eval(stringtoeval)
         str = w_list_vars(i,:);
         arraydata(1:w_nz,t_index) = eval([str,'(1:w_nz)']);
         eval([strtrim(str),'_array = arraydata;']);
	 t_index = t_index + 1;
     end
     disp(i);
end

% Reading the header from the sfc file
[sfc_filename,sfc_nz,sfc_z,sfc_ntimesteps,sfc_numvars,sfc_list_vars] = header_read([scm_path,sfcfile]);
 
% Read in sfc file's variables into MATLAB.
% Variables will be usable in the form <GrADS Variable Name>_array
for i=1:sfc_numvars
    t_index = 1;
    for timestep = t_start:sizet
        stringtoeval = [sfc_list_vars(i,:), ' = read_grads_clubb_endian([scm_path,sfc_filename],''ieee-le'',sfc_nz,t(timestep),t(timestep),i,sfc_numvars);'];
        eval(stringtoeval)
        str = sfc_list_vars(i,:);
        arraydata(1:sfc_nz,t_index) = eval([str,'(1:sfc_nz)']);
        eval([strtrim(str),'_array = arraydata(1:sfc_nz,:);']);
	t_index = t_index + 1;
    end
    disp(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Conversion Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform Necessary conversions
qtm_array = convert_units.total_water_mixing_ratio_to_specific_humidity( rtm_array );
T_forcing_array = convert_units.thlm_f_to_T_f( thlm_f_array, radht_array, exner_array );
ome_array = convert_units.w_wind_in_ms_to_Pas( wm_array, rho_array );
wt_array = convert_units.potential_temperature_to_temperature( wpthlp_array, exner_array );

ppt_array = rho_zm_array(1,:) .* morr_rain_rate_array ./ (sec_per_day * mm_per_m);

wq_array = wprtp_array ./ (1 + rtm_array);

icewaterpath_array = iwp_array + swp_array;

time_out = 1:(sizet - t_start + 1);
for i=1:(sizet - t_start + 1)
    time_out(i) =  i;
end

full_z  = convert_units.create_time_height_series( z, sizet - t_start + 1 );
full_w_z = convert_units.create_time_height_series( w_z, sizet - t_start + 1 );
full_sfc_z = convert_units.create_time_height_series( sfc_z, sizet - t_start + 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Definition Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the new file. By default it is in definition mode.
ncid = netcdf.create(output_path,'NC_WRITE');

% Define Global Attributes

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % General
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% % Reference to the model
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Reference_to_the_model','Golaz et. al 2002');

% % contact person
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'contact_person','Brandon Nielsen (nielsenb@uwm.edu)');
 
% % Type of model where the SCM is derived from (climate model, mesoscale
% weather prediction model, regional model) ?
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'type_of_model_where_the_SCM_is_derived_from','Standalone SCM');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Surface Scheme
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% % Is it a force-restore type or a multi-layer type?
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ... 
    'Is_it_a_force-restore_type_or_a_multi-layer_type','force-restore');

% %Does it have skin layer?
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'Does_it_have_skin_layer','No');

% %Is there a tile approach?
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'Is_there_a_tile_approach','No');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Turbulence Scheme
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Turbulence scheme  (e.g., K profile, TKE-l, ...)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ... 
    'Turbulence_scheme',...
    'Higher order closure');

% %Formulation of eddy diffusivity K.
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'Formulation_of_eddy_diffusivity_K.', ...
    'No eddy diffusivitiy; fluxes are prognosed');

% % For E-l and Louis-type scheme: give formulation length scale.
% % For K-profile: how is this  profile determined ? (e.g., based on
% Richardson, Brunt-Vaisala frequency (N^2),  Parcel method, other.
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), ...
    'How_is_this_profile_determined','Parcel method');

% Define dimensions

% Output Time
tdimid = netcdf.defdim(ncid,'time', (sizet - t_start + 1));

% Half Levels (zm)
levhdimid = netcdf.defdim(ncid,'levh', w_nz);

% Full Levels(zt)
levfdimid = netcdf.defdim(ncid,'levf', nzmax );

% Soil Levels(sfc)
levsdimid = netcdf.defdim(ncid,'levs', sfc_nz);

% Define variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timeseries Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timevarid = define_variable( 'time' ,'hours since 0000Z 20 Jan 2006', 'h', tdimid, ncid );
fsvarid = define_variable( 'Fs' ,'surface turbulent flux of sensible heat', 'W m^-2', [levsdimid tdimid], ncid );
fqvarid = define_variable( 'Fq' ,'surface turbulent flux of latent heat', 'W m^-2', [levsdimid tdimid], ncid );
dTFsw0varid = define_variable( 'dTFsw0' ,'total surface downwelling shortwave radiative flux', 'W m^-2', [levsdimid tdimid], ncid );
dTFlw0varid = define_variable( 'dTFlw0' ,'total surface downwelling longwave radiative flux', 'W m^-2', [levsdimid tdimid], ncid );
uTFsw0varid = define_variable( 'uTFsw0' ,'total surface upwelling shortwave radiative flux', 'W m^-2', [levsdimid tdimid], ncid );
uTFlw0varid = define_variable( 'uTFlw0' ,'total surface upwelling longwave radiative flux', 'W m^-2', [levsdimid tdimid], ncid );
pptvarid = define_variable( 'ppt' ,'surface precipitation', 'kg m^-2 s^-1', [levsdimid tdimid], ncid );
pwvarid = define_variable( 'PW' ,'preciptable water', 'kg m^-2', [levsdimid tdimid], ncid );
lwpvarid = define_variable( 'LWP' ,'liquid water path', 'kg m^-2', [levsdimid tdimid], ncid );
iwpvarid = define_variable( 'IWP' ,'ice water path', 'kg m^-2', [levsdimid tdimid], ncid );

netcdf.setFill(ncid,'NC_FILL');

% End definition
netcdf.endDef(ncid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Output File Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean State Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

netcdf.putVar( ncid, timevarid, time_out);
netcdf.putVar(ncid,fsvarid,sh_array);
netcdf.putVar(ncid,fqvarid,lh_array);
netcdf.putVar(ncid,dTFsw0varid,Frad_SW_down_array(1,:));
netcdf.putVar(ncid,dTFlw0varid,Frad_LW_down_array(1,:));
netcdf.putVar(ncid,uTFsw0varid,Frad_SW_up_array(1,:));
netcdf.putVar(ncid,uTFlw0varid,Frad_LW_up_array(1,:));
netcdf.putVar(ncid,pptvarid,ppt_array);
netcdf.putVar(ncid,pwvarid,vwp_array);
netcdf.putVar(ncid,lwpvarid,lwp_array);
netcdf.putVar(ncid,iwpvarid,icewaterpath_array);

% Close file
netcdf.close(ncid);
end

function varid = define_variable( shrt_name, long_name, units, dim_ids, file_id )

varid = netcdf.defVar(file_id, shrt_name, 'NC_FLOAT',dim_ids);
netcdf.putAtt(file_id, varid,'unit',units);
netcdf.putAtt(file_id, varid,'long_name',long_name);

end









