function[] = netcdf_gabls3_writer();

% Necessary paths to use Mexnc
addpath '/home/mjfalk/mexnc2/mexnc/'
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/mjfalk/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

% Source files of the GABLS3 case
scm_path = ['/home/faschinj/projects/experimental/hoc_v2.2_tuner/standalone/'];
smfile   = 'gabls3_20081201_zt.ctl';
swfile   = 'gabls3_20081201_zm.ctl';
sfcfile  = 'gabls3_20081201_sfc.ctl';

% Range of output data by time
t = 1:144;
sizet = size(t);
sizet = max(sizet);

% Reading the header from zt file
[filename,nz,z,ntimesteps,numvars,list_vars] = header_read([scm_path,smfile]);

% instantaneous, for initial time
for i=1:numvars
    for timestep = 1:sizet
    	stringtoeval = [list_vars(i,:), ' = read_grads_hoc_endian([scm_path,filename],''ieee-le'',nz,t(timestep),t(timestep),i,numvars);'];
    	eval(stringtoeval);
    	str = list_vars(i,:);
    	for j=1:nz
        	arraydata(j,timestep) = eval([str,'(j)']);
    	end
    	eval([strtrim(str),'_array = arraydata;']);
    end
    disp(i);
end

[w_filename,w_nz,w_z,w_ntimesteps,w_numvars,w_list_vars] = header_read([scm_path,swfile]);
 
% every 10 minutes
for i=1:w_numvars
     for timestep = 1:sizet
         stringtoeval = [w_list_vars(i,:), ' = read_grads_hoc_endian([scm_path,w_filename],''ieee-le'',w_nz,t(timestep),t(timestep),i,w_numvars);'];
         eval(stringtoeval)
         str = w_list_vars(i,:);
         for j=1:w_nz
             arraydata(j,timestep) = eval([str,'(j)']);
         end
         eval([strtrim(str),'_array = arraydata;']);
     end
       disp(i);
end

[sfc_filename,sfc_nz,sfc_z,sfc_ntimesteps,sfc_numvars,sfc_list_vars] = header_read([scm_path,sfcfile]);
 
% every 10 minutes
for i=1:sfc_numvars
    for timestep = 1:sizet
        stringtoeval = [sfc_list_vars(i,:), ' = read_grads_hoc_endian([scm_path,sfc_filename],''ieee-le'',sfc_nz,t(timestep),t(timestep),i,sfc_numvars);'];
        eval(stringtoeval)
        str = sfc_list_vars(i,:);
        for j=1:sfc_nz
             arraydata(j,timestep) = eval([str,'(j)']);
        end
        eval([strtrim(str),'_array = arraydata;']);
    end
       disp(i);
end


% Perform Necessary conversions
qtm_array = rtm_array ./ (1 + rtm_array);

g0 = 9.8;
p0 = 1e5;
R  = 287.04;
Cp = 1004.67;
Lv = 2.5e6;

%theta_array = thlm_array + (Lv/Cp).*rcm_array;
%T_array = theta_array .* exner_array;
%T_in_K_p_array = thlp_array .* exner_array; 
T_forcing_array = (thlm_f_array - radht_array) .* exner_array;
%ome_array = -(wm_array .* p_in_Pa_array)./ (R * g0 .* thvm_array);
ome_array = -(wm_array .* g0 .*rho_array);
wt_array = wpthlp_array .* exner_array;
wq_array = wprtp_array ./ (1 + rtm_array);

[ncid,status] = mexnc('create', 'gabls3_scm_UWM_CLUBB_v2.nc', nc_clobber_mode)
status = mexnc('close', ncid)

[ncid,status] = mexnc('open', 'gabls3_scm_UWM_CLUBB_v2.nc', nc_write_mode)
% Redefine the file
status = mexnc('redef',ncid)
% Define Global Attributes

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % General
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Reference to the model
% input = 'Golaz et. al 2002';
% [refmodid,status] = mexnc('put_att_text',ncid,NC_GLOBAL,'Reference to the model','string',length(input),input)
% 
% % contact person
% input = '';
% [contactPid,status] = mexnc('put_att_text',ncid,'NC_GLOBAL','contact person',length(input),input)
% 
% %type of model where the SCM is derived from (climate model, mesoscale weather prediction model, regional model) ?
% input = 'Standalone SCM';
% [modelTypeid,status] = mexnc('put_att_text',ncid,'NC_GLOBAL','type of model were the SCM is derived from',strlen(input),input)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Surface Scheme
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Is it a force-restore type or a multi-layer type?
% input = 'force-restore';
% [surfaceType,status] = mexnc('put_att_text',ncid,'NC_GLOBAL','Is it a force-restore type or a multi-layer type?',strlen(input),input)
% 
% %Does it have skin layer?
% input = 'No';
% [skinboolid,status] = mexnc('put_att_text',ncid,'NC_GLOBAL','Does it have skin layer',strlen(input),input)
% 
% %Is there a tile approach?
% input = 'No';
% [tileboolid,status] = mexnc('put_att_text',ncid,'NC_GLOBAL','Is there a time approach?',strlen(input),input)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Turbulence Scheme
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Turbulence scheme  (e.g., K profile, TKE-l, ...)
% input = 'Higher order closure';
% [turbulenceid,status] = mexnc('put_att_text',ncid,'NC_GLOBAL','Turbulence scheme  (e.g., K profile, TKE-l, ...)',strlen(input), input)
% 
% %Formulation of eddy diffusivity K.
% input = 'No eddy diffusivitiy; fluxes are prognosed';
% [diffuseid,status] = mexnc('put_att_text',ncid,'NC_GLOBAL','Formulation of eddy diffusivity K.',strlen(input),input)
% 
% % For E-l and Louis-type scheme: give formulation length scale.
% % For K-profile: how is this  profile determined ? (e.g., based on Richardson, Brunt-Vaisala frequency (N^2),  Parcel method, other.
% input = 'Parcel method';
% [lengthscaleid,status] = mexnc('put_att_text',ncid,'NC_GLOBAL','How is this profile determined?',strlen(input),input)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Other
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Any other model specific aspects you find relevent for this intercomparison.
% 
% % Any deviation from the prescribed set up that you had to make because of the specific structure of your model
% input = 'We needed to set the temperature of the top soil layer and vegetation to match the surface air at the initial time.';
% [deviationid,status] = mexnc('put_att_text',ncid,'NC_GLOBAL','time step',strlen(input),input)



% Define dimensions
% Output Time
[tdimid,status] = mexnc('def_dim',ncid,'time',sizet);
% Full Levels (zt)
[levhdimid,status] = mexnc('def_dim',ncid,'levh',w_nz);
% Half Levels(zm)
[levfdimid,status] = mexnc('def_dim',ncid,'levf',nz);
% Soil Levels(sfc?)
[levsdimid,status] = mexnc('def_dim',ncid,'levs',sfc_nz);

% Define variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time series output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tvarid,status] = define_variable_1d( ncid,'time' ,'seconds since 2006-07-01 12:00:00', 's', tdimid, ncid );
[ldwvarid,status] = define_variable_1d( ncid,'ldw' ,'long wave downward radiation at surface', 'W/m^2', tdimid, ncid );
[lupvarid,status] = define_variable_1d( ncid,'lup' ,'long wave upward radiation at surface', 'W/m^2', tdimid, ncid );
[qdwvarid,status] = define_variable_1d( ncid,'qdw' ,'short wave downward radiation at surface', 'W/m^2', tdimid, ncid );
[qupvarid,status] = define_variable_1d( ncid,'qup' ,'short wave upward radiation at surface', 'W/m^2', tdimid, ncid );
[tskvarid,status] = define_variable_1d( ncid,'tsk' ,'temperature skin layer', 'W/m^2', tdimid, ncid );
[qskvarid,status] = define_variable_1d( ncid,'qsk' ,'specific humiditiy top of vegetation layer', 'kg/kg', tdimid, ncid );
[gvarid,status] = define_variable_1d( ncid,'g' ,'soil heat flux', 'W/m^2', tdimid, ncid );
[shfvarid,status] = define_variable_1d( ncid,'shf' ,'sensible heat flux', 'W/m^2', tdimid, ncid );
[lhfvarid, status] = define_variable_1d( ncid,'lhf' ,'latent heat flux', 'W/m^2', tdimid, ncid );
[ustarvarid,status] = define_variable_1d( ncid,'ustar' ,'surface velocity', 'm/s', tdimid, ncid );
[hpblvarid, status] = define_variable_1d( ncid,'hpbl' ,'boundary layer height', 'm', tdimid, ncid );
[t2mvarid, status] = define_variable_1d( ncid,'t2m' ,'2m temperature', 'K', tdimid, ncid );
[q2mvarid, status] =  define_variable_1d( ncid,'q2m' ,'2m specific humidity', 'kg/kg', tdimid, ncid );
[u10mvarid, status] = define_variable_1d( ncid,'u10m' ,'10m u-component wind', 'm/s', tdimid, ncid );
[v10mvarid, status] =  define_variable_1d( ncid,'v10m' ,'10m v-component wind', 'm/s', tdimid, ncid );
[ccvarid, status] = define_variable_1d( ncid,'cc' ,'cloudcover fration', '0 1', tdimid, ncid );
[zfvarid, status] = define_variable_2d( ncid,'zf' ,'height of full level', 'm', tdimid, levfdimid, ncid );
[pfvarid, status] = define_variable_2d( ncid,'pf' ,'pressure at full level', 'Pa', tdimid, levfdimid, ncid );
[tkvarid, status] = define_variable_2d( ncid,'t' ,'temperature', 'K', tdimid, levfdimid, ncid );
[thvarid, status] = define_variable_2d( ncid,'th' ,'potential temperature', 'K', tdimid, levfdimid, ncid );
[qvarid, status] = define_variable_2d( ncid,'q' ,'specific humidity', 'kg/kg', tdimid, levfdimid, ncid );
[uvarid, status] = define_variable_2d( ncid,'u' ,'zonal component wind', 'm/s', tdimid, levfdimid, ncid );
[vvarid, status] = define_variable_2d( ncid,'v' ,'meridonal component wind', 'm/s', tdimid, levfdimid, ncid );
[ugeovarid, status] = define_variable_2d( ncid,'ugeo' ,'u-component geostrophic wind', 'm/s', tdimid, levfdimid, ncid );
[vgeovarid, status] = define_variable_2d( ncid,'vgeo' ,'v-component geostrophic wind', 'm/s', tdimid, levfdimid, ncid );
[dudt_lsvarid, status] = define_variable_2d( ncid,'dudt_ls' ,'u-component momentum advection', 'm/s/s', tdimid, levhdimid, ncid );
[dvdt_lsvarid, status] = define_variable_2d( ncid,'dvdt_ls' ,'v-component momentum advection', 'm/s/s', tdimid, levhdimid, ncid );
[dtdt_lsvarid, status] = define_variable_2d( ncid,'dtdt_ls' ,'temperature advection', 'K/s', tdimid, levhdimid, ncid );
[dqdt_lsvarid, status] = define_variable_2d( ncid,'dqdt_ls' ,'specific humidity advection', 'kg/kg/s', tdimid, levhdimid, ncid );
[omevarid, status] = define_variable_2d( ncid,'ome' ,'vertical movement', 'Pa/s', tdimid, levfdimid, ncid );
[zhvarid, status] = define_variable_2d( ncid,'zh' ,'height of half level', 'm', tdimid, levhdimid, ncid );
[phvarid, status] = define_variable_2d( ncid,'ph' ,'pressure at half level', 'Pa', tdimid, levhdimid, ncid );
[wtvarid, status] = define_variable_2d( ncid,'wt' ,'vertical temperature flux', 'Km/s', tdimid, levhdimid, ncid );
[wqvarid, status] = define_variable_2d( ncid,'wq' ,'vertical moisture flux', 'kg/kg m/s', tdimid, levhdimid, ncid );
[uwvarid, status] = define_variable_2d( ncid,'uw' ,'vertical flux u-component', 'm2/s2', tdimid, levhdimid, ncid );
[vwvarid, status] = define_variable_2d( ncid,'vw' ,'vertical flux v-component', 'm2/s2', tdimid, levhdimid, ncid );
[Kmvarid, status] = define_variable_2d( ncid,'Km' ,'eddy diffusivity momentum', 'm2/s', tdimid, levhdimid, ncid );
[Khvarid, status] = define_variable_2d( ncid,'Kh' ,'eddy diffusivity heat', 'm2/s', tdimid, levhdimid, ncid );
[TKEvarid, status] = define_variable_2d( ncid,'TKE' ,'turbulent kinetic energy', 'm^2/s^2', tdimid, levhdimid, ncid );
[shearvarid, status] = define_variable_2d( ncid,'shear' ,'shear production', 'm2/s3', tdimid, levhdimid, ncid );
[buoyvarid, status] = define_variable_2d( ncid,'buoy' ,'buoyancy production', 'm2/s3', tdimid, levhdimid, ncid );
[transvarid, status] = define_variable_2d( ncid,'trans' ,'total transport', 'm2/s3', tdimid, levhdimid, ncid );
[dissivarid, status] = define_variable_2d( ncid,'dissi' ,'dissipation', 'm2/s3', tdimid, levhdimid, ncid );
[zsvarid, status] = define_variable_2d( ncid,'zs' ,'height of soil level', 'm', tdimid, levsdimid, ncid );
[tsvarid, status] = define_variable_2d( ncid,'ts' ,'soil temperature', 'K', tdimid, levsdimid, ncid );
[thsvarid, status] = define_variable_2d( ncid,'ths' ,'soil water content', 'm2/s3', tdimid, levsdimid, ncid );%zeros



% Indicate which variables are not included
status = mexnc('put_att_float',ncid,tskvarid,'_FillValue','float',1,NaN);
status = mexnc('put_att_float',ncid ,gvarid,'_FillValue','float',1,NaN);
status = mexnc('put_att_float',ncid ,Kmvarid,'_FillValue','float',1,NaN);
%[Khvarid, status] = define_variable_2d( 'Kh' ,'eddy diffusivity heat', 'm2/s', tdimid, levf, ncid );
status = mexnc('put_att_float',ncid ,Khvarid,'_FillValue','float',1,NaN);

	%[transvarid, status] = define_variable_2d( 'trans' ,'total transport', 'm2/s3', tdimid, levf, ncid );
	% Zeros
status = mexnc('put_att_float',ncid ,transvarid,'_FillValue','float',1,NaN);

	%[dissivarid, status] = define_variable_2d( 'dissi' ,'dissipation', 'm2/s3', tdimid, levf, ncid );
	% Zeros
status = mexnc('put_att_float',ncid ,dissivarid,'_FillValue','float',1,NaN);

	%[thsvarid, status] = define_variable_2d( 'ths' ,'soil water content', 'm2/s3', tdimid, levf, ncid );
status = mexnc('put_att_float',ncid ,thsvarid,'_FillValue','float',1,NaN);

status = mexnc('put_att_float',ncid ,hpblvarid,'_FillValue','float',1,NaN);

status = mexnc('end_def',ncid);

% Apply values to file

for i=1:sizet
	% Example %
	% status = mexnc('VARPUT',ncid,thetalvarid,[i-1 k-1],[1 1],thlm_array(k,i),0);

	status = mexnc('VARPUT', ncid, tvarid, i-1,1,i*10.0*60.0,0);
	status = mexnc('VARPUT', ncid, ldwvarid, i-1,1,Frad_LW_down_array(1,i),0);
	status = mexnc('VARPUT', ncid, lupvarid, i-1,1,Frad_LW_up_array(1,i),0);
	status = mexnc('VARPUT', ncid, qdwvarid, i-1,1,Frad_SW_down_array(1,i),0);
	status = mexnc('VARPUT', ncid, qupvarid, i-1,1,Frad_SW_up_array(1,i),0);
	status = mexnc('VARPUT', ncid, shfvarid, i-1,1,sh_array(1,i),0);
	status = mexnc('VARPUT', ncid, lhfvarid, i-1,1,lh_array(1,i),0);
	status = mexnc('VARPUT', ncid, ustarvarid, i-1,1,ustar_array(1,i),0);
	%status = mexnc('VARPUT', ncid, hpblvarid, i-1,1,w_z(nz),0);
	status = mexnc('VARPUT', ncid, t2mvarid, i-1,1,T_in_K_array(1,i),0);
	status = mexnc('VARPUT', ncid, q2mvarid, i-1,1,qtm_array(1,i),0);
	status = mexnc('VARPUT', ncid, u10mvarid, i-1,1,um_array(2,i),0);
	status = mexnc('VARPUT', ncid, v10mvarid, i-1,1, vm_array(2,i),0);
	status = mexnc('VARPUT', ncid, ccvarid, i-1,1, cc_array(1,i),0);
end


for k=1:nz
    % k-1 comes from NetCDF starting variables at 0 and MATLAB starting
    % them at 1.
    
    % start, count, value, autoscale
%    status = mexnc('VARPUT',ncid,rhovarid,k-1,1,rho_array(k,1),0)
	for i=1:sizet
	
	%Mean State
	%Height at Full
	status = mexnc('VARPUT',ncid,zfvarid,[i-1 k-1],[1 1],z(k),0);
	%pressure at full level
	status = mexnc('VARPUT',ncid,pfvarid,[i-1 k-1],[1 1],p_in_Pa_array(k,i),0); % Grid switch
	%temperature
	status = mexnc('VARPUT',ncid,tkvarid,[i-1 k-1],[1 1],T_in_K_array(k,i),0); % Grid switch
	%potential temperature
	status = mexnc('VARPUT',ncid,thvarid,[i-1 k-1],[1 1],thlm_array(k,i),0); % Grid switch
	%specific humidity
	status = mexnc('VARPUT',ncid,qvarid,[i-1 k-1],[1 1],qtm_array(k,i),0); % Grid switch
	%zonal component wind;
	status = mexnc('VARPUT',ncid,uvarid,[i-1 k-1],[1 1],um_array(k,i),0); % Grid switch
	%meridonal component wind
	status = mexnc('VARPUT',ncid,vvarid,[i-1 k-1],[1 1],vm_array(k,i),0); % Grid switch

	%Fluxes
	%u-component geostrophic wind
	status = mexnc('VARPUT',ncid,ugeovarid,[i-1 k-1],[1 1],ug_array(k,i),0);
	%v-component geostrophic wind
	status = mexnc('VARPUT',ncid,vgeovarid,[i-1 k-1],[1 1],vg_array(k,i),0);	
	%u-component momentum advection;
	status = mexnc('VARPUT',ncid,dudt_lsvarid,[i-1 k-1],[1 1],um_f_array(k,i),0); % Grid switch
	%v-component momentum advection
	status = mexnc('VARPUT',ncid,dvdt_lsvarid,[i-1 k-1],[1 1],vm_f_array(k,i),0); % Grid switch
	%temperature advection
	status = mexnc('VARPUT',ncid,dtdt_lsvarid,[i-1 k-1],[1 1],T_forcing_array(k,i),0); % Grid switch
	%specific humidity advection
	status = mexnc('VARPUT',ncid,dqdt_lsvarid,[i-1 k-1],[1 1],rtm_f_array(k,i),0); % Grid switch
	%vertical movement
	status = mexnc('VARPUT',ncid,omevarid,[i-1 k-1],[1 1],ome_array(k,i),0); % Grid switch
	end
	
end

for k=1:w_nz
    % k-1 comes from NetCDF starting variables at 0 and MATLAB starting
    % them at 1.
    
    % start, count, value, autoscale
%    status = mexnc('VARPUT',ncid,rhovarid,k-1,1,rho_array(k,1),0)
	for i=1:sizet

	% Height at half level
	status = mexnc('VARPUT',ncid,zhvarid,[i-1 k-1],[1 1],w_z(k),0);
	
	% Pressure at half level ???
	
	%vertical temperature flux
	status = mexnc('VARPUT',ncid,wtvarid,[i-1 k-1],[1 1],wt_array(k,i),0);

	%vertical moisture flux ????
    status = mexnc('VARPUT',ncid,wqvarid,[i-1 k-1],[1 1],wq_array(k,i),0);
    %[wqvarid, status] = define_variable_2d( ncid,'wq' ,'vertical moisture flux', 'kg/kg m/s', tdimid, levhdimid, ncid );
	%vertical flux u-component
	status = mexnc('VARPUT',ncid,uwvarid,[i-1 k-1],[1 1],upwp_array(k,i),0);
	%vertical flux v-component
	status = mexnc('VARPUT',ncid,vwvarid,[i-1 k-1],[1 1],vpwp_array(k,i),0);
	%eddy diffusivity momentum
	% Zeros
	
	% Eddy diffusivity heat
	% Zeros
	
        % TKE
	status = mexnc('VARPUT',ncid, TKEvarid,[i-1 k-1],[1 1],em_array(k,i),0);

	% shear production
	status = mexnc('VARPUT',ncid, shearvarid,[i-1 k-1],[1 1],shear_array(k,i),0);
	% buoyancy production;
	status = mexnc('VARPUT',ncid, buoyvarid,[i-1 k-1],[1 1],wp2_bp_array(k,i),0);

	%total transport
	% Zeros

	% dissipation
	% Zeros

	end
	
end

for k=1:sfc_nz
	for i=1:sizet
		%[zsvarid, status] = define_variable_2d( 'zs' ,'height of soil level', 'm', tdimid, levf, ncid );
		status = mexnc('VARPUT',ncid,zsvarid,[i-1 k-1],[1 1],sfc_z(k),0); % Grid switch
		%[tsvarid, status] = define_variable_2d( 'ts' ,'soil temperature', 'K', tdimid, levf, ncid );
		status = mexnc('VARPUT',ncid,tsvarid,[i-1 k-1],[1 1],t_sfc_array(k,i),0); % Grid switch
		%[thsvarid, status] = define_variable_2d( 'ths' ,'soil water content', 'm2/s3', tdimid, levf, ncid );
	end
end

%[ndims,nvars,ngatts,unlimdim,status] = mexnc('inq',ncid)
status = mexnc('close',ncid);

function [varid,status] = define_variable_1d( ncid, shrt_name, long_name, units, dim_id, file_id );
	[varid,status] = mexnc('def_var',file_id,shrt_name,'float',1,[dim_id]);
	d_status = mexnc('put_att_text',ncid, varid,'unit','char',length(units), units);
	d_status = mexnc('put_att_text',ncid, varid,'long_name','char',length(long_name),long_name);


function [varid,status] = define_variable_2d( ncid, shrt_name, long_name, units, dim1_id, dim2_id, file_id );
	[varid,status] = mexnc('def_var',file_id,shrt_name,'float', 2, [dim1_id dim2_id]);
	d_status = mexnc('put_att_text',ncid, varid,'unit','char',length(units), units);
	d_status = mexnc('put_att_text',ncid, varid,'long_name','char',length(long_name),long_name);


