function[] = gabls3_output_creator_matlab()
% Source files of the GABLS3 case
scm_path = ['/home/faschinj/projects/experimental/clubb/output/'];
smfile   = 'gabls3_zt.ctl';
swfile   = 'gabls3_zm.ctl';
sfcfile  = 'gabls3_sfc.ctl';

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
ncid = netcdf.create('gabls3_scm_UWM_CLUBB_v3','NC_WRITE');
%[ncid,status] = mexnc('create', 'gabls3_scm_UWM_CLUBB_v2.nc', nc_clobber_mode)
%status = mexnc('close', ncid)

%[ncid,status] = mexnc('open', 'gabls3_scm_UWM_CLUBB_v2.nc', nc_write_mode)
% Redefine the file
%status = mexnc('redef',ncid)
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
tdimid = netcdf.defdim(ncid,'time', sizet);
%[tdimid,status] = mexnc('def_dim',ncid,'time',sizet);
% Full Levels (zt)
levhdimid = netcdf.defdim(ncid,'levh', w_nz);
%[levhdimid,status] = mexnc('def_dim',ncid,'levh',w_nz);
% Half Levels(zm)
levfdimid = netcdf.defdim(ncid,'levf', nz );
%[levfdimid,status] = mexnc('def_dim',ncid,'levf',nz);
% Soil Levels(sfc?)
levsdimid = netcdf.defdim(ncid,'levs', sfc_nz);
%[levsdimid,status] = mexnc('def_dim',ncid,'levs',sfc_nz);

% Define variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time series output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tvarid = define_variable( 'time' ,'seconds since 2006-07-01 12:00:00', 's', tdimid, ncid );
ldwvarid = define_variable( 'ldw' ,'long wave downward radiation at surface', 'W/m^2', tdimid, ncid );
lupvarid = define_variable( 'lup' ,'long wave upward radiation at surface', 'W/m^2', tdimid, ncid );
qdwvarid = define_variable( 'qdw' ,'short wave downward radiation at surface', 'W/m^2', tdimid, ncid );
qupvarid = define_variable( 'qup' ,'short wave upward radiation at surface', 'W/m^2', tdimid, ncid );
tskvarid = define_variable( 'tsk' ,'temperature skin layer', 'W/m^2', tdimid, ncid );
qskvarid = define_variable( 'qsk' ,'specific humiditiy top of vegetation layer', 'kg/kg', tdimid, ncid );
gvarid = define_variable( 'g' ,'soil heat flux', 'W/m^2', tdimid, ncid );
shfvarid = define_variable( 'shf' ,'sensible heat flux', 'W/m^2', tdimid, ncid );
lhfvarid = define_variable( 'lhf' ,'latent heat flux', 'W/m^2', tdimid, ncid );
ustarvarid = define_variable( 'ustar' ,'surface velocity', 'm/s', tdimid, ncid );
hpblvarid = define_variable( 'hpbl' ,'boundary layer height', 'm', tdimid, ncid );
t2mvarid = define_variable( 't2m' ,'2m temperature', 'K', tdimid, ncid );
q2mvarid =  define_variable( 'q2m' ,'2m specific humidity', 'kg/kg', tdimid, ncid );
u10mvarid = define_variable( 'u10m' ,'10m u-component wind', 'm/s', tdimid, ncid );
v10mvarid =  define_variable( 'v10m' ,'10m v-component wind', 'm/s', tdimid, ncid );
ccvarid = define_variable( 'cc' ,'cloudcover fration', '0 1', tdimid, ncid );
zfvarid = define_variable( 'zf' ,'height of full level', 'm', [tdimid levfdimid], ncid );
pfvarid = define_variable( 'pf' ,'pressure at full level', 'Pa', [tdimid levfdimid], ncid );
tkvarid = define_variable( 't' ,'temperature', 'K', [tdimid levfdimid], ncid );
thvarid = define_variable( 'th' ,'potential temperature', 'K', [tdimid levfdimid], ncid );
qvarid = define_variable( 'q' ,'specific humidity', 'kg/kg', [tdimid levfdimid], ncid );
uvarid = define_variable( 'u' ,'zonal component wind', 'm/s', [tdimid levfdimid], ncid );
vvarid = define_variable( 'v' ,'meridonal component wind', 'm/s', [tdimid levfdimid], ncid );
ugeovarid = define_variable( 'ugeo' ,'u-component geostrophic wind', 'm/s', [tdimid levfdimid], ncid );
vgeovarid = define_variable( 'vgeo' ,'v-component geostrophic wind', 'm/s', [tdimid levfdimid], ncid );
dudt_lsvarid = define_variable( 'dudt_ls' ,'u-component momentum advection', 'm/s/s', [tdimid levhdimid], ncid );
dvdt_lsvarid = define_variable( 'dvdt_ls' ,'v-component momentum advection', 'm/s/s', [tdimid levhdimid], ncid );
dtdt_lsvarid = define_variable( 'dtdt_ls' ,'temperature advection', 'K/s', [tdimid levhdimid], ncid );
dqdt_lsvarid = define_variable( 'dqdt_ls' ,'specific humidity advection', 'kg/kg/s', [tdimid levhdimid], ncid );
omevarid = define_variable( 'ome' ,'vertical movement', 'Pa/s', [tdimid levfdimid], ncid );
zhvarid = define_variable( 'zh' ,'height of half level', 'm', [tdimid levhdimid], ncid );
phvarid = define_variable( 'ph' ,'pressure at half level', 'Pa', [tdimid levhdimid], ncid );
wtvarid = define_variable( 'wt' ,'vertical temperature flux', 'Km/s', [tdimid levhdimid], ncid );
wqvarid = define_variable( 'wq' ,'vertical moisture flux', 'kg/kg m/s', [tdimid levhdimid], ncid );
uwvarid = define_variable( 'uw' ,'vertical flux u-component', 'm2/s2', [tdimid levhdimid], ncid );
vwvarid = define_variable( 'vw' ,'vertical flux v-component', 'm2/s2', [tdimid levhdimid], ncid );
Kmvarid = define_variable( 'Km' ,'eddy diffusivity momentum', 'm2/s', [tdimid levhdimid], ncid );
Khvarid = define_variable( 'Kh' ,'eddy diffusivity heat', 'm2/s', [tdimid levhdimid], ncid );
TKEvarid = define_variable( 'TKE' ,'turbulent kinetic energy', 'm^2/s^2', [tdimid levhdimid], ncid );
shearvarid = define_variable( 'shear' ,'shear production', 'm2/s3', [tdimid levhdimid], ncid );
buoyvarid = define_variable( 'buoy' ,'buoyancy production', 'm2/s3', [tdimid levhdimid], ncid );
transvarid = define_variable( 'trans' ,'total transport', 'm2/s3', [tdimid levhdimid], ncid );
dissivarid = define_variable( 'dissi' ,'dissipation', 'm2/s3', [tdimid levhdimid], ncid );
zsvarid = define_variable( 'zs' ,'height of soil level', 'm', [tdimid levsdimid], ncid );
tsvarid = define_variable( 'ts' ,'soil temperature', 'K', [tdimid levsdimid], ncid );
thsvarid = define_variable( 'ths' ,'soil water content', 'm2/s3', [tdimid levsdimid], ncid );%zeros


% 
% % Indicate which variables are not included
% status = mexnc('put_att_float',ncid,tskvarid,'_FillValue','float',1,NaN);
% status = mexnc('put_att_float',ncid ,gvarid,'_FillValue','float',1,NaN);
% status = mexnc('put_att_float',ncid ,Kmvarid,'_FillValue','float',1,NaN);
% %[Khvarid, status] = define_variable_2d( 'Kh' ,'eddy diffusivity heat', 'm2/s', tdimid, levf, ncid );
% status = mexnc('put_att_float',ncid ,Khvarid,'_FillValue','float',1,NaN);
% 
% 	%[transvarid, status] = define_variable_2d( 'trans' ,'total transport', 'm2/s3', tdimid, levf, ncid );
% 	% Zeros
% status = mexnc('put_att_float',ncid ,transvarid,'_FillValue','float',1,NaN);
% 
% 	%[dissivarid, status] = define_variable_2d( 'dissi' ,'dissipation', 'm2/s3', tdimid, levf, ncid );
% 	% Zeros
% status = mexnc('put_att_float',ncid ,dissivarid,'_FillValue','float',1,NaN);
% 
% 	%[thsvarid, status] = define_variable_2d( 'ths' ,'soil water content', 'm2/s3', tdimid, levf, ncid );
% status = mexnc('put_att_float',ncid ,thsvarid,'_FillValue','float',1,NaN);
% 
% status = mexnc('put_att_float',ncid ,hpblvarid,'_FillValue','float',1,NaN);
netcdf.endDef(ncid);
netcdf.setFill(ncid,'NC_FILL');
%status = mexnc('end_def',ncid);

% Apply values to file

for i=1:sizet
	% Example %
	% status = mexnc('VARPUT',ncid,thetalvarid,[i-1 k-1],[1 1],thlm_array(k,i),0);
    netcdf.putVar(ncid, tvarid, i-1, 1, i*10.0*60.0,0);
	%status = mexnc('VARPUT', ncid, tvarid, i-1,1,i*10.0*60.0,0);
	netcdf.putVar(ncid, ldwvarid, i-1,1,Frad_LW_down_array(1,i),0);
	netcdf.putVar( ncid, lupvarid, i-1,1,Frad_LW_up_array(1,i),0);
	netcdf.putVar( ncid, qdwvarid, i-1,1,Frad_SW_down_array(1,i),0);
	netcdf.putVar( ncid, qupvarid, i-1,1,Frad_SW_up_array(1,i),0);
	netcdf.putVar( ncid, shfvarid, i-1,1,sh_array(1,i),0);
	netcdf.putVar( ncid, lhfvarid, i-1,1,lh_array(1,i),0);
	netcdf.putVar( ncid, ustarvarid, i-1,1,ustar_array(1,i),0);
	%status = mexnc('VARPUT', ncid, hpblvarid, i-1,1,w_z(nz),0);
	netcdf.putVar( ncid, t2mvarid, i-1,1,T_in_K_array(1,i),0);
	netcdf.putVar( ncid, q2mvarid, i-1,1,qtm_array(1,i),0);
	netcdf.putVar( ncid, u10mvarid, i-1,1,um_array(2,i),0);
	netcdf.putVar( ncid, v10mvarid, i-1,1, vm_array(2,i),0);
	netcdf.putVar( ncid, ccvarid, i-1,1, cc_array(1,i),0);
end


for k=1:nz
    % k-1 comes from NetCDF starting variables at 0 and MATLAB starting
    % them at 1.
    
    % start, count, value, autoscale
%    status = mexnc('VARPUT',ncid,rhovarid,k-1,1,rho_array(k,1),0)
	for i=1:sizet
	
	%Mean State
	%Height at Full
    netcdf.putVar(ncid,zfvarid,[i-1 k-1],[1 1],z(k),0);
	%pressure at full level
	netcdf.putVar(ncid,pfvarid,[i-1 k-1],[1 1],p_in_Pa_array(k,i),0); % Grid switch
	%temperature
	netcdf.putVar(ncid,tkvarid,[i-1 k-1],[1 1],T_in_K_array(k,i),0); % Grid switch
	%potential temperature
	netcdf.putVar(ncid,thvarid,[i-1 k-1],[1 1],thlm_array(k,i),0); % Grid switch
	%specific humidity
	netcdf.putVar(ncid,qvarid,[i-1 k-1],[1 1],qtm_array(k,i),0); % Grid switch
	%zonal component wind;
	netcdf.putVar(ncid,uvarid,[i-1 k-1],[1 1],um_array(k,i),0); % Grid switch
	%meridonal component wind
	netcdf.putVar(ncid,vvarid,[i-1 k-1],[1 1],vm_array(k,i),0); % Grid switch

	%Fluxes
	%u-component geostrophic wind
	netcdf.putVar(ncid,ugeovarid,[i-1 k-1],[1 1],ug_array(k,i),0);
	%v-component geostrophic wind
	netcdf.putVar(ncid,vgeovarid,[i-1 k-1],[1 1],vg_array(k,i),0);	
	%u-component momentum advection;
	netcdf.putVar(ncid,dudt_lsvarid,[i-1 k-1],[1 1],um_f_array(k,i),0); % Grid switch
	%v-component momentum advection
	netcdf.putVar(ncid,dvdt_lsvarid,[i-1 k-1],[1 1],vm_f_array(k,i),0); % Grid switch
	%temperature advection
	netcdf.putVar(ncid,dtdt_lsvarid,[i-1 k-1],[1 1],T_forcing_array(k,i),0); % Grid switch
	%specific humidity advection
	netcdf.putVar(ncid,dqdt_lsvarid,[i-1 k-1],[1 1],rtm_f_array(k,i),0); % Grid switch
	%vertical movement
	netcdf.putVar(ncid,omevarid,[i-1 k-1],[1 1],ome_array(k,i),0); % Grid switch
	end
	
end

for k=1:w_nz
    % k-1 comes from NetCDF starting variables at 0 and MATLAB starting
    % them at 1.
    
    % start, count, value, autoscale
%    status = mexnc('VARPUT',ncid,rhovarid,k-1,1,rho_array(k,1),0)
	for i=1:sizet

	% Height at half level
	netcdf.putVar(ncid,zhvarid,[i-1 k-1],[1 1],w_z(k),0);
	
	% Pressure at half level ???
	
	%vertical temperature flux
	netcdf.putVar(ncid,wtvarid,[i-1 k-1],[1 1],wt_array(k,i),0);

	%vertical moisture flux ????
    netcdf.putVar(ncid,wqvarid,[i-1 k-1],[1 1],wq_array(k,i),0);
    %[wqvarid, status] = define_variable_2d( ncid,'wq' ,'vertical moisture flux', 'kg/kg m/s', tdimid, levhdimid, ncid );
	%vertical flux u-component
	netcdf.putVar(ncid,uwvarid,[i-1 k-1],[1 1],upwp_array(k,i),0);
	%vertical flux v-component
	netcdf.putVar(ncid,vwvarid,[i-1 k-1],[1 1],vpwp_array(k,i),0);
	%eddy diffusivity momentum
	% Zeros
	
	% Eddy diffusivity heat
	% Zeros
	
        % TKE
	netcdf.putVar(ncid, TKEvarid,[i-1 k-1],[1 1],em_array(k,i),0);

	% shear production
	netcdf.putVar(ncid, shearvarid,[i-1 k-1],[1 1],shear_array(k,i),0);
	% buoyancy production;
	netcdf.putVar(ncid, buoyvarid,[i-1 k-1],[1 1],wp2_bp_array(k,i),0);

	%total transport
	% Zeros

	% dissipation
	% Zeros

	end
	
end

for k=1:sfc_nz
	for i=1:sizet
		%[zsvarid, status] = define_variable_2d( 'zs' ,'height of soil level', 'm', tdimid, levf, ncid );
		netcdf.putVar(ncid,zsvarid,[i-1 k-1],[1 1],sfc_z(k),0); % Grid switch
		%[tsvarid, status] = define_variable_2d( 'ts' ,'soil temperature', 'K', tdimid, levf, ncid );
		netcdf.putVar(ncid,tsvarid,[i-1 k-1],[1 1],t_sfc_array(k,i),0); % Grid switch
		%[thsvarid, status] = define_variable_2d( 'ths' ,'soil water content', 'm2/s3', tdimid, levf, ncid );
	end
end

%[ndims,nvars,ngatts,unlimdim,status] = mexnc('inq',ncid)
%status = mexnc('close',ncid);
netcdf.close(ncid);
end

function varid = define_variable( shrt_name, long_name, units, dim_ids, file_id )
    varid = netcdf.defVar(file_id, shrt_name, 'NC_FLOAT',dim_ids);
	%[varid,status] = mexnc('def_var',file_id,shrt_name,'float',1,[dim_id]);
    netcdf.putAtt(file_id, varid,'unit',units);
	%d_status = mexnc('put_att_text',ncid, varid,'unit','char',length(units), units);
    netcdf.putAtt(file_id, varid,'long_name',long_name);
	%d_status = mexnc('put_att_text',ncid, varid,'long_name','char',length(long_name),long_name);
end 


