function[] = voca_output_creator( )
% VOCA_OUTPUT_CREATOR This function creates netCDF files required by the VOCA
% 3-d intercomparison. It uses WRF-CLUBB output files as source information.
%
%   This file is also meant to be an example for future MATLAB scripts to
%   convert WRF-CLUBB output for data submission in netCDF format.
%   Essentially a script for such a conversion will break down into these
%   sections:
%
%       File Input Section -
%           This is where the input files that are to be converted are
%           specified and read into MATLAB.
%
%       Interpolation, Calculation Section -
%           Since the input information produced by WRF-CLUBB may not match
%           one for one with the specifications of the output file, this
%           section is needed. Here variable interpolation to the correct
%           domain will occur and other necessary variable calculations
%           will be made such as density, temperature, etc.
%
%       Definition Section -
%           This is where netCDF definitions will be written to the output
%           file. This includes information such as variable names,
%           variable attributes, and global attributes.
%
%       Output File Section -
%           This section is responsible for writing variable data to the
%           output file and converting, calculating, etc. variables using
%           functions.
%


% Necessary include
addpath '../../matlab_include/' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File Input Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

curr_case = 'VOCA';

% set number of files
numfiles=3;

% define filenames
% These files are stored on hd1 in 
%   /share6/LS/ldgrant_a/VOCA_output/VOCA_runs_from_bluefire/VOCA_78/
% each file is in a subdirectory within VOCA_78/.  I copied these files
% over to steele in order to run this script.
WRFout_filename1 = 'VOCA_78/wrfout_d01_2008-10-12_00:00:00';
WRFout_filename2 = 'VOCA_78/wrfout_d01_2008-10-23_00:00:00';
WRFout_filename3 = 'VOCA_78/wrfout_d01_2008-11-05_12:00:00';
% Created file on 2011-7-15

% set up domain for output - do this here outside of the loop
% Define latitude, longitude range
% For inner domain: 0 to 40S by 1 degree intervals, 110W to 65W by 1-degree
voca_lats=[-40:1:0];
voca_lons=[-110:1:-65];

% create grid
[lats,lons]=meshgrid(voca_lats,voca_lons);

% define total number of times
% 8 3-hr intervals per day, and Oct 15 2008 00Z - Nov 16 2008 00Z = 32 days
totaltimes= 8.0*32.0+1; % = 257 times; for the output file, this will be 256 3-hr avg time intervals


% LOOP THROUGH FILES -- this is necessary for WRFout files from Bluefire
% runs, since they are split up into separate files completed with restart
% runs.  This loop ends after the interpolation to the voca domain.
for ifilenum=1:numfiles

if ifilenum==1
  % Open the WRF netCDF file and read information
  ncid = netcdf.open(WRFout_filename1,'NOWRITE'); % use ncid to reference WRF file
elseif ifilenum==2
  ncid = netcdf.open(WRFout_filename2,'NOWRITE');
elseif ifilenum==3
  ncid = netcdf.open(WRFout_filename3,'NOWRITE');
end
        

[numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);

% This loops through all dimensions and sets their values
% e.g. bottom_top = 44, the number of unstaggered vertical levels 
for i=1:numdims
    [dimname,dimlen] = netcdf.inqDim(ncid,i-1);
    eval([dimname,' = dimlen']) % print out dimensions
end

% Rename Time dimension name since it needs to be used later
numTimes=Time; clear Time

% Array of WRF variables to get
% Only get variables necessary since WRF files are so large
WRFvars = [
    'Times   ' % 19-length time characters: yyyy-mm-dd_hh:mm:ss
    'ZNU     ' % eta values unstaggered in z, half (mass) levels
    'ZNW     ' % eta values staggered in z, full (w) levels
    'P_TOP   ' % Pressure of model top, Pa
    'U       ' % u-wind, m/s (staggered in x)
    'V       ' % v-wind, m/s (staggerd in y)
    'W       ' % vertical velocity, m/s (staggered in z)
    'PH      ' % perturbation geopotential, m2/s2 (staggered in z)
    'PHB     ' % base-state geopotential, m2/s2 (staggered in z)
    'T       ' % perturbation potential temp, K (base state always 300 K)
    'P       ' % perturbation pressure, Pa
    'PB      ' % base state pressure, Pa
    'PSFC    ' % surface pressure, Pa
    'Q2      ' % 2-m Qv, kg/kg
    'T2      ' % 2-m Temp, K
    'U10     ' % 10-m u-wind, m/s
    'V10     ' % 10-m v-wind, m/s
    'QVAPOR  ' % water vapor mix ratio, kg/kg
    'QCLOUD  ' % cloud water mix ratio, kg/kg
    'QRAIN   ' % rain water mix ratio, kg/kg
    'QICE    ' % ice water mix ratio, kg/kg
    'RAINNC  ' % accumulated total grid scale precip, mm
    'LANDMASK' % 1=land, 0=water
    'CF_CLUBB' % clubb cloud frac
    'HGT     ' % terrain height, m
    'TSK     ' % surface skin temp, K
    'XLAT    ' % latitude
    'XLONG   ' % longitude
    'XLAT_U  ' % latitude staggered in x-dir
    'XLONG_U ' % longitude staggered in x-dir
    'XLAT_V  ' % latitude staggered in y-dir
    'XLONG_V ' % longitude staggered in y-dir
    'HFX     ' % Upward heat flux at sfc, W/m2
    'LH      ' % Latent heat flux at sfc, W/m2
    'OLR     ' % TOA outgoing LW, W/m2
    'GLW     ' % Sfc downward LW flux, W/m2
    'SWDOWN  ' % Sfc downward SW flux, W/m2
    ];

for i=1:numvars
    varname = netcdf.inqVar(ncid,i-1);
    getvar='F';
    for j=1:size(WRFvars,1)
        if strcmp( varname, strtrim(WRFvars(j,:)) ) == 1;
            getvar='T';
        end
    end
    
    if strcmp(getvar,'T') == 1
        data = netcdf.getVar(ncid,i-1);
        eval([varname,' = data;']);
    end
end
clear getvar varname data;

% Close the WRF output file which was read
netcdf.close(ncid)
clear ncid
disp('Done reading in WRF variables for file:'); disp(ifilenum);

% compress variables that don't vary in time
XLAT1(:,:)=XLAT(:,:,1); clear XLAT;
XLONG1(:,:)=XLONG(:,:,1); clear XLONG;
XLAT1_U(:,:)=XLAT_U(:,:,1); clear XLAT_U;
XLONG1_U(:,:)=XLONG_U(:,:,1); clear XLONG_U;
XLAT1_V(:,:)=XLAT_V(:,:,1); clear XLAT_V;
XLONG1_V(:,:)=XLONG_V(:,:,1); clear XLONG_V;
HGT1(:,:)=HGT(:,:,1); clear HGT;
LANDMASK1(:,:)=LANDMASK(:,:,1); clear LANDMASK
eta(:,1)=ZNU(:,1); clear ZNU
eta_stag(:,1)=ZNW(:,1); clear ZNW;
Ptop=P_TOP(1); clear P_TOP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Interpolation, Calculation Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ifilenum==1
  % Find index of initial time (Oct 15)
  % This checks the 10th row of Times, which is the second
  % digit of the date, and returns the first index for which 
  % the 2nd digit of the date is > 4, i.e. finds the index of
  % Oct 15 00Z.  This works if the model start time is not before Oct 06
  % Only do this if reading the first file
  tWRFstart = find(Times(10,:)=='5',1,'first');
  tWRFend = numTimes;
  
  tstart1=1; tend1=numTimes-tWRFstart+1;
  
  tstart=tstart1; tend=tend1; % set for loop
elseif ifilenum==2
  tWRFstart = 2; % because time 1 in this wrfout file is identical to the end time from the previous file
  tWRFend = numTimes;
  
  tstart2=tend1+1; tend2=tstart2+numTimes-2;
  
  tstart=tstart2; tend=tend2; % set for loop
elseif ifilenum==3
  tWRFstart = 2; % because time 1 in this wrfout file is identical to the end time from the previous file
  tWRFend = numTimes;
  
  tstart3=tend2+1; tend3=tstart3+numTimes-2;
  
  tstart=tstart3; tend=tend3; % set for loop
end


% allocate space for speed, but allocate for the entire duration
% (totaltimes); therefore, only do this on first pass through loop
if ifilenum==1
  u=zeros(length(voca_lons),length(voca_lats),bottom_top,totaltimes);
  v=u; theta=u; press=u; qvapor=u; qcloud=u; qrain=u; qice=u; cf=u;
end
% Interpolate data to lats,lons domain and clear other arrays for memory
for n = tstart:tend
  % unstaggered variables
  nwrf=n-tstart+tWRFstart;
  for k = 1:bottom_top % # of unstagerred vertical levels
      % U and V are staggered in x- and y- directions as well
      u(:,:,k,n) = interp2(XLAT1_U,XLONG1_U,U(:,:,k,nwrf),lats,lons); 
      v(:,:,k,n) = interp2(XLAT1_V,XLONG1_V,V(:,:,k,nwrf),lats,lons); 
      % Potential temp = base state theta (=300K) + perturbation theta
      theta(:,:,k,n) = interp2(XLAT1,XLONG1,T(:,:,k,nwrf)+300,lats,lons); 
      % pressure = base state + perturbation
      press(:,:,k,n) = interp2(XLAT1,XLONG1,P(:,:,k,nwrf)+PB(:,:,k,nwrf),lats,lons); 
      qvapor(:,:,k,n) = interp2(XLAT1,XLONG1,QVAPOR(:,:,k,nwrf),lats,lons); 
      qcloud(:,:,k,n) = interp2(XLAT1,XLONG1,QCLOUD(:,:,k,nwrf),lats,lons); 
      qrain(:,:,k,n) = interp2(XLAT1,XLONG1,QRAIN(:,:,k,nwrf),lats,lons); 
      qice(:,:,k,n) = interp2(XLAT1,XLONG1,QICE(:,:,k,nwrf),lats,lons);
      cf(:,:,k,n) = interp2(XLAT1,XLONG1,CF_CLUBB(:,:,k,nwrf),lats,lons);
  end
end
clear U V T P PB QVAPOR QCLOUD QRAIN QICE CF_CLUBB % for space
disp('Done interpolating unstaggered 3-D variables to VOCA domain for file:'); disp(ifilenum);

if ifilenum==1
  q2=zeros(length(voca_lons),length(voca_lats),totaltimes);
  psfc=q2; t2=q2; u10=q2; v10=q2; glw=q2; hfx=q2; hgt=q2; lh=q2; olr=q2;
    swdown=q2; tsk=q2; landmask=q2; rainnc=q2;
end
% 2-D unstaggered variables
for n = tstart:tend
  nwrf=n-tstart+tWRFstart;
    psfc(:,:,n) = interp2(XLAT1,XLONG1,PSFC(:,:,nwrf),lats,lons);
    q2(:,:,n) = interp2(XLAT1,XLONG1,Q2(:,:,nwrf),lats,lons); 
    t2(:,:,n) = interp2(XLAT1,XLONG1,T2(:,:,nwrf),lats,lons); 
    u10(:,:,n) = interp2(XLAT1,XLONG1,U10(:,:,nwrf),lats,lons);
    v10(:,:,n) = interp2(XLAT1,XLONG1,V10(:,:,nwrf),lats,lons);
    glw(:,:,n) = interp2(XLAT1,XLONG1,GLW(:,:,nwrf),lats,lons);
    hfx(:,:,n) = interp2(XLAT1,XLONG1,HFX(:,:,nwrf),lats,lons);
    hgt(:,:,n) = interp2(XLAT1,XLONG1,HGT1,lats,lons);
    lh(:,:,n) = interp2(XLAT1,XLONG1,LH(:,:,nwrf),lats,lons);
    olr(:,:,n) = interp2(XLAT1,XLONG1,OLR(:,:,nwrf),lats,lons);
    swdown(:,:,n) = interp2(XLAT1,XLONG1,SWDOWN(:,:,nwrf),lats,lons);
    tsk(:,:,n) = interp2(XLAT1,XLONG1,TSK(:,:,nwrf),lats,lons);
    landmask(:,:,n) = interp2(XLAT1,XLONG1,LANDMASK1,lats,lons);
    rainnc(:,:,n) = interp2(XLAT1,XLONG1,RAINNC(:,:,nwrf),lats,lons);
end
clear PSFC Q2 T2 U10 V10 GLW HFX HGT1 LH OLR SWDOWN TSK LANDMASK1 RAINNC % for space

if ifilenum==1
  w=zeros(length(voca_lons),length(voca_lats),bottom_top_stag,totaltimes);
  phi=w;
end
for n = tstart:tend
  % staggered variables
  nwrf=n-tstart+tWRFstart;
  for k = 1:bottom_top_stag % # of staggered vertical levels
      w(:,:,k,n) = interp2(XLAT1,XLONG1,W(:,:,k,nwrf),lats,lons); 
      % geopotential, phi, = base state geopotential, PH, + perturbation geopotential, PHB
      phi(:,:,k,n) = interp2(XLAT1,XLONG1,PH(:,:,k,nwrf)+PHB(:,:,k,nwrf),lats,lons); 
  end
end
clear W PH PHB

disp('Done interpolating variables to VOCA domain for file:'); disp(ifilenum);
%keyboard
end % for loop around file names!


% calculate a few necessary variables
Rd=287.; % J/kg/K
Cp=1005.; % J/kg/K
epsilon=0.622;
gravity=9.81; % m/s2, value used in WRF

temp=theta.*(press/100000).^(Rd/Cp); % equation for potential temp
tvirt = temp .* ( (1.+qvapor/epsilon)./(1.+qvapor) ); % equation for virtual temperature
rho = press./(Rd*tvirt); % equation of state

geopot_ht = phi/gravity; % m; ~= height

% calculate pressure on staggered grid
press_stag = calc_pressure_staggered( psfc, Ptop, eta_stag );

% calculate low cloud fraction below 700mb
cldlow=calc_low_cloud_frac( cf, press );

% VOCA spec wants average over 0-3Z, 3-6Z, etc. so start at hour 1.5 and
% add to days vector in 3-hour increments
daysarray = [ 1.5/24.0 : 3.0/24.0 : 3.0/24.0*(totaltimes-1) ];

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Definition Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the new files. By default they are in definition mode.
ncid = netcdf.create( 'WRF-CLUBB_UWM_outer_MET.nc','NC_WRITE' );

% Define Global Attributes
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','UWM WRF-CLUBB MET simulation for the VOCA Modeling Experiment');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Contact','Leah Grant (ldgrant@atmos.colostate.edu) and Vincent Larson (vlarson@uwm.edu)');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'WRF_model_version','3.1.1');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Model_time_step','one minute');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Horizontal_grid_spacing','25 km');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'PBL_scheme','CLUBB, a 1-D cloud and turbulence parameterization');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Microphysics_scheme','Morrison double-moment');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Radiation_schemes','RRTM Longwave and Dudhia Shortwave');

% Define dimensions
% time, zlevels, lat, lon
% Time has one less value than vars defined above because VOCA spec wants
% output averaged over 00-3Z, 3-6Z, etc.
timedimid = netcdf.defdim(ncid,'time', totaltimes-1 );
latdimid = netcdf.defdim(ncid,'lat', length(voca_lats) );
londimid = netcdf.defdim(ncid,'lon', length(voca_lons) );
zlevdimid = netcdf.defdim(ncid,'zlevels', bottom_top);
zlev_stagdimid = netcdf.defdim(ncid,'zlevels_stag',bottom_top_stag);

% Define variables:  function define_var( varname, dscrptn, unit, dim_ids, file_id )

% Variables not included because they aren't available outputs from the
% setup of our VOCA simulation (e.g. not available because we'd need WRF-Chem,
% not output from the radiation/microphysics scheme we used, etc.):
%
% emission_profile, surface momentum fluxes in x- and y- directions?, 
% DMS flux, SSFlux, BCEmiss, OCEmiss, S02Emiss, COEmiss,
% SWUPTOA, LWCF, SWCF, PFLUX, AOD, REFFTOP,
% reff, CCN, SO4, SS, BC, OC, SUMAERO, SUMAERO_PM1, SUMAERO_PM25,
% SUMAERO_N, DMS, S02, CO, O3

% 1-D variables
timevarid = define_var('time','Days since 00Z 15 Oct. 2008','Days since 2008-10-15 00:00:00Z',timedimid,ncid);
latvarid = define_var('lat','Latitude','Degrees north (south is negative)',latdimid,ncid);
lonvarid = define_var('lon','Longitude','Degrees east (West is negative)', londimid, ncid);
sigmavarid = define_var('sigma','sigma: = (P-Ptop)/(Psfc-Ptop) = eta in WRF, where Ptop = 5000 Pa', '[0-1]', zlev_stagdimid,ncid);

% 2-D variables
tlatlon_ids=[londimid latdimid timedimid];
% Define missing value
miss_val=-9999.;

TSvarid = define_var('TS','Surface radiative skin temp (TSK in WRF)','K',tlatlon_ids,ncid);
ZSvarid = define_var('ZS','Surface height (HGT in WRF)','m',tlatlon_ids,ncid);
landfracvarid = define_var('land_fraction','Land fraction (LANDMASK in WRF)','[0-1]',tlatlon_ids,ncid);
PSvarid = define_var('PS','Model surface pressure (PSFC in WRF)','Pa',tlatlon_ids,ncid);
%SLPvarid = define_var('SLP','Sea level pressure (computed from .. in WRF)','Pa',tlatlon_ids,ncid);
T2varid = define_var('T2','Temperature at 2m (T2 in WRF)','K',tlatlon_ids,ncid);
q2varid = define_var('q2','Specific humidity at 2m (computed from Q2 in WRF)','kg/kg',tlatlon_ids,ncid);
u10varid = define_var('u10','u-wind at 10 m (U10 in WRF)','m/s',tlatlon_ids,ncid);
v10varid = define_var('v10','v-wind at 10 m (V10 in WRF)','m/s',tlatlon_ids,ncid);
u850varid = define_var('u850','u-wind at 850 mb (interpolated from U in WRF)','m/s',tlatlon_ids,ncid);
  netcdf.putAtt(ncid, u850varid,'missing_value',miss_val);
v850varid = define_var('v850','v-wind at 850 mb (interpolated from V in WRF)','m/s',tlatlon_ids,ncid);
  netcdf.putAtt(ncid, v850varid,'missing_value',miss_val);
w850varid = define_var('w850','w-wind at 850 mb (interpolated from W in WRF)','m/s',tlatlon_ids,ncid);
  netcdf.putAtt(ncid, w850varid,'missing_value',miss_val);
SHFvarid = define_var('SHF','Sensible heat flux (HFX in WRF)','W/m2',tlatlon_ids,ncid);
LHFvarid = define_var('LHF','Latent heat flux (LH in WRF)','W/m2',tlatlon_ids,ncid);
LWUPTOAvarid = define_var('LWUPTOA','TOA upward longwave radiative flux (OLR in WRF)','W/m2',tlatlon_ids,ncid);
% SWUPTOAvarid = define_var('SWUPTOA','TOA upward shortwave radiative flux (not avail output from Dudhia scheme)','W/m2',tlatlon_ids,ncid);
LWDNSFCvarid = define_var('LWDNSFC','Surface downward longwave radiative flux (GLW in WRF)','W/m2',tlatlon_ids,ncid);
SWDNSFCvarid = define_var('SWDNSFC','Surface downward shortwave radiative flux (SWDOWN in WRF)','W/m2',tlatlon_ids, ncid);
WVPvarid = define_var('WVP','Water vapor path (computed from QVAPOR, T, P+PB, and PH+PHB in WRF)','kg/m2',tlatlon_ids,ncid);
LWPvarid = define_var('LWP','Liquid water path (computed from QCLOUD, QRAIN, T, P+PB, PH+PHB, and QVAPOR in WRF)','kg/m2',tlatlon_ids,ncid);
CLDLOWvarid = define_var('CLDLOW','Low cloud fraction below 700 mb (derived from CF_CLUBB in WRF-CLUBB)','[0-1]',tlatlon_ids,ncid);
PRECTvarid = define_var('PRECT','Surface precipitation (computed from RAINNC in WRF)','kg H20 / (m2 s)',tlatlon_ids,ncid);
% PFLUXvarid = define_var('PFLUX','Precipitation flux at 500m','kg H20 (m2 s)',tlatlon_ids,ncid); % not available output from this simulation
TCLDTOPvarid = define_var('TCLDTOP','Low cloud IR brightness temperature proxy (computed from CF_CLUBB, QCLOUD, and Temperature in WRF)','K',tlatlon_ids,ncid);
  netcdf.putAtt(ncid, TCLDTOPvarid,'missing_value',miss_val);
PCLDTOPvarid = define_var('PCLDTOP','Low cloud top height proxy (computed from CF_CLUBB, QCLOUD, and P+PB in WRF)','Pa',tlatlon_ids,ncid);
  netcdf.putAtt(ncid, PCLDTOPvarid,'missing_value',miss_val);
NDTOPvarid = define_var('NDTOP','Cloud droplet concentration at low cloud top (determined from CF_CLUBB and QCLOUD in WRF)','cm-3',tlatlon_ids,ncid);
  netcdf.putAtt(ncid, NDTOPvarid,'missing_value',miss_val);

% 3-D variables
alldim_ids = [londimid latdimid zlevdimid timedimid];
alldim_stag_ids = [londimid latdimid zlev_stagdimid timedimid];
Tvarid = define_var('T','Temperature (computed from T and P+PB in WRF)','K',alldim_ids,ncid);
qvvarid = define_var('qv','Specific humidity (computed from QVAPOR in WRF)','kg/kg',alldim_ids,ncid);
qlvarid = define_var('ql','Cloud liquid water content (QCLOUD in WRF)','kg/kg',alldim_ids,ncid);
qivarid = define_var('qi','Cloud ice water content (QICE in WRF)','kg/kg',alldim_ids,ncid);
cfvarid = define_var('cf','Cloud fraction (CF_CLUBB in WRF-CLUBB, diagnosed in CLUBB PBL scheme)','[0-1]',alldim_ids,ncid);
Nvarid = define_var('N','Cloud droplet number concentration (constant, = 150 cm-3 in WRFv3.1.1 version of Morrison microphysics)','cm-3',alldim_ids,ncid);
  netcdf.putAtt(ncid, Nvarid,'missing_value',miss_val);
Uvarid = define_var('u','U-wind (U in WRF)','m/s',alldim_ids,ncid);
Vvarid = define_var('v','V-wind (V in WRF)','m/s',alldim_ids,ncid);
Wvarid = define_var('w','W-wind (W in WRF)','m/s',alldim_stag_ids,ncid);
Zvarid = define_var('Z','Geopotential height ((PH+PHB)/g in WRF)','m',alldim_stag_ids,ncid);
Pvarid = define_var('P','Pressure (P+PB in WRF)','Pa',alldim_ids,ncid);

% Set the file to add values
netcdf.setFill(ncid,'NC_FILL');
netcdf.endDef(ncid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Output File Section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output data to the file
% 1-D vars
netcdf.putVar(ncid,timevarid,daysarray);
netcdf.putVar(ncid,latvarid,voca_lats);
netcdf.putVar(ncid,lonvarid,voca_lons);
netcdf.putVar(ncid,sigmavarid,eta_stag);
disp('Done writing 1-D variables')

% 2-D vars
netcdf.putVar( ncid,TSvarid,time_avg(tsk) );
netcdf.putVar( ncid,ZSvarid,time_avg(hgt) );
netcdf.putVar( ncid,landfracvarid,time_avg(landmask) );
netcdf.putVar( ncid,PSvarid,time_avg(psfc) );
%netcdf.putVar( ncid,SLPvarid,time_avg( calc_slp(  ) ) );
netcdf.putVar( ncid,T2varid,time_avg(t2) );
netcdf.putVar( ncid,q2varid,mix_rat_to_spec_humid(time_avg(q2)) );
netcdf.putVar( ncid,u10varid,time_avg(u10) );
netcdf.putVar( ncid,v10varid,time_avg(v10) );
netcdf.putVar( ncid,u850varid,time_avg( interp_to_850mb( u, press, miss_val ) ) );
netcdf.putVar( ncid,v850varid,time_avg( interp_to_850mb( v, press, miss_val ) ) );
netcdf.putVar( ncid,w850varid,time_avg( interp_to_850mb( w, press_stag, miss_val ) ) );
netcdf.putVar( ncid,SHFvarid,time_avg(hfx) );
netcdf.putVar( ncid,LHFvarid,time_avg(lh) );
netcdf.putVar( ncid,LWUPTOAvarid,time_avg(olr) );
netcdf.putVar( ncid,LWDNSFCvarid,time_avg(glw) );
netcdf.putVar( ncid,SWDNSFCvarid,time_avg(swdown) );
netcdf.putVar( ncid,WVPvarid,time_avg( calc_water_path( rho,qvapor,geopot_ht ) ) );
netcdf.putVar( ncid,LWPvarid,time_avg( calc_water_path( rho,qcloud+qrain,geopot_ht ) ) );
netcdf.putVar( ncid,CLDLOWvarid,time_avg( cldlow ) );
netcdf.putVar( ncid,PRECTvarid,convert_precip_flux( rainnc ) ); % function convert_precip_flux does time-avg
netcdf.putVar( ncid,TCLDTOPvarid, calc_t_p_nd_cldtop( temp, cldlow, qcloud, miss_val ) ); % this function does the time-avg
netcdf.putVar( ncid,PCLDTOPvarid, calc_t_p_nd_cldtop( press, cldlow, qcloud, miss_val ) );
  % create ND array, lat x lon x unstaggered height x time; just use size of cloud fraction array
  ND_array=zeros(size(cf));
  ND_array(:,:,:,:)=150.0; % 150 cm^-3, constant value
netcdf.putVar( ncid,NDTOPvarid, calc_t_p_nd_cldtop( ND_array, cldlow, qcloud, miss_val ) );
disp('Done writing 2-D variables')

% 3-D vars
netcdf.putVar( ncid,Tvarid,time_avg(temp) );
netcdf.putVar( ncid,qvvarid,mix_rat_to_spec_humid(time_avg(qvapor)) );
netcdf.putVar( ncid,qlvarid,time_avg(qcloud) );
netcdf.putVar( ncid,qivarid,time_avg(qice) );
netcdf.putVar( ncid,cfvarid,time_avg(cf) );
netcdf.putVar( ncid,Nvarid,fill_miss_values_ND_array( time_avg(ND_array),time_avg(qcloud),miss_val ) );
netcdf.putVar( ncid,Uvarid,time_avg(u) );
netcdf.putVar( ncid,Vvarid,time_avg(v) );
netcdf.putVar( ncid,Wvarid,time_avg(w) );
netcdf.putVar( ncid,Zvarid,time_avg(geopot_ht) );
netcdf.putVar( ncid,Pvarid,time_avg(press) );
disp('Done writing 3-D variables')

% Close file so it is written correctly
netcdf.close(ncid);

keyboard % pause if you want to plot anything
% clear; clc; close all; % clear memory
end % end function voca_output_creater

% Functions used to convert units, calculate variables, etc.

function varid = define_var( varname, dscrptn, unit, dim_ids, file_id )

varid = netcdf.defVar(file_id, varname, 'NC_FLOAT',dim_ids);
netcdf.putAtt(file_id, varid,'description',dscrptn);
netcdf.putAtt(file_id, varid,'units',unit);

end

function voca_3hr_avg_var = time_avg( var )
  miss_val=-9999.; % same as above
  if size(size(var),2)==4 % 4-D array
    % preallocate array for speed
    voca_3hr_avg_var = zeros(size(var,1),size(var,2),size(var,3),size(var,4)-1);
    for i=1:size(var,1)
    for j=1:size(var,2)
    for k=1:size(var,3);
      for itime=1:size(var,4)-1 % 4th dimension: time
        % first check if either the value of var at this time or the next
        % time is the missing value; if so, we don't want to average it
        if var(i,j,k,itime)==miss_val || var(i,j,k,itime+1)==miss_val
            voca_3hr_avg_var(i,j,k,itime)=miss_val;
        else
        % Average over only two points in time, since output every three
        % hours at 00Z, 3Z, etc. and we want time avg. over 0-3Z, 3-6Z,
        % etc.  This reduces the array size by one in the time-dimension
        voca_3hr_avg_var(i,j,k,itime) = 0.5*( var(i,j,k,itime) + var(i,j,k,itime+1) );
        end
      end % time loop
    end
    end
    end
  else % 3-D array
    voca_3hr_avg_var = zeros(size(var,1),size(var,2),size(var,3)-1);
    for i=1:size(var,1)
    for j=1:size(var,2);
      for itime=1:size(var,3)-1
        if var(i,j,itime)==miss_val || var(i,j,itime+1)==miss_val
            voca_3hr_avg_var(i,j,itime)=miss_val;
        else
            voca_3hr_avg_var(i,j,itime) = 0.5*( var(i,j,itime) + var(i,j,itime+1) );
        end
      end % time loop
    end
    end
  end % if statement to determine if var is 3-d or 4-d
end % function

function spec_humid = mix_rat_to_spec_humid( mix_rat )
    spec_humid = mix_rat ./ (1. + mix_rat );
end

function water_path = calc_water_path( rho,mix_rat,heights )
    % allocate space first
    water_path = zeros(size(rho,1),size(rho,2),size(rho,4)); % lat x lon x time
    % note: heights is approximated with geopotential heights, which are
    % staggered in z - so values of rho and mix_rat
    % are between height levels.  Therefore approx dz in the equation with
    % the difference in height at grid levels surrounding rho and rliq.
    % This might not be as accurate higher aloft, where the grid spacing is
    % more stretched, but VOCA is more concerned with Sc clouds in lower
    % atm
    for k=1:size(rho,3) % height dimension
        integral(:,:,:) = ...
            rho(:,:,k,:) .* mix_rat(:,:,k,:) .* ( heights(:,:,k+1,:)-heights(:,:,k,:) );
        water_path(:,:,:) = water_path(:,:,:) + integral(:,:,:);
        clear integral;
    end
end

function var850mb = interp_to_850mb( var, pressure, missing_value )
    % used for u, v, and w, which are 4-D arrays
    
    % first allocate space; lat x lon x time
    var850mb=zeros(size(var,1),size(var,2),size(var,4));
    for i=1:size(var,1)
    for j=1:size(var,2)
    for n=1:size(var,4)
        if pressure(i,j,1,n) < 85000. % Pa
            var850mb(i,j,n) = missing_value;
        elseif pressure(i,j,1,n) == 85000.
            var850mb(i,j,n) = var(i,j,1,n);
        else
            for k=2:size(var,3)
                if pressure(i,j,k,n) == 85000.;
                    var850mb(i,j,n) = var(i,j,k,n);
                elseif pressure(i,j,k,n) < 85000. 
                    klower=k-1;
                    kupper=k;
                    % interpolation formula for arbitrary pressure level
                    % 850mb somewhere between two other pressure levels
                    var850mb(i,j,n) = var(i,j,klower,n) + ...
                        ( pressure(i,j,klower,n)-85000. ) / ( pressure(i,j,klower,n)-pressure(i,j,kupper,n) ) * ...
                        ( var(i,j,kupper,n)-var(i,j,klower,n) );
                    break
                end
            end
        end
    end
    end
    end
end

function pressure_staggered = calc_pressure_staggered( psfc, ptop, eta_staggered )
    % eta = (p-ptop)/(psfc-ptop), in WRF, so 
    % p = eta*(psfc-ptop) + ptop
    
    % ptop a constant
    % eta varies in z
    % pressure varies in lat, lon, z, and time
    % psfc varies in lat, lon, and time
    
    % allocate space
    pressure_staggered = zeros(size(psfc,1),size(psfc,2),length(eta_staggered),size(psfc,3));
    
    for k=1:length(eta_staggered)
        pressure_staggered(:,:,k,:) = eta_staggered(k) * ( psfc(:,:,:)-ptop ) + ptop;
    end
end

function low_cloud_frac = calc_low_cloud_frac( cloud_frac, pressure )
    % calculate cloud fraction below 700mb
    % note: this calculates the column-maximum cloud fraction below 700mb
    
    % allocate space
    low_cloud_frac = zeros( size(cloud_frac,1),size(cloud_frac,2),size(cloud_frac,4) );
    for i=1:size(cloud_frac,1)
    for j=1:size(cloud_frac,2)
    for n=1:size(cloud_frac,4)
        if pressure(i,j,1,n) > 70000. % Pa ; skip this column if psfc<700mb since already set low_cloud_frac to zeros
            for k=1:size(cloud_frac,3)
                if low_cloud_frac(i,j,n) < cloud_frac(i,j,k,n)
                    low_cloud_frac(i,j,n) = cloud_frac(i,j,k,n); % this sets low_cloud_frac to the max value below 700mb
                    if pressure(i,j,k+1,n) < 70000.
                        break % leave loop if p<700mb
                    end
                end
            end
        end
    end
    end
    end
end

function sfc_precip_flux = convert_precip_flux( rainnc_mm )
    % sfc_precip, rainnc, in mm, but we need it in units of flux:
    % kg H20/(m2 s)
    rainnc_m = rainnc_mm / 1000.; % 1000 mm/m
    % rainnc_m units = m = m^3/m^2
    % rainnc_m * density_of_water units = kg / m^2
    rhowater = 1000.; % kg/m^3
    rainnc_kgm2 = rainnc_m * rhowater;
    
    % Loop over time to estimate flux, quantity/unit area/unit time
    % first allocate space
    sfc_precip_flux=zeros(size(rainnc_mm,1),size(rainnc_mm,2),size(rainnc_mm,3)-1);
    for n=2:size(rainnc_mm,3) % 3 is time dimension
        sfc_precip_flux(:,:,n-1) = ( rainnc_kgm2(:,:,n)-rainnc_kgm2(:,:,n-1) ) / (3*3600); % time output interval is 3 hours, 3600 s per hour
        % final units are kg/m2/s
    end
end

function final_nd_avg = fill_miss_values_ND_array( nd_avg, qcloud_avg, miss_val )
    % fill missing value if 3-hr time mean qcloud in grid cell is < 0.1
    % g/kg
    final_nd_avg=nd_avg;
    for i=1:size(nd_avg,1)
    for j=1:size(nd_avg,2)
    for k=1:size(nd_avg,3)
    for n=1:size(nd_avg,4)
        if qcloud_avg(i,j,k,n) < 0.1/1000. % < 0.1 g/kg
            final_nd_avg(i,j,k,n) = miss_val;
        end
    end
    end
    end
    end
end

function t_p_nd_cldtop_timeavg = calc_t_p_nd_cldtop( temp_press_nd, cldlow, qcloud, miss_val )
    cldtoplev=find_cldtop_lev( qcloud );
    t_p_nd_cldtop = find_cldtop_var( temp_press_nd, cldtoplev, miss_val );
    % temperatures or pressures or cloud droplet number concentrations
    % at corresponding grid points of cloud top levels
    % missing value if entire column has qcloud < 0.1 g/kg
    
    cldlow_avg=time_avg(cldlow); t_p_nd_cldtop_avg=time_avg(t_p_nd_cldtop);
    % also fill in missing value if cldlow_avg < 0.8
    t_p_nd_cldtop_timeavg = fill_miss_values_from_lowcld_frac( t_p_nd_cldtop_avg, cldlow_avg, miss_val );
end

function cldtop_lev = find_cldtop_lev( qcloud )
    % highest grid point where qcloud exceeds 0.1 g/kg
    % zero if qcloud < 0.1 g/kg for the whole column
    cldtop_lev=zeros(size(qcloud,1),size(qcloud,2),size(qcloud,4));
    for i=1:size(qcloud,1)
    for j=1:size(qcloud,2)
    for n=1:size(qcloud,4)
        for k=1:size(qcloud,3)
            if qcloud(i,j,k,n)>0.1/1000. % qcloud > .1 g/kg
                cldtop_lev(i,j,n)=k;
            end
        end
    end
    end
    end
end

function cldtop_var = find_cldtop_var( var, cldtoplev, miss_val )
    cldtop_var=zeros(size(var,1),size(var,2),size(var,4)); % allocate space
    for i=1:size(var,1)
    for j=1:size(var,2)
    for n=1:size(var,4)
        if cldtoplev(i,j,n)==0.
            cldtop_var(i,j,n)=miss_val; % missing value if entire column has qcloud < 0.1 g/kg
        else
            cldtop_var(i,j,n)=var(i,j,cldtoplev(i,j,n),n);
        end
    end
    end
    end
end

function var_avg_final = fill_miss_values_from_lowcld_frac( var_avg, lowcld_avg, miss_val )
    % write a missing value if 3-hr low cloud fraction < 0.8
    var_avg_final=var_avg;
    for i=1:size(var_avg,1)
    for j=1:size(var_avg,2)
    for n=1:size(var_avg,3)
        if lowcld_avg(i,j,n) < 0.8
            var_avg_final(i,j,n)=miss_val;
        end
    end
    end
    end
end


% plotting function

function contourf_plot( twod_var )
    % flip variable around so it plots correctly according to latitudes and
    % longitudes
    twodvarflip=flip_var_for_plot( twod_var );
    %landmaskflip=flip_var_for_plot( landmask );

    %figure
    contourf(twodvarflip); colorbar;
    %hold on
    %contour(landmaskflip,'w');

end

function flipped_var = flip_var_for_plot( var )
    flipped_var=zeros(size(var,2),size(var,1));
    for i=1:size(var,1)
        for j=1:size(var,2)
            flipped_var(j,i)=var(i,j);
        end
    end
end

function cfmax = calc_cfmax( cf )
    cfmax = squeeze(max(cf,[],3));
end