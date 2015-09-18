% $Id$
function [ varname_sam, units_corrector_type, ...
           num_var_sam, num_tot_var_sam ] = output_vars_3D( casename )

idx_count = 0;

%============= Add SAM Model Field Variables Here =========================

% w (m/s).
idx_count = idx_count + 1;
global idx_3D_w
idx_3D_w = idx_count;
units_corrector_type(idx_3D_w) = 0;
varname_sam(idx_3D_w,1:1) = 'W';

% rt (kg/kg).
idx_count = idx_count + 1;
global idx_3D_rt
idx_3D_rt = idx_count;
units_corrector_type(idx_3D_rt) = 0;
varname_sam(idx_3D_rt,1:2) = 'RT';

% theta_l (K).
idx_count = idx_count + 1;
global idx_3D_thl
idx_3D_thl = idx_count;
units_corrector_type(idx_3D_thl) = 0;
varname_sam(idx_3D_thl,1:3) = 'THL';

% chi (kg/kg).
idx_count = idx_count + 1;
global idx_3D_chi
idx_3D_chi = idx_count;
units_corrector_type(idx_3D_chi) = 0;
varname_sam(idx_3D_chi,1:3) = 'CHI';

% eta (kg/kg).
idx_count = idx_count + 1;
global idx_3D_eta
idx_3D_eta = idx_count;
units_corrector_type(idx_3D_eta) = 0;
varname_sam(idx_3D_eta,1:3) = 'ETA';

% rr (kg/kg)
idx_count = idx_count + 1;
global idx_3D_rr
idx_3D_rr = idx_count;
units_corrector_type(idx_3D_rr) = 0;
varname_sam(idx_3D_rr,1:2) = 'RR';

% Nr (num/kg)
idx_count = idx_count + 1;
global idx_3D_Nr
idx_3D_Nr = idx_count;
units_corrector_type(idx_3D_Nr) = 0;
varname_sam(idx_3D_Nr,1:2) = 'NR';

if ( strcmp( casename, 'LBA' ) )

   % ri (kg/kg)
   idx_count = idx_count + 1;
   global idx_3D_ri
   idx_3D_ri = idx_count;
   units_corrector_type(idx_3D_ri) = 0;
   varname_sam(idx_3D_ri,1:2) = 'RI';

   % Ni (num/kg)
   idx_count = idx_count + 1;
   global idx_3D_Ni
   idx_3D_Ni = idx_count;
   units_corrector_type(idx_3D_Ni) = 0;
   varname_sam(idx_3D_Ni,1:2) = 'NI';

   % rs (kg/kg)
   idx_count = idx_count + 1;
   global idx_3D_rs
   idx_3D_rs = idx_count;
   units_corrector_type(idx_3D_rs) = 0;
   varname_sam(idx_3D_rs,1:2) = 'RS';

   % Ns (num/kg)
   idx_count = idx_count + 1;
   global idx_3D_Ns
   idx_3D_Ns = idx_count;
   units_corrector_type(idx_3D_Ns) = 0;
   varname_sam(idx_3D_Ns,1:2) = 'NS';

   % rg (kg/kg)
   idx_count = idx_count + 1;
   global idx_3D_rg
   idx_3D_rg = idx_count;
   units_corrector_type(idx_3D_rg) = 0;
   varname_sam(idx_3D_rg,1:2) = 'RG';

   % Ng (num/kg)
   idx_count = idx_count + 1;
   global idx_3D_Ng
   idx_3D_Ng = idx_count;
   units_corrector_type(idx_3D_Ng) = 0;
   varname_sam(idx_3D_Ng,1:2) = 'NG';

end % strcmp( casename, 'LBA' )

num_var_sam = idx_count;

%============= Add SAM Model 1-D Variables Here ===========================

% West-East Distance (degrees)
idx_count = idx_count + 1;
global idx_3D_x
idx_3D_x = idx_count;
% Need to add a degree-to-meter unit corrector.
units_corrector_type(idx_3D_x) = 0;
varname_sam(idx_3D_x,1:1) = 'x';

% South-North Distance (degrees)
idx_count = idx_count + 1;
global idx_3D_y
idx_3D_y = idx_count;
% Need to add a degree-to-meter unit corrector.
units_corrector_type(idx_3D_y) = 0;
varname_sam(idx_3D_y,1:1) = 'y';

% Altitude (meters)
idx_count = idx_count + 1;
global idx_3D_z
idx_3D_z = idx_count;
units_corrector_type(idx_3D_z) = 0;
varname_sam(idx_3D_z,1:1) = 'z';

% Elapsed time (minutes)
idx_count = idx_count + 1;
global idx_3D_time
units_corrector_type(idx_3D_time) = 0;
idx_3D_time = idx_count;
varname_sam(idx_3D_time,1:4) = 'time';

num_tot_var_sam = idx_count;
