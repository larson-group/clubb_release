% $Id$
function [ varname_sam, units_corrector_type, ...
           num_var_sam, num_tot_var_sam ] = output_vars_3D

idx_count = 0;

%============= Add SAM Model Field Variables Here =========================

% rr (kg/kg)
idx_count = idx_count + 1;
global idx_3D_rr
idx_3D_rr = idx_count;
units_corrector_type(idx_3D_rr) = 3;
varname_sam(idx_3D_rr,1:2) = 'QP';

% w (m/s).
idx_count = idx_count + 1;
global idx_3D_w
idx_3D_w = idx_count;
units_corrector_type(idx_3D_w) = 0;
varname_sam(idx_3D_w,1:1) = 'W';

num_var_sam = idx_count;

%============= Add SAM Model 1-D Variables Here ===========================

% West-East Distance (meters)
idx_count = idx_count + 1;
global idx_3D_x
idx_3D_x = idx_count;
units_corrector_type(idx_3D_x) = 0;
varname_sam(idx_3D_x,1:1) = 'x';

% South-North Distance (meters)
idx_count = idx_count + 1;
global idx_3D_y
idx_3D_y = idx_count;
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
