% $Id$
function [ varid ] = nc_variable_id( ncid, varname )

% Paths to MEXNC and NetCDF Toolbox utilities.
addpath ( '/usr/share/mexcdf/mexnc' )
addpath ( '/usr/local/netcdf_toolbox/netcdf', '-end' )
addpath ( '/usr/local/netcdf_toolbox/netcdf/nctype', '-end' )
addpath ( '/usr/local/netcdf_toolbox/netcdf/ncutility', '-end' )

% Get the appropriate variable ID.
[varid, status] = mexnc ( 'inq_varid', ncid, varname );
