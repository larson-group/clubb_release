%-------------------------------------------------------------------------------
function[] = plot_spaghetti( title_str, str_xlabel, fdir, curr_case, mean_var, sub_var, time, xlimits, output_name );

% Description:
%   Creates a "spaghetti" plot of the ensemble created by the latin hypercube
%   sampling code.
% Assumptions:
%   The file names and dimensions of files are assumed to as they currently are 
%   the CLUBB code (i.e. COARDS conventions for the file with the mean variables
%   and the subcolumns stored in the 'y' of _nl_lh_sample_points_2D.
% References:
%   None
%-------------------------------------------------------------------------------

set(gca, 'LineStyleOrder',{'-','--'},'NextPlot','ReplaceChildren');
subplot('Position',[0.2 0.2 0.5 0.3])

% Open the subcolumn file
subfile = netcdf.open( [ fdir, '/', curr_case, '_nl_lh_sample_points_2D.nc'], 'NC_NOWRITE' );

% Open the zt file
meanfile = netcdf.open( [ fdir, '/', curr_case, '_zt.nc'], 'NC_NOWRITE' );

% Get the altitudes.  This is assumed to be the same between the mean and the subcolumns
z = netcdf.getVar( meanfile, netcdf.inqVarID( meanfile, 'altitude' ) );
z(1) = []; % Delete the ghost point

% Get all the subcolumns of the variable
id_sub = netcdf.inqVarId( subfile, sub_var );
var_all_subcolumns = netcdf.getVar( subfile, id_sub );

% Pick the time we want
var_subcolumn = squeeze( var_all_subcolumns(:,:,:,time) );
var_subcolumn(:,1) = [];

%% Plot all subcolumns in blue %%
% Determine the number of subcolumns present
%size_ens = size( squeeze( var_subcolumn(:,:) ) );
% Plot individual subcolumns
%for i=1:size_ens
%	plot( var_subcolumn(:,:), z );
%	hold on;
%end

%% Plot all subcolumns in different colors %%
plot( var_subcolumn(:,:), z );
hold on;

% Get the mean of the variable
id_mean = netcdf.inqVarId( meanfile, mean_var );
var_m = squeeze( netcdf.getVar( meanfile, id_mean ) );
var_m(1,:) = [];
% Plot mean
mean_plot = plot( var_m(:,time), z );
set(mean_plot,'Color','red','LineWidth',2)

title( title_str );
xlabel( str_xlabel );
ylabel( 'Altitude [m]' );

xlim( xlimits );

% Close the netCDF files
netcdf.close( subfile );
netcdf.close( meanfile );

% Write an encapsulated postscript file
output_file_name = [ 'eps/', output_name,'.eps' ];
print( '-depsc2', output_file_name );
close;

end
