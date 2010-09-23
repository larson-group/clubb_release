function[] = astex_a209_plot(nc_path, verify_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gca, 'ColorOrder', [0 0 1; 0 0.6 0; 1 0 0],'LineStyleOrder',{'-','--','o'},'NextPlot','ReplaceChildren');
set(gca, 'LineStyleOrder',{'-','--'},'NextPlot','ReplaceChildren');

curr_case = 'astex_a209';
sec_per_hour = 3600;
mm_per_m = 1000;
nz = 42;

t_start = 1;
t_end = t_start + 39;

% Array of the variables to plot for each file
file_1_vars = ['u     '; 'v     '; 'thetal'; 'qt    '; 'rho   '];
file_2_vars = ['time        '; 'zcb         '; 'ztop        '; 'zmaxcfrac   '; 'LWP         '; 'precw       '; 'cc          '; 'shf         '; 'lhf         '; 'smf         '; 'tke         '; 'v_lowlevel  '; 'qv_lowlevel '; 'th_lowlevel '; 'prec_300    '; 'prec_150    '; 'prec_srf    '; 'Kh_150      '; 'Kh_300      '; 'Kh_500      '; 'Kh_1250     '; 'tsair       '; 'ps          '; 'fsntc       '; 'fsnt        '; 'flntc       '; 'flnt        '; 'fsnsc       '; 'fsns        '; 'flnsc       '; 'flns        '];
file_3_vars_zf = ['u     '; 'v     '; 'thetal'; 'qt    '; 'rho   '; 'qs    '; 'ql    '; 'qr    '; 'cf    '];
file_3_vars_zh = ['wthl  '; 'wqt   '; 'uw    '; 'vw    '; 'prec  '; 'TKE   '];
%file_4_vars = [];

sfcfilepath1 = [nc_path, 'larson_profiles_ini.nc'];
sfcfilepath2 = [nc_path, 'larson_scalars.nc'];
sfcfilepath3 = [nc_path, 'larson_profiles_avg.nc'];
sfcfilepath4 = [nc_path, 'larson_profiles_forc.nc'];

% Plots for file 1
if ( exist(sfcfilepath1) )
	% Open the file for reading
	file1 = netcdf.open( sfcfilepath1, 'NC_NOWRITE' );
	
	% Read in zf to use as the y axis
	zfvarid = netcdf.inqVarID( file1, 'zf' );
	height = netcdf.getVar( file1, zfvarid );

	% Plot the timeseries variables
	for i=1:size(file_1_vars,1);

		var_to_plot = strtrim(file_1_vars(i, 1:size(file_1_vars,2)))

		% Read in variable
		varid = netcdf.inqVarID(file1, var_to_plot);
		file_var = netcdf.getVar(file1, varid);

		% Plot zf vs. the variable
		plot ( file_var, height );

		%Determine the unit
		attname = netcdf.inqAttName(file1, varid, 0);
		units = netcdf.getAtt( file1, varid, attname );

		% Set labels for axes and plot title
		ylabel('Altitude [m]');
		xlabel([var_to_plot, ' ', '[', units, ']']);
		title(['Astex a209 ', var_to_plot, ' verification ']);
		grid;

		%PDF
		output_file_name = [ verify_path, curr_case, '_', var_to_plot, '_ini_verify.pdf' ];
		print( '-dpdf', '-append', output_file_name );
	end
end

% Plots for file 2
if ( exist(sfcfilepath2) )
	% Open the file for reading
	file2 = netcdf.open( sfcfilepath2, 'NC_NOWRITE' );
	
	% Read in zf to use as the y axis
	timevarid = netcdf.inqVarID( file2, 'time' );
	time = netcdf.getVar( file2, timevarid );

	% Plot the timeseries variables
	for i=1:size(file_2_vars,1);

		var_to_plot = strtrim(file_2_vars(i, 1:size(file_2_vars,2)))

		% Read in variable
		varid = netcdf.inqVarID(file2, var_to_plot);
		file_var = netcdf.getVar(file2, varid);

		% Plot zf vs. the variable
		plot ( time, file_var );

		%Determine the unit
		attname = netcdf.inqAttName(file2, varid, 0);
		units = netcdf.getAtt( file2, varid, attname );

		% Set labels for axes and plot title
		ylabel([var_to_plot, ' ', '[', units, ']']);
		xlabel('Time [s]');
		title(['Astex a209 ', var_to_plot, ' verification ']);
		grid;

		%PDF
		output_file_name = [ verify_path, curr_case, '_', var_to_plot, '_scalars_verify.pdf' ];
		print( '-dpdf', '-append', output_file_name );
	end
end


% Plots for file 3
if ( exist(sfcfilepath2) )
	% Open the file for reading
	file3 = netcdf.open( sfcfilepath3, 'NC_NOWRITE' );

	% Read in zf and zh for the y axis
	zfvarid = netcdf.inqVarID( file3, 'zf' );
	zf = netcdf.getVar( file3, zfvarid );
	zhvarid = netcdf.inqVarID( file3, 'zh' );
	zh = netcdf.getVar( file3, zfvarid );

	% Plot the zf profiles as an average of the variable over the entire run time
	for i=1:size(file_3_vars_zf,1);

		var_to_plot = strtrim(file_3_vars_zf(i, 1:size(file_3_vars_zf,2)))

		varid = netcdf.inqVarID(file3, var_to_plot);
		file_var = netcdf.getVar(file3, varid);

		%Average
		var_avg = zeros(nz, 1);
		for t=t_start:t_end
			var_avg = var_avg + file_var(:,t);
		end
		var_avg = var_avg ./ (t_end - t_start);

		plot(var_avg, zf);

		%Determine the unit
		attname = netcdf.inqAttName(file3, varid, 0);
		units = netcdf.getAtt( file3, varid, attname );

		% Set labels for axes and plot title
		xlabel([var_to_plot, ' [', units, ']']); 
		ylabel('Height [m]');
		title([curr_case, ' ', var_to_plot, ' verification']);
		grid;

		%PDF
		output_file_name = [ verify_path, curr_case, '_', var_to_plot, '_avg_verify.pdf' ];
		print( '-dpdf', '-append', output_file_name );
	end

	% Plot the zh profiles as an average of the variable over the entire run time
	for i=1:size(file_3_vars_zh,1);

		var_to_plot = strtrim(file_3_vars_zh(i, 1:size(file_3_vars_zh,2)))

		varid = netcdf.inqVarID(file3, var_to_plot);
		file_var = netcdf.getVar(file3, varid);

		%Average
		var_avg = zeros(nz, 1);
		for t=t_start:t_end
			var_avg = var_avg + file_var(:,t);
		end
		var_avg = var_avg ./ (t_end - t_start);

		plot(var_avg, zh);

		%Determine the unit
		attname = netcdf.inqAttName(file3, varid, 0);
		units = netcdf.getAtt( file3, varid, attname );

		% Set labels for axes and plot title
		xlabel([var_to_plot, ' [', units, ']']); 
		ylabel('Height [m]');
		title([curr_case, ' ', var_to_plot, ' verification']);
		grid;

		%PDF
		output_file_name = [ verify_path, curr_case, '_', var_to_plot, '_avg_verify.pdf' ];
		print( '-dpdf', '-append', output_file_name );
	end
end
