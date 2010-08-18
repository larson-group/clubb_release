function[] = cloud_feedback_timeseries_plot(nc_path, verify_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gca, 'ColorOrder', [0 0 1; 0 0.6 0; 1 0 0],'LineStyleOrder',{'-','--','o'},'NextPlot','ReplaceChildren');
set(gca, 'LineStyleOrder',{'-','--'},'NextPlot','ReplaceChildren');

curr_case = 'cloud_feedback_s12_p2k';
sec_per_hour = 3600;
mm_per_m = 1000;
nz = 41;

t_start = 1;
t_end = t_start + 719;

timeseries_vars_to_plot = ['cldtot'; 'tglwp '; 'precw '; 'tsair '; 'ps    '; 'prect '; 'lh    '; 'sh    '; 'fsns  '; 'flns  '; 'fsnt  '; 'flnt  '; 'flntc '; 'fsntc '; 'fsnsc '; 'flnsc '];
sfcfilepath = [nc_path, curr_case, '_scm_UWM_CLUBB_v2.nc'];
press = ['p     '];

profile_vars_to_plot = [ 'T     '; 'ql    '; 'qv    '; 'cloud '; 'tdt_lw'; 'tdt_sw'; 'tdt_ls'; 'qdt_ls' ];

if ( exist(sfcfilepath) )
	% Plot the timeseries variables
	for i=1:size(timeseries_vars_to_plot,1);
		minVal = 0;
		maxVal = 0;

		var_to_plot = strtrim(timeseries_vars_to_plot(i, 1:size(timeseries_vars_to_plot,2)))	


		sfcfile = netcdf.open(sfcfilepath,'NC_NOWRITE');

		varid = netcdf.inqVarID(sfcfile,'time');
		file_time = netcdf.getVar(sfcfile,varid);
		%file_time = file_time ./ 5; % Timestep (output time in minutes)
		%file_time = file_time ./ 60;

		% Read in variable
		varid = netcdf.inqVarID(sfcfile, var_to_plot);
		file_var = netcdf.getVar(sfcfile, varid);

		plot (file_time, file_var);

		hold all

		%Retain values for axis scaling
		if (min(min(file_var)) < minVal)
			minVal = min(min(file_var));
		end

		if (max(max(file_var)) > maxVal)
			maxVal = max(max(file_var));
		end


		hold off

		%Determine the unit
		attname = netcdf.inqAttName(sfcfile, varid, 0);
		units = netcdf.getAtt( sfcfile, varid, attname );

		% Set labels for axes and plot title
		ylabel([var_to_plot, ' ', '[', units, ']']);
		xlabel('Time   [h]');
		title(['Cloud Feedback ', var_to_plot, ' verification ']);
		grid;

		%Scale the axis
		percentage = 95;
		
		% Normalize units to a factor of 1.
		percentage = percentage / 100;
		
		equiv_space = 0.01;
					
		if ( maxVal == minVal )
   			ymin = minVal - equiv_space;
   			ymax = maxVal + equiv_space;
		else
   			ydiff = maxVal - minVal;
   			yrange = ydiff / percentage;
   			ymedian = ( minVal + maxVal ) / 2;
   			ymin = ymedian - yrange / 2;
   			ymax = ymedian + yrange / 2;
		end

		axis([ min(min(file_time)) max(max(file_time)) ymin ymax ]);

		%PDF
		output_file_name = [ verify_path, curr_case, '_', var_to_plot, '_verify.pdf' ];
		print( '-dpdf', '-append', output_file_name );
	end

	% Read in the pressure to use as the y-axis in the profile plots
	pressure = zeros(1,nz);

	var_to_plot = strtrim(press(1))
	profilefile = netcdf.open(sfcfilepath,'NC_NOWRITE');
	varid = netcdf.inqVarID(sfcfile, 'p');
	pressure = netcdf.getVar(sfcfile, varid);

	% Plot the profiles as an average of the variable over the entire run time
	for i=1:size(profile_vars_to_plot,1);
		plot_index = 1;

		minVal = 0;
		maxVal = 0;

		var_to_plot = strtrim(profile_vars_to_plot(i, 1:size(profile_vars_to_plot,2)))

		profilefile = netcdf.open(sfcfilepath,'NC_NOWRITE');
		varid = netcdf.inqVarID(sfcfile,'time');
		file_time = netcdf.getVar(sfcfile,varid);
		%file_time = file_time ./ 5; % Timestep (output time in minutes)
		%file_time = file_time ./ 60;

		varid = netcdf.inqVarID(sfcfile, var_to_plot);
		file_var = netcdf.getVar(sfcfile, varid);

		%Average
		var_avg = zeros(nz, 1);
		for t=t_start:t_end
			var_avg = var_avg + file_var(:,t);
		end
		var_avg = var_avg ./ (t_end - t_start);

		plot(var_avg, pressure);

		hold all

		%Retain values for axis scaling
		if (min(min(var_avg)) < minVal)
			minVal = min(min(var_avg));
		end

		if (max(max(var_avg)) > maxVal)
			maxVal = max(max(var_avg));
		end

		plot_index = plot_index + 1;

		hold off

		%Determine the unit
		attname = netcdf.inqAttName(sfcfile, varid, 0);
		units = netcdf.getAtt( sfcfile, varid, attname );

		% Set labels for axes and plot title
		xlabel([var_to_plot, ' [', units, ']']); 
		ylabel('Pressure [mb]');
		title([curr_case, ' ', var_to_plot, ' verification']);
		grid;

		%Scale the axis
		percentage = 95;
		
		% Normalize units to a factor of 1.
		percentage = percentage / 100;
		
		equiv_space = 0.01;
				
		if ( maxVal == minVal )
   			xmin = minVal - equiv_space;
   			xmax = maxVal + equiv_space;
		else
   			xdiff = maxVal - minVal;
   			xrange = xdiff / percentage;
   			xmedian = ( minVal + maxVal ) / 2;
   			xmin = xmedian - xrange / 2;
   			xmax = xmedian + xrange / 2;
		end
		% Reverse the scale for pressure on the y-axis and set axes
		set(gca, 'YDir', 'reverse');
		axis([ xmin xmax min(min(pressure)) max(max(pressure)) ]);

		%PDF
		output_file_name = [ verify_path, curr_case, '_', var_to_plot, '_verify.pdf' ];
		print( '-dpdf', '-append', output_file_name );
	end
end
