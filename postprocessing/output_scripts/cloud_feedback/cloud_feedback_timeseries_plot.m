function[] = cloud_feedback_timeseries_plot(nc_path, verify_path);

addpath '/home/meyernr/mexnc2/mexnc/'
addpath '/home/meyernr/snctools'
addpath '/home/meyernr/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/meyernr/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/meyernr/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gca, 'ColorOrder', [0 0 1; 0 0.6 0; 1 0 0],'LineStyleOrder',{'-','--','o'},'NextPlot','ReplaceChildren');
set(gca, 'LineStyleOrder',{'-','--'},'NextPlot','ReplaceChildren');

curr_case = 'cloud_feedback_s12_p2k';
sec_per_hour = 3600;
mm_per_m = 1000;
nz = 41;
domain_top = 5950;

t_start = 1;
t_end = t_start + 719;

timeseries_vars_to_plot = ['cldtot'; 'tglwp '; 'precw '; 'tsair '; 'ps    '; 'prect '; 'lh    '; 'sh    '; 'fsns  '; 'flns  '; 'fsnt  '; 'flnt  '; 'flntc '; 'fsntc '; 'fsnsc '; 'flnsc ']
press = ['p     ']

profile_vars_to_plot = [ 'T     '; 'ql    '; 'qv    '; 'cloud '; 'tdt_lw'; 'tdt_sw'; 'tdt_ls'; 'qdt_ls' ]
% Plot the timeseries variables
for i=1:size(timeseries_vars_to_plot,1);
	minVal = 0;
	maxVal = 0;

	var_to_plot = strtrim(timeseries_vars_to_plot(i, 1:size(timeseries_vars_to_plot,2)))

	sfcfilepath = [nc_path, curr_case, '_scm_UWM_CLUBB_v2.nc'];	

	if ( exist(sfcfilepath) )
		sfcfile = netcdf(sfcfilepath,'nowrite');

		file_time = sfcfile{'time'}(:);
		%file_time = file_time ./ 5; % Timestep (output time in minutes)
		%file_time = file_time ./ 60;

		file_var = sfcfile{var_to_plot}(:);

		plot (file_time, file_var);

		hold all

		%Retain values for axis scaling
		if (min(min(file_var)) < minVal)
			minVal = min(min(file_var));
		end

		if (max(max(file_var)) > maxVal)
			maxVal = max(max(file_var));
		end
	end 

	hold off

	%Determine the unit
	units = nc_attget( sfcfilepath, var_to_plot, 'unit' );
	
	ylabel([var_to_plot, ' ', '[', units, ']'])
	xlabel('Time   [h]')
	title(['Cloud Feedback ', var_to_plot, ' verification '])
	grid

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

	axis([ min(min(file_time)) max(max(file_time)) ymin ymax ])

	%PDF
	output_file_name = [ verify_path, curr_case, '_', var_to_plot, '_verify.pdf' ];
	print( '-dpdf', '-append', output_file_name )
end

% Read in the pressure to use as the y-axis in the profile plots
profilefilepath = [nc_path, curr_case, '_scm_UWM_CLUBB_v2.nc'];
pressure = zeros(1,nz)
if( exist(profilefilepath) )
	var_to_plot = strtrim(press(1))
	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);
	pressure = profilefile{var_to_plot}(:)
	% Invert the pressure array?
end
% Plot the profiles as an average of the variable over the entire run time
for i=1:size(profile_vars_to_plot,1);
	plot_index = 1;

	minVal = 0;
	maxVal = 0;

	var_to_plot = strtrim(profile_vars_to_plot(i, 1:size(profile_vars_to_plot,2)))

	profilefilepath = [nc_path, curr_case, '_scm_UWM_CLUBB_v2.nc'];

	if ( exist(profilefilepath) )
		profilefile = netcdf(profilefilepath,'nowrite');
		file_time = profilefile{'time'}(:);
		%file_time = file_time ./ 5; % Timestep (output time in minutes)
		%file_time = file_time ./ 60;
		
		file_var = profilefile{var_to_plot}(:);

		for h=1:nz
			file_height(h) = (domain_top / nz) * (h - 1);
		end

		%Average
		var_avg = zeros(1, nz);
		for t=t_start:t_end
			var_avg = var_avg + file_var(t,:);
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
	end 

	hold off

	%Determine the unit
	units = nc_attget( profilefilepath, var_to_plot, 'unit' );

	xlabel([var_to_plot, ' [', units, ']']); 
	ylabel('Pressure [mb]')
	title([curr_case, ' ', var_to_plot, ' verification'])
	grid

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
	set(gca, 'YDir', 'reverse')
	axis([ xmin xmax min(min(pressure)) max(max(pressure)) ])

	%PDF
	output_file_name = [ verify_path, curr_case, '_', var_to_plot, '_verify.pdf' ];
	print( '-dpdf', '-append', output_file_name )
end
