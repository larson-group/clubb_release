function[] = twp_ice_profiles_plot();

addpath '/home/nielsenb/mexnc2/mexnc/'
addpath '/home/nielsenb/snctools'
addpath '/home/nielsenb/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/nielsenb/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/nielsenb/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gca, 'ColorOrder', [0 0 1; 0 0.6 0; 1 0 0],'LineStyleOrder',{'-','--','o'},'NextPlot','ReplaceChildren');
set(gca, 'LineStyleOrder',{'-','--'},'NextPlot','ReplaceChildren');

sec_per_hour = 3600;
mm_per_m = 1000;
nz = 128;

t_start = 1;
t_end = t_start + 120;

vars_to_plot = ['p    '; 'T    '; 'theta'; 'rho  '; 'w    '; 'qv   '; 'RH   '; 'qc   '; 'qi   '; 'qr   '; 'CF   '; 'Q1   '; 'Q2   '; 'TQsw '; 'TQlw ']

for i=1:size(vars_to_plot,1);
	plot_index = 1;

	minVal = 0;
	maxVal = 0;

	var_to_plot = strtrim(vars_to_plot(i, 1:size(vars_to_plot,2)))

	for j=0:99
		profilefilepath = ['submission/', 'profiles.CLUBB_p', sprintf('%02d',j), '.nc'];

		if ( exist(profilefilepath) )
			profilefile = netcdf(profilefilepath,'nowrite');
			file_time = profilefile{'time'}(:);

			%Convert time to  julian days
			file_time = file_time ./ 24;
			file_time = file_time + 20;

			file_var = profilefile{var_to_plot}(:);
			file_height = profilefile{'h'}(:);

			%Make sure RH is between 0 and 1
			if (strcmp('RH', var_to_plot))
				if (any(file_var(:) < 0.0))
					disp(['Invalid RH in ensemble member ', sprintf('%02d',j), '!']);
				end

				if (any(file_var(:) > 1.0))
					disp(['Invalid RH in ensemble member ', sprintf('%02d',j), '!']);
				end
			end

			%Average
			var_avg = zeros(1, nz);
			for t=t_start:t_end
				var_avg = var_avg + file_var(t,:);
			end
			var_avg = var_avg ./ (t_end - t_start);

			h(plot_index) = plot(var_avg, file_height(1,:));

			axis_label =  sprintf('%02d',j);
			legend_text(plot_index,1:length(axis_label)) = axis_label;

			hold all

			%Determine the unit
			units = nc_attget( profilefilepath, var_to_plot, 'unit' );

			%Retain values for axis scaling
			if (min(min(var_avg)) < minVal)
				minVal = min(min(var_avg));
			end

			if (max(max(var_avg)) > maxVal)
				maxVal = max(max(var_avg));
			end

			plot_index = plot_index + 1;

			close(profilefile);
		end 
	end

	%Now plot the best estimate	
	profilefilepath = ['submission/', 'profiles.CLUBB_best_estimate', '.nc'];

	if ( exist(profilefilepath) )
		profilefile = netcdf(profilefilepath,'nowrite');
		file_time = profilefile{'time'}(:);

		%Convert time to  julian days
		file_time = file_time ./ 24;
		file_time = file_time + 20;

		file_var = profilefile{var_to_plot}(:);
		file_height = profilefile{'h'}(:);

		%Make sure RH is between 0 and 1
		if (strcmp('RH', var_to_plot))
			if (any(file_var(:) < 0.0))
				disp('Invalid RH in best estimate!');
			end

			if (any(file_var(:) > 1.0))
				disp('Invalid RH in best estimate!');
			end
		end

		%Average
		var_avg = zeros(1, nz);
		for t=t_start:t_end
			var_avg = var_avg + file_var(t,:);
		end
		var_avg = var_avg ./ (t_end - t_start);

		h(plot_index) = plot(var_avg, file_height(1,:));

		axis_label =  'best estimate';
		legend_text(plot_index,1:length(axis_label)) = axis_label;

		%Determine the unit
		units = nc_attget( profilefilepath, var_to_plot, 'unit' );

		%Retain values for axis scaling
		if (min(min(var_avg)) < minVal)
			minVal = min(min(var_avg));
		end

		if (max(max(var_avg)) > maxVal)
			maxVal = max(max(var_avg));
		end

		hold all

		close(profilefile);
	end 

	hold off

	legend( h, legend_text, 'Location', 'NorthEast' )
	xlabel([var_to_plot, '     ', '[', units, ']'])
	ylabel('Height    [m]')
	title(['Ensemble ', var_to_plot, ' verification ',  '(5-min timestep)'])
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

	axis([ xmin xmax min(min(file_height)) max(max(file_height)) ])

	%PDF
	output_file_name = [ 'twp_ice_', var_to_plot, '_verify.pdf' ];
	print( '-dpdf', '-append', output_file_name )
	close;
end	
