function[] = twp_ice_timeseries_plot();

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

vars_to_plot = ['Fs    '; 'Fq    '; 'dTFsw0'; 'dTFlw0'; 'uTFsw0'; 'uTFlw0'; 'ppt   '; 'PW    '; 'LWP   '; 'IWP   ']

for i=1:size(vars_to_plot,1);
	plot_index = 1;

	minVal = 0;
	maxVal = 0;

	var_to_plot = strtrim(vars_to_plot(i, 1:size(vars_to_plot,2)))

	for j=0:99
		sfcfilepath = ['submission/', 'timeseries.CLUBB_p', sprintf('%02d',j), '.nc'];	

		if ( exist(sfcfilepath) )
			sfcfile = netcdf(sfcfilepath,'nowrite');

			file_time = sfcfile{'time'}(:);

			%Convert time to  julian days
			file_time = file_time ./ 24;
			file_time = file_time + 20;

			file_var = sfcfile{var_to_plot}(:);

			h(plot_index) = plot (file_time,file_var);

			axis_label =  sprintf('%02d',j);
			legend_text(plot_index,1:length(axis_label)) = axis_label;

			hold all

			%Determine the unit
			units = nc_attget( sfcfilepath, var_to_plot, 'unit' );

			%Retain values for axis scaling
			if (min(min(file_var)) < minVal)
				minVal = min(min(file_var));
			end

			if (max(max(file_var)) > maxVal)
				maxVal = max(max(file_var));
			end

			plot_index = plot_index + 1;

			close(sfcfile);
		end 
	end

	%Now plot the best estimate	
	sfcfilepath = ['submission/', 'timeseries.CLUBB_best_estimate', '.nc'];	

	if ( exist(sfcfilepath) )
		sfcfile = netcdf(sfcfilepath,'nowrite');
		file_time = sfcfile{'time'}(:);

		%Convert time to  julian days
		file_time = file_time ./ 24;
		file_time = file_time + 20;

		file_var = sfcfile{var_to_plot}(:);

		h(plot_index) = plot (file_time,file_var);

		axis_label = 'best estimate';
		legend_text(plot_index,1:length(axis_label)) = axis_label;

		hold all

		%Determine the unit
		units = nc_attget( sfcfilepath, var_to_plot, 'unit' );

		%Retain values for axis scaling
		if (min(min(file_var)) < minVal)
			minVal = min(min(file_var));
		end

		if (max(max(file_var)) > maxVal)
			maxVal = max(max(file_var));
		end

		close(sfcfile);
	end 


	hold off

	legend( h, legend_text, 'Location', 'NorthEast' )
	ylabel([var_to_plot, ' ', '[', units, ']'])
	xlabel('Time    [Julian day]')
	title(['Ensemble ', var_to_plot, ' verification ',  '(5-min timestep)'])
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
	output_file_name = [ 'twp_ice_', var_to_plot, '_verify.pdf' ];
	print( '-dpdf', '-append', output_file_name )
	close;
end	
