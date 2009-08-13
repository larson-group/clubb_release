function[] = cloud_feedback_timeseries_plot();

addpath '/home/senkbeir/mexnc2/mexnc/'
addpath '/home/senkbeir/snctools'
addpath '/home/senkbeir/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/senkbeir/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/senkbeir/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gca, 'ColorOrder', [0 0 1; 0 0.6 0; 1 0 0],'LineStyleOrder',{'-','--','o'},'NextPlot','ReplaceChildren');
set(gca, 'LineStyleOrder',{'-','--'},'NextPlot','ReplaceChildren');

curr_case = 's12_p2k';
sec_per_hour = 3600;
mm_per_m = 1000;
nz = 41;

t_start = 1;
t_end = t_start + 719;

vars_to_plot = ['cldtot'; 'tglwp '; 'precw '; 'tsair '; 'ps    '; 'prect '; 'lh    '; 'sh    '; 'fsnt  '; 'flnt  ']

for i=1:size(vars_to_plot,1);
	minVal = 0;
	maxVal = 0;

	var_to_plot = strtrim(vars_to_plot(i, 1:size(vars_to_plot,2)))

	sfcfilepath = ['/home/senkbeir/nc_output/', 'cloud_feedback_', curr_case, '_scm_UWM_CLUBB_v1.nc'];	

	if ( exist(sfcfilepath) )
		sfcfile = netcdf(sfcfilepath,'nowrite');

		file_time = sfcfile{'time'}(:);
		file_time = file_time ./ 5; % Timestep (output time in minutes)
		file_time = file_time ./ 60;

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
	output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_', curr_case, '_', var_to_plot, '_verify.pdf' ];
	print( '-dpdf', '-append', output_file_name )
end	
