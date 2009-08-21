function[] = cloud_feedback_profiles_plot();

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
domain_top = 5950;

t_start = 1;
t_end = t_start + 719;

vars_to_plot = [ 'T     '; 'ql    '; 'qv    '; 'cloud '; 'tdt_lw'; 'tdt_sw'; 'tdt_ls'; 'qdt_ls' ]

for i=1:size(vars_to_plot,1);
	plot_index = 1;

	minVal = 0;
	maxVal = 0;

	var_to_plot = strtrim(vars_to_plot(i, 1:size(vars_to_plot,2)))

	profilefilepath = ['/home/senkbeir/nc_output/', 'cloud_feedback_', curr_case, '_scm_UWM_CLUBB_v2.nc'];

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

		plot(var_avg, file_height(1,:));

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
	ylabel('Height    [m]')
	title(['Cloud Feedback ', curr_case, ' ', var_to_plot, ' verification'])
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
	output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_', curr_case, '_', var_to_plot, '_verify.pdf' ];
	print( '-dpdf', '-append', output_file_name )
end	
