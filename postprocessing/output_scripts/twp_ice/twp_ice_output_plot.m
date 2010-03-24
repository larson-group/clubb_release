function[] = twp_ice_output_plot();

addpath '/home/nielsenb/mexnc2/mexnc/'
addpath '/home/nielsenb/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/nielsenb/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/nielsenb/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gca, 'ColorOrder', [0 0 1; 0 0.6 0; 1 0 0],'LineStyleOrder',{'-','--','o'},'NextPlot','ReplaceChildren');
set(gca, 'LineStyleOrder',{'-','--'},'NextPlot','ReplaceChildren');

sec_per_hour = 3600;
mm_per_m = 1000;
nz = 128;

plot_index = 1;
for i=0:99
	sfcfilepath = ['submission/', 'timeseries.CLUBB_p', sprintf('%02d',i), '.nc'];	
	profilefilepath = ['submission/', 'profiles.CLUBB_p', sprintf('%02d',i), '.nc'];

	if ( exist(sfcfilepath) )
		sfcfile = netcdf(sfcfilepath,'nowrite');
		profilefile = netcdf(profilefilepath,'nowrite');
		file_time = sfcfile{'time'}(:);

		%Convert time to  julian days
		file_time = file_time ./ 24;
		file_time = file_time + 20;

		file_ppt = sfcfile{'ppt'}(:);
		file_rho = profilefile{'rho'}(:);

		%Convert ppt to mm/hr
		file_ppt = file_ppt ./ file_rho(:,1) .* sec_per_hour .* mm_per_m;

		h(plot_index) = plot (file_time,file_ppt);

		axis_label =  sprintf('%02d',i);
		legend_text(plot_index,1:length(axis_label)) = axis_label;

		hold all

		plot_index = plot_index + 1;
	end 
end

%Now plot the best estimate
sfcfilepath = ['submission/', 'timeseries.CLUBB_best_estimate', '.nc'];	
profilefilepath = ['submission/', 'profiles.CLUBB_best_estimate', '.nc'];

if ( exist(sfcfilepath) )
	sfcfile = netcdf(sfcfilepath,'nowrite');
	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = sfcfile{'time'}(:);

	%Convert time to  julian days
	file_time = file_time ./ 24;
	file_time = file_time + 20;

	file_ppt = sfcfile{'ppt'}(:);
	file_rho = profilefile{'rho'}(:);

	%Convert ppt to mm/hr
	file_ppt = file_ppt ./ file_rho(:,1) .* sec_per_hour .* mm_per_m;

	h(plot_index) = plot (file_time,file_ppt);

	axis_label = 'best estimate';
	legend_text(plot_index,1:length(axis_label)) = axis_label;

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time	[Julian Day]')
ylabel('Surface Precipitation    [mm/hr]')
title('Ensemble Surface Precip. Verification (5-min timestep)')
grid

%PDF
output_file_name = [ 'twp_ice_precip_verify.pdf' ];
print( '-dpdf', '-append', output_file_name )


plot_index = 1;
t_start = 1;
t_end = t_start + 120;
for i=0:99
	profilefilepath = ['submission/', 'profiles.CLUBB_p', sprintf('%02d',i), '.nc'];

	if ( exist(profilefilepath) )
		profilefile = netcdf(profilefilepath,'nowrite');
		file_time = profilefile{'time'}(:);

		%Convert time to  julian days
		file_time = file_time ./ 24;
		file_time = file_time + 20;

		file_w = profilefile{'w'}(:);
		file_height = profilefile{'h'}(:);

		%Average
		w_avg = zeros(1, nz);
		for t=t_start:t_end
			w_avg = w_avg + file_w(t,:);
		end
		w_avg = w_avg ./ (t_end - t_start);

		h(plot_index) = plot(w_avg, file_height(1,:));

		axis_label =  sprintf('%02d',i);
		legend_text(plot_index,1:length(axis_label)) = axis_label;

		hold all

		plot_index = plot_index + 1;
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

	file_w = profilefile{'w'}(:);
	file_height = profilefile{'h'}(:);

	%Average
	w_avg = zeros(1, nz);
	for t=t_start:t_end
		w_avg = w_avg + file_w(t,:);
	end
	w_avg = w_avg ./ (t_end - t_start);

	h(plot_index) = plot(w_avg, file_height(1,:));

	axis_label =  'best estimate';
	legend_text(plot_index,1:length(axis_label)) = axis_label;

	hold all
end 

%Scale the axis
ylim([100 20000])

hold off

legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Vertical Velocity    [m/s]')
ylabel('Height    [m]')
title('Ensemble Vertical Velocity Verification [Monsoon] (5-min timestep)')
grid

%PDF
output_file_name = [ 'twp_ice_monsoon_verify.pdf' ];
print( '-dpdf', '-append', output_file_name )


plot_index = 1;
t_start = 144;
t_end = t_start + 120;
for i=0:99
	profilefilepath = ['submission/', 'profiles.CLUBB_p', sprintf('%02d',i), '.nc'];

	if ( exist(profilefilepath) )
		profilefile = netcdf(profilefilepath,'nowrite');
		file_time = profilefile{'time'}(:);

		%Convert time to  julian days
		file_time = file_time ./ 24;
		file_time = file_time + 20;

		file_w = profilefile{'w'}(:);
		file_height = profilefile{'h'}(:);

		%Average
		w_avg = zeros(1, nz);
		for t=t_start:t_end
			w_avg = w_avg + file_w(t,:);
		end
		w_avg = w_avg ./ (t_end - t_start);

		h(plot_index) = plot(w_avg, file_height(1,:));

		axis_label =  sprintf('%02d',i);
		legend_text(plot_index,1:length(axis_label)) = axis_label;

		hold all

		plot_index = plot_index + 1;
	end 
end

%Now plot the best esimate
profilefilepath = ['submission/', 'profiles.CLUBB_best_estimate', '.nc'];

if ( exist(profilefilepath) )
	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	%Convert time to  julian days
	file_time = file_time ./ 24;
	file_time = file_time + 20;

	file_w = profilefile{'w'}(:);
	file_height = profilefile{'h'}(:);

	%Average
	w_avg = zeros(1, nz);
	for t=t_start:t_end
		w_avg = w_avg + file_w(t,:);
	end
	w_avg = w_avg ./ (t_end - t_start);

	h(plot_index) = plot(w_avg, file_height(1,:));

	axis_label =  'best estimate';
	legend_text(plot_index,1:length(axis_label)) = axis_label;

	hold all
end 

%Scale the axis
ylim([100 20000])

hold off

legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Vertical Velocity    [m/s]')
ylabel('Height    [m]')
title('Ensemble Vertical Velocity Verification [Suppressed Monsoon] (5-min timestep)')
grid

%PDF
output_file_name = [ 'twp_ice_suppressed_monsoon_verify.pdf' ];
print( '-dpdf', '-append', output_file_name )

plot_index = 1;
t_start = 1;
t_end = t_start + 120;
for i=0:99
	profilefilepath = ['submission/', 'profiles.CLUBB_p', sprintf('%02d',i), '.nc'];

	if ( exist(profilefilepath) )
		profilefile = netcdf(profilefilepath,'nowrite');
		file_time = profilefile{'time'}(:);

		%Convert time to  julian days
		file_time = file_time ./ 24;
		file_time = file_time + 20;

		file_w = profilefile{'Q1'}(:);
		file_height = profilefile{'h'}(:);

		%Average
		w_avg = zeros(1, nz);
		for t=t_start:t_end
			w_avg = w_avg + file_w(t,:);
		end
		w_avg = w_avg ./ (t_end - t_start);

		h(plot_index) = plot(w_avg, file_height(1,:));

		axis_label =  sprintf('%02d',i);
		legend_text(plot_index,1:length(axis_label)) = axis_label;

		hold all

		plot_index = plot_index + 1;
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

	file_w = profilefile{'Q1'}(:);
	file_height = profilefile{'h'}(:);

	%Average
	w_avg = zeros(1, nz);
	for t=t_start:t_end
		w_avg = w_avg + file_w(t,:);
	end
	w_avg = w_avg ./ (t_end - t_start);

	h(plot_index) = plot(w_avg, file_height(1,:));

	axis_label =  'best estimate';
	legend_text(plot_index,1:length(axis_label)) = axis_label;

	hold all
end 

%Scale the axis
ylim([100 15000])

hold off

legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Q1    [k/day]')
ylabel('Height    [m]')
title('Q1 Verification [Monsoon] (5-min timestep)')
grid

%PDF
output_file_name = [ 'twp_ice_q1_active_verify.pdf' ];
print( '-dpdf', '-append', output_file_name )

plot_index = 1;
t_start = 144;
t_end = t_start + 120;
for i=0:99
	profilefilepath = ['submission/', 'profiles.CLUBB_p', sprintf('%02d',i), '.nc'];

	if ( exist(profilefilepath) )
		profilefile = netcdf(profilefilepath,'nowrite');
		file_time = profilefile{'time'}(:);

		%Convert time to  julian days
		file_time = file_time ./ 24;
		file_time = file_time + 20;

		file_w = profilefile{'Q2'}(:);
		file_height = profilefile{'h'}(:);

		%Average
		w_avg = zeros(1, nz);
		for t=t_start:t_end
			w_avg = w_avg + file_w(t,:);
		end
		w_avg = w_avg ./ (t_end - t_start);

		h(plot_index) = plot(w_avg, file_height(1,:));

		axis_label =  sprintf('%02d',i);
		legend_text(plot_index,1:length(axis_label)) = axis_label;

		hold all

		plot_index = plot_index + 1;
	end 
end

%Now plot the best esimate
profilefilepath = ['submission/', 'profiles.CLUBB_best_estimate', '.nc'];

if ( exist(profilefilepath) )
	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	%Convert time to  julian days
	file_time = file_time ./ 24;
	file_time = file_time + 20;

	file_w = profilefile{'Q2'}(:);
	file_height = profilefile{'h'}(:);

	%Average
	w_avg = zeros(1, nz);
	for t=t_start:t_end
		w_avg = w_avg + file_w(t,:);
	end
	w_avg = w_avg ./ (t_end - t_start);

	h(plot_index) = plot(w_avg, file_height(1,:));

	axis_label =  'best estimate';
	legend_text(plot_index,1:length(axis_label)) = axis_label;

	hold all
end 

%Scale the axis
ylim([100 15000]) 

hold off

legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Q2    [K/day]')
ylabel('Height    [m]')
title('Q2 Verification [Suppressed Monsoon] (5-min timestep)')
grid

%PDF
output_file_name = [ 'twp_ice_q2_suppressed_verify.pdf' ];
print( '-dpdf', '-append', output_file_name )

%Now do profile verification
twp_ice_profiles_plot
twp_ice_timeseries_plot
