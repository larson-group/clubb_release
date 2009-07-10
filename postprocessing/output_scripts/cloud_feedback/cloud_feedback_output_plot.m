function[] = cloud_feedback_output_plot();

addpath '/home/senkbeir/mexnc2/mexnc/'
addpath '/home/senkbeir/netcdf_toolbox/netcdf_toolbox/netcdf/' -end
addpath '/home/senkbeir/netcdf_toolbox/netcdf_toolbox/netcdf/nctype/' -end
addpath '/home/senkbeir/netcdf_toolbox/netcdf_toolbox/netcdf/ncutility/' -end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gca, 'ColorOrder', [0 0 1; 0 0.6 0; 1 0 0],'LineStyleOrder',{'-','--','o'},'NextPlot','ReplaceChildren');
set(gca, 'LineStyleOrder',{'-','--'},'NextPlot','ReplaceChildren');

sec_per_hour = 3600;
mm_per_m = 1000;
nz = 128;

profilefilepath = ['/home/senkbeir/nc_output/', 'cloud_feedback_s6_scm_UWM_CLUBB_v1.nc'];

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_cldtot = profilefile{'cldtot'}(:);

	plot (file_time,file_cldtot);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('Total Cloud Cover')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_cldtot.pdf' ];
print( '-dpdf', '-append', output_file_name )

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_tglwp = profilefile{'tglwp'}(:);

	plot (file_time,file_tglwp);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('Vertically Integrated Liquid Water [kg/m^2]')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_tglwp.pdf' ];
print( '-dpdf', '-append', output_file_name )

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_precw = profilefile{'precw'}(:);

	plot (file_time,file_precw);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('Precipitable Water [kg/m^2]')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_precw.pdf' ];
print( '-dpdf', '-append', output_file_name )

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_tsair = profilefile{'tsair'}(:);

	plot (file_time,file_tsair);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('Surface Air Temperature [K]')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_tsair.pdf' ];
print( '-dpdf', '-append', output_file_name )

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_ps = profilefile{'ps'}(:);

	plot (file_time,file_ps);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('Surface Pressure [mb]')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_ps.pdf' ];
print( '-dpdf', '-append', output_file_name )

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_prect = profilefile{'prect'}(:);

	plot (file_time,file_prect);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('Total Precipitation [mm/day]')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_prect.pdf' ];
print( '-dpdf', '-append', output_file_name )

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_lh = profilefile{'lh'}(:);

	plot (file_time,file_lh);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('Surface latent heat flux [W/m^2]')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_lh.pdf' ];
print( '-dpdf', '-append', output_file_name )

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_sh = profilefile{'sh'}(:);

	plot (file_time,file_sh);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('Surface sensible heat flux [W/m^2]')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_sh.pdf' ];
print( '-dpdf', '-append', output_file_name )

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_fsnt = profilefile{'fsnt'}(:);

	plot (file_time,file_fsnt);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('TOA SW net downward total-sky radiation [W/m^2]')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_fsnt.pdf' ];
print( '-dpdf', '-append', output_file_name )

if ( exist(profilefilepath) )

	profilefile = netcdf(profilefilepath,'nowrite');
	file_time = profilefile{'time'}(:);

	file_flnt = profilefile{'flnt'}(:);

	plot (file_time,file_flnt);

	hold all
end 

hold off

%legend( h, legend_text, 'Location', 'NorthEast' )
xlabel('Time [h]')
ylabel('TOA LW total-sky upward radiation [W/m^2]')
title('Cloud Feedback')
grid

%PDF
output_file_name = [ '/home/matlabuser/cloud_feedback/cloud_feedback_verify_flnt.pdf' ];
print( '-dpdf', '-append', output_file_name )

end
