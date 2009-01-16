function plot_qt = check_balance(les_path,smfile,swfile,scm_path,ztfile,zmfile,simname,t1,t2,ztop)

%Program check_balance
%Checks that the _bt budget term actually is equal to the sums of all of
%the other budget terms.

%_bt is plotted as a dashed blue line, and the sum of the other terms is
%plotted as a solid red line.  If the stats files are consistent, the two
%lines should coincide, and you should see a line changing from blue to red
%in a dashed pattern.

%Checks seven fields for LES (wp2, wp3, wpthlp, wpqtp, thlp2, qtp2,
%qtpthlp) and four for SCM (wp2, wp3, wpthlp, wpqtp).
%By Michael J. Falk (mjfalk@uwm.edu), 1 August 2005
%
%INPUT PARAMETERS: (1) FULL pathname to the directory of the LES statistics
%                  file
%                  (2) Filename for the mass grid (stats_sm.ctl) control
%                  file
%                  (3) Filename for the w grid (stats_sw.ctl) control
%                  file
%                  (4) FULL pathname to the SCM output files
%                  (5) Filename for the mass grid (zt.ctl) control
%                  file
%                  (6) Filename for the w grid (zm.ctl) control
%                  file
%                  (7) Simulation name (to identify plots)
%                  (8) Time (in # of timesteps from the beginning) to start
%                  averaging
%                  (9) Time (in # of timesteps from the beginning) to stop
%                  averaging
%                  (10) The height (in meters) of the top of your plots.
%(this input is identical to plot_les_scm input)
%
% SAMPLE CALLS: check_balance('/home/mjfalk/KAUS/COAMPS/COAMPS_fire/','stats_sm.ctl','stats_sw.ctl','/home/mjfalk/KAUS/HOC/fire/','zt.ctl','zm.ctl','FIRE',61,181,1000)
%               check_balance('/home/mjfalk/KAUS/COAMPS/COAMPS_bomex/','stats_sm.ctl','stats_sw.ctl','/home/mjfalk/KAUS/HOC/bomex/','zt.ctl','zm.ctl','BOMEX',181,360,2500)


[les_mgrid_filename,les_mgrid_nz,les_mgrid_z,les_mgrid_timesteps,les_mgrid_numvars,les_mgrid_list_vars] = header_read([les_path,smfile]);
[les_wgrid_filename,les_wgrid_nz,les_wgrid_z,les_wgrid_timesteps,les_wgrid_numvars,les_wgrid_list_vars] = header_read([les_path,swfile]);

[scm_mgrid_filename,scm_mgrid_nz,scm_mgrid_z,scm_mgrid_timesteps,scm_mgrid_numvars,scm_mgrid_list_vars] = header_read([scm_path,ztfile]);
[scm_wgrid_filename,scm_wgrid_nz,scm_wgrid_z,scm_wgrid_timesteps,scm_wgrid_numvars,scm_wgrid_list_vars] = header_read([scm_path,zmfile]);

 disp ( ['t1 = ', num2str(t1)] )
 disp ( ['t2 = ', num2str(t2)] )

les_mgrid_filename = [les_path,les_mgrid_filename];
les_wgrid_filename = [les_path,les_wgrid_filename];
scm_mgrid_filename = [scm_path,scm_mgrid_filename];
scm_wgrid_filename = [scm_path,scm_wgrid_filename];
simname = ['COAMPS/HOC ',simname];

for i=1:les_wgrid_numvars
    varname=les_wgrid_list_vars(i,:);
    varname=deblank(varname);
    stringtoeval = ['les_wgrid_',varname,' = read_grads_hoc_endian(les_wgrid_filename,''ieee-be'',les_wgrid_nz,t1,t2,i,les_wgrid_numvars);'];
    eval(stringtoeval);
end

for i=1:les_mgrid_numvars
    varname=les_mgrid_list_vars(i,:);
    varname=deblank(varname);
    stringtoeval = ['les_mgrid_',varname,' = read_grads_hoc_endian(les_mgrid_filename,''ieee-be'',les_mgrid_nz,t1,t2,i,les_mgrid_numvars);'];
    eval(stringtoeval);
end

for i=1:scm_wgrid_numvars
    varname=scm_wgrid_list_vars(i,:);
    varname=deblank(varname);
    stringtoeval = ['scm_wgrid_',varname,' = read_grads_hoc_endian(scm_wgrid_filename,''ieee-le'',scm_wgrid_nz,t1,t2,i,scm_wgrid_numvars);'];
    eval(stringtoeval);
end

for i=1:scm_mgrid_numvars
    varname=scm_mgrid_list_vars(i,:);
    varname=deblank(varname);
    stringtoeval = ['scm_mgrid_',varname,' = read_grads_hoc_endian(scm_mgrid_filename,''ieee-le'',scm_mgrid_nz,t1,t2,i,scm_mgrid_numvars);'];
    eval(stringtoeval);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0 = 100000.;
Rd = 287;
cp = 1005;
L  = 2.5*10^6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COAMPS BUDGETS COMPARISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wp2
wp2sum = les_wgrid_wp2_dp + les_wgrid_wp2_tp + les_wgrid_wp2_pr + les_wgrid_wp2_bp;
wp2bt  = les_wgrid_wp2_bt;
%wp3
wp3sum = les_wgrid_wp3_dp + les_wgrid_wp3_tp + les_wgrid_wp3_pr + les_wgrid_wp3_bp + les_wgrid_wp3_rs;
wp3bt  = les_wgrid_wp3_bt;
%wpthlp
newpm = les_mgrid_pm;
newpm(les_wgrid_nz) = les_mgrid_pm(les_wgrid_nz-1) - (les_mgrid_pm(les_wgrid_nz-2) - les_mgrid_pm(les_wgrid_nz-1));
for i=1:les_wgrid_nz-1
    pavg = 0.5*(newpm(i)+newpm(i+1));
    wpthlpsum(i) = (les_wgrid_wpthp_dp(i) - (les_wgrid_wpqcp_dp(i) * (p0/(pavg))^(Rd/cp) * (L/cp))) + ...
    (les_wgrid_wpthp_tp(i) - (les_wgrid_wpqcp_tp(i) * (p0/(pavg))^(Rd/cp) * (L/cp))) + ...
    (les_wgrid_wpthp_mc(i) - (les_wgrid_wpqcp_mc(i) * (p0/(pavg))^(Rd/cp) * (L/cp))) + ...
    (les_wgrid_wpthp_rad(i))                                               + ...
    (les_wgrid_wpthp_pr(i) - (les_wgrid_wpqcp_pr(i) * (p0/(pavg))^(Rd/cp) * (L/cp))) + ...
    (les_wgrid_wpthp_bp(i) - (les_wgrid_wpqcp_bp(i) * (p0/(pavg))^(Rd/cp) * (L/cp))) + ...
    (les_wgrid_wpthp_rs(i) - (les_wgrid_wpqcp_rs(i) * (p0/(pavg))^(Rd/cp) * (L/cp)));

    wpthlpbt(i)  = les_wgrid_wpthp_bt(i) - (les_wgrid_wpqcp_bt(i) * (p0/(pavg))^(Rd/cp) * (L/cp));
end
wpthlpsum(les_wgrid_nz) = wpthlpsum(les_wgrid_nz-1);
wpthlpbt(les_wgrid_nz) = wpthlpbt(les_wgrid_nz-1);

%wpqtp
wpqtpsum = les_wgrid_wpqcp_dp + les_wgrid_wpqvp_dp + ...
    les_wgrid_wpqcp_tp + les_wgrid_wpqvp_tp + ...
    les_wgrid_wpqcp_mc + les_wgrid_wpqvp_mc + ...
    les_wgrid_wpqcp_pr + les_wgrid_wpqvp_pr + ...
    les_wgrid_wpqcp_bp + les_wgrid_wpqvp_bp + ...
    les_wgrid_wpqcp_rs + les_wgrid_wpqvp_rs;
wpqtpbt  = les_wgrid_wpqcp_bt + les_wgrid_wpqvp_bt;

%thlp2
thlp2sum = les_mgrid_thp2_dp - 2.*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_thpqcp_dp + ... 
    ((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_qcp2_dp + ...
    les_mgrid_thp2_tp - 2.*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_thpqcp_tp + ... 
    ((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_qcp2_tp + ...
    les_mgrid_thp2_mc - 2.*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_thpqcp_mc + ... 
    ((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_qcp2_mc + ...
    les_mgrid_thp2_rad - 2.*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_thpqcp_rad + ... 
    les_mgrid_thp2_rs - 2.*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_thpqcp_rs + ... 
    ((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_qcp2_rs;
thlp2bt  = les_mgrid_thp2_bt - 2.*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_thpqcp_bt + ... 
    ((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*((p0./les_mgrid_pm).^(Rd/cp).*(L/cp)).*les_mgrid_qcp2_bt;

%qtp2
qtp2sum = les_mgrid_qcp2_dp + les_mgrid_qvp2_dp + 2*les_mgrid_qvpqcp_dp + ...
    les_mgrid_qcp2_tp + les_mgrid_qvp2_tp + 2*les_mgrid_qvpqcp_tp + ...
    les_mgrid_qcp2_mc + les_mgrid_qvp2_mc + 2*les_mgrid_qvpqcp_mc + ...
    les_mgrid_qcp2_rs + les_mgrid_qvp2_rs + 2*les_mgrid_qvpqcp_rs;
qtp2bt  = les_mgrid_qcp2_bt + les_mgrid_qvp2_bt + 2*les_mgrid_qvpqcp_bt;

%qtpthlp
qtpthlpsum = les_mgrid_thpqcp_dp + les_mgrid_thpqvp_dp - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qcp2_dp ...
                                         - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qvpqcp_dp ...
           + les_mgrid_thpqcp_tp + les_mgrid_thpqvp_tp - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qcp2_tp ...
                                         - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qvpqcp_tp ...
           + les_mgrid_thpqcp_mc + les_mgrid_thpqvp_mc - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qcp2_mc ...
                                         - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qvpqcp_mc ...
           + les_mgrid_thpqcp_rad + les_mgrid_thpqvp_rad ...
           + les_mgrid_thpqcp_rs + les_mgrid_thpqvp_rs - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qcp2_rs ...
                                         - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qvpqcp_rs;

qtpthlpbt  = les_mgrid_thpqcp_bt + les_mgrid_thpqvp_bt - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qcp2_bt ...
                                         - (p0./les_mgrid_pm).^(Rd/cp).*(L/cp).*les_mgrid_qvpqcp_bt;

clf
overall = subplot(1,1,1);
plot(overall,wp2sum,les_wgrid_z,'-r',wp2bt,les_wgrid_z,'--b');
title('wp2 LES: Sum of terms vs. bt');
xlabel('m^2 s^-^3');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard

clf
overall = subplot(1,1,1);
plot(overall,wp3sum,les_wgrid_z,'-r',wp3bt,les_wgrid_z,'--b');
title('wp3 LES: Sum of terms vs. bt');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard

clf
overall = subplot(1,1,1);
plot(overall,wpthlpsum,les_wgrid_z,'-r',wpthlpbt,les_wgrid_z,'--b');
title('wpthlp LES: Sum of terms vs. bt');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard

clf
overall = subplot(1,1,1);
plot(overall,wpqtpsum,les_wgrid_z,'-r',wpqtpbt,les_wgrid_z,'--b');
title('wpqtp LES: Sum of terms vs. bt');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard

clf
overall = subplot(1,1,1);
plot(overall,thlp2sum,les_mgrid_z,'-r',thlp2bt,les_mgrid_z,'--b');
title('thlp2 LES: Sum of terms vs. bt');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard

clf
overall = subplot(1,1,1);
plot(overall,qtp2sum,les_mgrid_z,'-r',qtp2bt,les_mgrid_z,'--b');
title('qtp2 LES: Sum of terms vs. bt');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard

clf
overall = subplot(1,1,1);
plot(overall,qtpthlpsum,les_mgrid_z,'-r',qtpthlpbt,les_mgrid_z,'--b');
title('qtpthlp LES: Sum of terms vs. bt');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOC BUDGETS COMPARISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wp2
wp2sum = scm_wgrid_wp2_ma + scm_wgrid_wp2_ta + scm_wgrid_wp2_ac + scm_wgrid_wp2_bp + ...
    scm_wgrid_wp2_pr1 + scm_wgrid_wp2_pr2 + scm_wgrid_wp2_pr3 + scm_wgrid_wp2_dp1 + ...
    scm_wgrid_wp2_dp2 + scm_wgrid_wp2_cl;
wp2bt  = scm_wgrid_wp2_bt;

%wp3
wp3sum = scm_mgrid_wp3_ma + scm_mgrid_wp3_ta + scm_mgrid_wp3_tp + scm_mgrid_wp3_ac + scm_mgrid_wp3_bp + ...
    scm_mgrid_wp3_pr1 + scm_mgrid_wp3_pr2 + scm_mgrid_wp3_dp1 + scm_mgrid_wp3_cl;
wp3bt  = scm_mgrid_wp3_bt;

%wpthlp
wpthlpsum = scm_wgrid_wpthlp_ma + scm_wgrid_wpthlp_ta + scm_wgrid_wpthlp_tp + scm_wgrid_wpthlp_ac ... 
    + scm_wgrid_wpthlp_bp + scm_wgrid_wpthlp_pr1 + scm_wgrid_wpthlp_pr2 + scm_wgrid_wpthlp_pr3 ...
    + scm_wgrid_wpthlp_dp1;
wpthlpbt = scm_wgrid_wpthlp_bt;

%wpqtp
wpqtpsum = scm_wgrid_wpqtp_ma + scm_wgrid_wpqtp_ta + scm_wgrid_wpqtp_tp + scm_wgrid_wpqtp_ac ...
    + scm_wgrid_wpqtp_bp + scm_wgrid_wpqtp_pr1 + scm_wgrid_wpqtp_pr2 + scm_wgrid_wpqtp_pr3 ...
    + scm_wgrid_wpqtp_dp1;
wpqtpbt = scm_wgrid_wpqtp_bt;

clf
overall = subplot(1,1,1);
plot(overall,wp2sum,scm_wgrid_z,'-r',wp2bt,scm_wgrid_z,'--b');
title('wp2 SCM: Sum of terms vs. bt');
xlabel('m^2 s^-^3');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard

clf
overall = subplot(1,1,1);
plot(overall,wp3sum,scm_mgrid_z,'-r',wp3bt,scm_mgrid_z,'--b');
title('wp3 SCM: Sum of terms vs. bt');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard

clf
overall = subplot(1,1,1);
plot(overall,wpthlpsum,scm_wgrid_z,'-r',wpthlpbt,scm_wgrid_z,'--b');
title('wpthlp SCM: Sum of terms vs. bt');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard

clf
overall = subplot(1,1,1);
plot(overall,wpqtpsum,scm_wgrid_z,'-r',wpqtpbt,scm_wgrid_z,'--b');
title('wpqtp SCM: Sum of terms vs. bt');
ylabel('Height (m)');
legend('Sum of terms','bt term');
keyboard