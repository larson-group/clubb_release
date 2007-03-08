function [] = make_gabls_files();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_gabls_files()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, January-March 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters: none
% Input files: GABLS2 GrADS files (.ctl and .dat)
% Output parameters: none
% Output files: GABLS2 text files: mean_uwm.dat, 
%
% Requires: two Michael Falk utility scripts:
%           header_read.m
%           read_grads_hoc_endian.m
%
% Reference: http://people.su.se/~gsven/gabls/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program has the ability to plot profiles of several fields at the
% output times of 24, 32, and 36 hours for comparison with the GABLS web
% page.  This is commented out by default, but you can uncomment it if
% you'd like to scroll through the plots during each execution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set up input files and timestep
scm_path = ['/home/mjfalk/hoc/gabls/']
smfile   = 'gabls_zt.ctl';

t = 0:6:354;
t(1) = 1;
sizet = size(t);
sizet = max(sizet);

% Load GABLS2 data from mass grid files
[filename,nz,z,ntimesteps,numvars,list_vars] = header_read([scm_path,smfile]);
kmax = max(size(z))

for i=1:numvars
    for timestep = 1:sizet-1
        stringtoeval = [list_vars(i,:), ' = read_grads_hoc_endian([scm_path,filename],''ieee-le'',nz,t(timestep),t(timestep+1),i,numvars);'];
        eval(stringtoeval)
        str = list_vars(i,:);
        for j=1:nz
            arraydata(j,timestep) = eval([str,'(j)']);
        end
    end
    eval([strtrim(str),'_array = arraydata;']);
end

% Set constants
g0   = 9.8;
p0   = 1e5;
R    = 287.04;
Cp   = 1004.67;
Lv   = 2.5e6;
vonk = 0.4;

% Set derived variables
exner_array = (p_array./p0).^(R/Cp);
theta_array = thlm_array + (Lv/Cp).*rcm_array;
T_array = theta_array .* exner_array;
rvm_array = rtm_array - rcm_array;
RH_array = rvm_array ./ rsm_array;

% Open output file and print required data to it
fid=fopen('mean_uwm.dat','w');
for ts=1:sizet-1
    for k=1:kmax
        fprintf(fid,'%07.2f %09.3f %09.5f %09.5f %09.7f %09.5f %09.5f %09.7f\n', ...
            z(k),p_array(k,ts),T_array(k,ts),theta_array(k,ts),rvm_array(k,ts), ...
            um_array(k,ts),vm_array(k,ts),wm_array(k,ts) ...
        );
    end
    fprintf(fid,'\n');

%% Make plots if desired
%    if ((ts == 24) || (ts == 32) || (ts == 36))
%        plot (p_array(:,ts),z(:))
%        title(['pressure, time ',num2str(ts)])
%        keyboard
%        plot (T_array(:,ts),z(:))
%        title(['Temperature, time ',num2str(ts)])
%        keyboard
%        plot (theta_array(:,ts),z(:))
%        title(['theta, time ',num2str(ts)])
%        keyboard
%        plot (rvm_array(:,ts),z(:))
%        title(['rvm, time ',num2str(ts)])
%        keyboard
%        plot (sqrt(um_array(:,ts).^2 + vm_array(:,ts).^2),z(:))
%        title(['windspeed, time ',num2str(ts)])
%        keyboard
%        plot (um_array(:,ts),z(:))
%        title(['um, time ',num2str(ts)])
%        keyboard
%        plot (vm_array(:,ts),z(:))
%        title(['vm, time ',num2str(ts)])
%        keyboard
%        plot (wm_array(:,ts),z(:))
%        title(['wm, time ',num2str(ts)])
%        keyboard
%    end
end
fclose('all');


% Load GABLS2 data from momentum grid files
zw = 0:10:3990;
swfile   = 'gabls_zm.ctl';

[wfilename,wnz,wz,wntimesteps,wnumvars,wlist_vars] = header_read([scm_path,swfile]);

for i=1:wnumvars
    for timestep = 1:sizet-1
        stringtoeval = [wlist_vars(i,:), ' = read_grads_hoc_endian([scm_path,wfilename],''ieee-le'',wnz,t(timestep),t(timestep+1),i,wnumvars);'];
        eval(stringtoeval)
        str = wlist_vars(i,:);
        for j=1:nz
            arraydata(j,timestep) = eval([str,'(j)']);
        end
    end
    eval(['w_',strtrim(str),'_array = arraydata;']);
end

% Open output file and print required data to it
fid=fopen('turb_uwm.dat','w');

w_tke_array = (0.5 * w_up2_plus_vp2_array) + (0.5 * w_wp2_array);
kmaxw = max(size(zw))

for ts=1:sizet-1
    for k=1:kmaxw-1
        wp3_over_wp2_threehalves(k,ts) = (0.5 * (wp3_array(k,ts) + wp3_array(k+1,ts))) / ((w_wp2_array(k,ts)) ^ 1.5);
    end
    wp3_over_wp2_threehalves(kmax,ts) = 0.;
end

for ts=1:sizet-1
    for k=1:kmaxw-1
        fprintf(fid,'%07.2f %09.6f %09.6f %09.6f %09.7f %09.7f %09.7f %09.5f\n', ...
            zw(k),0.5*(w_wp2_array(k,ts)) + 0.5*(w_up2_plus_vp2_array(k,ts)),zeros(1,1),zeros(1,1),w_wp2_array(k,ts), ...
            w_thlp2_array(k,ts),w_wpthlp_array(k,ts),wp3_over_wp2_threehalves(k,ts) ...
        );
    end
    fprintf(fid,'\n');
%% Make plots if desired
%    if ((ts == 24) || (ts == 32) || (ts == 36))
%        plot (0.5*(w_wp2_array(:,ts)) + 0.5*(w_up2_plus_vp2_array(:,ts)),z(:))
%        title(['tke, time ',num2str(ts)])
%        keyboard
%        plot (w_up2_plus_vp2_array(:,ts),z(:))
%        title(['up2+vp2, time ',num2str(ts)])
%        keyboard
%        plot (w_wp2_array(:,ts),z(:))
%        title(['wp2, time ',num2str(ts)])
%        keyboard
%        plot (w_thlp2_array(:,ts),z(:))
%        title(['thlp2, time ',num2str(ts)])
%        keyboard
%        plot (w_wpthlp_array(:,ts),z(:))
%        title(['wpthlp, time ',num2str(ts)])
%        keyboard
%        plot (wp3_over_wp2_threehalves(:,ts),z(:))
%        title(['1.5*(wp3/wp2), time ',num2str(ts)])
%        keyboard
%    end
end


% Open output file and print required data to it
fid=fopen('flux_uwm.dat','w');
for ts=1:sizet-1
    for k=1:kmaxw-1
        fprintf(fid,'%07.2f %09.7f %09.7f %09.7f %09.7f %09.7f %09.7f\n', ...
            zw(k),w_khm_array(k,ts),zeros(1,1),w_upwp_array(k,ts), ...
            w_vpwp_array(k,ts),w_wpthlp_array(k,ts),w_wprtp_array(k,ts) ...
        );
    end
    fprintf(fid,'\n');
%% Make plots if desired
%    if ((ts == 24) || (ts == 32) || (ts == 36))
%        plot (w_khm_array(:,ts),z(:))
%        title(['khm, time ',num2str(ts)])
%        keyboard
%        plot (w_upwp_array(:,ts),z(:))
%        title(['upwp, time ',num2str(ts)])
%        keyboard
%        plot (w_vpwp_array(:,ts),z(:))
%        title(['vpwp, time ',num2str(ts)])
%        keyboard
%        plot (w_wpthlp_array(:,ts),z(:))
%        title(['wpthlp, time ',num2str(ts)])
%        keyboard
%        plot (w_wprtp_array(:,ts),z(:))
%        title(['wprtp, time ',num2str(ts)])
%        keyboard
%    end
end


% Open output file and print required data to it
fid=fopen('surf_uwm.dat','w');
for ts=1:sizet-1
    % Compute boundary layer height
    sfc_blh_term = sqrt(w_upwp_array(1,ts)^2 + w_vpwp_array(1,ts)^2);
    lowest_blh_level = kmax;
    for k=2:kmax-1
        blh_term = sqrt(w_upwp_array(k,ts)^2 + w_vpwp_array(k,ts)^2);
        if ((blh_term < (0.05 * sfc_blh_term)) && (k < lowest_blh_level))
            lowest_blh_level = k;
        end
    end

    blh(ts) = zw(lowest_blh_level) / 0.95;

    % Compute ustar and Monin-Obukhov length
    ustar(ts)       = (w_upwp_array(1,ts)^2 + w_vpwp_array(1,ts)^2) ^ 0.25;
    mo_length(ts)   = -thvm_array(1,ts) * (ustar(ts)^3) / (vonk * g0 * (w_wpthvp_array(1,ts)));
        fprintf(fid,'%07.4f %11.7f %11.5f %09.7f %09.7f %09.7f\n', ...
        t(ts)/6,mo_length(ts),blh(ts),ustar(ts), ...
        w_wpthlp_array(1,ts),w_wprtp_array(1,ts) ...
    );

    fprintf(fid,'\n');
%% Make plots if desired
%    if ((ts == 24) || (ts == 32) || (ts == 36))
%        plot (w_khm_array(:,ts),z(:))
%        title(['khm, time ',num2str(ts)])
%        keyboard
%        plot (w_upwp_array(:,ts),z(:))
%        title(['upwp, time ',num2str(ts)])
%        keyboard
%        plot (w_vpwp_array(:,ts),z(:))
%        title(['vpwp, time ',num2str(ts)])
%        keyboard
%        plot (w_wpthlp_array(:,ts),z(:))
%        title(['wpthlp, time ',num2str(ts)])
%        keyboard
%        plot (w_wprtp_array(:,ts),z(:))
%        title(['wprtp, time ',num2str(ts)])
%        keyboard
%    end
end

%% Make plots if desired
%plot (1:sizet-1,1./mo_length)
%keyboard
%plot (1:sizet-1,blh)
%keyboard
%plot (1:sizet-1,ustar)
%keyboard



%% Open budget file, but we are unsure about these formulations.

%fid=fopen('budget_uwm.dat','w');
%for ts=1:sizet-1
%    for k=1:kmax
%        fprintf(fid,'%09.1f %09.6f %09.6f %09.6f %09.6f\n', ...
%            zw(k),w_wp2_ma_array(1,1),w_wp2_bp_array(k,ts), ...
%            w_wp2_ta_array(k,ts),w_wp2_dp1_array(k,ts) + w_wp2_dp2_array(k,ts) ...
%        );
%    end
%    fprintf(fid,'\n');
%
%% Make plots if desired
%    if ((ts == 24) || (ts == 32) || (ts == 36))
%        plot (w_wp2_ma_array(:,ts),z(:))
%        title(['wp2_ma, time ',num2str(ts)])
%        keyboard
%        plot (w_wp2_bp_array(:,ts),z(:))
%        title(['wp2_bp, time ',num2str(ts)])
%        keyboard
%        plot (w_wp2_ta_array(:,ts),z(:))
%        title(['wp2_ta, time ',num2str(ts)])
%        keyboard
%        plot (w_wp2_dp1_array(:,ts) + w_wp2_dp2_array(:,ts),z(:))
%        title(['wp2_dp, time ',num2str(ts)])
%        keyboard
%    end
%end