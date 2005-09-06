function plot_qt = budget_plotter(top_path,smfile,swfile,simname1,type1,bot_path,ztfile,zmfile,simname2, ...
                                type2,t1_les,t2_les,t1_scm,t2_scm,zbottom,ztop)

%Program budget_plotter, Version 1.17
%Compares the turbulent budgets for two simulations, LES and/or SCM.
%By Michael J. Falk (mjfalk@uwm.edu), 6 September 2005
%Thanks for help, testing, and advice to Dave Schanen and Vince Larson
%
%INPUT PARAMETERS: (1) top_path: (FULL pathname to the directory of the
%                      first sim's files)
%                  (2) smfile: Filename for the mass grid (stats_sm.ctl or zt.ctl)
%                      control file
%                  (3) swfile: Filename for the w grid (stats_sw.ctl or zm.ctl)
%                      control file
%                  (4) simname1: Simulation Name for first sim (for plot legends)
%                      (enclosed in single quotes)
%                  (5) type1: 1 if first simulation is COAMPS LES; 2 if it is
%                      HOC SCM
%                  (6) bot_path: (FULL pathname to the directory of the
%                      second sim's files)
%                  (7) ztfile: Filename for the mass grid (stats_sm.ctl or zt.ctl)
%                      control file
%                  (8) zmfile: Filename for the w grid (stats_sw.ctl or zm.ctl)
%                      control file
%                  (9) simname2: Simulation Name for second sim (for plot legends)
%                      (enclosed in single quotes)
%                  (10) type2: 1 if second simulation is COAMPS LES; 2 if it is
%                       HOC SCM
%                  (11) t1_les: Time (in # of timesteps from the beginning) to start
%                       averaging for the first sim
%                  (12) t2_les: Time (in # of timesteps from the beginning) to stop
%                       averaging for the first sim
%                  (13) t1_scm: Time (in # of timesteps from the beginning) to start
%                       averaging for the second sim
%                  (14) t2_scm: Time (in # of timesteps from the beginning) to stop
%                       averaging for the second sim
%                  (15) zbottom: The height (in meters) of the bottom of your plots.
%                  (16) ztop: The height (in meters) of the top of your plots.
%
% REQUIRED PROGRAMS:
%  header_read.m
%  make_multi_plot.m
%  make_plot.m
%  read_grads_hoc_endian.m
%
% VERY HELPFUL DIFFERENT PROGRAM:
%  check_balance.m
%
% SAMPLE CALLS:
% Convention dictates that if you are comparing an LES sim and a SCM sim,
% you should put the LES on the top and the SCM on the bottom.
% budget_plotter('/home/mjfalk/KAUS/COAMPS/COAMPS_bomex/','stats_sm.ctl','stats_sw.ctl','COAMPS BOMEX',1, ...
%                 '/home/mjfalk/KAUS/HOC/bomex/','zt.ctl','zm.ctl','HOC BOMEX',2, ...
%                 180,360,180,360,0,2500)
% budget_plotter('/home/mjfalk/KAUS/COAMPS/COAMPS_fire/','stats_sm.ctl','stats_sw.ctl','COAMPS FIRE',1, ...
%                '/home/mjfalk/KAUS/COAMPS/COAMPS_atex/','stats_sm.ctl','stats_sw.ctl','COAMPS ATEX',1, ...
%                 60,180,60,120,250,1250)
%
% OUTPUT: This program makes .eps files for each of the seven parameters
% (wp2, wp3, wpthlp, wpqtp, thlp2, qtp2, qtpthlp).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read header information to learn about data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [top_mgrid_filename,top_mgrid_nz,top_mgrid_z,top_mgrid_ntimesteps,top_mgrid_numvars,top_mgrid_list_vars] = header_read([top_path,smfile]);
    [top_wgrid_filename,top_wgrid_nz,top_wgrid_z,top_wgrid_ntimesteps,top_wgrid_numvars,top_wgrid_list_vars] = header_read([top_path,swfile]);

    [bot_mgrid_filename,bot_mgrid_nz,bot_mgrid_z,bot_mgrid_ntimesteps,bot_mgrid_numvars,bot_mgrid_list_vars] = header_read([bot_path,ztfile]);
    [bot_wgrid_filename,bot_wgrid_nz,bot_wgrid_z,bot_wgrid_ntimesteps,bot_wgrid_numvars,bot_wgrid_list_vars] = header_read([bot_path,zmfile]);

top_mgrid_bottom_z_point = 2;
top_mgrid_top_z_point    = top_mgrid_nz-1;
for k=1:top_mgrid_nz-1
    if ((top_mgrid_z(k) < zbottom) & (top_mgrid_z(k+1) >= zbottom))
        top_mgrid_bottom_z_point = k+1;
    end
    
    if ((top_mgrid_z(k) <= ztop) & (top_mgrid_z(k+1) > ztop))
        top_mgrid_top_z_point    = k;
    end
end

top_wgrid_bottom_z_point = 2;
top_wgrid_top_z_point    = top_wgrid_nz-1;
for k=1:top_wgrid_nz-1
    if ((top_wgrid_z(k) < zbottom) & (top_wgrid_z(k+1) >= zbottom))
        top_wgrid_bottom_z_point = k+1;
    end
    
    if ((top_wgrid_z(k) <= ztop) & (top_wgrid_z(k+1) > ztop))
        top_wgrid_top_z_point    = k;
    end
end

bot_mgrid_bottom_z_point = 2;
bot_mgrid_top_z_point    = bot_mgrid_nz-1;
for k=1:bot_mgrid_nz-1
    if ((bot_mgrid_z(k) < zbottom) & (bot_mgrid_z(k+1) >= zbottom))
        bot_mgrid_bottom_z_point = k+1;
    end
    
    if ((bot_mgrid_z(k) <= ztop) & (bot_mgrid_z(k+1) > ztop))
        bot_mgrid_top_z_point    = k;
    end
end

bot_wgrid_bottom_z_point = 2
bot_wgrid_top_z_point    = bot_wgrid_nz-1
for k=1:bot_wgrid_nz-1
    if ((bot_wgrid_z(k) < zbottom) & (bot_wgrid_z(k+1) >= zbottom))
        bot_wgrid_bottom_z_point = k+1;
    end
    
    if ((bot_wgrid_z(k) <= ztop) & (bot_wgrid_z(k+1) > ztop))
        bot_wgrid_top_z_point    = k;
    end
end

top_mgrid_filename = [top_path,top_mgrid_filename];
top_wgrid_filename = [top_path,top_wgrid_filename];
bot_mgrid_filename = [bot_path,bot_mgrid_filename];
bot_wgrid_filename = [bot_path,bot_wgrid_filename];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data fields into MATLAB from .dat files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:top_wgrid_numvars
    varname=top_wgrid_list_vars(i,:);
    varname=deblank(varname);
    if (type1 == 1)
        stringtoeval = ['top_wgrid_',varname,' = read_grads_hoc_endian(top_wgrid_filename,''ieee-be'',top_wgrid_nz,t1_les,t2_les,i,top_wgrid_numvars);'];
    else
        stringtoeval = ['top_wgrid_',varname,' = read_grads_hoc_endian(top_wgrid_filename,''ieee-le'',top_wgrid_nz,t1_scm,t2_scm,i,top_wgrid_numvars);'];
    end
    eval(stringtoeval);
end

for i=1:top_mgrid_numvars
    varname=top_mgrid_list_vars(i,:);
    varname=deblank(varname);
    if (type1 == 1)
        stringtoeval = ['top_mgrid_',varname,' = read_grads_hoc_endian(top_mgrid_filename,''ieee-be'',top_mgrid_nz,t1_les,t2_les,i,top_mgrid_numvars);'];
    else
        stringtoeval = ['top_mgrid_',varname,' = read_grads_hoc_endian(top_mgrid_filename,''ieee-le'',top_mgrid_nz,t1_scm,t2_scm,i,top_mgrid_numvars);'];
    end
    eval(stringtoeval);
end

for i=1:bot_wgrid_numvars
    varname=bot_wgrid_list_vars(i,:);
    varname=deblank(varname);
    if (type2 == 1)
        stringtoeval = ['bot_wgrid_',varname,' = read_grads_hoc_endian(bot_wgrid_filename,''ieee-be'',bot_wgrid_nz,t1_les,t2_les,i,bot_wgrid_numvars);'];
    else
        stringtoeval = ['bot_wgrid_',varname,' = read_grads_hoc_endian(bot_wgrid_filename,''ieee-le'',bot_wgrid_nz,t1_scm,t2_scm,i,bot_wgrid_numvars);'];
    end
    eval(stringtoeval);
end

for i=1:bot_mgrid_numvars
    varname=bot_mgrid_list_vars(i,:);
    varname=deblank(varname);
    if (type2 == 1)
        stringtoeval = ['bot_mgrid_',varname,' = read_grads_hoc_endian(bot_mgrid_filename,''ieee-be'',bot_mgrid_nz,t1_les,t2_les,i,bot_mgrid_numvars);'];
    else
        stringtoeval = ['bot_mgrid_',varname,' = read_grads_hoc_endian(bot_mgrid_filename,''ieee-le'',bot_mgrid_nz,t1_scm,t2_scm,i,bot_mgrid_numvars);'];
    end
    eval(stringtoeval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0 = 100000.;
Rd = 287;
cp = 1005;
L  = 2.5*10^6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATIONS FOR wp2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (type1 == 1)
    for i=2:top_wgrid_nz-1
        top_wgrid_wp2_ta_eqn(i) = -(top_wgrid_wp3(i+1) - top_wgrid_wp3(i-1))/(top_wgrid_z(i+1)-top_wgrid_z(i-1));
    end

    top_wgrid_wp2_ta_eqn(1) = -(top_wgrid_wp3(2) - top_wgrid_wp3(1))/(top_wgrid_z(2)-top_wgrid_z(1));
    top_wgrid_wp2_ta_eqn(top_wgrid_nz) = -(top_wgrid_wp3(top_wgrid_nz) - top_wgrid_wp3(top_wgrid_nz-1))/(top_wgrid_z(top_wgrid_nz)-top_wgrid_z(top_wgrid_nz-1));
else
    
    for i=1:top_wgrid_nz-1
        top_wgrid_wp2_ta_eqn(i) = -(top_mgrid_wp3(i+1) - top_mgrid_wp3(i))/(top_mgrid_z(i+1)-top_mgrid_z(i));
    end

    top_wgrid_wp2_ta_eqn(top_wgrid_nz) = top_wgrid_wp2_ta_eqn(top_wgrid_nz-1);

    top_wgrid_wp2_tp = top_wgrid_wp2_ma + top_wgrid_wp2_ta + top_wgrid_wp2_ac;
    top_wgrid_wp2_pr = top_wgrid_wp2_pr1 + top_wgrid_wp2_pr2 + top_wgrid_wp2_pr3;
    top_wgrid_wp2_dp = top_wgrid_wp2_dp1 + top_wgrid_wp2_dp2 + top_wgrid_wp2_cl;
end

if (type2 == 1)
    for i=2:bot_wgrid_nz-1
        bot_wgrid_wp2_ta_eqn(i) = -(bot_wgrid_wp3(i+1) - bot_wgrid_wp3(i-1))/(bot_wgrid_z(i+1)-bot_wgrid_z(i-1));
    end

    bot_wgrid_wp2_ta_eqn(1) = -(bot_wgrid_wp3(2) - bot_wgrid_wp3(1))/(bot_wgrid_z(2)-bot_wgrid_z(1));
    bot_wgrid_wp2_ta_eqn(bot_wgrid_nz) = -(bot_wgrid_wp3(bot_wgrid_nz) - bot_wgrid_wp3(bot_wgrid_nz-1))/(bot_wgrid_z(bot_wgrid_nz)-bot_wgrid_z(bot_wgrid_nz-1));
else
    for i=1:bot_wgrid_nz-1
        bot_wgrid_wp2_ta_eqn(i) = -(bot_mgrid_wp3(i+1) - bot_mgrid_wp3(i))/(bot_mgrid_z(i+1)-bot_mgrid_z(i));
    end

%    bot_wgrid_wp2_ta_eqn(bot_wgrid_nz-1) = bot_wgrid_wp2_ta_eqn(bot_wgrid_nz-2);
    bot_wgrid_wp2_ta_eqn(bot_wgrid_nz) = bot_wgrid_wp2_ta_eqn(bot_wgrid_nz-1);
    bot_wgrid_wp2_tp = bot_wgrid_wp2_ma + bot_wgrid_wp2_ta + bot_wgrid_wp2_ac;
    bot_wgrid_wp2_pr = bot_wgrid_wp2_pr1 + bot_wgrid_wp2_pr2 + bot_wgrid_wp2_pr3;
    bot_wgrid_wp2_dp = bot_wgrid_wp2_dp1 + bot_wgrid_wp2_dp2 + bot_wgrid_wp2_cl;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATIONS FOR wp3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (type1 == 1)
    for i=2:top_wgrid_nz-1
        top_wgrid_wp3_ta_eqn(i) = -(top_wgrid_wp4(i+1) - top_wgrid_wp4(i-1))/(top_wgrid_z(i+1)-top_wgrid_z(i-1));
        top_wgrid_wp3_tp_eqn(i) = 3*top_wgrid_wp2(i)*((top_wgrid_wp2(i+1) - top_wgrid_wp2(i-1))/(top_wgrid_z(i+1)-top_wgrid_z(i-1)));
    end

    top_wgrid_wp3_ta_eqn(1) = -(top_wgrid_wp4(2) - top_wgrid_wp4(1))/(top_wgrid_z(2)-top_wgrid_z(1));
    top_wgrid_wp3_tp_eqn(1) = 3*top_wgrid_wp2(1)*((top_wgrid_wp2(2) - top_wgrid_wp2(1))/(top_wgrid_z(2)-top_wgrid_z(1)));

    top_wgrid_wp3_ta_eqn(top_wgrid_nz) = -(top_wgrid_wp4(top_wgrid_nz) - top_wgrid_wp4(top_wgrid_nz-1))/(top_wgrid_z(top_wgrid_nz)-top_wgrid_z(top_wgrid_nz-1));
    top_wgrid_wp3_tp_eqn(top_wgrid_nz) = 3*top_wgrid_wp2(top_wgrid_nz)*((top_wgrid_wp2(top_wgrid_nz) - top_wgrid_wp2(top_wgrid_nz-1))/(top_wgrid_z(top_wgrid_nz)-top_wgrid_z(top_wgrid_nz-1)));
else
    top_mgrid_wp3_tp = top_mgrid_wp3_ma + top_mgrid_wp3_ta + top_mgrid_wp3_tp + top_mgrid_wp3_ac;
    top_mgrid_wp3_pr = top_mgrid_wp3_pr1 + top_mgrid_wp3_pr2;
    top_mgrid_wp3_dp = top_mgrid_wp3_dp1 + top_mgrid_wp3_cl;

    for i=1:top_mgrid_nz-1
        top_mgrid_wp3_ta_eqn(i) = -(top_wgrid_wp4(i+1) - top_wgrid_wp4(i)) / (top_wgrid_z(i+1) - top_wgrid_z(i));
        top_mgrid_wp3_tp_eqn(i) = 3 * (top_wgrid_wp2(i) + top_wgrid_wp2(i+1)) * (top_wgrid_wp2(i+1) - top_wgrid_wp2(i)) / (top_wgrid_z(i+1) - top_wgrid_z(i));
    end

    top_mgrid_wp3_ta_eqn(top_mgrid_nz) = top_mgrid_wp3_ta_eqn(top_mgrid_nz-1);
    top_mgrid_wp3_tp_eqn(top_mgrid_nz) = top_mgrid_wp3_tp_eqn(top_mgrid_nz-1);
    top_mgrid_wp3_sum_eqn = top_mgrid_wp3_ta_eqn + top_mgrid_wp3_tp_eqn;
end

if (type2 == 1)
    for i=2:bot_wgrid_nz-1
        bot_wgrid_wp3_ta_eqn(i) = -(bot_wgrid_wp4(i+1) - bot_wgrid_wp4(i-1))/(bot_wgrid_z(i+1)-bot_wgrid_z(i-1));
        bot_wgrid_wp3_tp_eqn(i) = 3*bot_wgrid_wp2(i)*((bot_wgrid_wp2(i+1) - bot_wgrid_wp2(i-1))/(bot_wgrid_z(i+1)-bot_wgrid_z(i-1)));
    end

    bot_wgrid_wp3_ta_eqn(1) = -(bot_wgrid_wp4(2) - bot_wgrid_wp4(1))/(bot_wgrid_z(2)-bot_wgrid_z(1));
    bot_wgrid_wp3_tp_eqn(1) = 3*bot_wgrid_wp2(1)*((bot_wgrid_wp2(2) - bot_wgrid_wp2(1))/(bot_wgrid_z(2)-bot_wgrid_z(1)));

    bot_wgrid_wp3_ta_eqn(bot_wgrid_nz) = -(bot_wgrid_wp4(bot_wgrid_nz) - bot_wgrid_wp4(bot_wgrid_nz-1))/(bot_wgrid_z(bot_wgrid_nz)-bot_wgrid_z(bot_wgrid_nz-1));
    bot_wgrid_wp3_tp_eqn(bot_wgrid_nz) = 3*bot_wgrid_wp2(bot_wgrid_nz)*((bot_wgrid_wp2(bot_wgrid_nz) - bot_wgrid_wp2(bot_wgrid_nz-1))/(bot_wgrid_z(bot_wgrid_nz)-bot_wgrid_z(bot_wgrid_nz-1)));
else
    bot_mgrid_wp3_tp = bot_mgrid_wp3_ma + bot_mgrid_wp3_ta + bot_mgrid_wp3_tp + bot_mgrid_wp3_ac;
    bot_mgrid_wp3_pr = bot_mgrid_wp3_pr1 + bot_mgrid_wp3_pr2;
    bot_mgrid_wp3_dp = bot_mgrid_wp3_dp1 + bot_mgrid_wp3_cl;

    for i=1:bot_mgrid_nz-1
        bot_mgrid_wp3_ta_eqn(i) = -(bot_wgrid_wp4(i+1) - bot_wgrid_wp4(i)) / (bot_wgrid_z(i+1) - bot_wgrid_z(i));
        bot_mgrid_wp3_tp_eqn(i) = 3 * (bot_wgrid_wp2(i) + bot_wgrid_wp2(i+1)) * (bot_wgrid_wp2(i+1) - bot_wgrid_wp2(i)) / (bot_wgrid_z(i+1) - bot_wgrid_z(i));
    end

    bot_mgrid_wp3_ta_eqn(bot_mgrid_nz) = bot_mgrid_wp3_ta_eqn(bot_mgrid_nz-1);
    bot_mgrid_wp3_tp_eqn(bot_mgrid_nz) = bot_mgrid_wp3_tp_eqn(bot_mgrid_nz-1);
    bot_mgrid_wp3_sum_eqn = bot_mgrid_wp3_ta_eqn + bot_mgrid_wp3_tp_eqn;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATIONS OF wpthlp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (type1 == 1)
    for i=2:top_mgrid_nz-1
        top_mgrid_wpthlp_ta_eqn(i) = -(top_mgrid_wp2thlp(i+1) - top_mgrid_wp2thlp(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
        top_mgrid_wpthlp_tp_eqn(i) = -top_mgrid_wp2(i) * (top_mgrid_thlm(i+1) - top_mgrid_thlm(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
    end
    
    top_mgrid_wpthlp_ta_eqn(1) = -(top_mgrid_wp2thlp(2) - top_mgrid_wp2thlp(1))/(top_mgrid_z(2)-top_mgrid_z(1));
    top_mgrid_wpthlp_tp_eqn(1) = -top_mgrid_wp2(1) * (top_mgrid_thlm(2) - top_mgrid_thlm(1))/(top_mgrid_z(2)-top_mgrid_z(1));

    top_mgrid_wpthlp_ta_eqn(top_mgrid_nz) = -(top_mgrid_wp2thlp(top_mgrid_nz) - top_mgrid_wp2thlp(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));
    top_mgrid_wpthlp_tp_eqn(top_mgrid_nz) = -top_mgrid_wp2(top_mgrid_nz) * (top_mgrid_thlm(top_mgrid_nz) - top_mgrid_thlm(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));

%% Factor of 0.5 because wpqcp and wpthp are on w grid while pm is on m
%% grid.  It is a simple arithmetic average.

    for i=1:top_mgrid_nz
        thlterm = (((p0/top_mgrid_pm(i))^(Rd/cp))*(L/cp));
        top_mgrid_wpthlp_tp(i) = (0.5*(top_wgrid_wpthp_tp(i)+top_wgrid_wpthp_tp(i+1))) - (0.5*(top_wgrid_wpqcp_tp(i)+top_wgrid_wpqcp_tp(i+1)))*thlterm;
        top_mgrid_wpthlp_pr(i) = (0.5*(top_wgrid_wpthp_pr(i)+top_wgrid_wpthp_pr(i+1))) - (0.5*(top_wgrid_wpqcp_pr(i)+top_wgrid_wpqcp_pr(i+1)))*thlterm;
        top_mgrid_wpthlp_dp(i) = (0.5*(top_wgrid_wpthp_dp(i)+top_wgrid_wpthp_dp(i+1))) - (0.5*(top_wgrid_wpqcp_dp(i)+top_wgrid_wpqcp_dp(i+1)))*thlterm;
        top_mgrid_wpthlp_bp(i) = (0.5*(top_wgrid_wpthp_bp(i)+top_wgrid_wpthp_bp(i+1))) - (0.5*(top_wgrid_wpqcp_bp(i)+top_wgrid_wpqcp_bp(i+1)))*thlterm;
    end
else
    for i=2:top_wgrid_nz-1
        top_wgrid_wpthlp_ta_eqn(i) = -(top_mgrid_wp2thlp(i+1) - top_mgrid_wp2thlp(i))/(top_mgrid_z(i+1)-top_mgrid_z(i));
        top_wgrid_wpthlp_tp_eqn(i) = -top_wgrid_wp2(i) * (top_mgrid_thlm(i+1) - top_mgrid_thlm(i))/(top_mgrid_z(i+1)-top_mgrid_z(i));
    end

    top_wgrid_wpthlp_ta_eqn(1)    = top_wgrid_wpthlp_ta_eqn(2);
    top_wgrid_wpthlp_ta_eqn(top_wgrid_nz) = top_wgrid_wpthlp_ta_eqn(top_wgrid_nz-1);
    top_wgrid_wpthlp_tp_eqn(1)    = top_wgrid_wpthlp_tp_eqn(2);
    top_wgrid_wpthlp_tp_eqn(top_wgrid_nz) = top_wgrid_wpthlp_tp_eqn(top_wgrid_nz-1);

    top_wgrid_wpthlp_sum_eqn      = top_wgrid_wpthlp_tp_eqn + top_wgrid_wpthlp_ta_eqn;

    top_wgrid_wpthlp_tp = top_wgrid_wpthlp_ma + top_wgrid_wpthlp_ta + top_wgrid_wpthlp_tp + top_wgrid_wpthlp_ac;
    top_wgrid_wpthlp_pr = top_wgrid_wpthlp_pr1 + top_wgrid_wpthlp_pr2 + top_wgrid_wpthlp_pr3;
    top_wgrid_wpthlp_dp = top_wgrid_wpthlp_dp1;
end

if (type2 == 1)
    for i=2:bot_mgrid_nz-1
        bot_mgrid_wpthlp_ta_eqn(i) = -(bot_mgrid_wp2thlp(i+1) - bot_mgrid_wp2thlp(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
        bot_mgrid_wpthlp_tp_eqn(i) = -bot_mgrid_wp2(i) * (bot_mgrid_thlm(i+1) - bot_mgrid_thlm(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
    end
    
    bot_mgrid_wpthlp_ta_eqn(1) = -(bot_mgrid_wp2thlp(2) - bot_mgrid_wp2thlp(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));
    bot_mgrid_wpthlp_tp_eqn(1) = -bot_mgrid_wp2(1) * (bot_mgrid_thlm(2) - bot_mgrid_thlm(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));

    bot_mgrid_wpthlp_ta_eqn(bot_mgrid_nz) = -(bot_mgrid_wp2thlp(bot_mgrid_nz) - bot_mgrid_wp2thlp(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));
    bot_mgrid_wpthlp_tp_eqn(bot_mgrid_nz) = -bot_mgrid_wp2(bot_mgrid_nz) * (bot_mgrid_thlm(bot_mgrid_nz) - bot_mgrid_thlm(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));

%% Factor of 0.5 because wpqcp and wpthp are on w grid while pm is on m
%% grid.  It is a simple arithmetic average.

    for i=1:bot_mgrid_nz
        thlterm = (((p0/bot_mgrid_pm(i))^(Rd/cp))*(L/cp));
        bot_mgrid_wpthlp_tp(i) = (0.5*(bot_wgrid_wpthp_tp(i)+bot_wgrid_wpthp_tp(i+1))) - (0.5*(bot_wgrid_wpqcp_tp(i)+bot_wgrid_wpqcp_tp(i+1)))*thlterm;
        bot_mgrid_wpthlp_pr(i) = (0.5*(bot_wgrid_wpthp_pr(i)+bot_wgrid_wpthp_pr(i+1))) - (0.5*(bot_wgrid_wpqcp_pr(i)+bot_wgrid_wpqcp_pr(i+1)))*thlterm;
        bot_mgrid_wpthlp_dp(i) = (0.5*(bot_wgrid_wpthp_dp(i)+bot_wgrid_wpthp_dp(i+1))) - (0.5*(bot_wgrid_wpqcp_dp(i)+bot_wgrid_wpqcp_dp(i+1)))*thlterm;
        bot_mgrid_wpthlp_bp(i) = (0.5*(bot_wgrid_wpthp_bp(i)+bot_wgrid_wpthp_bp(i+1))) - (0.5*(bot_wgrid_wpqcp_bp(i)+bot_wgrid_wpqcp_bp(i+1)))*thlterm;
    end
else
    for i=2:bot_wgrid_nz-1
        bot_wgrid_wpthlp_ta_eqn(i) = -(bot_mgrid_wp2thlp(i+1) - bot_mgrid_wp2thlp(i))/(bot_mgrid_z(i+1)-bot_mgrid_z(i));
        bot_wgrid_wpthlp_tp_eqn(i) = -bot_wgrid_wp2(i) * (bot_mgrid_thlm(i+1) - bot_mgrid_thlm(i))/(bot_mgrid_z(i+1)-bot_mgrid_z(i));
    end

    bot_wgrid_wpthlp_ta_eqn(1)    = bot_wgrid_wpthlp_ta_eqn(2);
    bot_wgrid_wpthlp_ta_eqn(bot_wgrid_nz) = bot_wgrid_wpthlp_ta_eqn(bot_wgrid_nz-1);
    bot_wgrid_wpthlp_tp_eqn(1)    = bot_wgrid_wpthlp_tp_eqn(2);
    bot_wgrid_wpthlp_tp_eqn(bot_wgrid_nz) = bot_wgrid_wpthlp_tp_eqn(bot_wgrid_nz-1);

    bot_wgrid_wpthlp_sum_eqn      = bot_wgrid_wpthlp_tp_eqn + bot_wgrid_wpthlp_ta_eqn;

    bot_wgrid_wpthlp_tp = bot_wgrid_wpthlp_ma + bot_wgrid_wpthlp_ta + bot_wgrid_wpthlp_tp + bot_wgrid_wpthlp_ac;
    bot_wgrid_wpthlp_pr = bot_wgrid_wpthlp_pr1 + bot_wgrid_wpthlp_pr2 + bot_wgrid_wpthlp_pr3;
    bot_wgrid_wpthlp_dp = bot_wgrid_wpthlp_dp1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATIONS OF wpqtp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wpqtp is on m-grid
%budget is on w-grid

if (type1 == 1)
    for i=2:top_mgrid_nz-1
        top_mgrid_wpqtp_ta_eqn(i) = -(top_mgrid_wp2qtp(i+1) - top_mgrid_wp2qtp(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
        top_mgrid_wpqtp_tp_eqn(i) = -top_mgrid_wp2(i) * (top_mgrid_qtm(i+1) - top_mgrid_qtm(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
    end
    top_mgrid_wpqtp_ta_eqn(1) = -(top_mgrid_wp2qtp(2) - top_mgrid_wp2qtp(1))/(top_mgrid_z(2)-top_mgrid_z(1));
    top_mgrid_wpqtp_tp_eqn(1) = -top_mgrid_wp2(1) * (top_mgrid_qtm(2) - top_mgrid_qtm(1))/(top_mgrid_z(2)-top_mgrid_z(1));
    
    top_mgrid_wpqtp_ta_eqn(top_mgrid_nz) = -(top_mgrid_wp2qtp(top_mgrid_nz) - top_mgrid_wp2qtp(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));
    top_mgrid_wpqtp_tp_eqn(top_mgrid_nz) = -top_mgrid_wp2(top_mgrid_nz) * (top_mgrid_qtm(top_mgrid_nz) - top_mgrid_qtm(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));

    top_wgrid_wpqtp_dp = top_wgrid_wpqvp_dp + top_wgrid_wpqcp_dp;
    top_wgrid_wpqtp_tp = top_wgrid_wpqvp_tp + top_wgrid_wpqcp_tp;
    top_wgrid_wpqtp_pr = top_wgrid_wpqvp_pr + top_wgrid_wpqcp_pr;
    top_wgrid_wpqtp_bp = top_wgrid_wpqvp_bp + top_wgrid_wpqcp_bp;
else
    for i=1:top_mgrid_nz-1
        top_wgrid_wprtp_ta_eqn(i) = -(top_mgrid_wp2rtp(i+1)-top_mgrid_wp2rtp(i)) / (top_mgrid_z(i+1)-top_mgrid_z(i));
        top_wgrid_wprtp_tp_eqn(i) = -top_wgrid_wp2(i) * (top_mgrid_rtm(i+1)-top_mgrid_rtm(i)) / (top_mgrid_z(i+1)-top_mgrid_z(i));
    end

    top_wgrid_wprtp_ta_eqn(top_mgrid_nz) = top_wgrid_wprtp_ta_eqn(top_mgrid_nz-1);
    top_wgrid_wprtp_tp_eqn(top_mgrid_nz) = top_wgrid_wprtp_tp_eqn(top_mgrid_nz-1);
    top_wgrid_wprtp_sum_eqn = top_wgrid_wprtp_tp_eqn + top_wgrid_wprtp_ta_eqn;   
    top_wgrid_wprtp_tp = top_wgrid_wprtp_ma + top_wgrid_wprtp_ta + top_wgrid_wprtp_tp + top_wgrid_wprtp_ac;
    top_wgrid_wprtp_pr = top_wgrid_wprtp_pr1 + top_wgrid_wprtp_pr2 + top_wgrid_wprtp_pr3;
    top_wgrid_wprtp_dp = top_wgrid_wprtp_dp1;
end

if (type2 == 1)
    for i=2:bot_mgrid_nz-1
        bot_mgrid_wpqtp_ta_eqn(i) = -(bot_mgrid_wp2qtp(i+1) - bot_mgrid_wp2qtp(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
        bot_mgrid_wpqtp_tp_eqn(i) = -bot_mgrid_wp2(i) * (bot_mgrid_qtm(i+1) - bot_mgrid_qtm(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
    end
    bot_mgrid_wpqtp_ta_eqn(1) = -(bot_mgrid_wp2qtp(2) - bot_mgrid_wp2qtp(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));
    bot_mgrid_wpqtp_tp_eqn(1) = -bot_mgrid_wp2(1) * (bot_mgrid_qtm(2) - bot_mgrid_qtm(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));
    
    bot_mgrid_wpqtp_ta_eqn(bot_mgrid_nz) = -(bot_mgrid_wp2qtp(bot_mgrid_nz) - bot_mgrid_wp2qtp(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));
    bot_mgrid_wpqtp_tp_eqn(bot_mgrid_nz) = -bot_mgrid_wp2(bot_mgrid_nz) * (bot_mgrid_qtm(bot_mgrid_nz) - bot_mgrid_qtm(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));

    bot_wgrid_wpqtp_dp = bot_wgrid_wpqvp_dp + bot_wgrid_wpqcp_dp;
    bot_wgrid_wpqtp_tp = bot_wgrid_wpqvp_tp + bot_wgrid_wpqcp_tp;
    bot_wgrid_wpqtp_pr = bot_wgrid_wpqvp_pr + bot_wgrid_wpqcp_pr;
    bot_wgrid_wpqtp_bp = bot_wgrid_wpqvp_bp + bot_wgrid_wpqcp_bp;
else
    for i=1:bot_mgrid_nz-1
        bot_wgrid_wprtp_ta_eqn(i) = -(bot_mgrid_wp2rtp(i+1)-bot_mgrid_wp2rtp(i)) / (bot_mgrid_z(i+1)-bot_mgrid_z(i));
        bot_wgrid_wprtp_tp_eqn(i) = -bot_wgrid_wp2(i) * (bot_mgrid_rtm(i+1)-bot_mgrid_rtm(i)) / (bot_mgrid_z(i+1)-bot_mgrid_z(i));
    end

    bot_wgrid_wprtp_ta_eqn(bot_mgrid_nz) = bot_wgrid_wprtp_ta_eqn(bot_mgrid_nz-1);
    bot_wgrid_wprtp_tp_eqn(bot_mgrid_nz) = bot_wgrid_wprtp_tp_eqn(bot_mgrid_nz-1);
    bot_wgrid_wprtp_sum_eqn = bot_wgrid_wprtp_tp_eqn + bot_wgrid_wprtp_ta_eqn;   
    bot_wgrid_wprtp_tp = bot_wgrid_wprtp_ma + bot_wgrid_wprtp_ta + bot_wgrid_wprtp_tp + bot_wgrid_wprtp_ac;
    bot_wgrid_wprtp_pr = bot_wgrid_wprtp_pr1 + bot_wgrid_wprtp_pr2 + bot_wgrid_wprtp_pr3;
    bot_wgrid_wprtp_dp = bot_wgrid_wprtp_dp1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATIONS OF thlp2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (type1 == 1)
    for i=2:top_mgrid_nz-1
        top_mgrid_thlp2_ta_eqn(i) = -(top_mgrid_wpthlp2(i+1) - top_mgrid_wpthlp2(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
        top_mgrid_thlp2_tp_eqn(i) = -2*top_mgrid_wpthlp(i)*(top_mgrid_thlm(i+1) - top_mgrid_thlm(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
    end

    top_mgrid_thlp2_ta_eqn(1) = -(top_mgrid_wpthlp2(2) - top_mgrid_wpthlp2(1))/(top_mgrid_z(2)-top_mgrid_z(1));
    top_mgrid_thlp2_tp_eqn(1) = -2*top_mgrid_wpthlp(1)*(top_mgrid_thlm(2) - top_mgrid_thlm(1))/(top_mgrid_z(2)-top_mgrid_z(1));

    top_mgrid_thlp2_ta_eqn(top_mgrid_nz) = -(top_mgrid_wpthlp2(top_mgrid_nz) - top_mgrid_wpthlp2(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));
    top_mgrid_thlp2_tp_eqn(top_mgrid_nz) = -2*top_mgrid_wpthlp(top_mgrid_nz)*(top_mgrid_thlm(top_mgrid_nz) - top_mgrid_thlm(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));

    for i=1:top_mgrid_nz
        thlterm = (((p0/top_mgrid_pm(i)).^(Rd/cp))*(L/cp));
        top_mgrid_thlp2_tp(i)  = top_mgrid_thp2_tp(i)  + top_mgrid_qcp2_tp(i)*thlterm*thlterm - 2*top_mgrid_thpqcp_tp(i)*thlterm;
        top_mgrid_thlp2_dp(i)  = top_mgrid_thp2_dp(i)  + top_mgrid_qcp2_dp(i)*thlterm*thlterm - 2*top_mgrid_thpqcp_dp(i)*thlterm;
        top_mgrid_thlp2_bt(i)  = top_mgrid_thp2_bt(i)  + top_mgrid_qcp2_bt(i)*thlterm*thlterm - 2*top_mgrid_thpqcp_bt(i)*thlterm;
        top_mgrid_thlp2_rad(i) = top_mgrid_thp2_rad(i) - 2*top_mgrid_thpqcp_rad(i)*thlterm;
    end

    top_mgrid_thlp2_dp_eqn = top_mgrid_thlp2_bt - top_mgrid_thlp2_tp_eqn - top_mgrid_thlp2_ta_eqn - top_mgrid_thlp2_rad;
end

if (type2 == 1)
    for i=2:bot_mgrid_nz-1
        bot_mgrid_thlp2_ta_eqn(i) = -(bot_mgrid_wpthlp2(i+1) - bot_mgrid_wpthlp2(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
        bot_mgrid_thlp2_tp_eqn(i) = -2*bot_mgrid_wpthlp(i)*(bot_mgrid_thlm(i+1) - bot_mgrid_thlm(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
    end

    bot_mgrid_thlp2_ta_eqn(1) = -(bot_mgrid_wpthlp2(2) - bot_mgrid_wpthlp2(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));
    bot_mgrid_thlp2_tp_eqn(1) = -2*bot_mgrid_wpthlp(1)*(bot_mgrid_thlm(2) - bot_mgrid_thlm(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));

    bot_mgrid_thlp2_ta_eqn(bot_mgrid_nz) = -(bot_mgrid_wpthlp2(bot_mgrid_nz) - bot_mgrid_wpthlp2(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));
    bot_mgrid_thlp2_tp_eqn(bot_mgrid_nz) = -2*bot_mgrid_wpthlp(bot_mgrid_nz)*(bot_mgrid_thlm(bot_mgrid_nz) - bot_mgrid_thlm(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));

    for i=1:bot_mgrid_nz
        thlterm = (((p0/bot_mgrid_pm(i)).^(Rd/cp))*(L/cp));
        bot_mgrid_thlp2_tp(i)  = bot_mgrid_thp2_tp(i)  + bot_mgrid_qcp2_tp(i)*thlterm*thlterm - 2*bot_mgrid_thpqcp_tp(i)*thlterm;
        bot_mgrid_thlp2_dp(i)  = bot_mgrid_thp2_dp(i)  + bot_mgrid_qcp2_dp(i)*thlterm*thlterm - 2*bot_mgrid_thpqcp_dp(i)*thlterm;
        bot_mgrid_thlp2_bt(i)  = bot_mgrid_thp2_bt(i)  + bot_mgrid_qcp2_bt(i)*thlterm*thlterm - 2*bot_mgrid_thpqcp_bt(i)*thlterm;
        bot_mgrid_thlp2_rad(i) = bot_mgrid_thp2_rad(i) - 2*bot_mgrid_thpqcp_rad(i)*thlterm;
    end

    bot_mgrid_thlp2_dp_eqn = bot_mgrid_thlp2_bt - bot_mgrid_thlp2_tp_eqn - bot_mgrid_thlp2_ta_eqn - bot_mgrid_thlp2_rad;
end


%%%%%%%
% SCM %
%%%%%%%

%for i=2:bot_mgrid_nz-1
%    bot_mgrid_thlp2_ma_eqn(i) = -bot_mgrid_wm(i) * (bot_wgrid_thlp2(i) - bot_wgrid_thlp2(i-1))/(bot_wgrid_z(i)-bot_wgrid_z(i-1));
%    bot_mgrid_thlp2_ta_eqn(i) = -(bot_mgrid_wpthlp2(i+1)-(bot_mgrid_wpthlp2(i-1)))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
%    bot_mgrid_thlp2_tp_eqn(i) = -2*0.5*(bot_wgrid_wpthlp(i)+bot_wgrid_wpthlp(i-1))*(bot_mgrid_thlm(i+1)-bot_mgrid_thlm(i-1))/(bot_mgrid_z(i+1)/bot_mgrid_z(i-1));
%end
%bot_mgrid_thlp2_ma_eqn(1) = bot_mgrid_thlp2_ma_eqn(2);
%bot_mgrid_thlp2_ta_eqn(1) = bot_mgrid_thlp2_ta_eqn(2);
%bot_mgrid_thlp2_tp_eqn(1) = bot_mgrid_thlp2_tp_eqn(2);

%bot_mgrid_thlp2_ma_eqn(bot_mgrid_nz) = bot_mgrid_thlp2_ma_eqn(bot_mgrid_nz-1);
%bot_mgrid_thlp2_ta_eqn(bot_mgrid_nz) = bot_mgrid_thlp2_ta_eqn(bot_mgrid_nz-1);
%bot_mgrid_thlp2_tp_eqn(bot_mgrid_nz) = bot_mgrid_thlp2_tp_eqn(bot_mgrid_nz-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATIONS OF qtp2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (type1 == 1)
    for i=2:top_mgrid_nz-1
        top_mgrid_qtp2_ta_eqn(i) = -(top_mgrid_wpqtp2(i+1) - top_mgrid_wpqtp2(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
        top_mgrid_qtp2_tp_eqn(i) = -2*top_mgrid_wpqtp(i)*(top_mgrid_qtm(i+1) - top_mgrid_qtm(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
    end
    top_mgrid_qtp2_ta_eqn(1) = -(top_mgrid_wpqtp2(2) - top_mgrid_wpqtp2(1))/(top_mgrid_z(2)-top_mgrid_z(1));
    top_mgrid_qtp2_tp_eqn(1) = -2*top_mgrid_wpqtp(1)*(top_mgrid_qtm(2) - top_mgrid_qtm(1))/(top_mgrid_z(2)-top_mgrid_z(1));

    top_mgrid_qtp2_ta_eqn(top_mgrid_nz) = -(top_mgrid_wpqtp2(top_mgrid_nz) - top_mgrid_wpqtp2(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));
    top_mgrid_qtp2_tp_eqn(top_mgrid_nz) = -2*top_mgrid_wpqtp(top_mgrid_nz)*(top_mgrid_qtm(top_mgrid_nz) - top_mgrid_qtm(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));

    for i=1:top_mgrid_nz
        thlterm = (((p0/top_mgrid_pm(i)).^(Rd/cp))*(L/cp));
        top_mgrid_qtp2_tp(i)  = top_mgrid_qcp2_tp(i) + top_mgrid_qvp2_tp(i) + 2*top_mgrid_qvpqcp_tp(i);
        top_mgrid_qtp2_dp(i)  = top_mgrid_qcp2_dp(i) + top_mgrid_qvp2_dp(i) + 2*top_mgrid_qvpqcp_dp(i);
        top_mgrid_qtp2_bt(i)  = top_mgrid_qcp2_bt(i) + top_mgrid_qvp2_bt(i) + 2*top_mgrid_qvpqcp_bt(i);
    end

    top_mgrid_qtp2_dp_eqn = top_mgrid_qtp2_bt - top_mgrid_qtp2_tp_eqn - top_mgrid_qtp2_ta_eqn;
end

if (type2 == 1)
    for i=2:bot_mgrid_nz-1
        bot_mgrid_qtp2_ta_eqn(i) = -(bot_mgrid_wpqtp2(i+1) - bot_mgrid_wpqtp2(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
        bot_mgrid_qtp2_tp_eqn(i) = -2*bot_mgrid_wpqtp(i)*(bot_mgrid_qtm(i+1) - bot_mgrid_qtm(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
    end
    bot_mgrid_qtp2_ta_eqn(1) = -(bot_mgrid_wpqtp2(2) - bot_mgrid_wpqtp2(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));
    bot_mgrid_qtp2_tp_eqn(1) = -2*bot_mgrid_wpqtp(1)*(bot_mgrid_qtm(2) - bot_mgrid_qtm(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));

    bot_mgrid_qtp2_ta_eqn(bot_mgrid_nz) = -(bot_mgrid_wpqtp2(bot_mgrid_nz) - bot_mgrid_wpqtp2(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));
    bot_mgrid_qtp2_tp_eqn(bot_mgrid_nz) = -2*bot_mgrid_wpqtp(bot_mgrid_nz)*(bot_mgrid_qtm(bot_mgrid_nz) - bot_mgrid_qtm(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));

    for i=1:bot_mgrid_nz
        thlterm = (((p0/bot_mgrid_pm(i)).^(Rd/cp))*(L/cp));
        bot_mgrid_qtp2_tp(i)  = bot_mgrid_qcp2_tp(i) + bot_mgrid_qvp2_tp(i) + 2*bot_mgrid_qvpqcp_tp(i);
        bot_mgrid_qtp2_dp(i)  = bot_mgrid_qcp2_dp(i) + bot_mgrid_qvp2_dp(i) + 2*bot_mgrid_qvpqcp_dp(i);
        bot_mgrid_qtp2_bt(i)  = bot_mgrid_qcp2_bt(i) + bot_mgrid_qvp2_bt(i) + 2*bot_mgrid_qvpqcp_bt(i);
    end

    bot_mgrid_qtp2_dp_eqn = bot_mgrid_qtp2_bt - bot_mgrid_qtp2_tp_eqn - bot_mgrid_qtp2_ta_eqn;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATIONS OF qtpthlp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (type1 == 1)
    for i=2:top_mgrid_nz-1
        top_mgrid_qtpthlp_ta_eqn(i)  = -(top_mgrid_wpqtpthlp(i+1) - top_mgrid_wpqtpthlp(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
        top_mgrid_qtpthlp_tp1_eqn(i) = -top_mgrid_wpqtp(i) * (top_mgrid_thlm(i+1) - top_mgrid_thlm(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
        top_mgrid_qtpthlp_tp2_eqn(i) = -top_mgrid_wpthlp(i) * (top_mgrid_qtm(i+1) - top_mgrid_qtm(i-1))/(top_mgrid_z(i+1)-top_mgrid_z(i-1));
    end

    top_mgrid_qtpthlp_ta_eqn(1)  = -(top_mgrid_wpqtpthlp(2) - top_mgrid_wpqtpthlp(1))/(top_mgrid_z(2)-top_mgrid_z(1));
    top_mgrid_qtpthlp_tp1_eqn(1) = -top_mgrid_wpqtp(1) * (top_mgrid_thlm(2) - top_mgrid_thlm(1))/(top_mgrid_z(2)-top_mgrid_z(1));
    top_mgrid_qtpthlp_tp2_eqn(1) = -top_mgrid_wpthlp(1) * (top_mgrid_qtm(2) - top_mgrid_qtm(1))/(top_mgrid_z(2)-top_mgrid_z(1));

    top_mgrid_qtpthlp_ta_eqn(top_mgrid_nz)  = -(top_mgrid_wpqtpthlp(top_mgrid_nz) - top_mgrid_wpqtpthlp(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));
    top_mgrid_qtpthlp_tp1_eqn(top_mgrid_nz) = -top_mgrid_wpqtp(top_mgrid_nz) * (top_mgrid_thlm(top_mgrid_nz) - top_mgrid_thlm(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));
    top_mgrid_qtpthlp_tp2_eqn(top_mgrid_nz) = -top_mgrid_wpthlp(top_mgrid_nz) * (top_mgrid_qtm(top_mgrid_nz) - top_mgrid_qtm(top_mgrid_nz-1))/(top_mgrid_z(top_mgrid_nz)-top_mgrid_z(top_mgrid_nz-1));

    for i=1:top_mgrid_nz
        thlterm = (((p0/top_mgrid_pm(i)).^(Rd/cp))*(L/cp));
        top_mgrid_qtpthlp_tp(i) = top_mgrid_thpqcp_tp(i)  + top_mgrid_thpqvp_tp(i) - thlterm*top_mgrid_qcp2_tp(i) - thlterm*top_mgrid_qvpqcp_tp(i);
        top_mgrid_qtpthlp_dp(i) = top_mgrid_thpqcp_dp(i)  + top_mgrid_thpqvp_dp(i) - thlterm*top_mgrid_qcp2_dp(i) - thlterm*top_mgrid_qvpqcp_dp(i);
        top_mgrid_qtpthlp_bt(i) = top_mgrid_thpqcp_bt(i)  + top_mgrid_thpqvp_bt(i) - thlterm*top_mgrid_qcp2_bt(i) - thlterm*top_mgrid_qvpqcp_bt(i);
        top_mgrid_qtpthlp_rad(i) = top_mgrid_thpqcp_rad(i) + top_mgrid_thpqvp_rad(i);
    end

    top_mgrid_qtpthlp_dp_eqn = top_mgrid_qtpthlp_bt - top_mgrid_qtpthlp_tp1_eqn - top_mgrid_qtpthlp_tp2_eqn - top_mgrid_qtpthlp_ta_eqn - top_mgrid_qtpthlp_rad;
end

if (type2 == 1)
    for i=2:bot_mgrid_nz-1
        bot_mgrid_qtpthlp_ta_eqn(i)  = -(bot_mgrid_wpqtpthlp(i+1) - bot_mgrid_wpqtpthlp(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
        bot_mgrid_qtpthlp_tp1_eqn(i) = -bot_mgrid_wpqtp(i) * (bot_mgrid_thlm(i+1) - bot_mgrid_thlm(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
        bot_mgrid_qtpthlp_tp2_eqn(i) = -bot_mgrid_wpthlp(i) * (bot_mgrid_qtm(i+1) - bot_mgrid_qtm(i-1))/(bot_mgrid_z(i+1)-bot_mgrid_z(i-1));
    end

    bot_mgrid_qtpthlp_ta_eqn(1)  = -(bot_mgrid_wpqtpthlp(2) - bot_mgrid_wpqtpthlp(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));
    bot_mgrid_qtpthlp_tp1_eqn(1) = -bot_mgrid_wpqtp(1) * (bot_mgrid_thlm(2) - bot_mgrid_thlm(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));
    bot_mgrid_qtpthlp_tp2_eqn(1) = -bot_mgrid_wpthlp(1) * (bot_mgrid_qtm(2) - bot_mgrid_qtm(1))/(bot_mgrid_z(2)-bot_mgrid_z(1));

    bot_mgrid_qtpthlp_ta_eqn(bot_mgrid_nz)  = -(bot_mgrid_wpqtpthlp(bot_mgrid_nz) - bot_mgrid_wpqtpthlp(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));
    bot_mgrid_qtpthlp_tp1_eqn(bot_mgrid_nz) = -bot_mgrid_wpqtp(bot_mgrid_nz) * (bot_mgrid_thlm(bot_mgrid_nz) - bot_mgrid_thlm(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));
    bot_mgrid_qtpthlp_tp2_eqn(bot_mgrid_nz) = -bot_mgrid_wpthlp(bot_mgrid_nz) * (bot_mgrid_qtm(bot_mgrid_nz) - bot_mgrid_qtm(bot_mgrid_nz-1))/(bot_mgrid_z(bot_mgrid_nz)-bot_mgrid_z(bot_mgrid_nz-1));

    for i=1:bot_mgrid_nz
        thlterm = (((p0/bot_mgrid_pm(i)).^(Rd/cp))*(L/cp));
        bot_mgrid_qtpthlp_tp(i) = bot_mgrid_thpqcp_tp(i)  + bot_mgrid_thpqvp_tp(i) - thlterm*bot_mgrid_qcp2_tp(i) - thlterm*bot_mgrid_qvpqcp_tp(i);
        bot_mgrid_qtpthlp_dp(i) = bot_mgrid_thpqcp_dp(i)  + bot_mgrid_thpqvp_dp(i) - thlterm*bot_mgrid_qcp2_dp(i) - thlterm*bot_mgrid_qvpqcp_dp(i);
        bot_mgrid_qtpthlp_bt(i) = bot_mgrid_thpqcp_bt(i)  + bot_mgrid_thpqvp_bt(i) - thlterm*bot_mgrid_qcp2_bt(i) - thlterm*bot_mgrid_qvpqcp_bt(i);
        bot_mgrid_qtpthlp_rad(i) = bot_mgrid_thpqcp_rad(i) + bot_mgrid_thpqvp_rad(i);
    end

    bot_mgrid_qtpthlp_dp_eqn = bot_mgrid_qtpthlp_bt - bot_mgrid_qtpthlp_tp1_eqn - bot_mgrid_qtpthlp_tp2_eqn - bot_mgrid_qtpthlp_ta_eqn - bot_mgrid_qtpthlp_rad;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT WP2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orient landscape;

clf
maxdim = 1.1*max([abs(max(top_wgrid_wp2(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
                  abs(max(bot_wgrid_wp2(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
                  ]);
make_plot('1',top_wgrid_wp2,top_wgrid_z,'-k','wp2','(m^2 s^-^2)',maxdim,zbottom,ztop);
make_plot('6',bot_wgrid_wp2,bot_wgrid_z,'-k','wp2','(m^2 s^-^2)',maxdim,zbottom,ztop); 
ylabel('Height (m)');

if (type1 == 1)
    max1 = max ([ abs(max(top_wgrid_wp2_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ... 
    abs(max(top_wgrid_wp2_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ... 
    abs(max(top_wgrid_wp2_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ... 
    abs(max(top_wgrid_wp2_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    ]);
else
    max1 = max ([ abs(max(top_wgrid_wp2_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_pr1(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_pr1(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_pr2(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_pr2(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_pr3(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_pr3(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_dp1(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_dp1(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_dp2(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_dp2(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp2_cl(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp2_cl(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    ]);
end

if (type2 == 1)
    max2 = max ([ abs(max(bot_wgrid_wp2_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ... 
    abs(max(bot_wgrid_wp2_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ... 
    abs(max(bot_wgrid_wp2_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ... 
    abs(max(bot_wgrid_wp2_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    ]);
else
    max2 = max ([ abs(max(bot_wgrid_wp2_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_pr1(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_pr1(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_pr2(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_pr2(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_pr3(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_pr3(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_dp1(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_dp1(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_dp2(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_dp2(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp2_cl(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp2_cl(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    ]);
end

maxdim = 1.1*max([ max1 max2 ]);

if (type1 == 1)
    multi_var = [top_wgrid_wp2_ta_eqn; top_wgrid_wp2_tp'];
    multi_grid = [top_wgrid_z; top_wgrid_z];
    multi_marker = ['--k'; '-bh'];
    multi_legend = ['ta'];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection(m^2 s^-^3)',maxdim,zbottom,ztop);
    make_plot('3',top_wgrid_wp2_pr,top_wgrid_z,'-gx',simname1,'Pressure (m^2 s^-^3)',maxdim,zbottom,ztop);
    make_plot('4',top_wgrid_wp2_dp,top_wgrid_z,'-+ ','','Dissipation (m^2 s^-^3)',maxdim,zbottom,ztop);
    make_plot('5',top_wgrid_wp2_bp,top_wgrid_z,'-rp','','Buoyancy (m^2 s^-^3)',maxdim,zbottom,ztop);
%    title([simname1]);
else
    multi_var = [top_wgrid_wp2_ta_eqn; top_wgrid_wp2_tp'];
    multi_grid = [top_wgrid_z; top_wgrid_z];
    multi_marker = ['--k'; '-bh'];
    multi_legend = ['ta'];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection(m^2 s^-^3)',maxdim,zbottom,ztop);

    multi_var = [top_wgrid_wp2_pr3'; top_wgrid_wp2_pr2'; top_wgrid_wp2_pr1'; top_wgrid_wp2_pr'];
    multi_grid = [top_wgrid_z; top_wgrid_z; top_wgrid_z; top_wgrid_z;];
    multi_marker = ['-k '; ':k '; '--k'; '-gx'];
    multi_legend = ['pr1';'pr2';'pr3'];
    make_multi_plot('3',multi_var,multi_grid,multi_marker,multi_legend,'NE',simname1,'Pressure(m^2 s^-^3)',maxdim,zbottom,ztop);

    multi_var = [top_wgrid_wp2_cl'; top_wgrid_wp2_dp2'; top_wgrid_wp2_dp1'];
    multi_grid = [top_wgrid_z; top_wgrid_z; top_wgrid_z];
    multi_marker = ['-k '; ':k '; '--k'];
    multi_legend = ['dp1';'dp2';'cl '];
    make_multi_plot('4',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Dissipation(m^2 s^-^3)',maxdim,zbottom,ztop);

    make_plot('4',top_wgrid_wp2_dp,top_wgrid_z,'-+ ','','Dissipation (m^2 s^-^3)',maxdim,zbottom,ztop);
    make_plot('5',top_wgrid_wp2_bp,top_wgrid_z,'-rp','','Buoyancy (m^2 s^-^3)',maxdim,zbottom,ztop);
%    title([simname1]);
end

if (type2 == 1)
    multi_var = [bot_wgrid_wp2_ta_eqn; bot_wgrid_wp2_tp'];
    multi_grid = [bot_wgrid_z; bot_wgrid_z];
    multi_marker = ['--k'; '-bh'];
    multi_legend = ['ta'];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection(m^2 s^-^3)',maxdim,zbottom,ztop);
    make_plot('8',bot_wgrid_wp2_pr,bot_wgrid_z,'-gx',simname2,'Pressure (m^2 s^-^3)',maxdim,zbottom,ztop);
    make_plot('9',bot_wgrid_wp2_dp,bot_wgrid_z,'-+ ','','Dissipation (m^2 s^-^3)',maxdim,zbottom,ztop);
    make_plot('10',bot_wgrid_wp2_bp,bot_wgrid_z,'-rp','','Buoyancy (m^2 s^-^3)',maxdim,zbottom,ztop);
%    title([simname2]);
else
    multi_var = [bot_wgrid_wp2_ta_eqn; bot_wgrid_wp2_tp'];
    multi_grid = [bot_wgrid_z; bot_wgrid_z];
    multi_marker = ['--k'; '-bh'];
    multi_legend = ['ta'];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection(m^2 s^-^3)',maxdim,zbottom,ztop);

    multi_var = [bot_wgrid_wp2_pr3'; bot_wgrid_wp2_pr2'; bot_wgrid_wp2_pr1'; bot_wgrid_wp2_pr'];
    multi_grid = [bot_wgrid_z; bot_wgrid_z; bot_wgrid_z; bot_wgrid_z;];
    multi_marker = ['-k '; ':k '; '--k'; '-gx'];
    multi_legend = ['pr1';'pr2';'pr3'];
    make_multi_plot('8',multi_var,multi_grid,multi_marker,multi_legend,'NE',simname2,'Pressure(m^2 s^-^3)',maxdim,zbottom,ztop);

    multi_var = [bot_wgrid_wp2_cl'; bot_wgrid_wp2_dp2'; bot_wgrid_wp2_dp1'];
    multi_grid = [bot_wgrid_z; bot_wgrid_z; bot_wgrid_z];
    multi_marker = ['-k '; ':k '; '--k'];
    multi_legend = ['dp1';'dp2';'cl '];
    make_multi_plot('9',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Dissipation(m^2 s^-^3)',maxdim,zbottom,ztop);

    make_plot('9',bot_wgrid_wp2_dp,bot_wgrid_z,'-+ ','','Dissipation (m^2 s^-^3)',maxdim,zbottom,ztop);
    make_plot('10',bot_wgrid_wp2_bp,bot_wgrid_z,'-rp','','Buoyancy (m^2 s^-^3)',maxdim,zbottom,ztop);
%    title([simname2]);
end

print -depsc wp2;
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT WP3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

if (type1 == 1)
    top_wp3 = top_wgrid_wp3;
    top_grid = top_wgrid_z;
else
    top_wp3 = top_mgrid_wp3;
    top_grid = top_mgrid_z;
end    

if (type2 == 1)
    bot_wp3 = bot_wgrid_wp3;
    bot_grid = bot_wgrid_z;
else
    bot_wp3 = bot_mgrid_wp3;
    bot_grid = bot_mgrid_z;
end    



maxdim = 1.1*max([ abs(max(top_wp3(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wp3(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
                   abs(max(bot_wp3(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wp3(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
                   ]);
make_plot('1',top_wp3,top_grid,'-k','wp3','(m^3 s^-^3)',maxdim,zbottom,ztop);
make_plot('6',bot_wp3,bot_grid,'-k','wp3','(m^3 s^-^3)',maxdim,zbottom,ztop);
ylabel('Height (m)');

if (type1 == 1)
    max1 = ([abs(max(top_wgrid_wp3_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp3_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp3_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp3_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ... 
    abs(max(top_wgrid_wp3_tp_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp3_tp_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp3_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point) + top_wgrid_wp3_tp_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(min(top_wgrid_wp3_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point) + top_wgrid_wp3_tp_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp3_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp3_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp3_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp3_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wp3_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wp3_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    ]);
else
    max1 = ([abs(max(top_mgrid_wp3_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wp3_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wp3_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wp3_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ... 
    abs(max(top_mgrid_wp3_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wp3_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wp3_sum_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wp3_sum_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wp3_pr(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wp3_pr(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wp3_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wp3_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wp3_bp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wp3_bp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    ]);
end

if (type2 == 1)
    max2 = ([abs(max(bot_wgrid_wp3_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp3_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp3_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp3_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ... 
    abs(max(bot_wgrid_wp3_tp_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp3_tp_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp3_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point) + bot_wgrid_wp3_tp_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(min(bot_wgrid_wp3_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point) + bot_wgrid_wp3_tp_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp3_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp3_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp3_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp3_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wp3_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wp3_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    ]);
else
    max2 = ([abs(max(bot_mgrid_wp3_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wp3_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wp3_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wp3_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ... 
    abs(max(bot_mgrid_wp3_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wp3_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wp3_sum_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wp3_sum_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wp3_pr(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wp3_pr(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wp3_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wp3_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wp3_bp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wp3_bp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    ]);
end


maxdim = 1.1*max([ max1 max2 ]);

if (type1 == 1)
    multi_var = [top_wgrid_wp3_tp_eqn+top_wgrid_wp3_ta_eqn; top_wgrid_wp3_ta_eqn; top_wgrid_wp3_tp_eqn; top_wgrid_wp3_tp'];
    multi_grid = [top_wgrid_z; top_wgrid_z; top_wgrid_z; top_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum'; 'ta ';'tp '];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection(m^3 s^-^4)',maxdim,zbottom,ztop);
    make_plot('3',top_wgrid_wp3_pr,top_wgrid_z,'-gx',simname1,'Pressure (m^3 s^-^4)',maxdim,zbottom,ztop);
    make_plot('4',top_wgrid_wp3_dp,top_wgrid_z,'-+ ','','Dissipation (m^3 s^-^4)',maxdim,zbottom,ztop);
    make_plot('5',top_wgrid_wp3_bp,top_wgrid_z,'-rp','','Buoyancy (m^3 s^-^4)',maxdim,zbottom,ztop);
%    title([simname1]);
else
    multi_var = [top_mgrid_wp3_sum_eqn; top_mgrid_wp3_ta_eqn; top_mgrid_wp3_tp_eqn; top_mgrid_wp3_tp'];
    multi_grid = [top_mgrid_z; top_mgrid_z; top_mgrid_z; top_mgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp '];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection(m^3 s^-^4)',maxdim,zbottom,ztop);

    multi_var = [top_mgrid_wp3_pr2'; top_mgrid_wp3_pr1'; top_mgrid_wp3_pr'];
    multi_grid = [top_mgrid_z; top_mgrid_z; top_mgrid_z;];
    multi_marker = ['--k'; ':k '; '-gx'];
    multi_legend = ['pr2';'pr1'];
    make_multi_plot('3',multi_var,multi_grid,multi_marker,multi_legend,'NE',simname1,'Pressure(m^3 s^-^4)',maxdim,zbottom,ztop);

    multi_var = [ top_mgrid_wp3_cl'; top_mgrid_wp3_dp1'];
    multi_grid = [top_mgrid_z; top_mgrid_z];
    multi_marker = ['--k'; ':k '];
    multi_legend = ['cl ';'dp1'];
    make_multi_plot('4',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Dissipation(m^3 s^-^4)',maxdim,zbottom,ztop);

    make_plot('4',top_mgrid_wp3_dp,top_mgrid_z,'-+ ','','Dissipation (m^3 s^-^4)',maxdim,zbottom,ztop);
    make_plot('5',top_mgrid_wp3_bp,top_mgrid_z,'-rp','','Buoyancy (m^3 s^-^4)',maxdim,zbottom,ztop);
%    title([simname1]);
end

if (type2 == 1)
    multi_var = [bot_wgrid_wp3_tp_eqn+bot_wgrid_wp3_ta_eqn; bot_wgrid_wp3_ta_eqn; bot_wgrid_wp3_tp_eqn; bot_wgrid_wp3_tp'];
    multi_grid = [bot_wgrid_z; bot_wgrid_z; bot_wgrid_z; bot_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum'; 'ta ';'tp '];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection(m^3 s^-^4)',maxdim,zbottom,ztop);
    make_plot('8',bot_wgrid_wp3_pr,bot_wgrid_z,'-gx',simname2,'Pressure (m^3 s^-^4)',maxdim,zbottom,ztop);
    make_plot('9',bot_wgrid_wp3_dp,bot_wgrid_z,'-+ ','','Dissipation (m^3 s^-^4)',maxdim,zbottom,ztop);
    make_plot('10',bot_wgrid_wp3_bp,bot_wgrid_z,'-rp','','Buoyancy (m^3 s^-^4)',maxdim,zbottom,ztop);
%    title([simname2]);
else
    multi_var = [bot_mgrid_wp3_sum_eqn; bot_mgrid_wp3_ta_eqn; bot_mgrid_wp3_tp_eqn; bot_mgrid_wp3_tp'];
    multi_grid = [bot_mgrid_z; bot_mgrid_z; bot_mgrid_z; bot_mgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp '];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection(m^3 s^-^4)',maxdim,zbottom,ztop);

    multi_var = [bot_mgrid_wp3_pr2'; bot_mgrid_wp3_pr1'; bot_mgrid_wp3_pr'];
    multi_grid = [bot_mgrid_z; bot_mgrid_z; bot_mgrid_z;];
    multi_marker = ['--k'; ':k '; '-gx'];
    multi_legend = ['pr2';'pr1'];
    make_multi_plot('8',multi_var,multi_grid,multi_marker,multi_legend,'NE',simname2,'Pressure(m^3 s^-^4)',maxdim,zbottom,ztop);

    multi_var = [ bot_mgrid_wp3_cl'; bot_mgrid_wp3_dp1'];
    multi_grid = [bot_mgrid_z; bot_mgrid_z];
    multi_marker = ['--k'; ':k '];
    multi_legend = ['cl ';'dp1'];
    make_multi_plot('9',multi_var,multi_grid,multi_marker,multi_legend','NE','','Dissipation(m^3 s^-^4)',maxdim,zbottom,ztop);

    make_plot('9',bot_mgrid_wp3_dp,bot_mgrid_z,'-+ ','','Dissipation (m^3 s^-^4)',maxdim,zbottom,ztop);
    make_plot('10',bot_mgrid_wp3_bp,bot_mgrid_z,'-rp','','Buoyancy (m^3 s^-^4)',maxdim,zbottom,ztop);
%    title([simname2]);
end

print -depsc wp3;
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT WPTHLP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

if (type1 == 1)
    top_wpthlp = top_mgrid_wpthlp;
    top_grid = top_mgrid_z;
else
    top_wpthlp = top_wgrid_wpthlp;
    top_grid = top_wgrid_z;
end    

if (type2 == 1)
    bot_wpthlp = bot_mgrid_wpthlp;
    bot_grid = bot_mgrid_z;
else
    bot_wpthlp = bot_wgrid_wpthlp;
    bot_grid = bot_wgrid_z;
end

maxdim = 1.1*max([abs(max(top_wpthlp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wpthlp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
                  abs(max(bot_wpthlp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wpthlp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
                ]);

make_plot('1',top_wpthlp,top_grid,'-k','wpthlp','(K m s^-^1)',maxdim,zbottom,ztop);
ylabel('Height (m)');
make_plot('6',bot_wpthlp,bot_grid,'-k','wpthlp','(K m s^-^1)',maxdim,zbottom,ztop);
ylabel('Height (m)');

if (type1 == 1)
    max1 = ([abs(max(top_mgrid_wpthlp_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wpthlp_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wpthlp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wpthlp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ... 
    abs(max(top_mgrid_wpthlp_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wpthlp_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wpthlp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_wpthlp_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(min(top_mgrid_wpthlp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_wpthlp_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wpthlp_pr(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wpthlp_pr(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wpthlp_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wpthlp_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wpthlp_bp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wpthlp_bp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    ]);
else
    max1 = ([abs(max(top_wgrid_wpthlp_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpthlp_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wpthlp_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpthlp_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ... 
    abs(max(top_wgrid_wpthlp_tp_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpthlp_tp_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wpthlp_sum_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpthlp_sum_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wpthlp_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpthlp_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wpthlp_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpthlp_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wpthlp_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpthlp_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    ]);

end

if (type2 == 1)
    max2 = ([abs(max(bot_mgrid_wpthlp_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wpthlp_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wpthlp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wpthlp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ... 
    abs(max(bot_mgrid_wpthlp_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wpthlp_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wpthlp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_wpthlp_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(min(bot_mgrid_wpthlp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_wpthlp_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wpthlp_pr(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wpthlp_pr(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wpthlp_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wpthlp_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wpthlp_bp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wpthlp_bp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    ]);
else
    max2 = ([abs(max(bot_wgrid_wpthlp_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpthlp_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wpthlp_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpthlp_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ... 
    abs(max(bot_wgrid_wpthlp_tp_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpthlp_tp_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wpthlp_sum_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpthlp_sum_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wpthlp_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpthlp_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wpthlp_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpthlp_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wpthlp_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpthlp_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    ]);

end

maxdim = 1.1*max([ max1 max2 ]);

if (type1 == 1)
    multi_var = [top_mgrid_wpthlp_tp_eqn+top_mgrid_wpthlp_ta_eqn; top_mgrid_wpthlp_ta_eqn; top_mgrid_wpthlp_tp_eqn; top_mgrid_wpthlp_tp];
    multi_grid = [top_mgrid_z; top_mgrid_z; top_mgrid_z; top_mgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp '];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (K m s^-^2)',maxdim,zbottom,ztop);
    make_plot('3',top_mgrid_wpthlp_pr,top_mgrid_z,'-gx',simname1,'Pressure (K m s^-^2)',maxdim,zbottom,ztop);
    make_plot('4',top_mgrid_wpthlp_dp,top_mgrid_z,'-+ ','','Dissipation (K m s^-^2)',maxdim,zbottom,ztop);
    make_plot('5',top_mgrid_wpthlp_bp,top_mgrid_z,'-rp','','Buoyancy (K m s^-^2)',maxdim,zbottom,ztop);
%    title([simname1]);
else
    multi_var = [top_wgrid_wpthlp_sum_eqn; top_wgrid_wpthlp_ta_eqn; top_wgrid_wpthlp_tp_eqn; top_wgrid_wpthlp_tp'];
    multi_grid = [top_wgrid_z; top_wgrid_z; top_wgrid_z; top_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp '];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (K m s^-^2)',maxdim,zbottom,ztop);

    multi_var = [top_wgrid_wpthlp_pr3'; top_wgrid_wpthlp_pr2'; top_wgrid_wpthlp_pr1'; top_wgrid_wpthlp_pr'];
    multi_grid = [top_wgrid_z; top_wgrid_z; top_wgrid_z; top_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-gx'];
    multi_legend = ['pr3';'pr2';'pr1'];
    make_multi_plot('3',multi_var,multi_grid,multi_marker,multi_legend,'NE',simname1,'Pressure(K m s^-^2)',maxdim,zbottom,ztop);

    make_plot('4',top_wgrid_wpthlp_dp,top_wgrid_z,'-+ ','','Dissipation (K m s^-^2)',maxdim,zbottom,ztop);
    make_plot('5',top_wgrid_wpthlp_bp,top_wgrid_z,'-rp','','Buoyancy (K m s^-^2)',maxdim,zbottom,ztop);
%    title([simname1]);
end

if (type2 == 1)
    multi_var = [bot_mgrid_wpthlp_tp_eqn+bot_mgrid_wpthlp_ta_eqn; bot_mgrid_wpthlp_ta_eqn; bot_mgrid_wpthlp_tp_eqn; bot_mgrid_wpthlp_tp];
    multi_grid = [bot_mgrid_z; bot_mgrid_z; bot_mgrid_z; bot_mgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp '];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (K m s^-^2)',maxdim,zbottom,ztop);
    make_plot('8',bot_mgrid_wpthlp_pr,bot_mgrid_z,'-gx',simname2,'Pressure (K m s^-^2)',maxdim,zbottom,ztop);
    make_plot('9',bot_mgrid_wpthlp_dp,bot_mgrid_z,'-+ ','','Dissipation (K m s^-^2)',maxdim,zbottom,ztop);
    make_plot('10',bot_mgrid_wpthlp_bp,bot_mgrid_z,'-rp','','Buoyancy (K m s^-^2)',maxdim,zbottom,ztop);
%    title([simname2]);
else
    multi_var = [bot_wgrid_wpthlp_sum_eqn; bot_wgrid_wpthlp_ta_eqn; bot_wgrid_wpthlp_tp_eqn; bot_wgrid_wpthlp_tp'];
    multi_grid = [bot_wgrid_z; bot_wgrid_z; bot_wgrid_z; bot_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp '];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (K m s^-^2)',maxdim,zbottom,ztop);

    multi_var = [bot_wgrid_wpthlp_pr3'; bot_wgrid_wpthlp_pr2'; bot_wgrid_wpthlp_pr1'; bot_wgrid_wpthlp_pr'];
    multi_grid = [bot_wgrid_z; bot_wgrid_z; bot_wgrid_z; bot_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-gx'];
    multi_legend = ['pr3';'pr2';'pr1'];
    make_multi_plot('8',multi_var,multi_grid,multi_marker,multi_legend,'NE',simname2,'Pressure(K m s^-^2)',maxdim,zbottom,ztop);

    make_plot('9',bot_wgrid_wpthlp_dp,bot_wgrid_z,'-+ ','','Dissipation (K m s^-^2)',maxdim,zbottom,ztop);
    make_plot('10',bot_wgrid_wpthlp_bp,bot_wgrid_z,'-rp','','Buoyancy (K m s^-^2)',maxdim,zbottom,ztop);
%    title([simname2]);
end

print -depsc wpthlp;
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT WPQTP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

if (type1 == 1)
    top_wpqtp = top_mgrid_wpqtp;
    top_grid = top_mgrid_z;
else
    top_wpqtp = top_wgrid_wprtp;
    top_grid = top_wgrid_z;
end    

if (type2 == 1)
    bot_wpqtp = bot_mgrid_wpqtp;
    bot_grid = bot_mgrid_z;
else
    bot_wpqtp = bot_wgrid_wprtp;
    bot_grid = bot_wgrid_z;
end

maxdim = 1.1*max([abs(max(top_wpqtp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wpqtp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
                  abs(max(bot_wpqtp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wpqtp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
                  ]);
make_plot('1',top_wpqtp,top_grid,'-k','wpqtp','(m s^-^1)',maxdim,zbottom,ztop);
ylabel('Height (m)');
make_plot('6',bot_wpqtp,bot_grid,'-k','wpqtp','(m s^-^1)',maxdim,zbottom,ztop);
ylabel('Height (m)');

if (type1 == 1)
    max1 = ([abs(max(top_wgrid_wpqtp_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpqtp_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_mgrid_wpqtp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wpqtp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ... 
    abs(max(top_mgrid_wpqtp_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_wpqtp_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_wpqtp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_wpqtp_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(min(top_mgrid_wpqtp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_wpqtp_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_wgrid_wpqtp_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpqtp_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wpqtp_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpqtp_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wpqtp_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wpqtp_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    ]);

    top_mgrid_wpqtp_ta_eqn(top_wgrid_nz)=0;
    top_mgrid_wpqtp_tp_eqn(top_wgrid_nz)=0;
    zexpanded_les=top_mgrid_z;
    zexpanded_les(top_wgrid_nz)=zexpanded_les(top_wgrid_nz-1)+(zexpanded_les(top_wgrid_nz-1)-zexpanded_les(top_wgrid_nz-2));
else
    max1 = ([abs(max(top_wgrid_wprtp_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wprtp_tp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wprtp_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wprtp_ta_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ... 
    abs(max(top_wgrid_wprtp_tp_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wprtp_tp_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wprtp_sum_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wprtp_sum_eqn(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wprtp_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wprtp_pr(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wprtp_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wprtp_dp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    abs(max(top_wgrid_wprtp_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_wgrid_wprtp_bp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) ...
    ]);
end

if (type2 == 1)
    max2 = ([abs(max(bot_wgrid_wpqtp_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpqtp_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_mgrid_wpqtp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wpqtp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ... 
    abs(max(bot_mgrid_wpqtp_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_wpqtp_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_wpqtp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_wpqtp_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(min(bot_mgrid_wpqtp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_wpqtp_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_wgrid_wpqtp_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpqtp_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wpqtp_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpqtp_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wpqtp_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wpqtp_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    ]);

    bot_mgrid_wpqtp_ta_eqn(bot_wgrid_nz)=0;
    bot_mgrid_wpqtp_tp_eqn(bot_wgrid_nz)=0;
    zexpanded_scm=bot_mgrid_z;
    zexpanded_scm(bot_wgrid_nz)=zexpanded_scm(bot_wgrid_nz-1)+(zexpanded_scm(bot_wgrid_nz-1)-zexpanded_scm(bot_wgrid_nz-2));
else
    max2 = ([abs(max(bot_wgrid_wprtp_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wprtp_tp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wprtp_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wprtp_ta_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ... 
    abs(max(bot_wgrid_wprtp_tp_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wprtp_tp_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wprtp_sum_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wprtp_sum_eqn(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wprtp_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wprtp_pr(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wprtp_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wprtp_dp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    abs(max(bot_wgrid_wprtp_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) abs(min(bot_wgrid_wprtp_bp(bot_wgrid_bottom_z_point:bot_wgrid_top_z_point))) ...
    ]);
end

maxdim = 1.1*max([ max1 max2 ]);

if (type1 == 1)
    multi_var = [top_mgrid_wpqtp_tp_eqn+top_mgrid_wpqtp_ta_eqn; top_mgrid_wpqtp_ta_eqn; top_mgrid_wpqtp_tp_eqn; top_wgrid_wpqtp_tp'];
    multi_grid = [zexpanded_les; zexpanded_les; zexpanded_les; top_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum'; 'ta ';'tp '];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (m s^-^2)',maxdim,zbottom,ztop);
    make_plot('3',top_wgrid_wpqtp_pr,top_wgrid_z,'-gx',simname1,'Pressure (m s^-^2)',maxdim,zbottom,ztop);
    make_plot('4',top_wgrid_wpqtp_dp,top_wgrid_z,'-+ ','','Dissipation (m s^-^2)',maxdim,zbottom,ztop);
    make_plot('5',top_wgrid_wpqtp_bp,top_wgrid_z,'-rp','','Buoyancy (m s^-^2)',maxdim,zbottom,ztop);
%    title([simname1]);
else
    multi_var = [top_wgrid_wprtp_sum_eqn; top_wgrid_wprtp_ta_eqn; top_wgrid_wprtp_tp_eqn; top_wgrid_wprtp_tp'];
    multi_grid = [top_wgrid_z; top_wgrid_z; top_wgrid_z; top_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['ta ';'tp ';'sum'];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (m s^-^2)',maxdim,zbottom,ztop);

    multi_var = [top_wgrid_wprtp_pr3'; top_wgrid_wprtp_pr2'; top_wgrid_wprtp_pr1'; top_wgrid_wprtp_pr'];
    multi_grid = [top_wgrid_z; top_wgrid_z; top_wgrid_z; top_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-gx'];
    multi_legend = ['pr3';'pr2';'pr1'];
    make_multi_plot('3',multi_var,multi_grid,multi_marker,multi_legend,'NE',simname1,'Pressure(m s^-^2)',maxdim,zbottom,ztop);

    make_plot('4',top_wgrid_wprtp_dp,top_wgrid_z,'-+ ','','Dissipation (m s^-^2)',maxdim,zbottom,ztop);
    make_plot('5',top_wgrid_wprtp_bp,top_wgrid_z,'-rp','','Buoyancy (m s^-^2)',maxdim,zbottom,ztop);
%    title([simname1]);
end

if (type2 == 1)
    multi_var = [bot_mgrid_wpqtp_tp_eqn+bot_mgrid_wpqtp_ta_eqn; bot_mgrid_wpqtp_ta_eqn; bot_mgrid_wpqtp_tp_eqn; bot_wgrid_wpqtp_tp'];
    multi_grid = [zexpanded_scm; zexpanded_scm; zexpanded_scm; bot_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum'; 'ta ';'tp '];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (m s^-^2)',maxdim,zbottom,ztop);
    make_plot('8',bot_wgrid_wpqtp_pr,bot_wgrid_z,'-gx',simname2,'Pressure (m s^-^2)',maxdim,zbottom,ztop);
    make_plot('9',bot_wgrid_wpqtp_dp,bot_wgrid_z,'-+ ','','Dissipation (m s^-^2)',maxdim,zbottom,ztop);
    make_plot('10',bot_wgrid_wpqtp_bp,bot_wgrid_z,'-rp','','Buoyancy (m s^-^2)',maxdim,zbottom,ztop);
%    title([simname2]);
else
    multi_var = [bot_wgrid_wprtp_sum_eqn; bot_wgrid_wprtp_ta_eqn; bot_wgrid_wprtp_tp_eqn; bot_wgrid_wprtp_tp'];
    multi_grid = [bot_wgrid_z; bot_wgrid_z; bot_wgrid_z; bot_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['ta ';'tp ';'sum'];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (m s^-^2)',maxdim,zbottom,ztop);

    multi_var = [bot_wgrid_wprtp_pr3'; bot_wgrid_wprtp_pr2'; bot_wgrid_wprtp_pr1'; bot_wgrid_wprtp_pr'];
    multi_grid = [bot_wgrid_z; bot_wgrid_z; bot_wgrid_z; bot_wgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-gx'];
    multi_legend = ['pr3';'pr2';'pr1'];
    make_multi_plot('8',multi_var,multi_grid,multi_marker,multi_legend,'NE',simname2,'Pressure(m s^-^2)',maxdim,zbottom,ztop);

    make_plot('9',bot_wgrid_wprtp_dp,bot_wgrid_z,'-+ ','','Dissipation (m s^-^2)',maxdim,zbottom,ztop);
    make_plot('10',bot_wgrid_wprtp_bp,bot_wgrid_z,'-rp','','Buoyancy (m s^-^2)',maxdim,zbottom,ztop);
%    title([simname2]);
end

print -depsc wpqtp;
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THLP2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

%maxdim = 1.1*max([abs(max(top_mgrid_thlp2)) abs(min(top_mgrid_thlp2)) abs(max(bot_wgrid_thlp2)) abs(min(bot_wgrid_thlp2))]);
graphics_displayed = 0;

if (type1 == 1)
    maxdim = 1.1*max([ abs(max(top_mgrid_thlp2(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_thlp2(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ]);
    make_plot('1',top_mgrid_thlp2,top_mgrid_z,'-k','thlp2','(K^2)',maxdim,zbottom,ztop);
    ylabel('Height (m)');

maxdim = 1.1*max([abs(max(top_mgrid_thlp2_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_thlp2_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_thlp2_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_thlp2_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ... 
    abs(max(top_mgrid_thlp2_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_thlp2_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_thlp2_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_thlp2_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(min(top_mgrid_thlp2_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_thlp2_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_thlp2_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_thlp2_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_thlp2_rad(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_thlp2_rad(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    abs(max(top_mgrid_thlp2_dp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_thlp2_dp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
    % this plots the SCM-- we don't want that yet
%    abs(max(bot_mgrid_thlp2_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(max(bot_mgrid_thlp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(top_mgrid_thlp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ... 
%    abs(max(bot_mgrid_thlp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(max(bot_mgrid_thlp2_sum_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_sum_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    ]);

multi_var = [top_mgrid_thlp2_tp_eqn+top_mgrid_thlp2_ta_eqn; top_mgrid_thlp2_ta_eqn; top_mgrid_thlp2_tp_eqn; top_mgrid_thlp2_tp];
multi_grid = [top_mgrid_z; top_mgrid_z; top_mgrid_z; top_mgrid_z];
multi_marker = ['-k '; '--k'; ':k '; '-bh'];
multi_legend = ['sum';'ta ';'tp '];
make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (K^2 s^-^1)',maxdim,zbottom,ztop);

multi_var = [top_mgrid_thlp2_dp_eqn; top_mgrid_thlp2_dp];
multi_grid = [top_mgrid_z; top_mgrid_z];
multi_marker = ['-k '; '-+m'];
multi_legend = ['corrected'];
make_multi_plot('4',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Dissipation (K^2 s^-^1)',maxdim,zbottom,ztop);

make_plot('4',top_mgrid_thlp2_dp,top_mgrid_z,'-+ ',simname1,'Dissipation (K^2 s^-^1)',maxdim,zbottom,ztop);
make_plot('5',top_mgrid_thlp2_rad,top_mgrid_z,'-rp','','Radiation (K^2 s^-^1)',maxdim,zbottom,ztop);

%title([simname1]);
graphics_displayed = graphics_displayed + 1;
end

if (type2 == 1)
    maxdim = 1.1*max([ abs(max(bot_mgrid_thlp2)) abs(min(bot_mgrid_thlp2)) ]);
    make_plot('6',bot_mgrid_thlp2,bot_mgrid_z,'-k','thlp2','(K^2)',maxdim,zbottom,ztop);
    ylabel('Height (m)');

maxdim = 1.1*max([abs(max(bot_mgrid_thlp2_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_thlp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ... 
    abs(max(bot_mgrid_thlp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_thlp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_thlp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(min(bot_mgrid_thlp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_thlp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_thlp2_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_thlp2_rad(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_rad(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    abs(max(bot_mgrid_thlp2_dp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_dp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    % this plots the SCM-- we don't want that yet
%    abs(max(bot_mgrid_thlp2_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(max(bot_mgrid_thlp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ... 
%    abs(max(bot_mgrid_thlp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(max(bot_mgrid_thlp2_sum_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_thlp2_sum_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
    ]);

multi_var = [bot_mgrid_thlp2_tp_eqn+bot_mgrid_thlp2_ta_eqn; bot_mgrid_thlp2_ta_eqn; bot_mgrid_thlp2_tp_eqn; bot_mgrid_thlp2_tp];
multi_grid = [bot_mgrid_z; bot_mgrid_z; bot_mgrid_z; bot_mgrid_z];
multi_marker = ['-k '; '--k'; ':k '; '-bh'];
multi_legend = ['sum';'ta ';'tp '];
make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (K^2 s^-^1)',maxdim,zbottom,ztop);

multi_var = [bot_mgrid_thlp2_dp_eqn; bot_mgrid_thlp2_dp];
multi_grid = [bot_mgrid_z; bot_mgrid_z];
multi_marker = ['-k '; '-+m'];
multi_legend = ['corrected'];
make_multi_plot('9',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Dissipation (K^2 s^-^1)',maxdim,zbottom,ztop);

make_plot('9',bot_mgrid_thlp2_dp,bot_mgrid_z,'-+ ',simname2,'Dissipation (K^2 s^-^1)',maxdim,zbottom,ztop);
make_plot('10',bot_mgrid_thlp2_rad,bot_mgrid_z,'-rp','','Radiation (K^2 s^-^1)',maxdim,zbottom,ztop);

%title([simname2]);
graphics_displayed = graphics_displayed + 1;
end

if (graphics_displayed > 0)
    print -depsc thlp2;
    keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT QTP2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
graphics_displayed = 0;

if (type1 == 1)
    maxdim = 1.1*max([abs(max(top_mgrid_qtp2(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_mgrid_qtp2(top_wgrid_bottom_z_point:top_wgrid_top_z_point)))]);
    make_plot('1',top_mgrid_qtp2,top_mgrid_z,'-k','qtp2','(dimensionless)',maxdim,zbottom,ztop);
    ylabel('Height (m)');

    maxdim = 1.1*max([abs(max(top_mgrid_qtp2_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtp2_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(max(top_mgrid_qtp2_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtp2_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ... 
        abs(max(top_mgrid_qtp2_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtp2_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(max(top_mgrid_qtp2_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_qtp2_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(min(top_mgrid_qtp2_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_qtp2_tp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(max(top_mgrid_qtp2_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtp2_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ... 
        abs(max(top_mgrid_qtp2_dp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtp2_dp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        ]);

    multi_var = [top_mgrid_qtp2_tp_eqn+top_mgrid_qtp2_ta_eqn; top_mgrid_qtp2_ta_eqn; top_mgrid_qtp2_tp_eqn; top_mgrid_qtp2_tp];
    multi_grid = [top_mgrid_z; top_mgrid_z; top_mgrid_z; top_mgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp '];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (s^-^1)',maxdim,zbottom,ztop);

    multi_var = [top_mgrid_qtp2_dp_eqn; top_mgrid_qtp2_dp];
    multi_grid = [top_mgrid_z; top_mgrid_z; top_mgrid_z; top_mgrid_z];
    multi_marker = ['-k '; '-+m'];
    multi_legend = ['corrected'];
    make_multi_plot('4',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Dissipation (s^-^1)',maxdim,zbottom,ztop);
    make_plot('4',top_mgrid_qtp2_dp,top_mgrid_z,'-+ ',simname1,'Dissipation (s^-^1)',maxdim,zbottom,ztop);

%    title([simname1]);
    graphics_displayed = graphics_displayed + 1;
end

if (type2 == 1)
    maxdim = 1.1*max([abs(max(bot_mgrid_qtp2)) abs(min(bot_mgrid_qtp2))]);
    make_plot('6',bot_mgrid_qtp2,bot_mgrid_z,'-k','qtp2','(dimensionless)',maxdim,zbottom,ztop);
    ylabel('Height (m)');

    maxdim = 1.1*max([abs(max(bot_mgrid_qtp2_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtp2_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(max(bot_mgrid_qtp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ... 
        abs(max(bot_mgrid_qtp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(max(bot_mgrid_qtp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_qtp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(min(bot_mgrid_qtp2_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_qtp2_tp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(max(bot_mgrid_qtp2_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtp2_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ... 
        abs(max(bot_mgrid_qtp2_dp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtp2_dp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        ]);

    multi_var = [bot_mgrid_qtp2_tp_eqn+bot_mgrid_qtp2_ta_eqn; bot_mgrid_qtp2_ta_eqn; bot_mgrid_qtp2_tp_eqn; bot_mgrid_qtp2_tp];
    multi_grid = [bot_mgrid_z; bot_mgrid_z; bot_mgrid_z; bot_mgrid_z];
    multi_marker = ['-k '; '--k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp '];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (s^-^1)',maxdim,zbottom,ztop);

    multi_var = [bot_mgrid_qtp2_dp_eqn; bot_mgrid_qtp2_dp];
    multi_grid = [bot_mgrid_z; bot_mgrid_z; bot_mgrid_z; bot_mgrid_z];
    multi_marker = ['-k '; '-+m'];
    multi_legend = ['corrected'];
    make_multi_plot('9',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Dissipation (s^-^1)',maxdim,zbottom,ztop);
    make_plot('9',bot_mgrid_qtp2_dp,bot_mgrid_z,'-+ ',simname2,'Dissipation (s^-^1)',maxdim,zbottom,ztop);

%    title([simname2]);
    graphics_displayed = graphics_displayed + 1;
end

if (graphics_displayed > 0)
    print -depsc qtp2;
    keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT QTPTHLP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

graphics_displayed = 0;

if (type1 == 1)
    maxdim = 1.1*max([abs(max(top_mgrid_qtpthlp(top_wgrid_bottom_z_point:top_wgrid_top_z_point))) abs(min(top_mgrid_qtpthlp(top_wgrid_bottom_z_point:top_wgrid_top_z_point)))]);
    make_plot('1',top_mgrid_qtpthlp,top_mgrid_z,'-k','qtpthlp','(K)',maxdim,zbottom,ztop);
    ylabel('Height (m)');

    maxdim = 1.1*max([abs(max(top_mgrid_qtpthlp_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtpthlp_tp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(max(top_mgrid_qtpthlp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtpthlp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ... 
        abs(max(top_mgrid_qtpthlp_tp1_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtpthlp_tp1_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(max(top_mgrid_qtpthlp_tp2_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtpthlp_tp2_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(max(top_mgrid_qtpthlp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_qtpthlp_tp1_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_qtpthlp_tp2_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(min(top_mgrid_qtpthlp_ta_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_qtpthlp_tp1_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point) + top_mgrid_qtpthlp_tp2_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(max(top_mgrid_qtpthlp_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtpthlp_dp(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(max(top_mgrid_qtpthlp_rad(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtpthlp_rad(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        abs(max(top_mgrid_qtpthlp_dp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) abs(min(top_mgrid_qtpthlp_dp_eqn(top_mgrid_bottom_z_point:top_mgrid_top_z_point))) ...
        ]);

    multi_var = [top_mgrid_qtpthlp_tp1_eqn+top_mgrid_qtpthlp_tp2_eqn+top_mgrid_qtpthlp_ta_eqn; top_mgrid_qtpthlp_ta_eqn; top_mgrid_qtpthlp_tp2_eqn; ...
        top_mgrid_qtpthlp_tp1_eqn; top_mgrid_qtpthlp_tp];
    multi_grid = [top_mgrid_z; top_mgrid_z; top_mgrid_z; top_mgrid_z; top_mgrid_z];
    multi_marker = ['-k '; '--k'; '-.k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp2';'tp1'];
    make_multi_plot('2',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (K s^-^1)',maxdim,zbottom,ztop);

    multi_var = [top_mgrid_qtpthlp_dp_eqn; top_mgrid_qtpthlp_dp];
    multi_grid = [top_mgrid_z; top_mgrid_z];
    multi_marker = ['-k '; '-+m'];
    multi_legend = ['corrected'];
    make_multi_plot('4',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Dissipation (K s^-^1)',maxdim,zbottom,ztop);

    make_plot('4',top_mgrid_qtpthlp_dp,top_mgrid_z,'-+ ',simname1,'Dissipation (K s^-^1)',maxdim,zbottom,ztop);
    make_plot('5',top_mgrid_qtpthlp_rad,top_mgrid_z,'-rp','','Radiation (K s^-^1)',maxdim,zbottom,ztop);

%    title([simname1]);
    graphics_displayed = graphics_displayed + 1;
end

if (type2 == 1)
    maxdim = 1.1*max([abs(max(bot_mgrid_qtpthlp)) abs(min(bot_mgrid_qtpthlp))]);
    make_plot('6',bot_mgrid_qtpthlp,bot_mgrid_z,'-k','qtpthlp','(K)',maxdim,zbottom,ztop);
    ylabel('Height (m)');

    maxdim = 1.1*max([abs(max(bot_mgrid_qtpthlp_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtpthlp_tp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(max(bot_mgrid_qtpthlp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtpthlp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ... 
        abs(max(bot_mgrid_qtpthlp_tp1_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtpthlp_tp1_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(max(bot_mgrid_qtpthlp_tp2_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtpthlp_tp2_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(max(bot_mgrid_qtpthlp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_qtpthlp_tp1_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_qtpthlp_tp2_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(min(bot_mgrid_qtpthlp_ta_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_qtpthlp_tp1_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point) + bot_mgrid_qtpthlp_tp2_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(max(bot_mgrid_qtpthlp_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtpthlp_dp(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(max(bot_mgrid_qtpthlp_rad(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtpthlp_rad(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        abs(max(bot_mgrid_qtpthlp_dp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) abs(min(bot_mgrid_qtpthlp_dp_eqn(bot_mgrid_bottom_z_point:bot_mgrid_top_z_point))) ...
        ]);

    multi_var = [bot_mgrid_qtpthlp_tp1_eqn+bot_mgrid_qtpthlp_tp2_eqn+bot_mgrid_qtpthlp_ta_eqn; bot_mgrid_qtpthlp_ta_eqn; bot_mgrid_qtpthlp_tp2_eqn; ...
        bot_mgrid_qtpthlp_tp1_eqn; bot_mgrid_qtpthlp_tp];
    multi_grid = [bot_mgrid_z; bot_mgrid_z; bot_mgrid_z; bot_mgrid_z; bot_mgrid_z];
    multi_marker = ['-k '; '--k'; '-.k'; ':k '; '-bh'];
    multi_legend = ['sum';'ta ';'tp2';'tp1'];
    make_multi_plot('7',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Advection (K s^-^1)',maxdim,zbottom,ztop);

    multi_var = [bot_mgrid_qtpthlp_dp_eqn; bot_mgrid_qtpthlp_dp];
    multi_grid = [bot_mgrid_z; bot_mgrid_z];
    multi_marker = ['-k '; '-+m'];
    multi_legend = ['corrected'];
    make_multi_plot('9',multi_var,multi_grid,multi_marker,multi_legend,'NE','','Dissipation (K s^-^1)',maxdim,zbottom,ztop);

    make_plot('9',bot_mgrid_qtpthlp_dp,bot_mgrid_z,'-+ ',simname2,'Dissipation (K s^-^1)',maxdim,zbottom,ztop);
    make_plot('10',bot_mgrid_qtpthlp_rad,bot_mgrid_z,'-rp','','Radiation (K s^-^1)',maxdim,zbottom,ztop);

%    title([simname2]);
    graphics_displayed = graphics_displayed + 1;
end

if (graphics_displayed > 0)
    print -depsc qtpthlp;
    keyboard
end