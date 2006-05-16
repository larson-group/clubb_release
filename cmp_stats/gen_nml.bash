#!/bin/bash
# $Id: gen_nml.bash,v 1.4 2006-05-16 18:39:47 dschanen Exp $
# I'm too lazy to modify these one by one

MODEL_ALL=(arm atex bomex dycoms2_rf01 dycoms2_rf02_do dycoms2_rf02_ds\
 dycoms2_rf02_nd dycoms2_rf02_so fire nov11_altocu wangara)

for RUN_CASE in "${MODEL_ALL[@]}"; do
	STATS_IN=$RUN_CASE'_stats.in'

	echo "&statsnl" > $STATS_IN
	echo "lstats = t" >> $STATS_IN
	echo "stats_fmt = 'netcdf'" >> $STATS_IN
	echo "stats_tsamp = 60.0" >> $STATS_IN
	echo "stats_tout  = 3600.0" >> $STATS_IN
	echo "fname_zt = '$RUN_CASE""_zt'," >> $STATS_IN
	echo "vars_zt =" >> $STATS_IN
	echo "'thlm', 'thvm', 'rtm', 'rcm', 'um', 'vm', 'wmt', 'ug', 'vg'," >> $STATS_IN
	echo "'cf', 'p', 'Lscale', 'thlm_forcing', 'rtm_forcing'," >> $STATS_IN
	echo "'wp3', 'wpthlp2', 'wp2thlp', 'wprtp2', 'wp2rtp'," >> $STATS_IN
	echo "'lup', 'ldown', 'taut', 'kht', 'wp2thvp', 'wp2rcp', 'wprtpthlp'," >> $STATS_IN
	echo "'sc', 'rhot', 'Ncm'," >> $STATS_IN
	echo "'rsm', 'rrm', 'wp2zt'," >> $STATS_IN
	echo "''" >> $STATS_IN

	echo "fname_zm = '$RUN_CASE""_zm'," >> $STATS_IN
	echo "vars_zm =" >> $STATS_IN
	echo "'wp2', 'rtp2', 'thlp2', 'rtpthlp', 'wprtp', 'wpthlp', 'wp4'," >> $STATS_IN
	echo "'wpthvp', 'rtpthvp', 'thlpthvp', 'taum ', 'khm', 'wprcp'," >> $STATS_IN
	echo "'thlprcp', 'rtprcp', 'upwp', 'vpwp', 'rhom', 'sc', 'em'," >> $STATS_IN
	echo "'shear', 'Frad', 'Fprec', 'Fcsed'," >> $STATS_IN
	echo "''" >> $STATS_IN
	
	echo "fname_sfc = '$RUN_CASE""_sfc'," >> $STATS_IN
	echo "vars_sfc =" >> $STATS_IN
	echo "'ustar', " >> $STATS_IN
	echo "''" >> $STATS_IN
	
	echo "/" >> $STATS_IN
done
