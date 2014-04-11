#!/bin/bash
echo "******browser.bash started******" >> ~/nightly_plots.log

#nielsenb - We need to clean up the directory for the code browser too
rm -f ~/hoc_v2.2_tuner/hocbrowser/call_from/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/call_to/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/html_code/bugs_rad/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/html_code/latin_hc/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/html_code/main/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/html_code/COAMPS_Microphysics/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/html_code/HOC_parameterization/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/html_code/LH/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/html_code/src/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/html_code/BUGSrad/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/html_code/GCSS/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/ind/*.html
rm -f ~/hoc_v2.2_tuner/hocbrowser/*.html
#end nielsenb's changes

cd ~/f90tohtml/nwp_codes/ && \
perl hoc_prepare.pl && \
../f90tohtml hoc.f2h && \
cd ~/hoc_v2.2_tuner/ && \

#nielsenb - There is no "main" folder, so instead make a Main_index.html that
#mirrors the src index. This prevents 404 errors on main browser page
cp ~/hoc_v2.2_tuner/hocbrowser/ind/src_index.html ~/hoc_v2.2_tuner/hocbrowser/ind/main_index.html
#end nielsenb's changes

#cvs ci -m"Automated Generation" hocbrowser/ && \
echo "******browser.bash ended successfully******" >> ~/nightly_plots.log
