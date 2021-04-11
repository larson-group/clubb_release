## analyzeSensMatrix:
### A tool for analyzing sensitivities to changes in parameter values.

Suppose that a numerical model is well tuned.  Then, if a new structural model change is introduced,
the model will need to be retuned.  This can be time consuming.  To expedite the work,
here we provide a script that 
1. analyzes the sensitivity of model outputted metrics to tunable parameters, and
2. recommends a combination of tuning parameter values that brings the model closer to observations.

The script requires that a human user choose metrics and tuning parameters.  
The intended workflow is as follows:
1.  Start with an untuned numerical model.  Run a single simulation with parameters 
at their default values.  Output the metrics and parameter values to a "default" netcdf file.
2.  For each tunable parameter, run a single simulation with the parameter value changed,
but all other parameters at their default values.  For each sensitivity simulation, 
output to a netcdf file the metrics and parameter values.
3.  Based on this set of netcdf input files, observed values of the metrics, 
and user-chosen settings, call the function analyzeSensMatrix.
5.  If the function's recommended parameter values seem likely to better fit 
the observations than the default values, then run a single final simulation 
with the recommended parameter values.
6.  But if the parameters need to be better constrained, then repeat the analysis with more metrics.  
Or if the parameters are unable to match the obs, then find more effective parameters.

The scripts are written entirely in python.  The results may be plotted on a web-based
dashboard created using python dash.  Unit tests may be performed by running pytest.

This folder contains several files:

*analyze_sensitivity_matrix.py* contains the main workhorse function, analyzeSensMatrix,
plus some sample input (that won't work out of the box).

*test_analyzeSensMatrix.py* contains some unit tests that call analyzeSensMatrix 
and are intended to work out of the box.

*pytest.ini* is an initialization file for the pytest unit testing framework.

*sens_matrix_dashboard.py* creates a dashboard of diagnostic plots that can be viewed 
on a web browser.
