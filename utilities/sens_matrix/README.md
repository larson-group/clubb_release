## sens_matrix_dashboard.py:
### A tool for analyzing sensitivities of a global atmospheric model to changes in parameter values.

Suppose that a numerical model is well tuned.  Then, if a new structural model change is introduced,
the model will need to be retuned.  This can be time consuming.  To expedite the work,
here we provide a set of scripts that 
1. analyze the sensitivity of model outputted metrics to tunable parameters, and
2. recommend a combination of tuning parameter values that brings the model closer to observations.

#### Software requirements:

RegTune requires python and dash.  To install dash, run

`pip install dash`

#### To use the tuning tool:

1.  Start with an untuned global numerical model.  Run a single global simulation with parameters 
at their default values.  Output the metrics and parameter values to a "default" netcdf file.
2.  Choose which parameters you'd like to tune.  For each parameter to tune, run two simulations, 
one with the parameter value increased, and one with it decreased, 
but keep all other parameters at their default values.  
For each sensitivity simulation, output to a netcdf file the metrics and parameter values.
3.  Specify regional metrics (e.g. SWCF near Hawaii), their weights, and parameters to tune
in set_up_dashboard_inputs.py.
5.  Based on this set of netcdf input files, observed values of the metrics, 
and user-chosen settings, run `$> python3 sens_matrix_dashboard.py`.  
6.  If the tuner finds a better fit for the metrics of interest, then run a global simulation 
with the recommended parameter values.

The scripts are written entirely in python.  The results may be plotted on a web-based
dashboard created using python dash.  Unit tests may be performed by running pytest.

#### This folder contains several files:

*set_up_dashboard_inputs.py* allows the user to state which regional metrics to match,
how strongly to weight each metric, what parameters to tune, and what observed values
to match for each region.

*sens_matrix_dashboard.py* is a driver script that creates a plotly dashboard of diagnostic plots 
that can be viewed on a web browser.

*analyze_sensitivity_matrix.py* contains the main workhorse function, analyzeSensMatrix,
plus some sample input (that won't work out of the box).

*test_analyzeSensMatrix.py* contains some unit tests that call analyzeSensMatrix 
and are intended to work out of the box.

*pytest.ini* is an initialization file for the pytest unit testing framework.
