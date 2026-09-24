## Utilities Directory Contents ## 

The [CLUBBStandardsCheck.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/CLUBBStandardsCheck.py) python script can be used to check Fortran 
source files to determine if they follow certain good software engineering 
practices which are meant to be enforced for CLUBB source files. 

[check_for_errors.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/check_for_errors.py) checks CLUBB source code files for common mistakes or
CLUBB standards violations. It checks for uninitialized output variables,
magic flags, and magic numbers. An output variable is only initialized if it
is set in every part of an if or select case statement or if it is set outside
of these statements. A magic number is any number that appears
in a subroutine or function call or appears in an equation on the right-hand
side of an equals sign. Integers from -6 to +6 are not considered magic numbers,
and nothing in an if statement is considered a magic number. If the argument -w
or --show-warnings is provided, warnings will be printed when a variable is set
by a subroutine or function call. If no files are provided in the argument list,
the script will check all .F90 files under the src directory (recursively).
USAGE: "python3 check_for_errors.py [-w or --show-warnings] [<filename>.F90 ...]"

[print_tunable_parameters_table.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/print_tunable_parameters_table.py) prints out a table of the tunable parameters
with their values in markdown format. The output can be copied and pasted directly
into a github comment and will be formatted as a table there. If multiple files
are listed as arguments, then the output will have a column for each file listing
the values of the parameters if they are present or an indicator if commented out.
USAGE: "python print_tunable_parameters_table.py [options]
               <list of tunable_parameters.in files>"

Other files in this directory:

[__init__.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/__init__.py) makes utilities/ a Python package, so other code can use
"from utilities.<module> import ...".

[benchmark_converter.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/benchmark_converter.py) converts SAM/COAMPS LES benchmark NetCDF files into a
small NetCDF that uses CLUBB stats names and CLUBB-style time. Used by the
tuner/loss driver (through [create_case_namelist.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/create_case_namelist.py)) and the Dash plot and
tune tabs.

[clubb_settings_validation.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/clubb_settings_validation.py) is a hand-maintained Python copy of CLUBB's
parameter bounds and flag/parameter compatibility rules, so Dash and the tuner
can check settings without a compiled CLUBB. [tests/run_clubb_settings_validation_test.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/tests/run_clubb_settings_validation_test.py)
checks it against the compiled Fortran.

[convert_acc_to_omp.bash](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/convert_acc_to_omp.bash) adds OpenMP directives next to the OpenACC directives
in [src/](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/src), using Intel's OpenACC-to-OpenMP migration tool (cloned on first use).
Used by the nvhpc GPU-vs-CPU Jenkins test.

[convert_to_async.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/convert_to_async.py) adds async(1)/wait clauses to selected OpenACC directives
in [src/](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/src). Used by the GPU async output-matching Jenkins test.

[create_case_namelist.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/create_case_namelist.py) combines a case's model file, tunable parameters,
flags, SILHS settings, and stats list into the single <case>.in namelist that
CLUBB reads, applying any overrides. Called by [run_scripts/run_scm.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_scm.py),
[run_scm_loss.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_scm_loss.py), and the tuner.

[create_multi_col_params.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/create_multi_col_params.py) writes a multicolumn tunable-parameter file
(duplicated, perturbed, mirrored, or hypergrid columns). Called by
[create_case_namelist.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/create_case_namelist.py), the tuner, and several GPU Jenkins tests.

[flag_sets.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/flag_sets.py) reads a flag-set JSON file (see [input/flag_sets/](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/input/flag_sets)) and turns each
set into a namelist override string. Shared by
[run_scripts/run_clubb_w_varying_flags.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_clubb_w_varying_flags.py) and [tests/run_jax_vs_fortran_cases.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/tests/run_jax_vs_fortran_cases.py).

[les_chi_moments.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/les_chi_moments.py) derives chi (extended liquid water) moments from LES
fields using CLUBB's thermodynamic constants and saturation formula. Used by
[benchmark_converter.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/benchmark_converter.py) and [sam_3d_reference.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/sam_3d_reference.py).

[loss_metrics.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/loss_metrics.py) is a Python copy of the profile loss and Taylor metrics
computed in [src/clubb_loss_driver.F90](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/src/clubb_loss_driver.F90). Used by the Dash plot tab and the
loss-driver tests.

[output_paths.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/output_paths.py) turns a user-given output name into a directory under output/
(and rejects ".."). Shared by [run_scm.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_scm.py), [run_scm_loss.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_scm_loss.py), and Dash.

[sam_3d_reference.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/sam_3d_reference.py) reads raw 3D SAM snapshots from the shared benchmark-run
directory for Dash diagnostics and experimental PDF analyses.

[save_tunable_config.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/save_tunable_config.py) copies a config directory under
[input/parameter_and_flag_configs/](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/input/parameter_and_flag_configs) and applies Run/Tune-tab overrides to the
copy. Used only by Dash (it imports [dash_app](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/dash_app)).

[split_stats_to_legacy.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/split_stats_to_legacy.py) splits one combined stats NetCDF file into
legacy-style per-grid files (e.g. <case>_zt.nc, <case>_zm.nc). Nothing in the
repository calls it.

[time_clubb.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/time_clubb.py) benchmarks CLUBB by running groups of concurrent [run_scm.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_scm.py)
processes and writes a timing-profile directory. It is also the backend of the
Dash Profile tab.

[timing_profiles.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/timing_profiles.py) reads and writes the timing-profile directories produced by
[time_clubb.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/time_clubb.py). Used by [time_clubb.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/time_clubb.py) and the Dash Profile tab.
