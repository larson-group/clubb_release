# Quickstart
Pyplotgen takes parameters in the form `python3 ./PyPlotGen.py [OPTIONS]`
Pyplotgen only supports input in the netcdf (.nc) format.

## Valid options

| *Option Flag* | *Description* |
| --- | --- |
| -r --replace | Overwrite the output folder if it already exists |
| -c --clubb [FOLDER PATHNAME(S)] | Adds lines from CLUBB model output to plot. To plot lines in the specified folder for a given case, update the `'clubb_file': None,` parameter in the case's definition (`config/Case_Definitions.py`)for the desired case to contain the location of the clubb file. | 
| -e --e3sm [FOLDER PATHNAME(S)] | Adds lines from E3SM model output to plot. The filenames in the folder must either match the casename intended to run it, e.g. `bomex.nc` (using the lowercase version of the casename as defined in the `Case_definitions.py` file), or rewrite the `'e3sm_file': None,` parameter for the desired case to contain the location of the e3sm file.| 
| -s --sam [FOLDER PATHNAME(S)] | Adds lines from SAM model output to plot. To plot lines in the specified folder for a given case, update the `'sam_file': None,` parameter in the case's definition (`config/Case_Definitions.py`)for the desired case to contain the location of the SAM file. For example, adding SAM output for the ARM case may look like `'sam_file': sam_output_root+ "/GCSSARM_96x96x110_67m_40m_1s.nc",`| 
| -w --wrf [FOLDER PATHNAME(S)] | Adds lines from WRF model output to plot. To plot lines in the specified folder for a given case, update the `'wrf_file': None,` parameter in the case's definition (`config/Case_Definitions.py`)for the desired case to contain the location of the WRF file.| 
| -g --plot-golaz-best | Plots clubb r408 'best ever' plots |
| -l --les | Overplot LES output data. The nc files for les output can be overwritten by changing the directory listed for a given case inside the `Case_Definitions.py` file. |
| -d --plot-hoc-2005 | Plots HOC benchmark output |
| -a --all-best | Identical to -lgd. Adds golaz best, les, and hoc benchmark output to the plots |
| --thin | Plot lines with a thin width |
| -b --plot-budgets | Includes budget panels in output |
| --diff [FOLDER PATHNAME] | (Experimental) Plots the difference between the input folder and the folder specified after --diff instead of plotting a regular profile |
| --no-legends | Panels are drawn without a line legend |
| -o --output | Manually specify an output folder. If not specified, will automatically output to `pyplotgen/output` |
| --show-alphabetic-id | Adds an alphanumeric ID to each plot on a perc-case basis (e.g. the first plot will be labeled "a")
| --nightly | Apply special parameters only relevant when running as part of a nightly test. This is currently limited to disabling case output if not all models have data for a given case. E.g. this prevents wrf plots from including cases that only have clubb plots and no wrf plots. Do not plot this with clubb-only plots, just plot clubb normally for clubb nightly tests.
| --disable-multithreading | Turns off multithreading support. Useful for debugging as it ensures text is printed sequentially. |
| --hq --high-quality | Outputs higher resolution images. The dpi used for hi resolution images can be customized in Style_definitions.py |

## Installing Dependencies
To install the dependencies necessary for PyPlotgen to run, run the command
```
pip3 install -r requirements.txt
```
from within the PyPlotGen directory

## Example Run Commands

To plot clubb output located in `/home/USERNAME/clubb_nc_output` and save the generated plots in `/home/USERNAME/clubb/postprocessing/pyplotgen/output`, run this command:

`python3 ./pyplotgen.py -c /home/USERNAME/clubb_nc_output -o /home/USERNAME/clubb/postprocessing/pyplotgen/output`

If, in addition, you'd like to overplot LES lines and also separately plot CLUBB budgets, run this command:

`python3 ./pyplotgen.py --plot-budgets -l -c /home/USERNAME/clubb_nc_output -o /home/USERNAME/clubb/postprocessing/pyplotgen/output`

# Advanced Usage
## Reference Documentation
Reference documentation is available in the `pyplotgen/docs/html` folder. Open `pyplotgen/docs/html/index.html` in a web browser for easy viewing. More information is available in the `pyplotgen/docs/README.md` file.

## Running subsets of cases
If someone wants to run only a few cases, reguardless of how many datasets were outputted to a folder, they can do so by editing the last line of the `pyplotgen/config/Case_definitions.py` file such that it defines `ALL_CASES` to contain only the cases they wish to plot.  
This is what the bottom of `Case_definitions.py` normally looks like:
~~~~python
# DO NOT EDIT THIS LIST UNLESS YOU ARE ADDING A NEW CASE. NEVER REMOVE CASES FROM THIS LIST. You may define a subset of cases at the end of this file.
ALL_CASES = [ARM, ARM_97, ASTEX_A209, ATEX,
             BOMEX,
             CGILS_S6, CGILS_S11, CGILS_S12, CLEX9_NOV02, CLEX9_OCT14,
             DYCOMS2_RF01, DYCOMS2_RF01_FIXED_SST, DYCOMS2_RF02_DO, DYCOMS2_RF02_DS, DYCOMS2_RF02_ND, DYCOMS2_RF02_SO,
             FIRE,
             GABLS2, GABLS3, GABLS3_NIGHT,
             JUN25_ALTOCU,
             LBA,
             MC3E, MPACE_A, MPACE_B, MPACE_B_SILHS,
             NOV11_ALTOCU,
             RICO,
             TWP_ICE,
             WANGARA
             ]

# If uncommented, this line will override the real ALL_CASES given above, forcing pyplotgen to only plot some cases.
# ALL_CASES = [BOMEX]
~~~~

If someone wanted to plot _just_ the bomex case, even if other cases are in the same input folder, they can uncomment the last line to do so.

~~~~python
# DO NOT EDIT THIS LIST UNLESS YOU ARE ADDING A NEW CASE. NEVER REMOVE CASES FROM THIS LIST. You may define a subset of cases at the end of this file.
ALL_CASES = [ARM, ARM_97, ASTEX_A209, ATEX,
             BOMEX,
             CGILS_S6, CGILS_S11, CGILS_S12, CLEX9_NOV02, CLEX9_OCT14,
             DYCOMS2_RF01, DYCOMS2_RF01_FIXED_SST, DYCOMS2_RF02_DO, DYCOMS2_RF02_DS, DYCOMS2_RF02_ND, DYCOMS2_RF02_SO,
             FIRE,
             GABLS2, GABLS3, GABLS3_NIGHT,
             JUN25_ALTOCU,
             LBA,
             MC3E, MPACE_A, MPACE_B, MPACE_B_SILHS,
             NOV11_ALTOCU,
             RICO,
             TWP_ICE,
             WANGARA
             ]

# If uncommented, this line will override the real ALL_CASES given above, forcing pyplotgen to only plot some cases.
ALL_CASES = [BOMEX]
~~~~
`BOMEX` can be replaced with the name of any valid case, so long as it is one of the variables defined within this file. For example, to plot the fire and wangara case, one would do
~~~~python
ALL_CASES = [FIRE, WANGARA]
~~~~

## Making Publish Ready Plots
There are a few aspects to making a plot publish ready.
1. Ensuring the desired colors are used
2. Using math text formatting
3. Cleaning legend labels

### Changing color pallet
In the `config/Style_definitions.py` file there are two relavent variables: `COLOR_ROTATION` and `STYLE_ROTATION`. 
There are suggested presets available for use defined in comments. Only 1 set of definitions should be active at a given
time. All others should be commented out. For publication, it is suggested to use the print-safe definitions. Using a set 
of color/style's is as easy as uncommenting the desired set and commenting out all others.

### Using Math Text Formatting
Matplotlib (the plotting engine used in the background of pyplotgen) allows for LaTeX style math text formatting. 
Pyplotgen allows for titles and dependent axis labels to be overriden. To create a custom title/axis label, simply 
include a 'title' or 'axis_title' parameter in the variables definition. Additionally, these custom titles can utilize 
TeX style formatting. To do so, the math text must be wrapped in $'s, then matplotlib math text formatting can be used.
For more information on matplotlib's text formatting, see here: https://matplotlib.org/1.3.1/users/mathtext.html  
Note that super/subscript text by default only incorporates 1 character. To script multiple characters, wrap it in {}'s.
Example:
~~~~python
{'var_names':
    {
        'clubb': ['Skrt_zt'],
        'sam': [self.getSkrtZtLesCalc,'Skrt_zt'],
        'coamps': [self.getSkrtZtLesCalc, 'Skrt_zt'],
        'r408': ['Skrt_zt'],
        'hoc': ['Skrt_zt'],
        'e3sm': ['Skrt_zt'],
        'cam': ['Skrt_zt'],
        'wrf': ['Skrt_zt'],
    },
    'title': "$Skw_{zt}$ Skewness of modern dynamic levels",
    'axis_title': "Skw_zt [-]",
    'sci_scale': 0,
},
~~~~
### Cleaning legend labels
Changing the text on a legend label is easy. To do so, simply change the foldername of the data being imported. Note 
that underscores (`_`) will be replaced with spaces automatically. This is done to make the legends more readable, as it 
allows python to know how to split up new lines for longer legend labels.

## Adding new Variables (needs updating)
Note: additional documentation is available in the "Reference Documentation" listed above
Each variable is included within a _VariableGroup_,and each variable group is included into a case. To add a new variable, it must be added to a variable group. Once a variable is in a variable group, all cases that plot that group will now also plot the new variable. Variables have many properties available to describe them, here they are (pulled from VariableGroup.py src code comments):

#### Valid dict keys/options:

| Parameter | Description |
| --- | --- |  
| *var_names* | A list of names various models refer to this variable as. E.g. ['wprtp', 'WPRTP', 'wpqtp']. This list is to include the variable name for any models that Pyplotgen is plotting. It does not matter what order the variables are in; pyplotgen will search the list from left -> right so more frequent names are _preferred_ to go on the left for efficiency.|
| *sam_calc*| (optional) A functional reference to a method that calculates a sam variable. This is given as the name of the function *without* the () after the name. E.g. `self.getThlmSamCalc` |
| *clubb_calc*| (optional) A functional reference to a method that calculates a clubb variable. This is given as the name of the function *without* the () after the name. E.g. `self.getThlmClubbCalc` |
| *coamps_calc* | (optional) A functional reference to a method that calculates a coamps variable. This is given as the name of the function *without* the () after the name. E.g. `self.getThlmCoampsCalc`|
| *cam_calc*| (optional) A functional reference to a method that calculates a cam variable. This is given as the name of the function *without* the () after the name. E.g. `self.getThlmCamCalc` |
| *hoc_calc*| (optional) A functional reference to a method that calculates an hoc-2005 variable. This is given as the name of the function *without* the () after the name. E.g. `self.getThlmHocCalc` |
| *r408_calc*| (optional) A functional reference to a method that calculates a r408 variable. This is given as the name of the function *without* the () after the name. E.g. `self.getThlmR408SamCalc` |
| *e3sm_calc*| (optional) A functional reference to a method that calculates a e3sm variable. This is given as the name of the function *without* the () after the name. E.g. `self.getThlmE3smCalc` |
| *sam_conv_factor* | (optional) Numeric value to scale a sam variable by. E.g. `1/1000`, or `100`|
| *coamps_conv_factor*| (optional) Numeric value to scale a coamps variable by. E.g. `1/1000`, or `100`|
| *r408_conv_factor* | (optional) Numeric value to scale a clubb r408 variable by. E.g. `1/1000`, or `100`|
| *e3sm_conv_factor* | (optional) Numeric value to scale an e3sm variable by. E.g. `1/1000`, or `100`|
| *type* | (optional) Override the default type 'profile' with either 'budget' or 'timeseries'|
| *fallback_func* | (optional) (DEPRECATED: Use model_calc functions instead) If a variable is not found within a dataset and a `fallback_func` is specified, this function will be called to attempt retrieving the variable (before filling zeros if `fill_zeros=True`). Like the `model_calc` option, this is a functional reference to a method that calculates the given variable. E.g. `self.getWpthlpFallback`|
| *title* | (optional) Override the default panel title, or provide one if it's not specified in the netcdf file.|
| *axis_title* | (optional) Override the default dependent axis title, or provide one if it's not specified in the netcdf file.|
| *fill_zeros* | (optional) (DEPRECATED: Data arrays by default get filled with zeros now) If a variable isn't found in netcdf output, and there is no `_calc` function (or it failed) setting this to True allows PyPlotgen to fill all data points with 0 and continue plotting.|
| *lines* | (budget variables only) Defines lines to plot for budget cases. Passed separately because it's a lot of text. This is given in the form of a list of lines, below's an example:|
| *type* | (optional) Override the default type 'profile' with either 'budget', 'timeseries', or (not implemented) 'timeheight'

~~~~python
            thlm_lines = [
            {'var_names': ['thlm_bt'], 'label': 'thlm_bt'},
            {'var_names': ['thlm_ma'], 'label': 'thlm_ma'},
            {'var_names': ['thlm_ta'], 'label': 'thlm_ta'},
            {'var_names': ['thlm_mc'], 'label': 'thlm_mc'},
            {'var_names': ['thlm_clipping'], 'label': 'thlm_bt', 'clubb_calc': self.getThlmClipping},
            {'var_names': ['radht'], 'label': 'radht'},
            {'var_names': ['lsforcing'], 'label': 'lsforcing', 'clubb_calc': self.getLsforcing},
            {'var_names': ['thlm_residual'], 'label': 'thlm_residual', 'clubb_calc': self.getThlmResidual},
            ]
            self.variable_definitions = [
            {'var_names': ['thlm'], 'lines': thlm_lines, 'type': Panel.TYPE_BUDGET}
            ]
~~~~
Here are a couple examples of other (non-budget) variable definitions:

~~~~python{
            'var_names': ['Skrt_zt'],
            'sam_calc': self.getSkrtZtLesLine,
            'coamps_calc': self.getSkrtZtLesLine,
            },
            {
            'var_names': ['lwp', 'CWP'],
            'type': Panel.TYPE_TIMESERIES,
            'sam_conv_factor': 1/1000
            }
~~~~

### Troubleshooting/Advanced variables
A lot of the time the only parameter that is needed will be the list of aliases for the variable as used by multiple models. However there are times when a more complicated definition is needed, such as when variables must be calculated. For some of these situations, please check out the appropriate sections below.

#### Variable not found in case X
In the unfortunate event that a variable is exported in most cases except Case X Y Z, there are a few options available for handling this issue. The easy (and dirty) method of resolving this is to simply blacklist the variable and not plot it at all for that case. To do so, simply add the variable's name to the `blacklisted_vars` list for that case. For example, to blacklist `thlm` from the Wangara case, the blacklist would go from `'blacklisted_vars': []` to `'blacklisted_vars': ['thlm']`. If the missing variable should be plotted for the given case, then you must create a model_calc function (please see below). If neither option is applicable, the data array for the missing variable will be filled with zeros by default and a warning is written to the console.

#### Creating a new calculated function (for calculated variables)
Creating a new calculator function is relatively simple, but must be done correctly. 

Steps:  
1. Create a function following the naming scheme and method signature. This function should be created in the same `VariableGroup<Group Name>.py` file as the variable that's being handled.
~~~~python
    def get<VariableName><Model>Calc(self, dataset_override = None)
~~~~
2. Include the equation used for the calculation in a pydoc method comment, see example below
~~~~python
    def getWpthlpSamCalc(self, dataset_override = None):
        """

        :param self: 
        :param dataset_override: 
        """
        dependent_data
        # code
~~~~
3. Retrieve the variables needed for calculation using `self.getVarForCalculations(['list', 'of', 'var_names'], self.model_file) # self.model_file refers to variables like self.sam_file and self.coamps_file`. See example below
~~~~python
        tlflux = self.getVarForCalculations(['TLFLUX'], self.sam_file)
        rho = self.getVarForCalculations(['RHO'], self.sam_file)
        z = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)
~~~~
4. Perform the calculation using python functions
~~~~python
        wpthlp = tlflux / (rho * 1004)
~~~~
5. Return a tuple containing (calculated_variable, z)
~~~~python
        return wpthlp, z
~~~~

Here is the full example:
~~~~python
    def getWpthlpSamCalc(self, dataset_override = None):
        """

        :param self: 
        :param dataset_override: 
        :return: 
        """
        dependent_data
        tlflux = self.getVarForCalculations(['TLFLUX'], self.sam_file)
        rho = self.getVarForCalculations(['RHO'], self.sam_file)
        z = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)
        wpthlp = tlflux / (rho * 1004)
        return wpthlp, z
~~~~

Once this method is created, simply add this function to the variable's entry in the self.variable_definitions list at the top of the file. 
Example:
~~~~python
            {'var_names': { 'clubb': ['thlm'], ...} , 'sam_calc': self.getWpthlpSamCalc},
~~~~
**Note**: do NOT include () when giving the functions name, here we are using functional programming to pass the method itself as a parameter, not a call to the method, so we must not use the ().  
 
The variable data will now be calculated by the new function.

## Adding new Cases
Adding a new case is similar to adding a variable in that the process is generally simple, but must be done correctly. All cases are defined in the `pyplotgen/config/Case_definitions.py` file in the form of a dictionary. Copied over from the `Case_definitions.py` code comment, here are the parameters used to describe a case:

| **Parameter** | **Description** |
|---|---|
|*name*| must be the same as the filename without the extention. Eg. to use lba_zt.nc and lba_zm.nc the case's name must be 'lba'. Extensions are determined by the last instance of `_`|
|*start_time*| An integer value representing which timestep to begin the time-averaging interval. Valid options are from 1 -> list minute value. Give in terms of clubb minutes.|
|*end_time*| An integer value representing which timestep to end the time-averaging interval. Valid options are from 1 -> list minute value. Give in terms of clubb minutes. Also used to determine where to stop timeseries plots|
|*height_min_value*| The elevation to begin height plots at|
|*height_max_value*| The elevation to end height plots at|
|*blacklisted_vars*| List of variables to avoid plotting for this case. Names must use the clubb-name version|
|*sam_file*| Path to the SAM .nc file for this case (please use the SAM_OUTPUT_ROOT variable as a base).|
|*coamps_file*| dict containing paths to the coamps .nc files for this case (please use the LES_OUTPUT_ROOT variable as a base).|
|*r408_file*| dict containing paths to the r408 .nc files for this case (please use the R408_OUTPUT_ROOT variable as a base).|
|*hoc_file* | dict containing paths to the hod .nc files for this case (please use the HOC_OUTPUT_ROOT variable as a base).|
|*var_groups*| list of python class names, where the classes use the naming scheme VariableGroup____.py and define a variable group.|

Here's an example definition:
~~~~python
BOMEX = {'name': 'bomex', 'start_time': 181, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 
         'blacklisted_vars': [],
         'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/BOMEX_64x64x75/BOMEX_64x64x75_100m_40m_1s.nc",
         coamps_dataset: None,
         'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_zm.nc',
                       'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_zt.nc',
                       'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_sfc.nc'},
         'hoc_file': {'zm': HOC_OUTPUT_ROOT + '/bomex_zm.nc',
                      'zt': HOC_OUTPUT_ROOT + '/bomex_zt.nc',
                      'sfc': HOC_OUTPUT_ROOT + '/bomex_sfc.nc'},
         'var_groups': [VariableGroupBase, VariableGroupWs]}
~~~~~

Creating a new case in pyplotgen is as simple as creating one of these defintions (e.g. example above) and adding it to the `ALL_CASES` list:
~~~~python
ALL_CASES = [ARM, ARM_97, ASTEX_A209, ATEX,
             BOMEX,
             CGILS_S6, CGILS_S11, CGILS_S12, CLEX9_NOV02, CLEX9_OCT14,
             DYCOMS2_RF01, DYCOMS2_RF01_FIXED_SST, DYCOMS2_RF02_DO, DYCOMS2_RF02_DS, DYCOMS2_RF02_ND, DYCOMS2_RF02_SO,
             FIRE,
             GABLS2, GABLS3, GABLS3_NIGHT,
             JUN25_ALTOCU,
             LBA,
             MC3E, MPACE_A, MPACE_B, MPACE_B_SILHS,
             NOV11_ALTOCU,
             RICO,
             TWP_ICE,
             WANGARA
             ]
~~~~
Please maintain the list's alphabetical ordering for ease of use.

## Plotting SAM-exclusive VariableGroups
The VariableGroups `VariableGroupSamBudgets`, `VariableGroupSamProfiles` and `VariableGroupSamMultilineProfiles` were added to recreate the SAM plots done with corplot (python_sam_budgets_plotter).  
Pyplotgen will automatically generate the `VariableGroupSamBudgets` plots when using the `--plot-budgets` flag with SAM input.  
To generate the SAM profile plots, add the VariableGroups to the `var_groups` list in the definition of the cases to be plotted and, if not wanted, remove or comment out the other VariableGroups.  
These specific VariableGroups will only work with SAM input!

With the '-s' option multiple SAM input folders containing NetCDF files can be specified.  
To plot SAM data for a specific case the name of the input NetCDF file needs to be the same as the string given under the key `sam_file` in the same case defintion.  
Example for a modified case definition for BOMEX:
~~~~python
BOMEX = {'name': 'bomex', 'start_time': 181, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 2500,
         
         'blacklisted_vars': [],
         'sam_file': "bomex.nc", # Specify SAM NetCDF file here!
         'coamps_file': None,
         'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_zm.nc',
                       'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_zt.nc',
                       'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_sfc.nc'},
         'hoc_file': {'zm': HOC_OUTPUT_ROOT + '/bomex_zm.nc',
                      'zt': HOC_OUTPUT_ROOT + '/bomex_zt.nc',
                      'sfc': HOC_OUTPUT_ROOT + '/bomex_sfc.nc'},
         'e3sm_file': e3sm_output_root + '/bomex.nc',
         'var_groups': [#VariableGroupBase, VariableGroupWs,
				VariableGroupSamProfiles]} # New VariableGroup added and rest commented out
~~~~
Example command line call for SAM-only plots with two input folders:
```bash
python3 ./pyplotgen.py -s first/path/to/SAM/folder second/path/to/SAM/folder
```

## Running pyplotgen on a Windows system
The easiest way to run pyplotgen on Windows at the moment is to install a Linux shell emulator.
This will take care of most problems concerning interoperability.  
Running pyplotgen on Windows was tested using Git Bash which is part of the git installation for windows and can be found [here](https://git-scm.com/).  
Using Git Bash one can simply follow the same procedure as for regular Linux bash to run pyplotgen.  
Paths can be specified using either the Windows (C:\User\testuser\...) or the Unix format (/home/testuser/...).

# Pyplotgen code convention
These conventions build from, and may modify, the existing CLUBB naming convention.

### Naming
* **`UPPER_SNAKE_CASE`**: This naing scheme is used to represent "final" variables. These are variables that **not** changed during the runtime of pyplotgen. Python does not enforce this, as it is only a naming convention, so _please_ do not introduce any code that modifies these variables. The one and only exception to this is the `ALL_CASES` variable in `Case_definitions.py`. The variable is rewritten there IMMEDIATELY after it was originally written only because it allows us to run a small sample of cases without having to rewrite the true ALL_CASES variable.
    *  **Naming case definitions**: When naming a case definition, the convention is to set the 'name' parameter equal to a lower case version of the variable name used.
* **`lower_snake_case`**: This is the accepted norm for python variable names. Please use this for any regular use variables. Python notation is infamous for its heavy use of acronyms and shorthand, however due to the nature of our work it is recommended to spell variable names out if possible.
* **`\_\_underScoredCamelCase\_\_()`**: Two underscores around a camelcased method name represents a "helper method". In other langauges this is refered to as a private method. These methods help write clean code, but are not intended for use outside of the code within the class they've been written. E.g. `Foo.__bar__()` should only be written within the `Foo.py` file/ within the `Foo` class if the file contains multiple classes.
* **`lowerCamelCase()`**: This is used for normal function names. For non-helper functions (see previous), use this naming convention. This deviates from the typical python convention. There is no strong reason behind this other than it's what's been done and it was the format for pyplotgen's legacy code. It may be desireable to switch to a `lower_snake_case()` in the future. This would be a simple yet large refactor.
    * `UpperCamelCase()` is to be used when declaring a class name.
