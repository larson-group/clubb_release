# Model inputs

The four input catalogs define the cases, parameter values, flags and output
statistics used by the standalone model, JAX and loss drivers. Run scripts and
Dash edit JSON. [create_case_namelist.py](scripts/create_case_namelist.py) resolves those
settings and writes one complete namelist immediately before a top-level model
entry point reads it.

The maintained catalogs use JSONC, which permits `//` and `/* ... */` comments.
Shared explanations appear beside `defaults`; notes about a particular case or
configuration appear beside its preset. Ordinary `.json` overrides remain
accepted. Both formats require double-quoted keys and disallow trailing commas.
Dash preserves comments when updating a catalog.

Preparation modules live under [scripts/](scripts/); catalogs and other model
input data stay directly in this directory or their existing data subdirectories.

From the repository root:

```bash
python3 run_scripts/run_scm.py -cases bomex -max_iters 10
python3 run_scripts/run_scm.py -cases arm -config e3sm -flags ~/custom_flags.json
python3 run_scripts/run_scm.py bomex -params ~/custom_parameters.json -stats ~/custom_stats.json
```

## Defaults, shared configs and explicit selections

The case, parameter and flag catalogs contain `defaults` and sparse `presets`.
Case and flag catalogs can also contain inherited `groups` and ordered `lists`.
Their resolution is ordered:

1. Start with the section's catalog defaults.
2. Apply the preset named by `-config`, wherever that name exists.
3. Apply the explicit section preset or JSON file.
4. Apply direct CLI settings, followed by repeated `-override` arguments in order.

An omitted section retains its defaults when the shared config has no matching
preset. `-config default` is a built-in choice that uses the top-level defaults
without requiring a preset. Dash offers this same baseline choice, and saved
experiments must use another name. Other shared configs must exist in at least
one catalog; an explicitly requested section preset must exist in that section.
Empty object presets can be removed when no explicit section selector uses them.
Missing files and malformed definitions fail before the runner starts.
Relative, absolute and `~` paths are accepted.

| Selector | Catalog | Preset/custom-file contents |
| --- | --- | --- |
| Positional case or `-cases` | [case_definitions.jsonc](case_definitions.jsonc) | Overrides grouped by model, microphysics, radiation and stats namelist names |
| `-params` (`-params_file` also remains accepted) | [clubb_params.jsonc](clubb_params.jsonc) | Parameter overrides within `model` and/or `silhs` |
| `-flags` | [clubb_config_flags.jsonc](clubb_config_flags.jsonc) | Boolean flags and integer selectors within `model` and/or `silhs` |
| `-stats` | [stats.jsonc](stats.jsonc) | A list of variable names; metadata comes from the catalog's `all` tree |

A custom case, parameter or flag file has exactly the same contents as one
preset entry. To promote a custom definition, copy its contents under a new name
in `presets`. Objects merge recursively and arrays replace whole arrays. Explicit
files and presets are sparse patches; omitted keys preserve values supplied by
the shared config.
`null` setting values are rejected, including inside arrays. Use an explicit
number, boolean, string or empty list for the intended exception.

Use explicit values for case exceptions, such as `"sfctype": 0`. Inactive
configuration groups can inherit their catalog defaults; individual disabled
input fields can be written as `false`. Saved configurations retain inherited
settings rather than deleting them to select a hidden model default.

Keep the experiment's date (`day`, `month`, `year`), location (`lat_vals`,
`lon_vals`) and simulation period (`time_initial`, `time_final`) in its case
definition. These fields have no shared catalog defaults. Latitude/longitude
are degrees north/east; start/end times are seconds since midnight UTC on the
start date. Dataset-driven cases can retain their intentional omissions.

Case defaults keep optional processes disabled while defining their settings.
Parameter defaults favor common values among cases that use each feature to
reduce repeated overrides. They can differ from Fortran initialization defaults;
each case retains the exceptions needed by its active processes. For example,
the BUGSrad diffuse visible albedo defaults to `0.07`, with `0.1` and `0.23`
retained as case overrides where needed. Explicit choices for inactive processes
also stay with their cases; otherwise inactive parameters inherit the common
catalog values.
Comments group settings shared across schemes and those specific to Morrison,
Khairoutdinov-Kogan, SILHS, simplified radiation and BUGSrad. Each named case
selects its schemes, enables the processes it uses and overrides its scientific
parameters. A custom case file inherits the same disabled processes unless it
enables them explicitly. For example, selecting `rad_scheme: "simplified"` enables
radiation; set `l_calc_thlp2_rad: true` separately to include its contribution to
liquid-water potential-temperature variance. Parameters for inactive processes
can stay defined, such as `eff_drop_radius: 1.0e-5` m for simplified shortwave
radiation. Time-dependent schedules and case-specific aerosol distributions stay
with their cases.

For example, a parameter file can contain only:

```json
{"model": {"C8": 0.4}, "silhs": {"vert_decorr_coef": 0.1}}
```

A custom case declares its experiment context:

```json
{
  "model_setting": {
    "runtype": "bomex",
    "day": 22, "month": 6, "year": 1969,
    "lat_vals": 15.0, "lon_vals": -56.5,
    "time_initial": 0.0, "time_final": 21600.0,
    "dt_main": 30.0
  }
}
```

A custom stats file can contain:

```json
["thlm", "rtm", "wp2"]
```

Stats use an independent hierarchy in `stats.jsonc`, without defaults or sparse
overrides. `-config` does not affect stats. Omitting `-stats` selects `standard`;
selecting a category includes all its descendants. The standard category is a
reduced selection of the legacy standard outputs and includes the fields required
by native restart reads. Some pyplotgen plots need optional diagnostics from
`non-standard`; for example, add scalar outputs with
`-stats standard+non-standard/scalars`.
`-stats standard/budgets` selects one subtree, while
`-stats core+radiation` combines matching categories across branches.
`var:NAME` selects one variable. `all` selects everything and `none` disables stats.

Stats presets and custom JSON lists are complete selections. A preset alone
selects exactly its members; `-stats tuning_stats+core` selects their union.
Every variable is emitted once. A single preset or custom list keeps its stored
order; unions follow catalog order. An empty custom list selects no variables.
List predefined selectors without creating a namelist:

```bash
python3 input/scripts/create_case_namelist.py -list parameters
python3 input/scripts/create_case_namelist.py -list configs
python3 input/scripts/create_case_namelist.py -list stats
```

## Inherited groups and run lists

Case definitions may be nested under `groups`. Every group name ends in
`_group`, using its family's common name. Every leaf case name remains globally
unique and can still be selected directly. Group `settings` apply to every
child; nested groups add another inheritance level:

```jsonc
"groups": {
  "new_cases_group": {
    "settings": {"model_setting": {"dt_main": 30}},
    "presets": {
      "new_case_default": {},
      "new_case_cold": {"model_setting": {"t_sfc_nl": 270}},
      "new_case_hot": {"model_setting": {"t_sfc_nl": 300}}
    }
  }
},
"lists": {
  "smoke_cases": ["bomex", "new_cases_group"]
}
```

A case resolves catalog defaults, a matching `-config` preset, outer-to-inner
group settings, its own patch, then direct run options and general overrides.
Objects merge recursively and arrays replace whole arrays at every level.
Selecting a leaf directly retains its ancestry. Saved case edits update only
that leaf relative to its inherited base; copying a case into a custom file
stores its settings relative to catalog defaults.

`lists` are ordered selections without settings. They can reference cases,
groups or other lists; cycles, unknown members and ambiguous names fail before
execution. Duplicates run once in first-use order. The five maintained lists are
`all_cases`, `standard_cases`, `priority_cases`, `min_cases` and `short_cases`.
Their membership preserves the previous batch-run selections, rather than
implicitly including every new case added to the catalog.

Flag collections use the same group structure with `model`/`silhs` settings.
Their member names are local: two groups can each define `flag1`. A single
member is selected as `GROUP/flag1`. A flag list can combine groups or qualified
members. Each selected member starts from the selected baseline independently.
For example:

```sh
python3 run_scripts/run_scm.py -cases clex_group -workers 2
python3 run_scripts/run_scm.py -flags run_bindiff_w_flags_config_core_flags_group/flag1 bomex
python3 run_scripts/run_clubb_w_varying_flags.py -flag_group run_bindiff_w_flags_config_core_flags_group -cases priority_cases
```

The migrated flag collections retain the original collection stems followed by
`_group`, all local member names and all sparse override values. The JAX-specific
collection keeps its unsupported-algorithm and roundoff exclusions; explanatory
notes live beside the groups in the flag catalog. Comparison runners retain the
optional unmodified baseline and apply sweep overrides after their curated case
setup. Flags are logical/integer algorithm switches; real-valued tunable
coefficients remain in the parameter catalog.

## Preparation tools

The central builder can also be called directly:

```bash
python3 input/scripts/create_case_namelist.py bomex -output_dir /tmp/bomex
python3 input/scripts/create_case_namelist.py -cases arm -config e3sm -stats core -max_iters 10
```

The colocated [case](scripts/case_json_to_namelist.py),
[parameter](scripts/parameters_json_to_namelist.py), [flag](scripts/flags_json_to_namelist.py) and
[stats](scripts/stats_json_to_namelist.py) converters are internal library functions.
[input_catalog.py](scripts/input_catalog.py) owns shared catalog selection and
precedence; [namelist.py](scripts/namelist.py) owns JSON reading/merging and
low-level serialization. The final writer serializes the assembled groups
without converting them back into section-specific settings.
Preparing or displaying settings does not generate intermediate namelist files.
Python and JAX retain their namelist readers at the model entry point.

[create_multi_col_params.py](scripts/create_multi_col_params.py) supplies the
in-memory column generation used by normal runs and loss preparation. These
paths no longer write intermediate parameter JSON snapshots. Its standalone
CLI exports JSON column arrays and optional batch dimensions:

```bash
python3 input/scripts/create_multi_col_params.py -multicol 4 -output_file /tmp/columns.json
python3 input/scripts/create_multi_col_params.py -params_file input/clubb_params.jsonc -multicol C8/0.2:0.8/4 -batch_size 2 -output_file /tmp/grid.json
```

[save_tunable_config.py](scripts/save_tunable_config.py) saves Dash parameter/flag edits
as sparse named presets with provenance. The case editor saves sparse
`*.case.json` overrides; stats lists use `*.stats.json`. These files remain JSON
when read by another preparation tool.

[benchmark_converter.py](scripts/benchmark_converter.py) prepares LES benchmark NetCDF
inputs for loss evaluation:

```bash
python3 input/scripts/benchmark_converter.py -list_fields
python3 input/scripts/benchmark_converter.py input/les_and_clubb_benchmark_runs/les_runs/bomex_coamps_sm.nc /tmp/bomex_benchmark.nc -source coamps -fields cloud_frac
```

[flag_sets.py](scripts/flag_sets.py) expands catalog flag selections for
comparison sweeps and renders their final overrides.

The tables in [case_data/](case_data/), grids, atmospheric profiles and
benchmark NetCDF files remain data inputs. Their `.in` or other suffixes do not
make them model/parameter/flag/stats namelist definitions.


## Case input files

Cases select their data files explicitly in `model_setting` in
[case_definitions.jsonc](case_definitions.jsonc). Several cases, or several input
roles within one case, can reference the same file. For example, BOMEX uses one
scalar profile for both passive and eddy-diffusivity scalars:

```json
"sounding_file": "../input/case_data/bomex_sounding.in",
"sclr_sounding_file": "../input/case_data/bomex_sclr_sounding.in",
"edsclr_sounding_file": "../input/case_data/bomex_sclr_sounding.in"
```

The filename settings are `sounding_file`, `sclr_sounding_file`,
`edsclr_sounding_file`, `ozone_sounding_file`, `forcing_file`, `sfc_file`,
`corr_file_cloud`, `corr_file_below`, and `iop_file`. Leave unused inputs empty.
`forcings_file_path` remains the directory setting for cases that read specialized
forcing tables, such as LBA and MPACE-A.

Paths are absolute or relative to `run_scripts/`, matching the standalone's run
directory. Python and JAX resolve these paths once during initialization,
independently of where the generated namelist is stored. Readers receive the
selected paths directly; `runtype` selects case behavior and does not select
input filenames. Hand-written namelists must supply the paths needed by their
enabled features.

To use a different sounding without creating a new case preset, run from the
repository root:

```bash
python3 run_scripts/run_scm.py -cases bomex -override 'model_setting.sounding_file=/path/to/my_profile.in'
```
