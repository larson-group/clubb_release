# P3 microphysics in CLUBB

This folder holds the P3 (Predicted Particle Properties) bulk microphysics
scheme of Morrison and Milbrandt (2015). CLUBB standalone can use it as the
microphysics scheme with `microphys_scheme = "p3"`. The connection to CLUBB is
[`p3_microphys_module.F90`](../p3_microphys_module.F90); its header explains how
CLUBB's fields are passed to P3 and how P3's results return as tendencies.

## Running a case with P3

Select an existing case from
[`input/case_definitions.jsonc`](../../../input/case_definitions.jsonc)
and override its microphysics settings with
[`run_scm.py`](../../../run_scripts/run_scm.py). Run these examples from the
repository root:

```bash
python3 run_scripts/run_scm.py dycoms2_rf02_ds \
  -override "microphysics_setting.microphys_scheme=p3,microphysics_setting.l_cloud_sed=false"
python3 run_scripts/run_scm.py rico \
  -override "microphysics_setting.microphys_scheme=p3,microphysics_setting.l_ice_microphys=true"
python3 run_scripts/run_scm.py mpace_b \
  -override "microphysics_setting.microphys_scheme=p3,microphysics_setting.l_predict_nc=false"
```

Ice is carried in CLUBB only with `l_ice_microphys = true`, already enabled
in `mpace_b`. Without it, a run stops if P3 forms ice. The RICO examples
enable ice to retain small amounts that P3 can form in cold SILHS samples.
RF02 disables CLUBB cloud sedimentation because
P3 already sediments cloud droplets. Other settings retain the selected
case's values, including droplet number and PDF parameters.

These examples run P3 on grid means and print a warning at startup: P3
uses no cloud fraction, so rain formation in partly cloudy layers is likely
wrong. For CLUBB's subgrid variability to reach P3, run it on SILHS samples
as described in [P3 on SILHS samples](#p3-on-silhs-samples).

The microphysics settings `p3_iparam` and `p3_lookup_table_dir` select the
warm-rain scheme and the lookup table directory. For example, to use the
Khairoutdinov and Kogan (2000) warm-rain formulas instead of the default
Seifert and Beheng (2001):

```bash
python3 run_scripts/run_scm.py dycoms2_rf02_ds \
  -override "microphysics_setting.microphys_scheme=p3,microphysics_setting.l_cloud_sed=false,microphysics_setting.p3_iparam=3"
```

### P3 on SILHS samples

P3 expects each grid box to be either fully cloudy or clear, as in CM1. A
CLUBB grid box is usually partly cloudy, so with grid means P3 sees a thin
cloud spread over the whole box and forms little rain. With
`lh_microphys_type = "interactive"`, SILHS draws sample subcolumns from
CLUBB's subgrid PDF instead. Each sample is fully cloudy or clear at each
level and carries its own cloud water, vapor and rain. P3 runs on every
sample, and the tendencies are averaged with the sample weights. SILHS also
samples the variation within the cloud, to which the Seifert and Beheng
(2001) autoconversion is sensitive.

Use `rico_silhs` to retain its existing SILHS and PDF settings with P3:

```bash
python3 run_scripts/run_scm.py rico_silhs \
  -override "microphysics_setting.microphys_scheme=p3,microphysics_setting.l_ice_microphys=true"
```

The number of samples is `lh_num_samples` (default 8). P3 runs once per
sample, so the microphysics costs that many times more. SILHS with P3
requires prescribed cloud-droplet number (`l_predict_Nc = false`). Ice can
be carried with `l_ice_microphys = true`, for example:

```bash
python3 run_scripts/run_scm.py mpace_b \
  -override "microphysics_setting.microphys_scheme=p3,microphysics_setting.l_predict_nc=false,microphysics_setting.lh_microphys_type=interactive"
```

With ice enabled, SILHS samples total ice mass and ice-particle number.
Each sample uses the grid-mean rime mass fraction and rime density: rime
mass scales with sampled total ice, and rime volume scales with rime mass.
These rime samples are constructed after sample clipping, before sample
statistics and microphysics read them.
This keeps rime mass within total ice without introducing a separate rime
distribution. Incoming mean rime mass must not exceed total ice mass;
sampling stops if that relationship is violated. For valid inputs, rime
fields are zero when mean ice mass, rime mass or rime volume is nonpositive.
The sample means inherit any sampling error in
total ice; the method does not normalize them back to the grid means.

After hydrometeors are advanced, CLUBB also bounds mean rime mass by mean
total ice and reduces rime volume proportionally. This preserves rime
density without changing total ice water when sampling error or transport
would otherwise leave inconsistent ice properties.

## What is supported

- One ice category with two-moment ice, the configuration of the CM1 RICO
  runs that this P3 version comes from.
- Rain mass and number (`rrm`, `Nrm`) are CLUBB hydrometeors. With
  `l_ice_microphys`, so are P3's ice mass, ice number, rime mass and rime
  volume (`rim`, `Nim`, `qirim`, `birim`). Rime mass is part of `rim`, so CLUBB
  does not count it as separate water. Without ice support, a run stops if
  P3 forms ice, including within a SILHS sample.
- P3 sediments cloud droplets, rain and ice itself, so `l_cloud_sed` must be
  false.
- P3 runs without its subgrid cloud and precipitation fraction option
  (SCPF), as in CM1. Without SILHS, it receives grid means and treats each
  grid box as fully cloudy, so the grid means serve as in-cloud values. In
  partly cloudy layers this underestimates autoconversion: for given grid
  means, the Seifert and Beheng (2001) rate grows about as 1 / cloud
  fraction. P3 does not use CLUBB's cloud fraction, subgrid PDF or
  precipitation fraction, except through SILHS. A run prints a warning at
  startup unless `lh_microphys_type = "interactive"`.
- SILHS supports warm rain and ice with prescribed cloud-droplet number.
  The predicted-cloud-droplet-number budget with SILHS has not been validated.
  Ice samples use the grid-mean rime fraction and density described above.
  P3 sediments cloud droplets, rain and ice within each sample subcolumn.
- The Seifert and Onishi (2016) turbulent collision option and the Python and
  JAX drivers are not supported.

## Source of the code and lookup table

[`module_mp_p3.F90`](module_mp_p3.F90) is `src/module_mp_p3.F` from the
`rico_turb` branch of [CM1R21-SDM](https://github.com/kkchandr/CM1R21-SDM)
at commit `8ecef10` (P3 version 4.0.7-HM-TurbCC, last updated 2026-11-02).
That fork of CM1 release 21 adds selectable warm-rain schemes, RICO aerosol
for droplet activation and the Seifert and Onishi (2016) option. Its RICO
runs use this P3 with one ice category and two-moment cloud water (CM1
`ptype = 51`), as CLUBB does.

P3 reads its ice lookup table at initialization even without ice, so every
P3 run needs it. The table is not stored in this repository. Download CM1
release 21.0 (`cm1r21.0.tar.gz`) from the
[CM1 website](https://www2.mmm.ucar.edu/people/bryan/cm1/getcode.html) and
copy `run/p3_lookupTable_1.dat-2momI_v5.1.6_oldDimax` into
`input/p3_lookup_tables/`. Run from the repository root:

```bash
tar -xzf cm1r21.0.tar.gz cm1r21.0/run/p3_lookupTable_1.dat-2momI_v5.1.6_oldDimax
mkdir -p input/p3_lookup_tables
mv cm1r21.0/run/p3_lookupTable_1.dat-2momI_v5.1.6_oldDimax input/p3_lookup_tables/
```

The same file is in the `run/` directory of the CM1R21-SDM commit above. P3
checks the table version, so later CM1 releases, which ship `v5.4` tables,
do not work with this P3. Use `p3_lookup_table_dir` to read the table from
another directory. The ice-ice table `p3_lookupTable_2.dat-4.1` is needed only
for more than one ice category, which CLUBB does not support.

## Changes from the CM1R21-SDM source

- The file is renamed from `.F` to `.F90` so that it is compiled as free-form
  source, as in CM1.
- The three-moment and multi-category lookup tables are allocatable and are
  allocated in `p3_init` only when those options are used. As static arrays
  they took about 80 MB, which debug builds (`-finit-real=snan`) stored in
  every executable that links P3. This change was a TODO in the source.

## References

- Morrison, H. and J. A. Milbrandt, 2015: Parameterization of cloud
  microphysics based on the prediction of bulk ice particle properties.
  Part I: Scheme description and idealized tests. *J. Atmos. Sci.*, 72,
  287-311.
- Seifert, A. and K. D. Beheng, 2001: A double-moment parameterization for
  simulating autoconversion, accretion and selfcollection. *Atmos. Res.*,
  59-60, 265-281.
- Khairoutdinov, M. and Y. Kogan, 2000: A new cloud physics parameterization
  in a large-eddy simulation model of marine stratocumulus. *Mon. Wea. Rev.*,
  128, 229-243.
