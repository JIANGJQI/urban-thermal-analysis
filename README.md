# Urban Thermal Analysis on the Qinghai-Tibet Plateau

This repository contains a public, reduced version of a research workflow for
analyzing season-specific thermal contrasts of hillside urbanization on the
Qinghai-Tibet Plateau (QTP). The workflow applies the elevation-constrained,
patch-specific framework of Xu et al. (2026) to summer and winter Landsat LST
composites for 2000 and 2020.

The associated manuscript is currently unpublished, so the full locked
analytical dataset is not included here.

## Scope

- Study region: Qinghai-Tibet Plateau.
- Endpoint years: 2000 and 2020.
- Seasons: summer and winter.
- Terrain structure: six slope classes and six 500 m elevation bands from
  0 to 3000 m.
- Thermal contrast: patch-footprint LST (`T1`) minus same-elevation
  non-impervious reference LST (`T0`).
- Main code path: patch-level thermal reproduction, locked-output assembly,
  and landscape/thermal association calculations.

## Repository Contents

```text
src/
  reproduce_xu_method.py
  lock_xu_outputs.py
  compute_xu_warming_area_correlations.py

data/
  locked_sample/
    analysis_lock.json
    DATA_DICTIONARY.json
    input_manifest.csv
    input_hash_manifest.csv
    quality_audit.csv
    patch_thermal_sample.csv
    summary_sample.csv
    SAMPLE_MANIFEST.json

figures/
  fig5_seasonal_asymmetry_reproduced.png
```

`data/locked_sample/` is a structured sample of the internal locked analytical
outputs. It preserves the schema and year-season-slope-elevation dimensions,
but it is not the full dataset used to produce the manuscript results.

`figures/` contains one representative result figure showing the
elevation-dependent seasonal asymmetry pattern.

## Data Availability

Raw raster/vector inputs and the full locked analytical outputs are not
publicly released at this stage because the manuscript is unpublished and some
source datasets have redistribution constraints.

The public sample includes metadata, audit files, input manifests, and a small
patch-level sample so that readers can inspect the data structure and workflow
logic without reproducing the full unpublished analysis.

Included sample files are stored under `data/locked_sample/`:

- `analysis_lock.json`
- `DATA_DICTIONARY.json`
- `input_manifest.csv`
- `input_hash_manifest.csv`
- `quality_audit.csv`
- `patch_thermal_sample.csv`
- `summary_sample.csv`
- `SAMPLE_MANIFEST.json`

The sample preserves the schema and the main analytical dimensions
(year, season, slope class, and elevation band), but it is not sufficient to
reproduce the full manuscript results.

## Reproduction Notes

The scripts in `src/` document the internal reproduction workflow:

1. `reproduce_xu_method.py` computes patch-level `T1`, `T0`, and `DeltaT`.
2. `lock_xu_outputs.py` assembles year-season outputs into locked summary
   tables and audit files.
3. `compute_xu_warming_area_correlations.py` computes positive-contrast
   footprint metrics and landscape associations.

Running the full workflow requires the internal raw rasters, patch geometries,
and full locked intermediate files, which are not included in this public
repository.
