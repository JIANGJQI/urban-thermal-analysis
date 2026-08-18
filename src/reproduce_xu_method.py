"""Independent reproduction of the Xu et al. (2026) thermal method."""

from __future__ import annotations

import argparse
import math
import pickle
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import shapely
from rasterio.enums import Resampling
from rasterio.features import geometry_mask, rasterize
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window, from_bounds

ROOT = Path(r"D:\hill")
PROJECT = ROOT / "paper2"
SOURCE = ROOT / "数据输出"
PATCH_DIR = SOURCE / "03_Results" / "patches"
LST_DIR = SOURCE / "04_LST"
IS_DIR = SOURCE / "02_LandCover"
DEM_PATH = SOURCE / "01_DEM_Slope" / "QTP_DEM_30m.tif"
OUT = PROJECT / "data" / "locked_xu2026"
OUT.mkdir(parents=True, exist_ok=True)

YEARS = (2000, 2020)
SEASONS = ("summer", "winter")
FOOTPRINT_MULTIPLIER = 6.0
ELEVATION_BIN_M = 1.0
SLOPE_LABELS = ("0-5", "5-10", "10-15", "15-20", "20-30", "30-plus")


def paths_for(year: int, season: str) -> dict[str, Path]:
    return {
        "patches": PATCH_DIR / f"is_patches_{year}.shp",
        "lst": LST_DIR / f"LST_{year}_{season}_albers.tif",
        "dem": DEM_PATH,
        "impervious": IS_DIR / f"QTP_IS_{year}.tif",
    }


def validate(paths: dict[str, Path]) -> None:
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing inputs:\n" + "\n".join(missing))
    with rasterio.open(paths["lst"]) as lst, rasterio.open(paths["dem"]) as dem:
        if lst.crs != dem.crs:
            raise ValueError("LST and DEM CRS differ")
        if not math.isclose(abs(lst.res[0]), 100.0, abs_tol=1e-6):
            raise ValueError(f"Expected delivered 100 m LST, found {lst.res}")
        if not math.isclose(abs(dem.res[0]), 30.0, abs_tol=1e-6):
            raise ValueError(f"Expected 30 m DEM, found {dem.res}")


def clipped_window(dataset, geometry):
    window = from_bounds(*geometry.bounds, transform=dataset.transform)
    return window.round_offsets().round_lengths().intersection(
        Window(0, 0, dataset.width, dataset.height)
    )


def overlap_weights(dataset, geometry, window):
    height, width = int(window.height), int(window.width)
    transform = dataset.window_transform(window)
    rows, cols = np.indices((height, width))
    left = transform.c + cols * transform.a
    right = left + transform.a
    top = transform.f + rows * transform.e
    bottom = top + transform.e
    cells = shapely.box(
        np.minimum(left, right), np.minimum(bottom, top),
        np.maximum(left, right), np.maximum(bottom, top),
    )
    intersections = shapely.intersection(cells, geometry)
    return shapely.area(intersections) / abs(transform.a * transform.e)


def weighted_mean(values, weights):
    valid = np.isfinite(values) & (weights > 0)
    if not valid.any():
        return np.nan, 0, 0.0
    selected = weights[valid]
    return (
        float(np.average(values[valid], weights=selected)),
        int(valid.sum()), float(selected.sum()),
    )


def footprint_for_multiplier(geometry, multiplier=FOOTPRINT_MULTIPLIER):
    clean = geometry if geometry.is_valid else geometry.buffer(0)
    target = clean.area * multiplier
    low, high = 0.0, max(1.0, clean.area ** 0.5)
    while clean.buffer(high).area < target:
        high *= 2
    for _ in range(70):
        middle = (low + high) / 2
        if clean.buffer(middle).area < target:
            low = middle
        else:
            high = middle
    return clean.buffer((low + high) / 2)


def footprint_elevation_range(dem, footprint):
    window = clipped_window(dem, footprint)
    values = dem.read(1, window=window, masked=True)
    inside = geometry_mask(
        [footprint], out_shape=values.shape,
        transform=dem.window_transform(window),
        invert=True, all_touched=True,
    )
    valid = inside & (~np.ma.getmaskarray(values)) & np.isfinite(values.data)
    if not valid.any():
        return np.nan, np.nan
    selected = values.data[valid]
    return float(selected.min()), float(selected.max())


def make_footprints(patches, year=None):
    cache = OUT / f"footprints_{year}.pkl" if year is not None else None
    if cache is not None and cache.exists():
        with cache.open("rb") as stream:
            footprints = pickle.load(stream)
        if len(footprints) != len(patches):
            raise ValueError("Cached footprint count does not match patches")
        print(f"loaded footprint cache: {cache}", flush=True)
        return footprints
    started = time.perf_counter()
    footprints = []
    for position, geometry in enumerate(patches.geometry, start=1):
        footprints.append(footprint_for_multiplier(geometry))
        if position % 1000 == 0:
            print(f"footprints: {position}/{len(patches)}", flush=True)
    print(f"footprints completed in {time.perf_counter()-started:.1f}s", flush=True)
    if cache is not None:
        with cache.open("wb") as stream:
            pickle.dump(footprints, stream, protocol=pickle.HIGHEST_PROTOCOL)
    return footprints


def reference_histogram(lst, dem_on_lst, is_on_lst, footprints):
    """Region-wide elevation lookup after excluding IS and all T1 footprints."""
    sums, counts = {}, {}
    index = gpd.GeoSeries(footprints, crs=lst.crs).sindex
    for block_number, (_, window) in enumerate(lst.block_windows(1), start=1):
        values = lst.read(1, window=window, masked=True).filled(np.nan)
        dem = dem_on_lst.read(1, window=window, masked=True).filled(np.nan)
        impervious = is_on_lst.read(1, window=window, masked=True).filled(1)
        bounds = rasterio.windows.bounds(window, lst.transform)
        positions = index.query(shapely.box(*bounds), predicate="intersects")
        excluded = np.zeros(values.shape, dtype=bool)
        if len(positions):
            excluded = rasterize(
                [(footprints[int(position)], 1) for position in positions],
                out_shape=values.shape,
                transform=lst.window_transform(window),
                fill=0, dtype="uint8", all_touched=True,
            ).astype(bool)
        valid = (
            np.isfinite(values) & np.isfinite(dem)
            & (impervious == 0) & (~excluded)
        )
        if valid.any():
            bins = np.floor(dem[valid] / ELEVATION_BIN_M).astype(np.int32)
            unique, inverse = np.unique(bins, return_inverse=True)
            block_counts = np.bincount(inverse)
            block_sums = np.bincount(inverse, weights=values[valid])
            for key, count, total in zip(unique, block_counts, block_sums):
                key = int(key)
                counts[key] = counts.get(key, 0) + int(count)
                sums[key] = sums.get(key, 0.0) + float(total)
        if block_number % 1000 == 0:
            print(f"reference blocks: {block_number}", flush=True)
    keys = np.array(sorted(counts), dtype=np.int32)
    table = pd.DataFrame({
        "elevation_m": keys * ELEVATION_BIN_M,
        "sum_lst": [sums[int(key)] for key in keys],
        "n_pixels": [counts[int(key)] for key in keys],
    })
    table["cumulative_sum_lst"] = table.sum_lst.cumsum()
    table["cumulative_n_pixels"] = table.n_pixels.cumsum()
    return table


def reference_mean(table, low, high):
    elevation = table.elevation_m.to_numpy()
    cumulative_sum = table.cumulative_sum_lst.to_numpy()
    cumulative_n = table.cumulative_n_pixels.to_numpy()
    left = int(np.searchsorted(elevation, low, side="left"))
    right = int(np.searchsorted(elevation, high, side="right")) - 1
    if right < left or right < 0 or left >= len(table):
        return np.nan, 0
    total = cumulative_sum[right] - (cumulative_sum[left - 1] if left else 0.0)
    count = cumulative_n[right] - (cumulative_n[left - 1] if left else 0)
    return (float(total / count), int(count)) if count else (np.nan, 0)


def extract_patches(patches, selected_indices, footprints, lst, dem, reference):
    rows, started = [], time.perf_counter()
    for position, source_index in enumerate(selected_indices, start=1):
        patch = patches.loc[source_index]
        footprint = footprints[source_index]
        elev_min, elev_max = footprint_elevation_range(dem, footprint)
        window = clipped_window(lst, footprint)
        values = lst.read(1, window=window, masked=True).filled(np.nan)
        weights = overlap_weights(lst, footprint, window)
        t1, t1_pixels, t1_weight = weighted_mean(values, weights)
        expected = float(weights[weights > 0].sum())
        t0, t0_pixels = reference_mean(reference, elev_min, elev_max)
        slope_code = int(patch.slope_cls)
        rows.append({
            "source_index": int(source_index), "pid": int(patch.pid),
            "area_km2": float(patch.area_km2), "slope_class": slope_code,
            "slope_label": SLOPE_LABELS[slope_code],
            "is_hillside": slope_code >= 1,
            "footprint_area_km2": float(footprint.area / 1_000_000),
            "footprint_elev_min_m": elev_min,
            "footprint_elev_max_m": elev_max,
            "T1_C": t1, "T1_pixels": t1_pixels,
            "T1_valid_weight": t1_weight, "T1_expected_weight": expected,
            "T1_valid_fraction": t1_weight / expected if expected else np.nan,
            "T0_C": t0, "T0_pixels": t0_pixels,
            "dT_C": t1 - t0 if np.isfinite(t1) and np.isfinite(t0) else np.nan,
        })
        if position % 500 == 0 or position == len(selected_indices):
            print(
                f"patches: {position}/{len(selected_indices)} "
                f"({time.perf_counter()-started:.1f}s)", flush=True
            )
    return pd.DataFrame(rows)


def summarize(data, year, season):
    valid = data[np.isfinite(data.dT_C)].copy()
    valid["elev_band"] = pd.cut(
        valid.footprint_elev_min_m,
        bins=[0, 500, 1000, 1500, 2000, 2500, 3000],
        labels=["0-500", "500-1000", "1000-1500", "1500-2000",
                "2000-2500", "2500-3000"],
        right=False,
    )
    rows = []
    for dimension, column in (
        ("slope", "slope_label"), ("elevation", "elev_band")
    ):
        for group, part in valid.dropna(subset=[column]).groupby(
            column, observed=True
        ):
            rows.append({
                "year": year, "season": season, "dimension": dimension,
                "group": str(group), "n": len(part),
                "T1_mean_C": part.T1_C.mean(), "T0_mean_C": part.T0_C.mean(),
                "dT_mean_C": part.dT_C.mean(),
                "dT_median_C": part.dT_C.median(),
                "dT_max_C": part.dT_C.max(),
                "positive_fraction": (part.dT_C > 0).mean(),
                "positive_patch_area_km2":
                    part.loc[part.dT_C > 0, "area_km2"].sum(),
                "mean_T1_valid_fraction": part.T1_valid_fraction.mean(),
            })
    return pd.DataFrame(rows)


def run(year, season, limit=None):
    paths = paths_for(year, season)
    validate(paths)
    patches = gpd.read_file(paths["patches"]).reset_index(drop=True)
    footprints = make_footprints(patches, year)
    if limit:
        selected = np.linspace(
            0, len(patches) - 1, min(limit, len(patches)), dtype=int
        ).tolist()
    else:
        selected = list(range(len(patches)))
    suffix = f"_sample{limit}" if limit else ""
    with rasterio.open(paths["lst"]) as lst, \
         rasterio.open(paths["dem"]) as dem, \
         rasterio.open(paths["impervious"]) as impervious:
        if patches.crs != lst.crs:
            patches = patches.to_crs(lst.crs)
            footprints = make_footprints(patches, year)
        with WarpedVRT(
            dem, crs=lst.crs, transform=lst.transform,
            width=lst.width, height=lst.height,
            resampling=Resampling.average,
        ) as dem_on_lst, WarpedVRT(
            impervious, crs=lst.crs, transform=lst.transform,
            width=lst.width, height=lst.height,
            resampling=Resampling.max,
        ) as is_on_lst:
            reference_path = OUT / f"reference_elevation_{year}_{season}.csv"
            if reference_path.exists():
                reference = pd.read_csv(reference_path)
                print(f"loaded reference cache: {reference_path}", flush=True)
            else:
                reference = reference_histogram(
                    lst, dem_on_lst, is_on_lst, footprints
                )
                reference.to_csv(reference_path, index=False)
            result = extract_patches(
                patches, selected, footprints, lst, dem, reference
            )
    result.insert(0, "season", season)
    result.insert(0, "year", year)
    result.to_csv(
        OUT / f"patch_thermal_{year}_{season}{suffix}.csv", index=False
    )
    summarize(result, year, season).to_csv(
        OUT / f"summary_{year}_{season}{suffix}.csv", index=False
    )
    print(result.dT_C.describe().to_string())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, choices=YEARS, required=True)
    parser.add_argument("--season", choices=SEASONS, required=True)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    run(args.year, args.season, args.limit)


if __name__ == "__main__":
    main()



