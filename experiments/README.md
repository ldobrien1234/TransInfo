# experiments

Cross-checks a much cheaper center estimate (a pixel-intensity-weighted
centroid) against the two TI-minimizing centers that `RotRef/TransInfo_images.m`
finds for each image:

1. **rotation center** — the point about which rotational TI is minimized.
2. **reflection center** — the point about which reflectional TI is minimized.
3. **weighted centroid** — `sum(w.*coord)/sum(w)` over pixels inside the
   flower domain (grayscale intensity > `Dthresh`, matching the TI
   function's own domain definition), where `w = weightFun(intensity)`.
   This is orders of magnitude cheaper than (1) and (2), since it doesn't
   require a grid search over candidate transforms — the question this
   folder answers is how good an approximation it is.

There are two variants of the center search, differing only in how TI is
computed:

| script | TI function | non-overlap handling | output folders |
|---|---|---|---|
| `compare_centers_centroid.m` | `transinfo.m` | averages over the **intersection** of the image domain and its transformed domain; rejects candidate centers whose domain coverage is too low (**coverage masking**, below) | `plotsAll/`, `resultsAll/` |
| `compare_centers_centroid_union.m` | `transinfoUnion.m` | pads both images out to the **union** of their domains and sets non-overlapping padding to 1, so there is no shrinking-intersection failure mode and no coverage output to mask on | `plotsUnion/`, `resultsUnion/` |

`summarize_center_comparison.m` is a third script that post-processes the
`center_comparison.csv` written by either of the above into per-metric
summary statistics and comparable histograms.

## `compare_centers_centroid.m` / `compare_centers_centroid_union.m`

```matlab
results = compare_centers_centroid();                    % defaults to the 'angiosperms' folder
results = compare_centers_centroid('angiosperms', 'MaxImages', 5);
results = compare_centers_centroid('angiosperms', 'MaxImages', 5, 'Sample', 'random', 'Seed', 1);
results = compare_centers_centroid({'angiosperms/early_angiosperms/4_nymphaeaceae.png', ...
    'angiosperms/eudicots/some_flower.png'});             % explicit, hand-picked subset
```

`compare_centers_centroid_union.m` takes the same arguments. Both scripts
also accept:

- `'OutputSuffix'` (default `'All'` for `compare_centers_centroid.m`,
  `'Union'` for `compare_centers_centroid_union.m`) names the output folders
  (`plots<suffix>/` and `results<suffix>/`, created if needed), so a run on a
  different image set doesn't overwrite an earlier one.
- `'MaxDim'` (default `Inf`, i.e. no downscaling) caps the larger dimension
  of the image the center searches run on. A bigger image is downscaled by
  `maxDim/max(width,height)` first, and the centers are mapped back to the
  original image's coordinates, so results stay comparable across runs (the
  plots are still drawn over the full-resolution image). The searches cost
  about one pass over the image per candidate transform — 21×21×11 of them
  for rotation — so this is what makes full-resolution photographs
  tractable, at the price of locating the centers on a coarser grid. The
  weighted centroid is computed on the same downscaled image, so all three
  centers see identical data.

Both scripts also resolve a
relative `target` against this folder if it isn't found under the repository
root, so an image set kept alongside the scripts can be named directly:

```matlab
results = compare_centers_centroid('other_images', 'OutputSuffix', 'Other');
% experiments/other_images -> experiments/plotsOther, experiments/resultsOther
results = compare_centers_centroid('other_images', 'MaxDim', 800, ...
    'OutputSuffix', 'Other');                 % search at <=800 px per side
results = compare_centers_centroid_union('other_images', 'MaxDim', 800, ...
    'OutputSuffix', 'OtherUnion');            % same, union-domain TI
results = summarize_center_comparison({'resultsOther', 'resultsOtherUnion'});
```

`other_images` holds full-resolution photographs (several thousand pixels on
a side, against a few hundred for the `angiosperms` set), so a run there
without `'MaxDim'` will take far longer per image.

`target` may be a single image, a folder (searched recursively for
`.png`/`.jpg`, optionally capped/subsampled via `'MaxImages'`/`'Sample'`/`'Seed'`),
or an explicit cell array of image paths. See the function's help
(`help compare_centers_centroid`) for the full option list, including
`'WeightFun'` (default `@(v) v`, i.e. weight by intensity).

For each image, the script:
- finds the rotation center and reflection center via the same kind of grid
  search as `TransInfo_images.m` (for `compare_centers_centroid.m`, also
  rejecting candidate centers whose transformed domain overlaps less than
  `COVthresh` of the original flower's domain — see **Coverage masking**
  below);
- computes the weighted centroid;
- plots all three centers over the image (saved to `plotsAll/` or
  `plotsUnion/`, mirroring the input folder structure, as
  `<image>_centers.png`);
- records the three centers, the image's width/height, and the pairwise
  Euclidean distances between the centers — in pixels, and as a percentage
  of the image's diagonal (`100 * dist / hypot(width,height)`) for a
  size-independent comparison across differently-sized images.

All rows are written to `center_comparison.csv` in the results folder, and a
histogram (pixel and percentage versions, 30 bins) of each of the three
pairwise-distance columns is saved alongside it as `<column>_hist.png`.

### Coverage masking (`compare_centers_centroid.m` only)

`transinfo.m`'s TI value is an *average* divergence over the intersection of
an image's domain and its transformed domain. A candidate center that pushes
most of the flower out of its own footprint can leave only a handful of
unrepresentative pixels behind, whose small average divergence mimics a true
symmetry center. `find_rotation_center`/`find_reflection_center` (local
functions in `compare_centers_centroid.m`) reject candidates below a
coverage threshold (`COVthresh`) using `transinfo.m`'s `coverage` output,
falling back to the unmasked search only if no candidate clears the bar.
`compare_centers_centroid_union.m` has no equivalent step because
`transinfoUnion.m` never shrinks the comparison region.

### `center_comparison.csv` columns

| column | meaning |
|---|---|
| `image` | file name (no extension) |
| `rot_x`, `rot_y` | rotation center |
| `ref_x`, `ref_y` | reflection center |
| `centroid_x`, `centroid_y` | weighted centroid |
| `img_width`, `img_height` | image dimensions (px) |
| `search_scale` | scale factor the center searches ran at (1 unless `'MaxDim'` downscaled the image); the centers above are in original-image coordinates either way |
| `dist_rot_centroid`, `dist_ref_centroid`, `dist_rot_ref` | pairwise distances (px) |
| `dist_rot_centroid_pct`, `dist_ref_centroid_pct`, `dist_rot_ref_pct` | same distances, as % of image diagonal |

## `summarize_center_comparison.m`

```matlab
summarize_center_comparison();              % default: both resultsAll and resultsUnion
summarize_center_comparison('resultsUnion'); % just one folder
summarize_center_comparison({'resultsAll', 'resultsUnion'});
```

`target` is a results folder (or cell array of them) that already contains a
`center_comparison.csv`; a bare name is resolved against this folder.
Default: `{'resultsAll', 'resultsUnion'}`, the two scripts' default output
folders.

For each folder, it derives the rest of that folder's output from
`center_comparison.csv`:

- **`center_comparison_pct_summary.csv`** — one row per `*_pct` metric: `n`,
  mean, median, standard deviation, 90th/95th percentile (linear-interpolation
  / numpy-R7 convention), and the percentage of images falling under 5% / 10%.
  Only the percentage-of-diagonal columns are summarized, since raw pixel
  distances aren't comparable across differently-sized images. Because the
  distances are right-skewed (a handful of images disagree much more than the
  rest), prefer the median and the percent-under-threshold columns over the
  mean when summarizing "how close" — the mean is pulled up by the tail.
- **`dist_rot_centroid_pct_hist.png`**, **`dist_ref_centroid_pct_hist.png`**,
  **`dist_rot_ref_pct_hist.png`** — one histogram per `*_pct` metric, all
  three sharing a common x-axis and y-axis range so their shapes are directly
  comparable (unlike the per-run histograms, which autoscale independently).

It also deletes any stale raw-pixel `*_hist.png` files left over from an
older run, since it no longer produces them.

## Output folders

- **`plotsAll/`**, **`resultsAll/`** — output of a full `compare_centers_centroid.m`
  run (171 images, the `angiosperms` default target). `results_all/` is the
  same run's output under the script's older default folder name.
- **`plotsUnion/`**, **`resultsUnion/`** — the same for `compare_centers_centroid_union.m`.
- **`plotsMasked/`**, **`resultsMasked/`** — output of an earlier, partial
  run (11 images) kept for reference; regenerate rather than relying on it
  as current.
- **`plotsOther/`**, **`resultsOther/`** — the suggested folders for a
  run of either script over `other_images/` (created on demand by
  `'OutputSuffix', 'Other'`; use distinct suffixes if you run both).
- **`other_images/`** — assorted non-`angiosperms` images (photographs and
  processed figures) to try the center comparison on.
- **`scratch_ti_curves/`** — ad-hoc TI-curve plots for individual images,
  not part of the pipeline above.
