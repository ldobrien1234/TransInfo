# experiments

Cross-checks a much cheaper center estimate (a pixel-intensity-weighted
centroid) against the two TI-minimizing centers that `RotRef/TransInfo_images.m`
finds for each image:

1. **rotation center** — the point about which rotational TI is minimized.
2. **reflection center** — the point about which reflectional TI is minimized.
3. **weighted centroid** — `sum(w.*coord)/sum(w)` over pixels inside the
   flower domain (grayscale intensity > `Dthresh`, matching `transinfo.m`'s
   own domain definition), where `w = weightFun(intensity)`. This is orders
   of magnitude cheaper than (1) and (2), since it doesn't require a grid
   search over candidate transforms — the question this folder answers is
   how good an approximation it is.

## `compare_centers_centroid.m`

```matlab
results = compare_centers_centroid();                    % defaults to the 'angiosperms' folder
results = compare_centers_centroid('angiosperms', 'MaxImages', 5);
results = compare_centers_centroid('angiosperms', 'MaxImages', 5, 'Sample', 'random', 'Seed', 1);
results = compare_centers_centroid({'angiosperms/early_angiosperms/4_nymphaeaceae.png', ...
    'angiosperms/eudicots/some_flower.png'});             % explicit, hand-picked subset
```

`target` may be a single image, a folder (searched recursively for
`.png`/`.jpg`, optionally capped/subsampled via `'MaxImages'`/`'Sample'`/`'Seed'`),
or an explicit cell array of image paths. See the function's help
(`help compare_centers_centroid`) for the full option list, including
`'WeightFun'` (default `@(v) v`, i.e. weight by intensity).

For each image, the script:
- finds the rotation center and reflection center via the same kind of grid
  search as `TransInfo_images.m`, rejecting candidate centers whose
  transformed domain overlaps less than `COVthresh` of the original
  flower's domain (see **Coverage masking** below);
- computes the weighted centroid;
- plots all three centers over the image (saved to `plotsAll/`, mirroring
  the input folder structure, as `<image>_centers.png`);
- records the three centers, the image's width/height, and the pairwise
  Euclidean distances between the centers — in pixels, and as a percentage
  of the image's diagonal (`100 * dist / hypot(width,height)`) for a
  size-independent comparison across differently-sized images.

All rows are written to `results_all/center_comparison.csv`, and a
histogram (pixel and percentage versions, 30 bins) of each of the three
pairwise-distance columns is saved alongside it as
`<column>_hist.png`.

### Coverage masking

`transinfo.m`'s TI value is an *average* divergence over the intersection of
an image's domain and its transformed domain. A candidate center that pushes
most of the flower out of its own footprint can leave only a handful of
unrepresentative pixels behind, whose small average divergence mimics a true
symmetry center. `find_rotation_center`/`find_reflection_center` (local
functions in `compare_centers_centroid.m`) reject candidates below a
coverage threshold (`COVthresh`) using `transinfo.m`'s `coverage` output,
falling back to the unmasked search only if no candidate clears the bar.

### `center_comparison.csv` columns

| column | meaning |
|---|---|
| `image` | file name (no extension) |
| `rot_x`, `rot_y` | rotation center |
| `ref_x`, `ref_y` | reflection center |
| `centroid_x`, `centroid_y` | weighted centroid |
| `img_width`, `img_height` | image dimensions (px) |
| `dist_rot_centroid`, `dist_ref_centroid`, `dist_rot_ref` | pairwise distances (px) |
| `dist_rot_centroid_pct`, `dist_ref_centroid_pct`, `dist_rot_ref_pct` | same distances, as % of image diagonal |

### `center_comparison_pct_summary.csv`

A one-row-per-metric summary of the three `*_pct` columns above (n, mean,
median, standard deviation, 90th/95th percentile, and the percentage of
images falling under 5%/10%) — useful for judging how closely the cheap
centroid tracks the TI-minimizing centers without eyeballing histograms.
Because the underlying distances are right-skewed (a handful of images have
much larger disagreement than the rest), prefer the median and the
percent-under-threshold columns over the mean when summarizing "how close":
the mean is pulled up disproportionately by the tail. This file is a
derived summary of `center_comparison.csv`, not regenerated automatically by
`compare_centers_centroid.m`.

## Output folders

- **`plotsAll/`**, **`results_all/`** — output of a full run (171 images,
  the `angiosperms` default target) with the current version of the script.
- **`plotsMasked/`**, **`resultsMasked/`** — output of an earlier, partial
  run (11 images) kept for reference; regenerate via
  `compare_centers_centroid.m` if needed rather than relying on it as
  current.
