# Measuring a root-tip angle over 0.200 mm

The upstream point lies 0.200 mm from the optical tip **along the fitted centerline**. The angle is the direction of the straight chord from that point to the tip. This definition fixes the sampled physical length even when root curvature and scanner resolution vary.

The fit starts with the reviewed mask component. Its medial path locates the corresponding region on the original scan through the recorded lane-to-source transform. Across the distal shaft, the fitter models each transverse grayscale profile as a Gaussian-blurred bright bar on a sloping background. A robust fit estimates the center between pixel centers. A cubic smoothing spline joins those centers, and a second pass samples profiles perpendicular to the updated curve. The optical tip lies at the steepest local intensity decrease while moving distally along the curve. The source's pixel size converts the fitted path to millimetres.

At 1200 dpi, 0.200 mm spans approximately 9.45 original pixels of arc. Interpret signed angles relative to the original scan border, after mirroring X for roots that originally grew left. Downward bending is positive. A fit that exceeds the usable centerline or lacks reliable tip evidence cannot supply an accepted angle. Review the overlays, especially for faint tissue, debris beside the cap, and curved or tapered tips. The profile parameters were tested on 1200 dpi flatbed scans. Review fit quality on other image scales.

## Using a different measurement length

Pass `--length-mm` to `measure_subpixel.py` for both the initial fit and the decision replay. Keep one physical length per reported batch so its accession means and time courses have a common definition. Record the length when comparing experiments. ACORBA's native vector setting counts steps on its internal skeleton and is recorded separately.

## Col-0 comparison

The original test compared 35 matched root-time measurements from six seedlings with two verification series. At 0.200 mm, squared Pearson correlations were 0.9742 and 0.9813. Direct agreement R² values were 0.9686 and 0.9784, with mean absolute differences of 4.69° and 3.73°. Direct agreement uses the measured angle without fitting an offset or scale. The first reference scored highest at 0.200 mm across the tested range of 0.100–0.400 mm, and the second at 0.325 mm. Both reference series have undocumented physical segment lengths. The skill uses 0.200 mm as a consistent measurement definition across roots.
