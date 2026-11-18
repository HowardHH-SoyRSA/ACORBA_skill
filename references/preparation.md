# Preparing scanner images

Run `python scripts/inspect_plates.py --help` for discovery options. Read the label panels against the master workbook and review the time ordering, plate boxes and registration. Handwriting can remain ambiguous even when one catalogue spelling is close. Use another time point or the sample record before assigning it. For a sample with a confirmed name outside the master list, set `confirmed_custom_name` to `true` on that plate. Keep a copy of the reviewed plan with the analysis. Run it with:

```text
python scripts/prepare.py reviewed-plan.json --out prepared
```

### Plan fields

The JSON contains a `plates` array and an optional `catalog` path to the .xlsx file. Give each plate a confirmed `accession` and a stable `plate_id`. Distinct replicate plates of one accession need separate IDs. Supply absolute paths or paths relative to the command's working directory. The physical time in hours is required for each frame. For multi-page TIFFs, `source_page` starts at zero.

| Input stage | Fields and behavior |
|---|---|
| `raw` | Each plate has a registered `size: [width, height]` excluding the scale strip. Each frame supplies `source`, `source_page`, `time_hours`, and a rigid `source_to_plate` matrix of two rows and three columns. For a simple crop starting at (x, y), use `[[1,0,-x],[0,1,-y]]`. A moved plate needs its own matrix at each time. The script crops and registers before adding a calibrated 5-pixel-wide vertical bar to a white strip on the lower right. |
| `split` | Supply the individual accession/time-point images. A calibrated 5-pixel black bar on a white side strip is recognized on either side. To identify a known bar explicitly, supply `existing_scale: [left, top, 5, length]`. A new bar is added if none is present. |
| `stack` | Supply the same stack path with explicit frame pages and times. Already scaled stacks with the desired page order are reused. Dimensions, time labels and calibration must remain consistent across the physical plate. |
| `lane` | Supply the lane number in `lane`, plus the current `direction` (`left` or `right`). For multiple lanes from one plate, repeat the plate ID and use distinct lane numbers. Existing right-facing, uncompressed, 8-bit grayscale stacks are copied byte-for-byte. To enable plant identification, supply per-frame `anchors` in supplied lane coordinates, e.g. `[{"plant":1,"x":95,"y":212,"observed":true}]`. Anchor x is the expected proximal root attachment after clearing the shoot, not the distal tip. |

Mixed starting stages are supported across plates. A partially processed raw plate can point any frame to an `existing_split`, and the plate to an `existing_stack`. The recovered registered pixels, scale bar, DPI and frame order are checked before reuse. For a processed source with raw provenance, add `original` per frame with `source`, `source_page`, `source_stage: "raw"` and `source_to_plate` (raw → supplied split/stack/lane). With only cropped files and no invertible provenance the available input is recorded as the overlay source.

### Orientation and lane review

Raw and split color inputs use the earliest frame for orientation. Green shoots locate seedling columns, and pale tissue extending from each shoot votes on root direction. Confirm the result on the source and prepared lanes, particularly on handwritten, scratched or weakly pigmented plates. The scanner backend needs roots facing right. The program mirrors left-facing plates horizontally and uses fixed lane boxes over the time series. It removes shoots on the proximal side of each planting band using visible green tissue. Photograph pairs that do not register cleanly, rotated layouts with lanes running horizontally, and differently oriented plants within one plate need reviewed crops before a batch.

The `direction` override describes the **supplied** image. `shoot_x` optionally sets original-plate x positions in increasing order. `lane_boxes` optionally sets the left-to-right lane crops **after orientation**, as `[left,top,right,bottom]` with right/bottom exclusive. Coordinates refer to the full split image including the scale strip. For a reviewed row layout, set `anchor_rows` to a dictionary keyed by lane number, e.g. `{"1":[216,810,1400]}`; y coordinates refer to the full plate. For monochrome plates lacking green shoot cues, specify direction, lane boxes, shoot positions and row anchors, or supply a verified prepared lane with manually placed anchors.

Resolution comes from calibrated source metadata or a measured plate-level `dpi` value. For a known ruler interval, multiply measured pixels per centimetre by 2.54. At 1200 dpi, 1 cm is 472 pixels. For JPEG/PNG photographs, verify the physical scale against the ruler before using the metadata. RGB data use Pillow's grayscale conversion; unsigned 16-bit grayscale uses a fixed 0–65535 to 0–255 mapping across frames. Signed or floating-point inputs need a reviewed, common intensity range. Verify plate orientation and identity using the early and late frames. An empty sowing position or a cluster of seedlings needs review during segmentation.

### Outputs

`Split plate photos` contains named frames when needed. `TIFF stacks` contains the time-ordered accession stacks with their scale bars. `ACORBA lanes` contains independently verified, uncompressed, 8-bit grayscale stacks. Prepared lane images omit the side strip. `manifest.json` records the acquisition times, source paths/pages, per-frame shoot anchors, mirroring, lane boxes, calibration bar and inverse transform for each frame. Coordinates in `lane_to_source` map zero-based **pixel centers** in the final prepared image to the corresponding source image. The scale bar is outside the specimen pixels and is never a root measurement.

Use the same output folder to resume. A plate-level cache checks the plan, source identity (size and modification time), preprocessing code and saved lane hashes. Use a different output folder to retain a previous analysis with different settings.
