---
name: acorba-root-angles
description: Measure Arabidopsis root-tip angles from plate scans using reviewed ACORBA masks and a grayscale subpixel centerline fit over a 0.200 mm tip segment. Prepare TIFF stacks, review fitted tips, and export original-image annotations and curvature statistics.
---

# ACORBA root angles

Use the scripts relative to this skill directory. If `references/workstation.local.json` exists, read the saved paths to the local scanner, accession workbook, Python interpreters, and optional GPU libraries. Otherwise use the supplied workbook and ACORBA location, or look for the master accession workbook in the experiment's parent directories. A single Python environment with the packages in `requirements.txt` can run the workflow. For separate environments, use scientific Python for preparation and QC and one with ReportLab and openpyxl for reporting. Keep source data local, use a separate output folder, and pass the requested input paths when selecting a subset of a larger folder.

## Prepare

1. Check for existing registered split photos, scale bars, accession stacks, lane stacks, and their manifests. Use the furthest verified stage that preserves the required original-image mapping. For raw scans without an existing plan, run `python scripts/inspect_plates.py <scans...> --catalog <master.xlsx> --out <inspection>`. For previously processed inputs use `--stage split`, `stack`, or `lane` as appropriate. Supply `--hours` in page order if a stack lacks a time index. For other scan layouts, supply each plate's crop with a repeated `--box X0 Y0 X1 Y1`.
2. Open the generated label views. Read handwriting in the rotated edge strips and cross-check the physical plate across time. Use `python scripts/inspect_plates.py --catalog <master.xlsx> --lookup <readings...>` for canonical spellings. Resolve ambiguous handwriting from another frame or the experiment record. Fill the reviewed plan's `accession` and `plate_id` (use the accession name for a single plate and distinct IDs for biological replicates). Verify per-frame registration, especially any box-only translation. A moved or swapped plate needs its own mapping. Preserve original paths, page indices, DPI, and affine transforms when resuming. If the original source or mapping is unavailable, identify the earliest available image to the user and label its provenance accurately.
3. Run `python scripts/prepare.py <reviewed-plan.json> --out <prepared>`. Review an early frame for right-facing roots, one planting column per lane, and retention of distal tips throughout the series. The 1 cm scale is vertical, 5 pixels wide and calibrated from Y resolution. Check scanner calibration if the DPI is unknown. Adjust lane boxes or planting rows if the detected geometry disagrees with visible plants. The preparation plan format and supported entry stages are in [Preparation](references/preparation.md).

## Measure and review

4. Run `python scripts/run_acorba.py --manifest <prepared/manifest.json> --scanner <scanner.exe> --output <native>`. This calls ACORBA Scanner with Deep Machine Learning, Super accuracy, vector length 5 pixels, first-angle normalization off, and smoothing disabled. It saves segmentation masks, analysis plots, and native logs. Keep native records alongside reports. For an existing compatible GPU installation, repeat `--gpu-libs <directory>` for the required runtime DLL folders. CPU execution is also supported.
   For Regular mode, pass `--accuracy regular`. The reviewer removes the known ACORBA v1.3 tile padding when mapping saved Regular masks to the scan. The native angle values retain the effects of ACORBA's internal resize.
5. Run `python scripts/review_acorba.py --manifest <prepared/manifest.json> --native <native> --output <qc>`. Inspect the generated frame sheets against the prepared lane and original when needed. Match each native skeleton to its mask, distal tip, and biological root. Use position across frames to disambiguate similar roots. Dust, scratches, plate rims, broken segments and merged seedlings can all produce native angle values. Enter a reason for exclusions and keep unclear assignments unresolved. A candidate match is a proposal until reviewed.
6. Save the completed mask decisions in a separate CSV using the generated template columns and evidence keys. Rerun the same review command with `--decisions <reviewed.csv>`. The resulting `<qc>/detections.csv` identifies accepted mask components. Keep ACORBA's final native angles and plots with the segmentation evidence.

7. **Measure the grayscale subpixel angle, using 0.200 mm by default.** Run the following command with the scientific Python environment:

   ```text
   python scripts/measure_subpixel.py --manifest <prepared/manifest.json> --native <native> --detections <qc/detections.csv> --output <subpixel>
   ```

   The fitter follows the reviewed mask to the original scan, fits the grayscale centerline and locates the optical tip. It measures the direction from a point 0.200 mm back **along the fitted curve** to the tip. Angles use the source-image horizontal, with root growth oriented rightwards and downward bending positive. Scanner vector construction pertains to the native audit angle. For an explicitly requested physical length pass `--length-mm <mm>` and use a distinct output folder.

8. Inspect the fit overlays against the unprocessed image. Verify the biological tip, centerline, and upstream endpoint. Reject debris and spurious edges. Keep weak, truncated or ambiguous tips unresolved. Record decisions in a copy of `<subpixel>/review-template.csv` and repeat step 7 with `--decisions <fit-decisions.csv>`. Only accepted fits contribute angles to the default statistics. Reuse compatible native masks and review decisions for previously processed scans. When starting from an existing stack or lane, use its recorded original-image mapping. If a calibrated processed source is all that is available, identify that limitation and use `--allow-processed-source` only for a measurement intended on those pixels. Supply `--source-dpi` only from verified calibration.

## Report

9. Run `python scripts/report.py --manifest <prepared/manifest.json> --detections <subpixel/detections.csv> --out <reports>`. Verify the individual worksheet against the accepted optical fits, and the per-accession/time summary against accepted root counts. The annotated PDF shows the fitted measurement chord on source pixels. The curvature graph shows means and sample SD. For an explicit request for historic ACORBA vector angles, pass `--angle-source native` to the report command and supply the native QC table from step 6.

10. Return the workbook and both PDFs. State the physical segment length, source-image stage, accepted count and exclusions. Report unassigned planting positions separately. When a native workbook fails to export, verify the final angle lists, masks and plots before using native evidence. Identify any missing frames. See the [measurement method](references/subpixel.md) and the [ACORBA user manual](references/ACORBA-user-manual-v1.2.pdf).
