# ACORBA root angles for Codex

Measure Arabidopsis root-tip bending over time from flatbed scans. ACORBA segments the roots. The skill then fits a subpixel centerline to the grayscale intensities in each original scan and measures the direction from a point **0.200 mm back along the centerline** to the optical tip. This gives every accepted root the same physical measurement length. The outputs include individual angles, accession summaries and annotated scans.

## Use in Codex

Install this repository as a user skill named `acorba-root-angles`. From a directory outside an existing checkout, a local Windows installation can use the following commands. Use a new folder for the first installation.

```powershell
git clone https://github.com/HowardHH-SoyRSA/ACORBA_skill.git "$env:USERPROFILE\.codex\skills\acorba-root-angles"
cd "$env:USERPROFILE\.codex\skills\acorba-root-angles"
python -m pip install -r requirements.txt
```

Use Python 3.10 or later. Install ACORBA v1.3 for Windows separately and locate its `scanner.exe`. If Codex uses more than one Python environment, the preparation and review commands need Pillow, NumPy, SciPy and scikit-image. The reporting command can use a separate interpreter with Pillow, ReportLab and openpyxl. The scanner runs its own packaged segmentation model. [Upstream software and models](https://sourceforge.net/projects/acorba/)

In Codex, give the skill the paths to the experiment images, the accession workbook when available, and the output folder. For example:

> Use $acorba-root-angles on the six scans in series 1. Check plate names against my master accession workbook, retain existing verified stacks, and export the reviewed angles and PDFs.

## Inputs and measurements

- Start from raw scans, split accession photos, time stacks or individual lanes. Existing registrations, verified images and calibrated scale bars can be reused. Raw label views expose all four plate edges so Codex can read handwriting. Accession lookup returns possible spellings for visual confirmation.
- Split photos carry a vertical, 5-pixel-wide 1 cm scale bar. Time stacks follow acquisition order. The preparation step orients root growth to the right, crops individual planting lanes and writes single-channel, 8-bit uncompressed TIFFs for the scanner.
- ACORBA's scanner defaults are **Deep Machine Learning**, **Super accuracy**, **5-pixel native vector construction** and **first-angle normalization off**. Either scanner accuracy mode can supply reviewed masks. The physical distance for the reported angle is set by the subpixel measurement step.
- Review each mask against its root, then review the fitted centerline, optical tip and measurement chord on the scan. Only accepted fits enter the statistics. Angles use the image border as horizontal and downward bending as positive after orienting the root to the right. Each measurement retains its native component, audit angle and review record.

Once a batch has a prepared manifest, native results and a reviewed mask table, run:

```text
python scripts/measure_subpixel.py --manifest prepared/manifest.json --native native --detections qc/detections.csv --output subpixel
```

Review the overlays, copy `subpixel/review-template.csv` to `fit-decisions.csv` and enter the decisions. Apply them and export the results:

```text
python scripts/measure_subpixel.py --manifest prepared/manifest.json --native native --detections qc/detections.csv --output subpixel --decisions fit-decisions.csv
python scripts/report.py --manifest prepared/manifest.json --detections subpixel/detections.csv --out reports
```

The default is **0.200 mm**, equivalent to **9.45 original pixels at 1200 dpi**. Use `--length-mm` on both measurement commands to study a different distance. The calibration, optical-tip definition and Col-0 verification results are described in the [measurement notes](references/subpixel.md).

## Output

| File | Contents |
| --- | --- |
| `root_angles.xlsx` | Individual subpixel angles, physical segment length, QC, source identifiers and native ACORBA audit angles. The second sheet gives mean, sample size, sample SD and SE by accession and acquisition time. |
| `annotated_angles.pdf` | Original-scan annotations for the measured accessions, including each accepted angle and its tip-to-upstream chord. |
| `curvature.pdf` | All accessions on one graph of mean curvature versus time, with sample SD error bars and a legend on the right. |

The experiment output also retains the prepared TIFFs, transformations, native analysis files and QC decisions for reinspection. The full procedure and command options are in the [skill instructions](SKILL.md) and the [preparation guide](references/preparation.md).

## Citation and original manual

Nelson B. C. Serre and Matyáš Fendrych (2022). *ACORBA: Automated workflow to measure Arabidopsis thaliana root tip angle dynamics*. **Quantitative Plant Biology** 3, e9. [Paper and supplementary material](https://doi.org/10.1017/qpb.2022.4)

The repository includes Nelson B. C. Serre’s *ACORBA User Manual v1.2* (September 2021). Its SHA-256 matches `S2632882822000042sup001.pdf` in the paper’s supplementary archive. Section 1c of the manual describes the upstream software, methods and scripts under CC BY-NC 2.0. The journal publishes the article under CC BY-NC-SA 4.0. [Included manual (PDF)](references/ACORBA-user-manual-v1.2.pdf) · [Original supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10095971/supplementaryFiles) · [Upstream license described in the manual](https://creativecommons.org/licenses/by-nc/2.0/) · [Article license](https://creativecommons.org/licenses/by-nc-sa/4.0/)

Carsten Steger (1996). *An Unbiased Detector of Curvilinear Structures*. Technical Report FGBV–96–03, TU München. The paper provides background on subpixel localization of line profiles. This skill implements the grayscale profile fit and optical-edge procedure described in its measurement notes. [Technical report (PDF)](https://mv.in.tum.de/_media/members/steger/publications/1996/fgbv-96-03-steger.pdf)
