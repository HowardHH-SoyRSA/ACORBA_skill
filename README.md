# ACORBA root angles for Codex

Measure Arabidopsis root-tip bending across time from flatbed plate scans. This skill helps Codex read handwritten plate labels, prepare calibrated image stacks, call a local ACORBA scanner installation, review the segmentations, and produce a workbook and publication-style PDFs. ACORBA supplies the signed angles. Each accepted value retains its native frame, component, and review record.

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
- The default run uses **Scanner**, **Deep Machine Learning**, **Super accuracy**, **5-pixel vector construction** and **first-angle normalization off**. Signed values refer to the right-facing measurement orientation. Annotations retain that angle sign when mapped back to a left-facing raw image.
- Review the native segmentation and final skeleton plots. Review sheets link native components to candidate roots. Accept, exclude or leave individual detections unresolved. The statistical unit is an accepted physical root at one time point. Scratches and dust can receive real numeric values from ACORBA and require explicit exclusion. For dense clusters, distinct root identities must be traceable.

## Output

| File | Contents |
| --- | --- |
| `root_angles.xlsx` | Individual native angles with QC and source identifiers in one sheet. The second sheet gives mean, sample size, sample SD and SE by accession and acquisition time. |
| `annotated_angles.pdf` | One combined document for the measured accessions with signed angles overlaid on source-image pixels. |
| `curvature.pdf` | All accessions on one graph of mean curvature versus time, with sample SD error bars and a legend on the right. |

The experiment output also retains the prepared TIFFs, transformations, native analysis files and QC decisions for reinspection. The full procedure and command options are in the [skill instructions](SKILL.md) and the [preparation guide](references/preparation.md).

## Citation and original manual

Nelson B. C. Serre and Matyáš Fendrych (2022). *ACORBA: Automated workflow to measure Arabidopsis thaliana root tip angle dynamics*. **Quantitative Plant Biology** 3, e9. [Paper and supplementary material](https://doi.org/10.1017/qpb.2022.4)

The repository includes Nelson B. C. Serre’s *ACORBA User Manual v1.2* (September 2021). Its SHA-256 matches `S2632882822000042sup001.pdf` in the paper’s supplementary archive. Section 1c of the manual describes the upstream software, methods and scripts under CC BY-NC 2.0. The journal publishes the article under CC BY-NC-SA 4.0. [Included manual (PDF)](references/ACORBA-user-manual-v1.2.pdf) · [Original supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10095971/supplementaryFiles) · [Upstream license described in the manual](https://creativecommons.org/licenses/by-nc/2.0/) · [Article license](https://creativecommons.org/licenses/by-nc-sa/4.0/)
