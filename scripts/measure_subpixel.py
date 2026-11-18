"""Measure signed tip chords on original grayscale scans using reviewed ACORBA mask identities."""
import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

from review_acorba import FIELDS, KEY_FIELDS, font, frame_evidence_key, mask_crop, ordered_path, row_key
from run_acorba import evidence, file_hash, final_lists, lane_key, write_json
from subpixel_fit import DEFAULT_LENGTH_MM, fit_tip


EXTRA_FIELDS = [
    "measurement_method", "length_mm", "source_stage", "dpi_x", "dpi_y", "dpi_source",
    "mirror_source_x", "native_angle_deg", "native_qc_status", "native_qc_reason",
    "fit_angle_deg", "fit_status", "fit_reason", "fit_evidence_key", "review_sheet",
    "tip_source_x", "tip_source_y", "upstream_source_x", "upstream_source_y",
    "upstream_x", "upstream_y", "fit_diagnostics",
]
DECISION_FIELDS = KEY_FIELDS + [
    "root_id", "candidate_id", "fit_evidence_key", "fit_status", "fit_reason",
    "fit_angle_deg", "length_mm", "review_sheet",
]


def load_rows(path, required):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        return list(reader)


def write_rows(path, rows, fields):
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False,
                                     encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(stream.name, path)


def digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def reviewed_decisions(path):
    decisions = {}
    if path is None:
        return decisions
    for decision in load_rows(path, DECISION_FIELDS[:10]):
        key = row_key(decision)
        if key in decisions:
            raise ValueError(f"Duplicate fit decision for {key}")
        decisions[key] = decision
    return decisions


def physical_resolution(source, frame, override):
    if override is not None:
        resolution, provenance = (override, override), "CLI source calibration"
    elif frame.get("source_dpi") is not None:
        resolution, provenance = frame["source_dpi"], "manifest source calibration"
    else:
        resolution, provenance = source.info.get("dpi"), "source image metadata"
    dpi = np.asarray(resolution, dtype=float)
    if dpi.ndim == 0:
        dpi = np.repeat(dpi, 2)
    if dpi.shape != (2,) or not np.isfinite(dpi).all() or min(dpi) <= 1:
        raise ValueError("Source requires calibrated X/Y DPI or a measured --source-dpi")
    return dpi.tolist(), provenance


def frame_transform(frame):
    affine = np.asarray(frame["lane_to_source"], dtype=float)
    if affine.shape != (2, 3) or not np.isfinite(affine).all():
        raise ValueError("Invalid lane-to-source transform")
    linear = affine[:, :2]
    if abs(np.linalg.det(linear)) < 1e-10:
        raise ValueError("Singular lane-to-source transform")
    if abs(linear[0, 0]) < 1e-8 * np.linalg.norm(linear[:, 0]):
        raise ValueError("Right-facing lane maps vertically in the source; verify source orientation")
    return affine, bool(linear[0, 0] < 0)


def source_patch(source, labels, boxes, candidate_id, affine):
    match = re.fullmatch(r"M([1-9][0-9]*)", candidate_id)
    if not match or int(match[1]) > len(boxes):
        raise ValueError(f"Reviewed component {candidate_id!r} is absent in this mask")
    index = int(match[1])
    box = boxes[index - 1]
    if box is None:
        raise ValueError(f"Reviewed component {candidate_id} is empty")
    seed = ordered_path(labels[box] == index)
    if seed is None or len(seed) < 4:
        raise ValueError("Reviewed component has insufficient continuous seed path")
    seed = seed + [box[1].start, box[0].start]
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(seed, axis=0), axis=1))]
    seed = seed[arc >= arc[-1] - 130]
    seed = seed @ affine[:, :2].T + affine[:, 2]
    origin = np.maximum(0, np.floor(seed.min(axis=0) - 45).astype(int))
    end = np.minimum(source.size, np.ceil(seed.max(axis=0) + 45).astype(int))
    if np.any(end <= origin):
        raise ValueError("Reviewed mask seed is outside the declared source image")
    crop = source.crop((*origin, *end))
    if crop.mode.startswith("I;16"):
        values = np.asarray(crop, dtype=np.uint32)
        gray = ((values + 128) // 257).astype(np.uint8)
    elif crop.mode in ("I", "F"):
        raise ValueError("Source needs a reviewed grayscale conversion for signed or floating intensities")
    else:
        gray = np.asarray(crop.convert("L"))
    return gray, seed - origin, origin


def native_evidence(lane, native, manifest_dir):
    meta_path = native / lane_key(lane) / "run.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if meta["status"] not in {"complete", "native_export_failed"}:
        raise ValueError(f"Reviewed detections refer to an incomplete native run: {meta_path}")
    attempt = meta_path.parent / meta["attempt"]
    files = evidence(attempt, len(lane["frames"]))
    signature = meta["signature"]
    prepared = manifest_dir / lane["path"]
    if signature["input_sha256"] != file_hash(prepared):
        raise ValueError(f"Prepared lane differs from reviewed native input: {prepared}")
    hashes = {}
    for path in files:
        hashes[path] = file_hash(path)
        if hashes[path] != meta["sha256"][path.relative_to(attempt).as_posix()]:
            raise ValueError(f"Changed native evidence: {path}")
    with Image.open(prepared) as image:
        size = image.size
    return files[1], final_lists(files[0]), attempt, hashes, signature, size


def cache_fit(cache_dir, key, gray, seed, length_mm, dpi, mirror):
    path = cache_dir / (key + ".json")
    if path.is_file():
        saved = json.loads(path.read_text(encoding="utf-8"))
        if saved["input_key"] != key or saved["fit_hash"] != digest(saved["fit"]):
            raise ValueError(f"Damaged fit cache: {path}")
        return saved["fit"], saved["fit_hash"]
    fitted = fit_tip(gray, seed, length_mm=length_mm, dpi=dpi, mirror_x=mirror)
    fit = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in fitted.items()}
    fit_hash = digest(fit)
    write_json(path, {"input_key": key, "fit_hash": fit_hash, "fit": fit})
    return fit, fit_hash


def review_panel(gray, fit, row):
    # Both views show identical original pixels at the same zoom. The plain
    # image keeps faint tips visible beside the marked centerline and chord.
    tip = np.array(fit["tip_xy"])
    upstream = np.array(fit["upstream_xy"])
    radius = max(32, math.ceil(np.linalg.norm(tip - upstream) * 1.35))
    start = np.floor(tip - radius).astype(int)
    side = 2 * radius + 1
    image = Image.fromarray(gray).crop((*start, *(start + side))).resize(
        (210, 210), Image.Resampling.NEAREST).convert("RGB")
    marked = image.copy()
    draw = ImageDraw.Draw(marked)

    def project(p):
        return tuple((np.asarray(p) - start + 0.5) * 210 / side)

    draw.line([project(p) for p in fit["curve_xy"]], fill="#1bebff", width=2)
    draw.line([project(tip), project(upstream)], fill="#ffe653", width=2)
    for point, color in [(tip, "#ff4949"), (upstream, "#ffe653")]:
        x, y = project(point)
        draw.ellipse((x-2, y-2, x+2, y+2), fill=color)
    panel = Image.new("RGB", (426, 240), "white")
    panel.paste(image, (0, 30))
    panel.paste(marked, (216, 30))
    title = f"Root {row['root_id']}   {row['candidate_id']}   #{row['native_index']}   {fit['angle_deg']:+.2f} deg"
    ImageDraw.Draw(panel).text((2, 6), title, fill="black", font=font(14))
    return panel


def store_sheet(path, panels, lane, frame, length):
    width = min(2, len(panels)) * 432
    height = math.ceil(len(panels) / 2) * 245 + 58
    sheet = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(sheet)
    draw.text((4, 4), f"{lane['accession']} / {lane['plate_id']} / lane {lane['lane']} / {frame['time_hours']:g} h",
              fill="black", font=font(15))
    draw.text((4, 28), f"{length:g} mm | Red: optical tip | Yellow: upstream | Cyan: fitted centerline",
              fill="black", font=font(12))
    for index, panel in enumerate(panels):
        sheet.paste(panel, ((index % 2)*432, 55+(index // 2)*245))
    sheet.save(path)


def initialize_row(row, length, frame):
    return dict(row, measurement_method="subpixel_centerline", length_mm=length,
                source_stage=frame["source_stage"], native_angle_deg=row["angle_deg"],
                native_qc_status=row["qc_status"], native_qc_reason=row["qc_reason"],
                angle_deg="", tip_x="", tip_y="", anchor_method="", fit_status="not_measured",
                qc_status="unresolved" if row["qc_status"] == "accepted" else row["qc_status"])


def apply_fit_decision(row, decision):
    if decision is None:
        return
    for field in ("fit_evidence_key", "root_id", "candidate_id"):
        if row[field] != decision[field]:
            raise ValueError(f"Fit decision belongs to different evidence or root: {row_key(row)}")
    status = decision["fit_status"].strip()
    if status not in {"accepted", "excluded", "unresolved", "pending"}:
        raise ValueError(f"Invalid fit_status {status!r} for {row_key(row)}")
    reason = decision["fit_reason"].strip()
    if status != "pending" and not reason:
        raise ValueError(f"A fit review reason is required for {row_key(row)}")
    row.update(fit_status=status, fit_reason=reason, qc_reason=reason,
               qc_status=status if status != "pending" else "unresolved")
    if status == "accepted":
        row["angle_deg"] = row["fit_angle_deg"]


def fail(row, message):
    row.update(fit_status="failed", qc_status="unresolved", fit_reason=str(message), qc_reason=str(message))


def measure_frame(lane, index, rows, masks, crop, native_data, options):
    _, records, attempt, hashes, signature, size = native_data
    mask = native_data[0]
    frame = lane["frames"][index]
    accepted = [r for r in rows if r["native_qc_status"] == "accepted"]
    if not accepted:
        return
    native = records[index]
    suffix = "_recalculated" if native["stage"] == "recalculated" else ""
    plot = attempt / f"lane.tif_plateplot_{index+1}{suffix}.png"
    native_key = frame_evidence_key(signature["input_sha256"], hashes[plot], hashes[mask],
                                    native["angles"], frame, crop if masks.size != size else None)
    for row in accepted:
        j = int(row["native_index"]) - 1
        if (j < 0 or j >= len(native["angles"]) or row["evidence_key"] != native_key
                or native["angles"][j] is None or float(row["native_angle_deg"]) != native["angles"][j]):
            raise ValueError(f"Reviewed native evidence does not match this manifest and run: {row_key(row)}")
    masks.seek(index)
    labels, _ = ndimage.label(np.asarray(masks.crop(crop)) != 0, structure=np.ones((3, 3), dtype="u1"))
    boxes = ndimage.find_objects(labels)
    panels = []
    folder = options.output / "review" / lane_key(lane)
    folder.mkdir(parents=True, exist_ok=True)
    sheet = folder / f"frame-{index:03}.png"
    with ExitStack() as opened:
        try:
            if frame["source_stage"] != "raw" and not options.allow_processed_source:
                raise ValueError(f"Source stage is {frame['source_stage']}; supply original mapping or --allow-processed-source")
            affine, mirror = frame_transform(frame)
            source_path = options.manifest.parent / frame["source"]
            source = opened.enter_context(Image.open(source_path))
            source.seek(int(frame["source_page"]))
            dpi, dpi_source = physical_resolution(source, frame, options.source_dpi)
        except (OSError, EOFError, ValueError) as exc:
            for row in accepted:
                fail(row, exc)
            return
        for row in accepted:
            row.update(dpi_x=dpi[0], dpi_y=dpi[1], dpi_source=dpi_source,
                       mirror_source_x=str(mirror).lower())
            try:
                gray, seed, origin = source_patch(source, labels, boxes, row["candidate_id"], affine)
                key = digest([row_key(row), row["candidate_id"], row["root_id"], native_key,
                              options.algorithm_hash, options.length_mm, dpi, mirror, origin.tolist(),
                              gray.shape, hashlib.sha256(gray.tobytes()).hexdigest(),
                              hashlib.sha256(seed.tobytes()).hexdigest()])
                fit, fit_hash = cache_fit(options.output / "fits", key, gray, seed,
                                          options.length_mm, dpi, mirror)
            except (ValueError, np.linalg.LinAlgError) as exc:
                fail(row, exc)
                continue
            tip = origin + np.array(fit["tip_xy"])
            upstream = origin + np.array(fit["upstream_xy"])
            lane_tip = np.linalg.solve(affine[:, :2], tip-affine[:, 2])
            lane_upstream = np.linalg.solve(affine[:, :2], upstream-affine[:, 2])
            row.update(fit_angle_deg=fit["angle_deg"], fit_status="pending",
                       fit_reason="Inspect original grayscale, fitted centerline and optical tip",
                       qc_reason="Subpixel fit awaits visual review", fit_evidence_key=digest([key, fit_hash]),
                       tip_x=lane_tip[0], tip_y=lane_tip[1], upstream_x=lane_upstream[0], upstream_y=lane_upstream[1],
                       tip_source_x=tip[0], tip_source_y=tip[1], upstream_source_x=upstream[0], upstream_source_y=upstream[1],
                       anchor_method="Original grayscale optical tip", review_sheet=str(sheet),
                       fit_diagnostics=json.dumps(fit["diagnostics"], sort_keys=True))
            panels.append(review_panel(gray, fit, row))
    if panels:
        store_sheet(sheet, panels, lane, frame, options.length_mm)


def run(options):
    options.manifest = options.manifest.resolve()
    options.output = options.output.resolve()
    options.output.mkdir(parents=True, exist_ok=True)
    (options.output / "fits").mkdir(exist_ok=True)
    destinations = {options.output / name for name in ("detections.csv", "review-template.csv")}
    if any(p and p.resolve() in destinations for p in (options.detections, options.decisions)):
        raise ValueError("Use a separate output folder and a separate reviewed copy of the decision template")
    if not math.isfinite(options.length_mm) or options.length_mm <= 0:
        raise ValueError("--length-mm must be finite and positive")
    if options.source_dpi is not None and (not math.isfinite(options.source_dpi) or options.source_dpi <= 1):
        raise ValueError("--source-dpi must be a calibrated source resolution above 1 DPI")
    scripts = Path(__file__).parent
    options.algorithm_hash = digest([file_hash(scripts / name) for name in
                                     ("subpixel_fit.py", "measure_subpixel.py", "review_acorba.py")])
    lanes = json.loads(options.manifest.read_text(encoding="utf-8-sig"))["lanes"]
    selections = defaultdict(list)
    accepted_roots, accepted_components, native_ids = set(), set(), set()
    for row in load_rows(options.detections, FIELDS):
        if "measurement_method" in row:
            raise ValueError("Supply native reviewed detections as the input to subpixel measurement")
        key = row_key(row)
        if row["qc_status"] not in {"accepted", "unresolved", "excluded"}:
            raise ValueError(f"Unknown native QC state: {key}")
        if row["native_index"]:
            if key in native_ids:
                raise ValueError(f"Duplicate native detection: {key}")
            native_ids.add(key)
        if row["qc_status"] == "accepted":
            root, component = (*key[:4], row["root_id"]), (*key[:4], row["candidate_id"])
            if (not row["root_id"].strip() or not row["native_index"] or not row["candidate_id"]
                    or root in accepted_roots or component in accepted_components):
                raise ValueError(f"Incomplete or duplicate reviewed plant/mask: {key}")
            accepted_roots.add(root)
            accepted_components.add(component)
        selections[key[:4]].append(row)
    decisions = reviewed_decisions(options.decisions)
    output_rows = []
    for lane in lanes:
        key = (str(lane["accession"]), str(lane["plate_id"]), str(lane["lane"]))
        entries = []
        for index, frame in enumerate(lane["frames"]):
            selected = selections.pop((*key, str(index)), [])
            if not selected or any(float(r["time_h"]) != float(frame["time_hours"]) for r in selected):
                raise ValueError(f"Reviewed frame absent or time differs from manifest: {key}, frame {index}")
            entries.append([initialize_row(r, options.length_mm, frame) for r in selected])
        if any(r["native_qc_status"] == "accepted" for batch in entries for r in batch):
            data = native_evidence(lane, options.native, options.manifest.parent)
            with Image.open(data[0]) as masks:
                crop = mask_crop(masks.size, data[-1], data[-2]["options"])
                for index, batch in enumerate(entries):
                    measure_frame(lane, index, batch, masks, crop, data, options)
        output_rows.extend(r for batch in entries for r in batch)
    if selections:
        raise ValueError(f"Detections outside manifest: {next(iter(selections))}")
    template = []
    for row in output_rows:
        if row["fit_status"] == "pending":
            template.append(row.copy())
            apply_fit_decision(row, decisions.pop(row_key(row), None))
    if decisions:
        raise ValueError(f"Fit decisions have no matching successful candidate: {next(iter(decisions))}")
    write_rows(options.output / "review-template.csv", template, DECISION_FIELDS)
    write_rows(options.output / "detections.csv", output_rows, FIELDS + EXTRA_FIELDS)
    summary = {"method": "subpixel_centerline", "length_mm": options.length_mm,
               "fit_status_counts": dict(Counter(r["fit_status"] for r in output_rows)),
               "qc_status_counts": dict(Counter(r["qc_status"] for r in output_rows))}
    write_json(options.output / "measurement-summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True, help="ACORBA native run directory")
    parser.add_argument("--detections", type=Path, required=True, help="Reviewed native detections.csv")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--length-mm", type=float, default=DEFAULT_LENGTH_MM, help="Physical arc distance from the optical tip (default 0.200 mm)")
    parser.add_argument("--decisions", type=Path, help="Reviewed copy of the subpixel review-template.csv")
    parser.add_argument("--source-dpi", type=float, help="Measured isotropic source DPI override")
    parser.add_argument("--allow-processed-source", action="store_true", help="Fit the manifest's earliest available intensity image when raw mapping is unavailable")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
