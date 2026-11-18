"""Prepare reviewed raw plates, split photos, accession stacks, or individual lanes for ACORBA."""

import argparse
import hashlib
import json
import os
import re
import shutil
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks


def catalog_names(path):
    """Read the name / accession column of the first worksheet without a second runtime."""
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(path) as book:
        shared = []
        if "xl/sharedStrings.xml" in book.namelist():
            root = ET.fromstring(book.read("xl/sharedStrings.xml"))
            shared = ["".join(el.itertext()) for el in root]
        sheet = ET.fromstring(book.read("xl/worksheets/sheet1.xml"))
    column = None
    names = set()
    for row in sheet.findall(".//s:row", ns):
        for cell in row.findall("s:c", ns):
            col = re.match(r"[A-Z]+", cell.get("r")).group()
            val = cell.findtext("s:v", "", ns)
            if cell.get("t") == "s":
                val = shared[int(val)]
            if cell.get("t") == "inlineStr":
                val = "".join(cell.find("s:is", ns).itertext())
            if column is None and val.strip().lower() in ("name", "accession", "accession name"):
                column = col
            elif col == column and val.strip():
                names.add(val.strip())
    if not names:
        raise ValueError("Could not locate accession names in the first worksheet")
    return names


def homogeneous(affine):
    m = np.vstack((affine, [0.0, 0.0, 1.0]))
    if m.shape != (3, 3) or not np.isfinite(m).all() or abs(np.linalg.det(m)) < 1e-8:
        raise ValueError("Expected an invertible 2 x 3 pixel-center affine transform")
    return m


def atomic_json(path, data):
    tmp = path.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf8")
    os.replace(tmp, path)


def save_tiff(path, frames, dpi, compression):
    tmp = path.with_suffix(".partial.tif")
    frames[0].save(tmp, save_all=True, append_images=frames[1:], compression=compression, dpi=dpi)
    os.replace(tmp, path)


def safe_name(name):
    if not isinstance(name, str) or not name.strip():
        raise ValueError(
            "Read the plate label views and set a confirmed accession/sample name in the preparation plan"
        )
    if name in (".", "..") or re.search(r'[<>:"/\\|?*\x00-\x1f]', name) or name.endswith((" ", ".")):
        raise ValueError(f"Accession/plate name contains characters unsuitable for filenames: {name!r}")
    return name


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(2**20), b""):
            h.update(part)
    return h.hexdigest()


def get_page(path, page):
    with Image.open(path) as im:
        im.seek(page)
        dpi = im.info.get("dpi")
        if im.mode.startswith("I;16"):
            # A fixed full-range conversion preserves relative intensity across time.
            values = np.asarray(im, dtype=np.uint32)
            return Image.fromarray(((values + 128) // 257).astype(np.uint8)), dpi
        if im.mode in ("I", "F"):
            raise ValueError(
                "Convert signed or floating-point microscopy data with a reviewed, consistent intensity range first"
            )
        return im.copy(), dpi


def measured_dpi(meta, supplied):
    dpi = supplied if supplied is not None else meta
    if isinstance(dpi, (int, float)):
        dpi = [dpi, dpi]
    if dpi:
        dpi = tuple(map(float, dpi))
    if not dpi or len(dpi) != 2 or min(dpi) <= 1 or not np.isfinite(dpi).all():
        raise ValueError("A calibrated DPI or measured pixels-per-cm × 2.54 is required")
    if abs(dpi[0] / dpi[1] - 1) > 0.005:
        raise ValueError("ACORBA needs isotropic pixels; review calibration before analysis")
    return tuple(map(float, dpi))


def add_bar(plate, dpi):
    length = round(dpi[1] / 2.54)
    margin = max(24, round(0.466 * length))
    pad = max(6, min(round(length / 4), plate.height // 16))
    if length + 2 * pad > plate.height:
        raise ValueError("The calibrated 1 cm bar does not fit in this plate image")
    result = Image.new(plate.mode, (plate.width + margin, plate.height), "white")
    result.paste(plate, (0, 0))
    x, y = result.width - margin // 3, result.height - pad - length
    ImageDraw.Draw(result).rectangle([x, y, x + 4, y + length - 1], fill="black")
    return result, [x, y, 5, length]


def verify_bar(im, box, dpi):
    x, y, w, h = box
    if not (0 < x < im.width - w and 0 < y < im.height - h and w == 5 and h == round(dpi[1] / 2.54)):
        raise ValueError("Scale bar geometry does not match the calibration")
    crop = np.array(im.crop((x - 1, y - 1, x + w + 1, y + h + 1)).convert("L"))
    side_white = max(np.min(crop[1:-1, 0]), np.min(crop[1:-1, -1])) >= 246
    if crop[1:-1, 1:-1].max() > 12 or np.min(crop[[0, -1], 1:-1]) < 246 or not side_white:
        raise ValueError("Existing scale bar does not match its declared position")


def find_bar(im, dpi):
    """Identify an exact calibrated vertical bar in either lower side margin."""
    length = round(dpi[1] / 2.54)
    if length < 2 or length >= im.height:
        return None
    top = max(0, im.height // 2 - length - 1)
    width = min(im.width, max(16, im.width // 6))
    candidates = set()
    for left in {0, im.width - width}:
        gray = np.asarray(im.crop((left, top, left + width, im.height)).convert("L"))
        dark = gray <= 12
        # A black surround or scanner edge isn't a 5-pixel bar. Only try
        # columns with a white→black transition and an exact-length black run.
        starts = (gray[: -length - 1] >= 246) & dark[1:-length] & (gray[length + 1 :] >= 246)
        starts &= dark[length:-1]
        y_candidates, x_candidates = np.nonzero(starts)
        for y, x in zip(y_candidates, x_candidates):
            if x + 6 >= width or not np.all(dark[y + 1 : y + 1 + length, x : x + 5]):
                continue
            box = (int(x + left), int(y + 1 + top), 5, length)
            if not dark[y + length // 2 + 1, x + 5] or (x and not dark[y + length // 2 + 1, x - 1]):
                try:
                    verify_bar(im, box, dpi)
                    candidates.add(box)
                except ValueError:
                    pass
    if len(candidates) > 1:
        raise ValueError("Multiple possible scale bars. Specify the reviewed existing_scale rectangle")
    return list(next(iter(candidates))) if candidates else None


def raw_crop(source, matrix, size):
    inv = np.linalg.inv(matrix)
    # Pillow evaluates affine transforms in edge coordinates, whereas registration
    # records use integer coordinates for pixel centers.
    edge = inv.copy()
    edge[:2, 2] += 0.5 - inv[:2, :2] @ [0.5, 0.5]
    return source.transform(
        tuple(size), Image.Transform.AFFINE, edge[:2].ravel(), Image.Resampling.BILINEAR
    )


def same_crop(raw, saved):
    if raw.width > saved.width or raw.height != saved.height:
        raise ValueError("Existing split photo has incompatible dimensions")
    a = np.asarray(raw.convert("L"))[16::37, 17::31].astype(float)
    b = np.asarray(saved.crop((0, 0, raw.width, raw.height)).convert("L"))[16::37, 17::31].astype(float)
    diff = np.abs(a - b)
    if np.mean(diff) > 1.5 or np.quantile(diff, 0.99) > 5:
        raise ValueError("Existing split does not agree with registered source pixels")


def copy_tiff(src, dst):
    temp = dst.with_suffix(".partial.tif")
    if temp.exists():
        temp.unlink()
    shutil.copyfile(src, temp)
    os.replace(temp, dst)


def single_tiff(path):
    with Image.open(path) as image:
        return image.format == "TIFF" and image.n_frames == 1


def geometry(im, plate):
    """Find colored shoots and vote on root direction in the earliest image."""
    scale = max(1, round(im.height / 1100))
    small = im.convert("RGB").resize(
        (max(1, im.width // scale), max(1, im.height // scale)), Image.Resampling.BOX
    )
    r, g, b = np.moveaxis(np.asarray(small, dtype=np.int16), 2, 0)
    green = (g > r + 5) & (g > b + 10) & (r > 45)
    h, w = green.shape
    green[: max(1, h // 85)] = False
    green[-max(1, h // 24) :] = False
    green[:, : max(1, w // 45)] = False
    green[:, -max(1, w // 45) :] = False
    # For a single narrow lane use its full width; here the target is spacing between lanes.
    distance = max(10, round(h * 0.10))
    hist = gaussian_filter1d(green.sum(axis=0).astype(float), max(2, h * 0.013))
    peaks = find_peaks(hist, height=max(0.15, hist.max() * 0.11), distance=distance)[0]
    if "shoot_x" in plate:
        peaks = np.array(plate["shoot_x"]) // scale
    pitch = np.median(np.diff(peaks)) if len(peaks) > 1 else w * 0.75
    lum = 0.299 * r + 0.587 * g + 0.114 * b
    rows, votes = [], []
    for x in peaks:
        rad = max(4, round(0.18 * pitch))
        strip = green[:, max(0, x - rad) : min(w, x + rad + 1)]
        density = gaussian_filter1d(strip.sum(1).astype(float), max(2, h * 0.008))
        ys = find_peaks(density, height=max(0.1, density.max() * 0.045), distance=max(8, h * 0.061))[0]
        rows.append(ys)
        for y in ys:
            sy, sx = np.nonzero(strip[max(0, y - rad) : min(h, y + rad + 1)])
            if not len(sx):
                continue
            cx = round(np.median(sx)) + max(0, x - rad)
            ry = max(5, round(pitch * 0.11))
            sub = lum[max(0, y - ry) : min(h, y + ry + 1)]
            outer, inner = round(0.43 * pitch), round(0.10 * pitch)
            left, right = (
                sub[:, max(0, cx - outer) : max(0, cx - inner)],
                sub[:, min(w, cx + inner) : min(w, cx + outer)],
            )
            if min(left.size, right.size) < 0.8 * max(left.size, right.size):
                continue
            background = np.median(sub[:, max(0, cx - outer) : min(w, cx + outer)])
            l, rr = np.mean(left > background + 16), np.mean(right > background + 16)
            if max(l, rr) > 0.008 and min(l, rr) < 0.7 * max(l, rr):
                votes.append(1 if rr > l else -1)
    score = np.mean(votes) if votes else 0
    direction = plate.get("direction")
    if direction is None and (len(votes) < 3 or abs(score) < 0.22):
        raise ValueError(
            "Uncertain tip orientation; inspect an early frame and set direction to left or right"
        )
    direction = direction or ("right" if score > 0 else "left")
    if direction not in ("left", "right"):
        raise ValueError("direction must be left or right in the supplied image")
    mirrored = direction == "left"
    full_w = im.width
    shoots = peaks[::-1] if mirrored else peaks
    shoots = [full_w - 1 - int(x * scale) if mirrored else int(x * scale) for x in shoots]
    ordered_rows = rows[::-1] if mirrored else rows
    row_template = []
    if ordered_rows:
        sizes = [len(v) for v in ordered_rows]
        # Prefer the modal count of observed planting rows; leave absent shoots unverified.
        count = max(set(sizes), key=lambda n: (sizes.count(n), n))
        if count:
            row_template = np.median([v for v in ordered_rows if len(v) == count], axis=0) * scale
    lane_boxes = plate.get("lane_boxes")
    if lane_boxes is None:
        if not shoots:
            raise ValueError(
                "No shoot columns found; specify lane_boxes and anchor_rows from visual review"
            )
        lane_boxes = [
            [
                max(0, round(x - 0.24 * pitch * scale)),
                round(im.height * 0.014),
                min(full_w, round(x + 0.70 * pitch * scale)),
                round(im.height * 0.948),
            ]
            for x in shoots
        ]
    if len(shoots) != len(lane_boxes) and "anchor_rows" not in plate:
        raise ValueError(
            "Detected and reviewed lane counts disagree; specify shoot_x in original plate coordinates"
        )
    return (
        mirrored,
        lane_boxes,
        shoots,
        list(map(float, row_template)),
        {"votes": len(votes), "right_minus_left_fraction": float(score)},
    )


def clear_shoots(frame, box, shoot_x, row_ys, frame_index):
    x0, y0, x1, y1 = box
    lane = frame.crop(box)
    gray = lane.convert("L")
    if not row_ys:
        return gray, []
    sample = 4
    rgb = np.asarray(
        lane.convert("RGB").resize(
            (max(1, lane.width // sample), max(1, lane.height // sample)), Image.Resampling.BOX
        ),
        dtype=np.int16,
    )
    r, g, b = np.moveaxis(rgb, 2, 0)
    green = (g > r + 5) & (g > b + 10) & (r > 45)
    yy, xx = np.nonzero(green)
    xs = xx * sample + (sample - 1) / 2
    ys = yy * sample + (sample - 1) / 2
    mid = np.median(np.diff(row_ys)) if len(row_ys) > 1 else lane.height
    boundaries = np.round(
        [0] + [(a + b) / 2 - y0 for a, b in zip(row_ys, row_ys[1:])] + [lane.height]
    ).astype(int)
    anchors = []
    for i, ry in enumerate(row_ys):
        use = (np.abs(xs - (shoot_x - x0)) < lane.width * 0.25) & (np.abs(ys - (ry - y0)) < mid * 0.39)
        if use.sum() >= 6:
            # Keep the distal root. Remove shoot and proximal attachment in this planting band.
            cutoff = round(float(np.quantile(xs[use], 0.98))) + max(4, round(lane.width * 0.013))
            center = float(np.median(ys[use]))
        else:
            cutoff = int(shoot_x - x0)
            center = float(ry - y0)
        cutoff = max(0, min(lane.width - 1, cutoff))
        if cutoff:
            ImageDraw.Draw(gray).rectangle(
                (0, int(boundaries[i]), cutoff - 1, int(boundaries[i + 1]) - 1), fill=0
            )
        anchors.append(dict(plant=i + 1, x=cutoff, y=center, observed=bool(use.sum() >= 6)))
    return gray, anchors


def prepare_plate(plate, out, stem):
    name = safe_name(plate["accession"])
    ident = safe_name(plate.get("plate_id", name))
    frames = sorted(plate["frames"], key=lambda f: float(f["time_hours"]))
    times = [float(f["time_hours"]) for f in frames]
    if len(times) != len(set(times)) or not np.isfinite(times).all() or not times:
        raise ValueError(f"Missing or duplicate acquisition times for {ident}")
    stage = plate["stage"]
    if stage not in ("raw", "split", "stack", "lane"):
        raise ValueError(f"Unknown stage {stage}")
    entries, split_paths = [], []
    stack_dir = out / "TIFF stacks"
    lane_dir = out / "ACORBA lanes" / ident
    split_dir = out / "Split plate photos"
    for directory in (split_dir, stack_dir, lane_dir):
        directory.mkdir(parents=True, exist_ok=True)
    first_size = None
    for f, t in zip(frames, times):
        src = Path(f["source"]).resolve()
        page = int(f.get("source_page", 0))
        im, meta = get_page(src, page)
        dpi = measured_dpi(meta, plate.get("dpi"))
        mat = homogeneous(f.get("source_to_plate", np.eye(3)[:2]))
        if not np.allclose(mat[:2, :2].T @ mat[:2, :2], np.eye(2), atol=0.001):
            raise ValueError(
                "Use a rigid raw-to-plate registration to preserve scanner calibration and angles"
            )
        read_origin = src
        original_page = page
        original_stage = stage
        if stage == "raw":
            if "source_to_plate" not in f or "size" not in plate:
                raise ValueError(
                    "Raw input needs a reviewed plate size and per-frame source_to_plate transform"
                )
            aligned = raw_crop(im, mat, plate["size"])
            im.close()
            if f.get("existing_split"):
                candidate = Path(f["existing_split"]).resolve()
                im, existing_dpi = get_page(candidate, 0)
                if not np.allclose(dpi, measured_dpi(existing_dpi, None)):
                    raise ValueError("Registered photo and source calibration differ")
                same_crop(aligned, im)
                bar = f.get("existing_scale", plate.get("existing_scale")) or find_bar(im, dpi)
                if not bar:
                    raise ValueError(
                        "No calibrated 5-pixel scale bar in existing_split; omit existing_split to regenerate it"
                    )
                verify_bar(im, bar, dpi)
            else:
                im, bar = add_bar(aligned, dpi)
                candidate = None
            aligned.close()
        else:
            candidate = None
            bar = f.get("existing_scale", plate.get("existing_scale"))
            if not bar and stage in ("split", "stack"):
                bar = find_bar(im, dpi)
            already_has_bar = bar is not None
            if bar:
                verify_bar(im, bar, dpi)
            elif stage in ("split", "stack"):
                im, bar = add_bar(im, dpi)
            origin = f.get("original")
            if origin:
                read_origin = Path(origin["source"]).resolve()
                original_page = int(origin.get("source_page", 0))
                original_stage = origin.get("source_stage", "raw")
                mat = homogeneous(origin["source_to_plate"])
        if first_size is None:
            first_size = im.size
            first_dpi = dpi
        if im.size != first_size or not np.allclose(first_dpi, dpi):
            raise ValueError(f"Register images to equal geometry and DPI across {ident}")
        dst = split_dir / f"{stem} {t:g}h.tif"
        if stage != "lane":
            if candidate and single_tiff(candidate):
                copy_tiff(candidate, dst)
            elif stage == "split" and already_has_bar and single_tiff(src):
                copy_tiff(src, dst)
            elif stage != "stack" or not already_has_bar:
                save_tiff(dst, [im], dpi, "tiff_lzw")
            else:
                dst = src
            split_paths.append((dst, page if dst == src else 0))
        else:
            split_paths.append((src, page))
        im.close()
        entries.append(
            dict(
                time_hours=t,
                source=str(read_origin),
                source_page=original_page,
                source_stage=original_stage,
                source_to_plate=mat[:2].tolist(),
                scale_bar=bar if stage != "lane" else None,
            )
        )
    if stage != "lane":
        stack = stack_dir / (stem + ".tif")
        existing = plate.get("existing_stack")
        if (
            stage == "stack"
            and len(set(p for p, n in split_paths)) == 1
            and all(n == i for i, (p, n) in enumerate(split_paths))
        ):
            existing = split_paths[0][0]
        if existing:
            with Image.open(existing) as sim:
                if sim.format != "TIFF" or sim.n_frames != len(times):
                    raise ValueError("Existing stack frame count mismatch")
                for i, (p, n) in enumerate(split_paths):
                    reference, _ = get_page(p, n)
                    actual, _ = get_page(existing, i)
                    if (
                        actual.mode != reference.mode
                        or actual.size != reference.size
                        or sha_image(actual) != sha_image(reference)
                    ):
                        raise ValueError(f"Existing stack differs from split photo at {times[i]} h")
                    reference.close()
                    actual.close()
            copy_tiff(existing, stack)
        else:
            pages = [get_page(p, n)[0] for p, n in split_paths]
            save_tiff(stack, pages, first_dpi, "tiff_lzw")
            for p in pages:
                p.close()
    first, _ = get_page(*split_paths[0])
    if stage == "lane":
        if plate.get("direction") not in ("left", "right"):
            raise ValueError("For a preprocessed lane, confirm direction in the supplied image")
        mirror = plate["direction"] == "left"
        boxes = [[0, 0, *first.size]]
        shoots = []
        rows = []
        vote = {}
    else:
        mirror, boxes, shoots, rows, vote = geometry(first, plate)
    first.close()
    lanes = []
    for j, box in enumerate(boxes):
        x0, y0, x1, y1 = box
        if not (0 <= x0 < x1 <= first_size[0] and 0 <= y0 < y1 <= first_size[1]):
            raise ValueError(f"Out-of-bounds lane box {box}")
        number = int(plate.get("lane", j + 1)) if stage == "lane" else j + 1
        dst = lane_dir / f"{ident}-lane-{number:02d}.tif"
        pages = []
        records = []
        for i, ((p, n), entry) in enumerate(zip(split_paths, entries)):
            im, _ = get_page(p, n)
            if mirror:
                im = im.transpose(Image.Transpose.FLIP_LEFT_RIGHT)
            if stage == "lane":
                gray = im.convert("L")
                anchors = frames[i].get("anchors", [])
                if mirror:
                    anchors = [dict(a, x=im.width - 1 - a["x"]) for a in anchors]
            else:
                ys = plate.get("anchor_rows", {}).get(str(number), rows)
                sx = shoots[j] if j < len(shoots) else x0 + (x1 - x0) * 0.24
                gray, anchors = clear_shoots(im, box, sx, ys, i)
            im.close()
            pages.append(gray)
            lane_to_plate = (
                np.array([[-1.0, 0.0, first_size[0] - 1 - x0], [0.0, 1.0, y0], [0.0, 0.0, 1.0]])
                if mirror
                else homogeneous([[1, 0, x0], [0, 1, y0]])
            )
            mapping = np.linalg.inv(homogeneous(entry["source_to_plate"])) @ lane_to_plate
            records.append(dict(entry, lane_to_source=mapping[:2].tolist(), anchors=anchors))
        compliant = False
        if stage == "lane" and not mirror and len(set(p for p, n in split_paths)) == 1:
            src = split_paths[0][0]
            with Image.open(src) as im:
                compliant = im.n_frames == len(times) and all(
                    i == n for i, (p, n) in enumerate(split_paths)
                )
                for i in range(im.n_frames):
                    im.seek(i)
                    compliant &= im.mode == "L" and im.tag_v2.get(259) == 1
        if compliant:
            copy_tiff(src, dst)
        else:
            save_tiff(dst, pages, first_dpi, "raw")
        with Image.open(dst) as saved:
            if saved.n_frames != len(times):
                raise ValueError("Prepared TIFF page count mismatch")
            for i, page in enumerate(pages):
                saved.seek(i)
                if (
                    saved.mode != "L"
                    or saved.tag_v2.get(259) != 1
                    or saved.size != page.size
                    or sha_image(saved) != sha_image(page)
                ):
                    raise ValueError("Prepared 8-bit grayscale TIFF validation failed")
                page.close()
        lane = dict(
            accession=name,
            plate_id=ident,
            lane=number,
            path=str(dst.resolve()),
            frames=records,
            mirrored=mirror,
            orientation_evidence=vote,
            crop_xyxy=box,
            lane_sha256=sha(dst),
        )
        lanes.append(lane)
        print(f"{ident} lane {number}: {dst.name}", flush=True)
    return lanes


def sha_image(image):
    return hashlib.sha256(image.tobytes()).digest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text(encoding="utf-8-sig"))
    names = catalog_names(plan["catalog"]) if plan.get("catalog") else None
    plates = plan["plates"]
    id_groups = {}
    accession_plates = {}
    for p in plates:
        name = safe_name(p["accession"])
        if names and name not in names and p.get("confirmed_custom_name") is not True:
            raise ValueError(
                f"{name!r} missing from master accession name column. Check the label, or confirm a custom sample name explicitly"
            )
        ident = safe_name(p.get("plate_id", name)).casefold()
        accession_plates.setdefault(name.casefold(), set()).add(ident)
        stream = int(p.get("lane", 1)) if p["stage"] == "lane" else None
        occupied = id_groups.setdefault(ident, set())
        if stream in occupied or None in occupied or (stream is None and occupied):
            raise ValueError(
                "Give replicate plates distinct plate_id values and each pre-split lane its own lane number"
            )
        occupied.add(stream)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    # Reuse an entire plate only when its plan, inputs, code, and outputs still match.
    lanes = []
    for p in plates:
        ident = safe_name(p.get("plate_id", p["accession"]))
        stem = p["accession"]
        if len(accession_plates[stem.casefold()]) > 1:
            stem = ident if ident.casefold().startswith(stem.casefold() + " ") else f"{stem} {ident}"
        inputs = {Path(f["source"]).resolve() for f in p["frames"]}
        inputs.update(
            Path(f["existing_split"]).resolve() for f in p["frames"] if f.get("existing_split")
        )
        if p.get("existing_stack"):
            inputs.add(Path(p["existing_stack"]).resolve())
        if any(path.is_relative_to(out) for path in inputs):
            raise ValueError("Choose an output folder separate from your input files")
        signature = hashlib.sha256(
            json.dumps(
                [
                    p,
                    stem,
                    sha(__file__),
                    [(str(f), f.stat().st_size, f.stat().st_mtime_ns) for f in sorted(inputs)],
                ],
                sort_keys=True,
            ).encode()
        ).hexdigest()
        key = f"{ident}.lane-{int(p.get('lane',1)):02d}" if p["stage"] == "lane" else ident
        rec = out / f"{key}.prepared.json"
        cached = json.loads(rec.read_text(encoding="utf8")) if rec.exists() else {}
        if (
            cached.get("signature") == signature
            and all(
                Path(s["path"]).is_file() and sha(s["path"]) == s["lane_sha256"]
                for s in cached["lanes"]
            )
            and all(
                Path(f).is_file() and Path(f).stat().st_size == size
                for f, size in cached.get("intermediates", [])
            )
        ):
            print(f"Reusing verified {ident} lane stacks", flush=True)
            group = cached["lanes"]
        else:
            group = prepare_plate(p, out, stem)
            mids = [] if p["stage"] == "lane" else [out / "TIFF stacks" / (stem + ".tif")]
            mids.extend(
                f
                for f in (out / "Split plate photos").glob("*.tif")
                if f.stem.rsplit(" ", 1)[0] == stem
            )
            atomic_json(
                rec,
                dict(
                    signature=signature,
                    lanes=group,
                    intermediates=[(str(f), f.stat().st_size) for f in mids],
                ),
            )
        lanes.extend(group)
    atomic_json(
        out / "manifest.json",
        dict(
            version=1,
            coordinates="integer pixel centers, zero-based; boxes right/bottom exclusive",
            lanes=lanes,
        ),
    )


if __name__ == "__main__":
    main()
