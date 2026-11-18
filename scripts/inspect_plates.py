"""Locate plates, expose their labels for reading, and draft a registered preparation plan."""

import argparse
import csv
from difflib import get_close_matches
import inspect
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageOps
from scipy import ndimage
from skimage.feature import ORB, match_descriptors
from skimage.measure import ransac
from skimage.transform import EuclideanTransform

from prepare import catalog_names, find_bar, measured_dpi


TIME = re.compile(r"(?<![^\W_])(?<!\.)(\d+(?:\.\d+)?)\s*(?:hours?|hrs?|h)(?![a-z])", re.I)
IMAGE_SUFFIXES = {".tif", ".tiff", ".jpg", ".jpeg", ".png"}


def list_images(inputs):
    files = []
    for item in inputs:
        paths = sorted(item.iterdir()) if item.is_dir() else [item]
        files.extend(p.resolve() for p in paths if p.is_file() and p.suffix.lower() in IMAGE_SUFFIXES)
    if not files or len(files) != len(set(files)):
        raise ValueError("Supply a nonempty set of distinct input images at one processing stage")
    return files


def preview(image, width, height):
    img = image.convert("RGB")
    img.thumbnail((width, height), Image.Resampling.LANCZOS)
    return img


def exposure_pages(image):
    pages = []
    for index in range(getattr(image, "n_frames", 1)):
        image.seek(index)
        pages.append((index, image.size, bool(getattr(image, "tag_v2", {}).get(254, 0) & 1)))
    width, height = max((size for _, size, _ in pages), key=lambda s: s[0] * s[1])
    real = [i for i, (w, h), reduced in pages if not reduced and (w > width / 8 or h > height / 8)]
    if not real:
        raise ValueError("The TIFF contains only reduced-resolution pages")
    image.seek(real[0])
    return real


def plate_boxes(image):
    """Candidate rectangular plate interiors on dark scanner beds, including touching plates."""
    small = preview(image, 900, 900)
    a = np.asarray(small, dtype="int16")
    r, g, b = a.transpose(2, 0, 1)
    # Agar in these transmitted-light scans is blue-grey. Other illumination can use --box.
    agar = (b > r + 8) & (g > r + 2) & (r > 18)
    clean = ndimage.binary_opening(agar, structure=np.ones((5, 5), bool))
    spans = ndimage.find_objects(ndimage.label(clean)[0])
    w, h = small.size
    sx, sy = image.width / w, image.height / h
    boxes = []
    for ys, xs in spans:
        bh, bw = ys.stop - ys.start, xs.stop - xs.start
        if bh < h * 0.21 or bw < w * 0.12 or bw * bh < w * h * 0.037:
            continue
        across = max(1, round(bw / bh))
        down = max(1, round(bh / bw))
        for row in range(down):
            for col in range(across):
                boxes.append(
                    [
                        max(0, round((xs.start + col * bw / across - 1) * sx)),
                        max(0, round((ys.start + row * bh / down - 1) * sy)),
                        min(image.width, round((xs.start + (col + 1) * bw / across + 1) * sx)),
                        min(image.height, round((ys.start + (row + 1) * bh / down + 1) * sy)),
                    ]
                )
    return sorted(boxes, key=lambda box: (round(box[1] / (h * sy * 0.12)), box[0]))


def label_sheet(image, box, path, title):
    """Keep a full context image and turn all four writing edges horizontal."""
    plate = image.crop(box)
    w, h = plate.size
    tiles = [
        ("Top", plate.crop((0, 0, w, max(1, h // 5)))),
        ("Right", plate.crop((w * 4 // 5, 0, w, h)).transpose(Image.Transpose.ROTATE_90)),
        ("Bottom", plate.crop((0, h * 4 // 5, w, h))),
        ("Left", plate.crop((0, 0, max(1, w // 5), h)).transpose(Image.Transpose.ROTATE_270)),
    ]
    sheet = Image.new("RGB", (1420, 740), "white")
    d = ImageDraw.Draw(sheet)
    d.text((8, 6), title.encode("ascii", "replace").decode(), fill="black")
    whole = preview(plate, 510, 650)
    sheet.paste(whole, (5, 56))
    y = 30
    for heading, tile in tiles:
        d.text((522, y), heading, fill="black")
        tile = ImageOps.autocontrast(preview(tile, 884, 149), cutoff=0.5)
        sheet.paste(tile, (522, y + 15))
        y += 173
    sheet.save(path, quality=91)


def features(im, box):
    crop = im.crop(box).convert("L")
    scale = min(1, 900 / max(crop.size))
    resized = crop.resize(tuple(max(1, round(n * scale)) for n in crop.size), Image.Resampling.LANCZOS)
    data = np.asarray(ImageOps.autocontrast(resized, cutoff=0.2), dtype=float) / 255
    orb = ORB(n_keypoints=850, fast_threshold=0.035)
    try:
        orb.detect_and_extract(data)
    except RuntimeError:
        return None
    points = (orb.keypoints[:, ::-1] + 0.5) * (np.array(crop.size) / resized.size) - 0.5
    return points, orb.descriptors


def register(reference, current, box):
    x, y = box[:2]
    translate = np.array([[1.0, 0.0, -x], [0.0, 1.0, -y], [0.0, 0.0, 1.0]])
    if reference is not None and current is not None:
        matches = match_descriptors(current[1], reference[1], max_ratio=0.80, cross_check=True)
        if len(matches) >= 8:
            source, target = current[0][matches[:, 0]], reference[0][matches[:, 1]]
            kw = {"rng" if "rng" in inspect.signature(ransac).parameters else "random_state": 0}
            fit, keep = ransac(
                (source, target),
                EuclideanTransform,
                min_samples=3,
                residual_threshold=8.0,
                max_trials=500,
                **kw,
            )
            if fit is not None and np.count_nonzero(keep) >= 8:
                span = np.ptp(source[keep], axis=0)
                well_spread = np.all(span > 0.25 * (np.array(box[2:]) - box[:2]))
                if well_spread and abs(fit.rotation) < math.radians(20):
                    return (fit.params @ translate)[:2].tolist(), {
                        "method": "rigid",
                        "inliers": int(keep.sum()),
                        "median_residual_px": round(
                            float(np.median(fit.residuals(source[keep], target[keep]))), 2
                        ),
                    }
    return translate[:2].tolist(), {
        "method": "box translation only; verify alignment or edit transform"
    }


def existing_bar(img):
    """Find the calibrated 5-pixel black scale on either margin (including mirrored stacks)."""
    try:
        dpi = measured_dpi(img.info.get("dpi"), None)
    except ValueError:
        return None
    return find_bar(img, dpi)


def lookup_labels(catalog, readings):
    names = sorted(catalog_names(catalog))

    def simplified(text):
        return "".join(c for c in text.casefold() if c.isalnum())

    for reading in readings:
        simple = simplified(reading)
        keys = {simplified(name) for name in names}
        near = get_close_matches(simple, sorted(keys), n=5, cutoff=0.52)
        ranked = [name for key in near for name in names if simplified(name) == key]
        print(
            json.dumps(
                {"reading": reading, "exact": reading in names, "candidates": ranked},
                ensure_ascii=False,
            )
        )


def hours_of_stack(path, count, provided):
    if provided:
        hours = provided
    else:
        index = path.parent / "stack_index.csv"
        hours = []
        if index.is_file():
            with index.open(encoding="utf-8-sig", newline="") as f:
                rows = sorted(
                    (r for r in csv.DictReader(f) if r["stack_file"] == path.name),
                    key=lambda r: int(r["page_number"]),
                )
            if [int(r["page_number"]) for r in rows] == list(range(1, count + 1)):
                hours = [float(r["time_hours"]) for r in rows]
    if len(hours) != count or len(set(hours)) != count:
        raise ValueError(f"Supply --hours in page order for {path.name} ({count} pages)")
    return hours


def drafts(args):
    files = list_images(args.inputs)
    names = catalog_names(args.catalog) if args.catalog else set()
    groups = defaultdict(list)
    ready = []
    for path in files:
        with Image.open(path) as im:
            pages = exposure_pages(im)
            match = TIME.search(path.stem)
            prefix = path.stem[: match.start()].rstrip(" _-") if match else path.stem
            lane_match = re.search(r"[-_ ]lane[-_ ]?(\d+)$", prefix, re.I)
            name = prefix[: lane_match.start()] if lane_match else prefix
            name = re.sub(r"[-_ ]mirrored$", "", name, flags=re.I)
            stage = args.stage
            if stage == "auto":
                if lane_match:
                    stage = "lane"
                elif len(pages) > 1:
                    stage = "stack"
                elif match and name in names:
                    stage = "split"
                else:
                    stage = "raw"
            if stage in ("stack", "lane"):
                hrs = hours_of_stack(path, len(pages), args.hours)
                bar = existing_bar(im) if stage == "stack" else None
                rec = {
                    "accession": name,
                    "plate_id": name,
                    "stage": stage,
                    "frames": [
                        {"source": str(path), "source_page": i, "time_hours": t}
                        for i, t in zip(pages, hrs)
                    ],
                }
                if bar:
                    rec["existing_scale"] = bar
                if stage == "lane":
                    rec["lane"] = int(lane_match[1]) if lane_match else 1
                label_sheet(im, (0, 0, *im.size), args.out / f"stack-{len(ready)+1}.jpg", path.name)
                ready.append(rec)
                continue
            if not match or len(pages) != 1:
                raise ValueError(
                    f"Raw/split frames need one page and a time in the filename: {path.name}"
                )
            groups[(stage, name)].append((float(match[1]), path, pages[0]))
    for group_idx, ((stage, series), sources) in enumerate(groups.items(), 1):
        sources.sort()
        if len(set(t for t, _, _ in sources)) != len(sources):
            raise ValueError(f"Duplicate time in series {series}. Separate biological replicates.")
        if stage == "split":
            frames = []
            for t, path, page in sources:
                frame = {"source": str(path), "source_page": page, "time_hours": t}
                with Image.open(path) as im:
                    im.seek(page)
                    bar = existing_bar(im)
                if bar:
                    frame["existing_scale"] = bar
                frames.append(frame)
            ready.append({"accession": series, "plate_id": series, "stage": stage, "frames": frames})
            continue
        plates, baseline = [], []
        for i, (t, path, page) in enumerate(sources):
            with Image.open(path) as im:
                im.seek(page)
                boxes = plate_boxes(im)
                if args.box:
                    boxes = args.box
                if not boxes or (i and len(boxes) != len(plates)):
                    # Save the scene so a new layout can be boxed in original pixels.
                    preview(im, 1400, 1000).save(args.out / f"group-{group_idx}-{t:g}h-scene.jpg")
                    raise ValueError(
                        f"Check plate boxes for {path.name}, size {im.size}; use --box per plate"
                    )
                for j, box in enumerate(boxes):
                    if not (0 <= box[0] < box[2] <= im.width and 0 <= box[1] < box[3] <= im.height):
                        raise ValueError(f"Out-of-bounds box {box} for {path.name}")
                    label = f"group-{group_idx}-plate-{j+1}-{t:g}h.jpg"
                    label_sheet(im, box, args.out / label, f"{path.name}  plate {j+1}  box {box}")
                    feats = features(im, box)
                    if i == 0:
                        baseline.append(feats)
                        plates.append(
                            {
                                "accession": None,
                                "plate_id": f"{series or group_idx}-p{j+1}",
                                "stage": "raw",
                                "size": [box[2] - box[0], box[3] - box[1]],
                                "frames": [],
                            }
                        )
                        affine = [[1, 0, -box[0]], [0, 1, -box[1]]]
                        qc = {"method": "reference"}
                    else:
                        affine, qc = register(baseline[j], feats, box)
                    plates[j]["frames"].append(
                        {
                            "time_hours": t,
                            "source": str(path),
                            "source_page": page,
                            "source_to_plate": affine,
                            "label_view": label,
                            "registration": qc,
                        }
                    )
        ready.extend(plates)
    return {"catalog": str(args.catalog.resolve()) if args.catalog else None, "plates": ready}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "inputs",
        nargs="*",
        type=Path,
        help="Raw scans, split photos, stacks, lanes, or a folder at one stage",
    )
    p.add_argument("--out", type=Path)
    p.add_argument("--stage", choices=("auto", "raw", "split", "stack", "lane"), default="auto")
    p.add_argument("--catalog", type=Path, help="Master accession workbook")
    p.add_argument(
        "--lookup",
        nargs="+",
        help="Check visual readings against the catalogue and return candidate spellings",
    )
    p.add_argument(
        "--box",
        nargs=4,
        type=int,
        action="append",
        metavar=("X0", "Y0", "X1", "Y1"),
        help="Override automatic plate location; repeat in physical plate order",
    )
    p.add_argument(
        "--hours",
        nargs="+",
        type=float,
        help="Stack page times unless an adjacent stack_index.csv provides them",
    )
    args = p.parse_args()
    if args.lookup:
        if not args.catalog:
            p.error("--lookup requires --catalog")
        lookup_labels(args.catalog, args.lookup)
        return
    if not args.inputs or not args.out:
        p.error("supply image paths and --out")
    args.out.mkdir(parents=True, exist_ok=True)
    plan = drafts(args)
    (args.out / "draft-plan.json").write_text(
        json.dumps(plan, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(
        f"{len(plan['plates'])} plate/lane records. Read the label views, verify registration and naming, then save a reviewed plan."
    )


if __name__ == "__main__":
    main()
