"""Match native angle lists to segmentation evidence; export only explicitly reviewed roots as accepted."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage, sparse
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import dijkstra
from skimage.morphology import skeletonize

from run_acorba import evidence, file_hash, final_lists, lane_key


FIELDS = ["accession", "plate_id", "lane_id", "frame_index", "time_h", "native_index",
          "angle_deg", "native_stage", "native_plot", "evidence_key", "candidate_id",
          "root_id", "tip_x", "tip_y", "anchor_method", "qc_status", "qc_reason"]
DECISION_FIELDS = ["accession", "plate_id", "lane_id", "frame_index", "native_index",
                   "evidence_key", "candidate_id", "root_id", "qc_status", "qc_reason"]
KEY_FIELDS = DECISION_FIELDS[:5]


def save_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        table = csv.DictWriter(stream, fields, extrasaction="ignore")
        table.writeheader()
        table.writerows(rows)


def ordered_path(mask):
    """Trace the most distant skeleton endpoint from the proximal (left) endpoint for display."""
    skel = skeletonize(mask)
    yy, xx = np.nonzero(skel)
    pixels = np.full(mask.shape, -1, dtype=int)
    pixels[yy, xx] = np.arange(len(xx))
    ypad, xpad = np.pad(pixels, 1), xx + 1
    neighbours = ndimage.convolve(skel.astype("u1"), np.ones((3, 3), dtype="u1"), mode="constant", cval=0)
    ends = np.flatnonzero(neighbours[yy, xx] == 2)
    if len(ends) < 2:
        return None
    start = int(ends[np.argmin(xx[ends])])
    edges, dest, weights = [], [], []
    for dy, dx in ((0, 1), (1, -1), (1, 0), (1, 1)):
        ids = ypad[yy + 1 + dy, xpad + dx]
        ok = ids >= 0
        edges.extend(np.flatnonzero(ok))
        dest.extend(ids[ok])
        weights.extend([math.hypot(dx, dy)] * int(ok.sum()))
    graph = sparse.coo_matrix((weights, (edges, dest)), shape=(len(xx), len(xx))).tocsr()
    distances, parents = dijkstra(graph, directed=False, indices=start, return_predecessors=True)
    reachable = ends[np.isfinite(distances[ends])]
    if len(reachable) < 2:
        return None
    tip = int(reachable[np.argmax(distances[reachable])])
    chain = [tip]
    while chain[-1] != start:
        chain.append(int(parents[chain[-1]]))
    return np.column_stack((xx[chain[::-1]], yy[chain[::-1]]))


def profile(x, y):
    lo, hi = x.min(), x.max()
    bins = np.clip(np.rint((x - lo) / max(hi - lo, 1) * 31).astype(int), 0, 31)
    columns = np.unique(bins)
    medians = [np.median(y[bins == i]) for i in columns]
    curve = np.interp(np.arange(32), columns, medians)
    return (curve - np.min(curve)) / max(float(np.ptp(curve)), 1)


def candidates(mask, anchors):
    labels, _ = ndimage.label(mask, structure=np.ones((3, 3), dtype="u1"))
    result = []
    for idx, box in enumerate(ndimage.find_objects(labels), 1):
        if box is None:
            continue
        yy, xx = box
        piece = labels[box] == idx
        # The native contour filter is stricter. Keep additional objects for visual adjudication.
        if max(piece.shape) < 65 or piece.sum() < 100:
            continue
        path = ordered_path(piece)
        if path is None:
            continue
        path = path + (xx.start, yy.start)
        px, py = path[0].tolist()
        anchor = min(anchors, key=lambda a: np.hypot(px - a["x"], py - a["y"])) if anchors else {}
        close = (bool(anchor) and abs(px - anchor["x"]) < max(70, .12 * mask.shape[1])
                 and abs(py - anchor["y"]) < max(65, .03 * mask.shape[0]))
        result.append(dict(id=f"M{idx}", box=(xx.start, yy.start, xx.stop, yy.stop), path=path,
                           profile=profile(*path.T), width=float(np.ptp(path[:, 0])), close=close,
                           plant=str(anchor["plant"]) if close else "", tip=path[-1].tolist(),
                           base=path[0].tolist()))
    return result, labels


def native_panels(plot, count):
    image = np.asarray(Image.open(plot).convert("RGB"))
    h, w = image.shape[:2]
    # Native v1.3 saves a 6.4 x 4.8-inch Matplotlib figure with 0.2 subplot spacing,
    # tight bounding box, shared x range and a separate y range for every skeleton.
    body = min(h, round(w * .77 * 4.8 / (.775 * 6.4)))
    slot = body / (count + .2 * (count - 1))
    result = []
    for top in range(count):
        y = round(h - body + top * slot * 1.2)
        panel = image[y:round(y + slot)]
        ink = (np.ptp(panel, axis=2) < 7) & (panel[:, :, 0] >= 105) & (panel[:, :, 0] < 220)
        colored = (panel[:, :, 0].astype(float) - panel[:, :, 1] > 65)
        yy, xx = np.nonzero(ink | colored)
        if not len(xx):
            result.append(dict(index=count-top, image=panel, width=0., profile=np.zeros(32)))
        else:
            result.append(dict(index=count-top, image=panel, width=float(np.ptp(xx)), profile=profile(xx, yy)))
    return result


def propose(panels, items):
    if not panels or not items:
        return {}
    # Native repairs may split/join a raw mask. This score suggests correspondences;
    # the reviewer compares the native shape, the mask and the source before acceptance.
    dist = np.array([[np.mean(np.abs(p["profile"] - obj["profile"])) for obj in items] for p in panels])
    widths = np.array([p["width"] for p in panels])[:, None] / np.maximum([o["width"] for o in items], 1)
    ratios = widths[np.arange(len(panels)), dist.argmin(axis=1)]
    scale = max(.01, float(np.median(ratios)))
    cost = dist + .42 * np.abs(np.log(np.maximum(widths, .01) / scale))
    # OpenCV contour enumeration typically goes from bottom to top. A small rank tie-breaker
    # helps with similar traces; rank alone is never sufficient for review.
    cost += .08 * np.abs(np.linspace(0, 1, len(panels))[:, None] - np.linspace(0, 1, len(items)))
    rows, columns = linear_sum_assignment(cost)
    return {panels[i]["index"]: items[j] for i, j in zip(rows, columns)}


def scaled(image, maxw, maxh):
    im = image.copy()
    im.thumbnail((maxw, maxh))
    return im


def font(size):
    for face in ("DejaVuSans.ttf", "arial.ttf"):
        try:
            return ImageFont.truetype(face, size)
        except OSError:
            pass
    return ImageFont.load_default()


def draw_review(path, title, gray, labels, panels, items, matches, anchors):
    w, h = 1490, max(765, 90 + max(len(panels), len(items), 1) * 76)
    sheet = Image.new("RGB", (w, h), "white")
    draw = ImageDraw.Draw(sheet)
    draw.text((12, 8), title, font=font(18), fill="black")
    draw.text((12, 37), "Prepared scan / mask components", fill="black", font=font(13))
    draw.text((355, 37), "Native final skeleton / proposed mask ID", fill="black", font=font(13))
    draw.text((790, 37), "Mask candidates in vertical order", fill="black", font=font(13))
    draw.text((1142, 37), "Same crop without mask", fill="black", font=font(13))
    raw = Image.fromarray(gray).convert("RGB")
    overlay = np.asarray(raw).copy()
    foreground = labels > 0
    overlay[foreground] = .55 * overlay[foreground] + .45 * np.array([53, 250, 120])
    full = scaled(Image.fromarray(overlay), 305, h - 74)
    sheet.paste(full, (26, 67))
    scale = full.height / gray.shape[0]
    for anchor in anchors:
        xx, yy = 26 + anchor["x"]*scale, 67+anchor["y"]*scale
        draw.text((max(1, xx-23), yy-7), str(anchor["plant"]), font=font(12), fill="yellow", stroke_width=1, stroke_fill="black")
    for obj in items:
        xx, yy = 26+obj["tip"][0]*scale, 67+obj["tip"][1]*scale
        draw.text((xx+1, yy-6), obj["id"], fill="#53ffb7", font=font(12), stroke_width=1, stroke_fill="black")
    for i, native in enumerate(panels):
        y = 64 + i*76
        obj = matches.get(native["index"], {})
        angle = native.get("angle")
        legend = f"#{native['index']}  ->  {obj.get('id', '?')}"
        if angle is not None:
            legend += f"   {angle:+.1f} deg"
        draw.text((355, y), legend, fill="#253850", font=font(12))
        shot = scaled(Image.fromarray(native["image"]), 390, 54)
        sheet.paste(shot, (367, y+16))
    for i, obj in enumerate(items):
        y = 64 + i*76
        x0,y0,x1,y1 = obj["box"]
        roi = tuple([max(0,x0-18), max(0,y0-9), min(gray.shape[1],x1+15), min(gray.shape[0],y1+9)])
        crop = Image.fromarray(overlay).crop(roi)
        shot = scaled(crop, 320, 54)
        draw.text((790, y), f"{obj['id']}   planting position {obj['plant'] or '?'}   base {obj['base']}", font=font(12), fill="black")
        sheet.paste(shot, (806, y+16))
        sheet.paste(scaled(raw.crop(roi), 320, 54), (1142, y+16))
    sheet.save(path)


def detection_row(lane, index):
    return dict(accession=lane["accession"], plate_id=lane["plate_id"], lane_id=lane["lane"],
                frame_index=index, time_h=lane["frames"][index]["time_hours"],
                qc_status="unresolved", qc_reason="No reviewed native measurement")


def row_key(row):
    return tuple(str(row[k]) for k in KEY_FIELDS)


def apply_decision(row, decision, items):
    if decision["evidence_key"] != row["evidence_key"]:
        raise ValueError(f"Review belongs to a different run: {row_key(row)}")
    status = decision["qc_status"]
    if status not in ("accepted", "excluded", "unresolved"):
        raise ValueError(f"Unrecognized QC decision: {status}")
    if not decision["qc_reason"].strip():
        raise ValueError(f"Review reason required for {row_key(row)}")
    selected = next((v for v in items if v["id"] == decision["candidate_id"]), None)
    if status == "accepted" and (selected is None or not decision["root_id"].strip() or row["angle_deg"] is None):
        raise ValueError(f"Accepted measurement needs a candidate, finite native angle and plant ID: {row_key(row)}")
    row.update(candidate_id=decision["candidate_id"], root_id=decision["root_id"],
               qc_reason=decision["qc_reason"], qc_status=status, tip_x=None, tip_y=None, anchor_method="")
    if selected:
        row.update(tip_x=selected["tip"][0], tip_y=selected["tip"][1],
                   anchor_method="Distal mask skeleton endpoint (display estimate)")


def process(lane, native, out, decisions, used, manifest_dir):
    key = lane_key(lane)
    meta_file = native / key / "run.json"
    frames = lane["frames"]
    if not meta_file.is_file():
        return [dict(detection_row(lane, i), qc_reason="Native run unavailable") for i in range(len(frames))]
    meta = json.loads(meta_file.read_text(encoding="utf-8"))
    if meta["status"] not in ("complete", "native_export_failed"):
        return [dict(detection_row(lane, i), qc_reason="Native run incomplete; see native log") for i in range(len(frames))]
    attempt = meta_file.parent / meta["attempt"]
    fingerprint = meta["signature"]["input_sha256"]
    if fingerprint != file_hash(manifest_dir / lane["path"]):
        raise ValueError(f"Prepared TIFF differs from the native run: {key}")
    files = evidence(attempt, len(frames))
    log, mask = files[:2]
    for file in files:
        if file_hash(file) != meta["sha256"][file.relative_to(attempt).as_posix()]:
            raise ValueError(f"Evidence changed: {file}")
    native_records = final_lists(log)
    records = []
    review_dir = out / key
    review_dir.mkdir(exist_ok=True)
    with Image.open(mask) as masks, Image.open(manifest_dir / lane["path"]) as source:
        if masks.size != source.size:
            raise ValueError(f"{key}: native mask {masks.size} differs from the prepared lane {source.size}. "
                             "Regular mode can add tile padding. Use Super accuracy or establish the exact "
                             "mask-to-source transform before reviewing these native results.")
        for i, frame in enumerate(frames):
            masks.seek(i)
            source.seek(i)
            anchors = frame.get("anchors", [])
            items, labels = candidates(np.asarray(masks) != 0, anchors)
            result = native_records[i]
            values = result["angles"]
            suffix = "_recalculated" if result["stage"] == "recalculated" else ""
            plot = attempt / f"lane.tif_plateplot_{i+1}{suffix}.png"
            panels = native_panels(plot, len(values)) if values else []
            for panel in panels:
                panel["angle"] = values[panel["index"]-1]
            matches = propose(panels, items)
            figure = review_dir / f"frame-{i:03}.png"
            draw_review(figure, f"{lane['accession']} / {lane['plate_id']} / lane {lane['lane']} / {frame['time_hours']} h", np.asarray(source), labels, panels, items, matches, anchors)
            evidence_key = hashlib.sha256(json.dumps([
                fingerprint, file_hash(plot), meta["sha256"][mask.relative_to(attempt).as_posix()],
                values, frame, "review-v1"], sort_keys=True).encode()).hexdigest()
            seen_plants, seen_masks, reviewed_positions = set(), set(), set()
            for j, angle in enumerate(values, 1):
                row = dict(detection_row(lane, i), native_index=j, angle_deg=angle,
                           native_stage=result["stage"], native_plot=str(plot), evidence_key=evidence_key)
                obj = matches.get(j)
                if obj:
                    row.update(candidate_id=obj["id"], root_id=obj["plant"], tip_x=obj["tip"][0], tip_y=obj["tip"][1],
                               anchor_method="Distal mask skeleton endpoint (display estimate)",
                               qc_reason="Proposed mask correspondence; requires review")
                key_row = row_key(row)
                if key_row in decisions:
                    apply_decision(row, decisions[key_row], items)
                    used.add(key_row)
                    if row["root_id"]:
                        reviewed_positions.add(str(row["root_id"]).split(".")[0])
                if row["qc_status"] == "accepted":
                    if row["root_id"] in seen_plants or row["candidate_id"] in seen_masks:
                        raise ValueError(f"Duplicate accepted plant or mask in {key_row}")
                    seen_plants.add(row["root_id"])
                    seen_masks.add(row["candidate_id"])
                records.append(row)
            for anchor in anchors:
                plant = str(anchor["plant"])
                if plant not in reviewed_positions:
                    records.append(dict(detection_row(lane, i), root_id=plant,
                                        qc_reason="Planting position has no accepted native measurement"))
            if not values and not anchors:
                records.append(detection_row(lane, i))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--decisions", type=Path, help="Reviewed copy of the generated review-template.csv")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    decisions = {}
    if args.decisions:
        with args.decisions.open(encoding="utf-8-sig", newline="") as stream:
            for row in csv.DictReader(stream):
                k = row_key(row)
                if k in decisions:
                    raise ValueError(f"Duplicate decision for {k}")
                decisions[k] = row
    manifest = json.loads(args.manifest.read_text(encoding="utf-8-sig"))
    rows, used = [], set()
    for lane in manifest["lanes"]:
        rows.extend(process(lane, args.native.resolve(), output, decisions, used, args.manifest.resolve().parent))
    if used != set(decisions):
        raise ValueError("Review includes detections absent from this run")
    save_csv(output / "detections.csv", rows, FIELDS)
    save_csv(output / "review-template.csv", [r for r in rows if r.get("native_index")], DECISION_FIELDS)
    print(json.dumps({state: sum(r["qc_status"] == state for r in rows)
                      for state in ("accepted", "excluded", "unresolved")}))


if __name__ == "__main__":
    main()
