"""Export a reviewed ACORBA batch to Excel, a time-course graph and an original-image atlas."""
import argparse
import colorsys
import csv
import io
import json
import math
import os
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import groupby
from pathlib import Path

from PIL import Image
from reportlab.lib.colors import HexColor, white
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from reportlab import __file__ as reportlab_path
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont


def report_fonts():
    bundled = Path(reportlab_path).resolve().parent / "fonts"
    win_fonts = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"
    if (win_fonts / "arial.ttf").is_file() and (win_fonts / "arialbd.ttf").is_file():
        normal, bold = win_fonts / "arial.ttf", win_fonts / "arialbd.ttf"
    else:
        normal, bold = bundled / "Vera.ttf", bundled / "VeraBd.ttf"
    pdfmetrics.registerFont(TTFont("Report-Regular", str(normal)))
    pdfmetrics.registerFont(TTFont("Report-Bold", str(bold)))
    return "Report-Regular", "Report-Bold"


BODY_FONT, BOLD_FONT = report_fonts()


def accession_key(name):
    return (name.casefold() != "col-0", name.casefold())


def number(value):
    if value is None or str(value).strip() == "":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite measurement: {value!r}")
    return result


def transformed(matrix, x, y):
    return tuple(a*x + b*y + c for a, b, c in matrix)


@dataclass(frozen=True)
class AngleMethod:
    source: str
    arc_mm: float | None = None

    @property
    def label(self):
        if self.source == "native":
            return "Native ACORBA"
        return f"Subpixel centerline, {self.arc_mm:.3f} mm distal arc"


def method_from_records(records, source):
    if source == "native":
        if any(r.get("measurement_method") == "subpixel_centerline" for r in records):
            raise ValueError("Use the original native QC detections.csv for --angle-source native.")
        return AngleMethod(source)
    if any(r.get("measurement_method") != "subpixel_centerline" for r in records):
        raise ValueError("The default report requires the reviewed subpixel detections.csv. "
                         "To report native ACORBA values, use --angle-source native with the native QC file.")
    lengths = {number(r["length_mm"]) for r in records if r.get("length_mm")}
    if len(lengths) != 1 or min(lengths) <= 0:
        raise ValueError("A subpixel report must have a single positive physical arc length.")
    return AngleMethod(source, lengths.pop())


def coordinate(record, x_field, y_field):
    x, y = number(record.get(x_field)), number(record.get(y_field))
    if (x is None) != (y is None):
        raise ValueError(f"Incomplete coordinate in {x_field} and {y_field}")
    return None if x is None else (x, y)


def subpixel_geometry(record, angle, status):
    fit_status = record["fit_status"]
    if (status == "accepted") != (fit_status == "accepted"):
        raise ValueError("The report status and subpixel fit review do not agree.")
    tip = coordinate(record, "tip_source_x", "tip_source_y")
    upstream = coordinate(record, "upstream_source_x", "upstream_source_y")
    if status == "accepted":
        proposed = number(record["fit_angle_deg"])
        if tip is None or upstream is None or proposed is None or angle is None or tip == upstream:
            raise ValueError("An accepted subpixel angle requires the fitted angle and both source endpoints.")
        if not math.isclose(angle, proposed, abs_tol=1e-6):
            raise ValueError("The reported angle differs from the reviewed subpixel fit.")
        mx = record["mirror_source_x"].strip().lower()
        dpi_x, dpi_y = number(record["dpi_x"]), number(record["dpi_y"])
        if mx not in ("true", "false") or not dpi_x or not dpi_y or min(dpi_x, dpi_y) <= 0:
            raise ValueError("An accepted fit requires its source orientation and positive X/Y resolution.")
        dx = (tip[0]-upstream[0]) / dpi_x
        dy = (tip[1]-upstream[1]) / dpi_y
        geometric_angle = math.degrees(math.atan2(dy, -dx if mx == "true" else dx))
        if abs(math.remainder(angle-geometric_angle, 360)) > 0.001:
            raise ValueError(f"The signed angle {angle:.4f}° disagrees with its calibrated source endpoints "
                             f"({geometric_angle:.4f}°).")
    elif angle is not None:
        raise ValueError("A provisional fit angle belongs in fit_angle_deg until the fit is accepted.")
    return tip, upstream


def read_inputs(manifest, detections, angle_source="subpixel"):
    frames = {}
    for lane in json.loads(manifest.read_text(encoding="utf-8-sig"))["lanes"]:
        for index, entry in enumerate(lane["frames"]):
            key = lane["accession"], str(lane["plate_id"]), int(lane["lane"]), index
            if key in frames:
                raise ValueError(f"Duplicate manifest frame {key}")
            path = Path(entry["source"])
            frames[key] = dict(entry, source=path if path.is_absolute() else manifest.parent/path)
    rows, seen, covered = [], set(), set()
    with detections.open(encoding="utf-8-sig", newline="") as handle:
        records = list(csv.DictReader(handle))
    if not records:
        raise ValueError("No detections or missing-frame records.")
    method = method_from_records(records, angle_source)
    for record in records:
        key = record["accession"], record["plate_id"], int(record["lane_id"]), int(record["frame_index"])
        frame = frames[key]
        covered.add(key)
        time_h = number(record["time_h"])
        if time_h != number(frame["time_hours"]):
            raise ValueError(f"Detection time differs from source manifest for {key}")
        status = record["qc_status"]
        if status not in ("accepted", "excluded", "unresolved"):
            raise ValueError(f"Unrecognized QC status for {key}: {status}")
        angle = number(record.get("angle_deg"))
        upstream = None
        if method.source == "subpixel":
            point, upstream = subpixel_geometry(record, angle, status)
        else:
            anchor = coordinate(record, "tip_x", "tip_y")
            point = transformed(frame["lane_to_source"], *anchor) if anchor is not None else None
        root = record.get("root_id", "").strip()
        if status == "accepted":
            identity = (*key[:3], time_h, root)
            if angle is None or point is None or not root or identity in seen:
                raise ValueError(f"Missing or duplicate accepted measurement for {identity}")
            seen.add(identity)
        rows.append(dict(record, key=key, time=time_h, angle=angle, point=point, upstream=upstream, root=root))
    absent = frames.keys() - covered
    if absent:
        raise ValueError(f"{len(absent)} manifest frames have no QC records. Complete QC, including missing-frame rows. First: {min(absent)}")
    return method, frames, sorted(rows, key=lambda r: (accession_key(r["accession"]), r["time"], r["plate_id"],
                                                r["key"][2], r["root"], r["native_index"]))


def summarize(rows):
    result = []
    for (accession, time_h), items in groupby(rows, key=lambda r: (r["accession"], r["time"])):
        records = list(items)
        counts = Counter(r["qc_status"] for r in records)
        unresolved_objects = sum(r["qc_status"] == "unresolved" and bool(r["native_index"]) for r in records)
        angles = [r["angle"] for r in records if r["qc_status"] == "accepted"]
        n = len(angles)
        sd = statistics.stdev(angles) if n > 1 else None
        result.append([accession, time_h, n, statistics.mean(angles) if n else None,
                       sd, sd/math.sqrt(n) if n > 1 else None, counts["excluded"], unresolved_objects,
                       counts["unresolved"] - unresolved_objects])
    return result


def excel_sheets(rows, frames, summary, method):
    subpixel = method.source == "subpixel"
    angle_header = "Accepted subpixel angle (°)" if subpixel else "Native ACORBA angle (°)"
    header = ["Accession", "Time (h)", "Plate", "Lane", "Root ID", angle_header, "QC status", "QC reason",
              "Native object", "Source X (px)", "Source Y (px)", "Anchor method", "Source image",
              "Source page (1-based)", "Source stage", "Native plot", "Angle method"]
    widths = [18,12,17,8,12,31,15,72,15,18,18,54,60,23,16,60,47]
    decimal_cols = [5,9,10]
    if subpixel:
        header += ["Distal arc (mm)", "Upstream source X (px)", "Upstream source Y (px)", "Fit review",
                   "Fitted candidate angle (°)", "Native ACORBA angle (°)", "Native QC", "Native QC reason",
                   "Fit review image", "Source X resolution (dpi)", "Source Y resolution (dpi)", "Source X mirrored"]
        widths += [19,26,26,18,30,28,17,64,60,28,28,22]
        decimal_cols += [18,19,21,22]
    data = []
    for row in rows:
        frame = frames[row["key"]]
        x, y = row["point"] or (None, None)
        values = [row["accession"], row["time"], row["plate_id"], row["key"][2], row["root"] or None,
                     row["angle"], row["qc_status"], row["qc_reason"], row["native_index"] or None,
                     x, y, ("Fitted optical tip" if x is not None else None) if subpixel else row.get("anchor_method"),
                     str(frame["source"]), int(frame["source_page"])+1, frame["source_stage"], row.get("native_plot"),
                     method.label]
        if subpixel:
            ux, uy = row["upstream"] or (None, None)
            values += [number(row["length_mm"]), ux, uy, row["fit_status"], number(row["fit_angle_deg"]),
                       number(row["native_angle_deg"]), row["native_qc_status"], row["native_qc_reason"],
                       row["review_sheet"], number(row["dpi_x"]), number(row["dpi_y"]), row["mirror_source_x"]]
        data.append(values)
    return [
        dict(name="Measurements", header=header, data=data, widths=widths, decimal_cols=decimal_cols),
        dict(name="Summary", header=["Accession", "Time (h)", "N roots", "Mean (°)", "Sample SD (°)",
                                     "SE (°)", "Excluded objects", "Unresolved objects", "Unassigned positions / unavailable frames",
                                     "Angle method"],
             data=[[*s, method.label] for s in summary], widths=[18,12,12,16,19,16,22,23,49,47], decimal_cols=[3,4,5])]


def workbook(path, sheets):
    # Use the desktop's bundled spreadsheet library when present. A public, Python-only
    # installation can create the same two-sheet workbook using openpyxl.
    runtime_roots = [Path(sys.executable).resolve().parents[1],
                     Path.home()/".cache/codex-runtimes/codex-primary-runtime/dependencies"]
    for root in runtime_roots:
        node = root/"node"/"bin"/("node.exe" if os.name == "nt" else "node")
        lib = root/"node"/"node_modules"/"@oai"/"artifact-tool"
        if node.is_file() and (lib/"package.json").is_file():
            script = Path(__file__).with_name("report_workbook.mjs")
            subprocess.run([str(node), str(script), str(lib), str(path)], input=json.dumps(sheets).encode(), check=True)
            return
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter
    book = Workbook()
    book.remove(book.active)
    for source in sheets:
        sheet = book.create_sheet(source["name"])
        sheet.append(source["header"])
        for data_row, values in enumerate(source["data"], 2):
            sheet.append(values)
            for col, value in enumerate(values, 1):
                if isinstance(value, str):
                    sheet.cell(data_row, col).data_type = "s"
            for i in source["decimal_cols"]:
                sheet.cell(data_row, i+1).number_format = '0.00'
        for cell in sheet[1]:
            cell.fill = PatternFill("solid", fgColor="254C58")
            cell.font = Font(color="FFFFFF", bold=True)
        for i, width in enumerate(source["widths"], 1):
            sheet.column_dimensions[get_column_letter(i)].width = width
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
    book.save(path)


def nice_ticks(lo, hi, target=9):
    raw = (hi-lo)/target if hi > lo else 1.0
    magnitude = 10**math.floor(math.log10(raw))
    step = next(n*magnitude for n in (1,2,5,10) if n*magnitude >= raw)
    return [i*step for i in range(math.floor(lo/step), math.ceil(hi/step)+1)]


def chart(path, summary, method):
    groups = {key: list(vals) for key, vals in groupby(summary, key=lambda s:s[0])}
    times = sorted({row[1] for row in summary})
    palette = ["#202020", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#7B3294", "#E69F00",
               "#56B4E9", "#C8102E", "#4D9221", "#31688E", "#AA4499", "#8C510A", "#00A6A6"]
    width, height = 1008, max(576, 100+16*len(groups))
    c = canvas.Canvas(str(path), pagesize=(width,height))
    c.setTitle(f"Root curvature by accession | {method.label}")
    legend_width = max(105, max(c.stringWidth(k, BODY_FONT, 10) for k in groups)+45)
    left, bottom, right, top = 75, 68, width-legend_width-22, height-22
    bounds = [v for _,_,_,mean,sd,*_ in summary if mean is not None for v in (mean-(sd or 0), mean+(sd or 0))]
    ymin, ymax = min([0,*bounds]), max([0,*bounds])
    if ymin == ymax:
        ymin -= 1
        ymax += 1
    ticks = nice_ticks(ymin, ymax)
    xpad = max(.04, (times[-1]-times[0])*.014)
    xmin, xmax = times[0]-xpad, times[-1]+xpad
    X = lambda t: left+(right-left)*(t-xmin)/(xmax-xmin)
    Y = lambda a: bottom+(top-bottom)*(a-ticks[0])/(ticks[-1]-ticks[0])
    c.setFillColor(HexColor("#EBEBEB"))
    c.rect(left,bottom,right-left,top-bottom,fill=1,stroke=0)
    xticks = times if len(times)<=16 else [t for t in nice_ticks(times[0],times[-1],14) if xmin <= t <= xmax]
    for t in xticks:
        c.setStrokeColor(white)
        c.setLineWidth(.65)
        c.line(X(t),bottom,X(t),top)
        c.setFillColor(HexColor("#333333"))
        c.setFont(BODY_FONT,10)
        c.drawCentredString(X(t), bottom-16, f"{t:g}")
    for tick in ticks:
        c.setStrokeColor(white)
        c.line(left,Y(tick),right,Y(tick))
        c.setFillColor(HexColor("#333333"))
        c.drawRightString(left-9,Y(tick)-3,f"{tick:g}")
    c.setLineWidth(.6)
    c.setStrokeColor(HexColor("#666666"))
    c.rect(left,bottom,right-left,top-bottom,fill=0)
    c.setFillColor(HexColor("#161616"))
    c.setFont(BODY_FONT,16)
    c.drawCentredString((right+left)/2, bottom-40, "time point")
    c.saveState()
    c.translate(24,(top+bottom)/2)
    c.rotate(90)
    c.drawCentredString(0,0,"curvature/°")
    c.restoreState()
    legend_x, legend_top = right+20,(bottom+top)/2 + 8*len(groups)
    c.setFont(BOLD_FONT,11)
    c.drawString(legend_x,legend_top+12,"Accession")
    for i,(name,series) in enumerate(groups.items()):
        if len(groups) <= len(palette):
            color = palette[i]
        else:
            color = "#202020" if name == "Col-0" else "#%02x%02x%02x"%tuple(
                round(255*v) for v in colorsys.hsv_to_rgb(i/len(groups),.78,.65))
        ink = HexColor(color)
        c.setStrokeColor(ink)
        c.setFillColor(ink)
        previous = None
        for _,time_h,_,mean,sd,*_ in series:
            if mean is None:
                previous = None
                continue
            x,y = X(time_h),Y(mean)
            if previous is not None:
                c.setLineWidth(1.25)
                c.line(*previous,x,y)
            if sd is not None:
                c.setLineWidth(.6)
                c.line(x,Y(mean-sd),x,Y(mean+sd))
                for edge in (Y(mean-sd),Y(mean+sd)):
                    c.line(x-2.3,edge,x+2.3,edge)
            c.circle(x,y,2,stroke=0,fill=1)
            previous = (x,y)
        y = legend_top-16*i-5
        c.setLineWidth(1.25)
        c.line(legend_x,y,legend_x+19,y)
        c.circle(legend_x+9.5,y,2.2,fill=1,stroke=0)
        c.setFillColor(HexColor("#161616"))
        c.setFont(BODY_FONT,10)
        c.drawString(legend_x+26,y-3.3,name)
    c.setFont(BODY_FONT,9)
    c.drawString(left,10,f"{method.label}. Time in hours. Points show accepted root means. Bars show ±1 sample SD for two or more roots.")
    c.showPage()
    c.save()


def draw_image(c, img, x, y, w, h):
    # Embed a compressed display copy. All cropping and overlay positions use original pixel coordinates.
    buf = io.BytesIO()
    if img.mode not in ("RGB", "L"):
        img = img.convert("RGB")
    img.save(buf, "JPEG", quality=94)
    buf.seek(0)
    c.drawImage(ImageReader(buf), x, y, w, h)


def draw_tip(c, im, frame, row, x, y, size):
    ax, ay = row["point"]
    origin = row["upstream"]
    for px, py in ([row["point"], origin] if origin is not None else [row["point"]]):
        if not (-.5 <= px < im.width-.5 and -.5 <= py < im.height-.5):
            raise ValueError(f"Measurement endpoint outside source image for {row['key']}: {(px,py)}")
    if origin is None:
        span, center_x, center_y = 220, ax, ay
    else:
        # Magnify the measured chord. Both ends stay in their original pixel positions.
        span = max(52, math.ceil(6 * math.dist(origin, row["point"])))
        center_x, center_y = (ax+origin[0])/2, (ay+origin[1])/2
    span = min(span, im.width, im.height)
    lx = min(max(round(center_x)-span//2,0), im.width-span)
    ly = min(max(round(center_y)-span//2,0), im.height-span)
    draw_image(c, im.crop((lx, ly, lx+span, ly+span)), x,y,size,size)
    cx,cy = x+size*(ax-lx+.5)/span, y+size*(ly+span-ay-.5)/span
    c.saveState()
    rect = c.beginPath()
    rect.rect(x,y,size,size)
    c.clipPath(rect,stroke=0)
    if origin is None:
        mat = frame["lane_to_source"]
        radians = math.radians(row["angle"])
        def direction(vx, vy):
            dx,dy = mat[0][0]*vx+mat[0][1]*vy, -(mat[1][0]*vx+mat[1][1]*vy)
            norm = math.hypot(dx,dy)
            return dx/norm,dy/norm
        bx,by = direction(1,0)
        dx,dy = direction(math.cos(radians), math.sin(radians))
        tx,ty = cx+dx*size*.31, cy+dy*size*.31
    else:
        tx,ty = cx,cy
        cx,cy = x+size*(origin[0]-lx+.5)/span, y+size*(ly+span-origin[1]-.5)/span
        bx,by = (-1 if row["mirror_source_x"].strip().lower()=="true" else 1), 0
    c.setStrokeColor(white)
    c.setDash(3,2)
    c.setLineWidth(1)
    c.line(cx,cy,cx+bx*size*.32,cy+by*size*.32)
    c.setDash()
    for shade,weight in (("#14282D",3.6),("#35F5CF",1.8)):
        c.setStrokeColor(HexColor(shade))
        c.setLineWidth(weight)
        c.line(cx,cy,tx,ty)
        if origin is None:
            for sign in (-1,1):
                c.line(tx,ty,tx-6*dx+sign*3*dy,ty-6*dy-sign*3*dx)
    c.setFillColor(white)
    c.circle(cx,cy,2.1,fill=1,stroke=0)
    if origin is not None:
        c.setFillColor(HexColor("#FFC857"))
        c.circle(tx,ty,2.1,fill=1,stroke=0)
    c.restoreState()
    caption = f"{row['angle']:+.2f}°"
    c.setFillColor(HexColor("#193944"))
    c.roundRect(x+5, y+size-20, c.stringWidth(caption, BOLD_FONT,10)+10, 15, 3, stroke=0, fill=1)
    c.setFont(BOLD_FONT, 10)
    c.setFillColor(white)
    c.drawString(x+10, y+size-15.6, caption)


def annotate_overview(c, points, bounds):
    """Keep numeric labels and their leader lines on the source image."""
    left, bottom, right, top = bounds
    occupied = [(px-3,py-3,px+3,py+3) for px,py,_ in points]
    def overlap(a, b):
        return max(0,min(a[2],b[2])-max(a[0],b[0]))*max(0,min(a[3],b[3])-max(a[1],b[1]))
    for px, py, value in sorted(points, key=lambda p:(-p[1],p[0])):
        label = f"{value:+.2f}°"
        width, height = c.stringWidth(label, BOLD_FONT, 6.8)+4, 9
        boxes = []
        for dy in (-4,5,-13,14,-22,23,-31,32):
            for dx in (5,-width-5):
                xx = max(left, min(px+dx, right-width))
                yy = max(bottom, min(py+dy, top-height))
                box = (xx, yy, xx+width, yy+height)
                cost = sum(overlap(box, other) for other in occupied)
                boxes.append((cost,abs(dy+4),box))
        _,_,best = min(boxes, key=lambda item:(item[0],item[1]))
        occupied.append(best)
        x,y,rx,ty = best
        c.setStrokeColor(HexColor("#35F5CF"))
        c.setLineWidth(.7)
        c.line(px,py,max(x,min(px,rx)),max(y,min(py,ty)))
        c.circle(px,py,1.8,fill=0)
        c.setFillColor(HexColor("#193944"))
        c.roundRect(x,y,width,height,1.5,fill=1,stroke=0)
        c.setFillColor(white)
        c.setFont(BOLD_FONT,6.8)
        c.drawString(x+2,y+2.1,label)


def atlas(path, frames, rows, method):
    by_key = defaultdict(list)
    for row in rows:
        by_key[row["key"]].append(row)
    by_source = defaultdict(list)
    for key, frame in frames.items():
        by_source[(frame["time_hours"], frame["source"], frame["source_page"])].append(key)
    w,h = 842,595
    c = canvas.Canvas(str(path),pagesize=(w,h),pageCompression=1)
    c.setTitle(f"Root-tip angles on source images | {method.label}")
    def heading(title, subtitle):
        c.setFillColor(HexColor("#193944"))
        c.setFont(BOLD_FONT,20)
        c.drawString(30,h-37,title)
        c.setFillColor(HexColor("#45616A"))
        c.setFont(BODY_FONT,10)
        c.drawString(30,h-56,subtitle)
    def footer(note):
        c.setStrokeColor(HexColor("#D0D9D9"))
        c.line(30,37,w-30,37)
        c.setFillColor(HexColor("#45616A"))
        c.setFont(BODY_FONT,9)
        c.drawString(30,21,note)
        c.drawRightString(w-30,21,str(c.getPageNumber()))
        c.showPage()
    for (time_h, source, page), keys in sorted(by_source.items(),key=lambda k:(k[0][0],str(k[0][1]),k[0][2])):
        with Image.open(source) as im:
            im.seek(page)
            rows_by_source = [r for key in keys for r in by_key[key] if r["qc_status"]=="accepted"]
            heading(f"Source overview  |  {time_h:g} h",f"{source.name}   /   source page {page+1}   /   {len(rows_by_source)} accepted angles   /   {method.label}")
            thumb=im.copy()
            thumb.thumbnail((1750,975),Image.Resampling.LANCZOS)
            factor=min((w-65)/im.width,(h-123)/im.height)
            iw,ih=im.width*factor,im.height*factor
            x,y=(w-iw)/2,48+(h-123-ih)/2
            draw_image(c,thumb,x,y,iw,ih)
            points = [(x+(r["point"][0]+.5)*factor, y+(im.height-r["point"][1]-.5)*factor, r["angle"])
                      for r in rows_by_source]
            if method.source == "subpixel":
                c.setLineWidth(.65)
                c.setStrokeColor(HexColor("#35F5CF"))
                for (px, py, _), row in zip(points, rows_by_source):
                    ux, uy = row["upstream"]
                    c.line(x+(ux+.5)*factor, y+(im.height-uy-.5)*factor, px, py)
            annotate_overview(c,points,(x,y,x+iw,y+ih))
            stages = ", ".join(sorted({frames[k]["source_stage"] for k in keys}))
            footer(f"Input stage: {stages}. Labels show {method.label} angles. Individual roots appear on the following pages.")
            for key in sorted(keys,key=lambda k:(accession_key(k[0]), k[1], k[2])):
                records = by_key[key]
                accepted = [r for r in records if r["qc_status"]=="accepted"]
                accepted.sort(key=lambda r:(not r["root"].isdigit(), int(r["root"]) if r["root"].isdigit() else r["root"]))
                counts = Counter(r["qc_status"] for r in records)
                no_assignment = sum(r["qc_status"] == "unresolved" and not r["native_index"] for r in records)
                for start in range(0,max(len(accepted),1),8):
                    label=f"{key[0]}   /   plate {key[1]}   /   lane {key[2]}   /   {time_h:g} h"
                    heading(label, f"{source.name}   /   {counts['accepted']} accepted     {counts['excluded']} excluded objects"
                                   f"     {counts['unresolved']-no_assignment} unresolved objects     {no_assignment} unassigned records")
                    if not accepted:
                        c.setFont(BODY_FONT,14)
                        c.drawCentredString(w/2,h/2,"No accepted angles for this lane and time point.")
                    for i,r in enumerate(accepted[start:start+8]):
                        x=30+(i%4)*198
                        y=h-267-(i//4)*229
                        size=173
                        draw_tip(c,im,frames[key],r,x,y,size)
                        c.setFillColor(HexColor("#193944"))
                        c.setFont(BOLD_FONT,12)
                        c.drawString(x,y-16,f"Root {r['root']}   {r['angle']:.2f}°")
                        c.setFont(BODY_FONT,8.2)
                        c.setFillColor(HexColor("#47636B"))
                        c.drawString(x,y-30,f"Source x {r['point'][0]:.1f} px, y {r['point'][1]:.1f} px")
                    note = ("Cyan chord connects the upstream point (white) to the optical tip (gold). Dashed line shows zero."
                            if method.source == "subpixel" else
                            "Arrow shows the signed angle at the mapped display anchor. Dashed line shows zero.")
                    footer(f"{method.label}. {note}")
    c.save()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("manifest","detections","out"):
        p.add_argument("--"+flag, type=Path, required=True)
    p.add_argument("--angle-source", choices=("subpixel","native"), default="subpixel",
                   help="Use the reviewed subpixel table (default) or the original native ACORBA QC table.")
    args=p.parse_args()
    method, frames, rows = read_inputs(args.manifest,args.detections,args.angle_source)
    summary=summarize(rows)
    names = ["root_angles.xlsx", "curvature.pdf", "annotated_angles.pdf"]
    args.out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.out) as folder:
        stage=Path(folder)
        workbook(stage/names[0],excel_sheets(rows,frames,summary,method))
        chart(stage/names[1],summary,method)
        atlas(stage/names[2],frames,rows,method)
        for name in names:
            os.replace(stage/name,args.out/name)
    print(json.dumps({"angle_method":method.label,"groups":len(summary),"qc":dict(Counter(r["qc_status"] for r in rows)),
                      "files":[str(args.out/name) for name in names]},ensure_ascii=True))


if __name__=="__main__":
    main()
