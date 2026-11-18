"""Run ACORBA v1.3 Scanner on prepared lanes and retain native evidence for review."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

from PIL import Image


OPTIONS = {
    "exp_type": "Scanner", "method": "Deep Machine Learning", "superaccuracy": "True",
    "vector": "5", "normalization": "False", "smooth": "True",  # native switch deactivates smoothing
    "rootplot": "True", "savesegmentation": "True", "prediction": "None",
    "binary_folder": "", "custom": "", "tradmethod": "Entropy", "circlepix": "40", "broken": "50",
}
ANGLE_LIST = re.compile(
    r"Angle list for plate_(?P<name>[^\r\n]+?)_timeframe_(?P<frame>\d+):"
    r"(?P<stage>[^\r\n]*)[\r\n]+\s*(?P<angles>\[[^\]\r\n]*\])"
)


def file_hash(path):
    sha = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            sha.update(block)
    return sha.hexdigest()


def write_json(path, data):
    path = Path(path)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, encoding="utf-8") as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2, allow_nan=False)
    os.replace(stream.name, path)


def final_lists(log):
    records = {}
    text = Path(log).read_bytes().decode("utf-8", errors="replace")
    for match in ANGLE_LIST.finditer(text):
        body = match["angles"][1:-1].strip()
        angles = [float(token) for token in body.split(",")] if body else []
        records[int(match["frame"]) - 1] = {
            "angles": [a if math.isfinite(a) else None for a in angles],
            "stage": "recalculated" if "recalculated" in match["stage"] else "initial",
            "log_line": text.count("\n", 0, match.start()) + 1,
        }
    return records


def lane_key(lane):
    identity = [lane["accession"], str(lane["plate_id"]), str(lane["lane"])]
    readable = re.sub(r"[^a-zA-Z0-9_-]", "_", "-".join(identity))[:75]
    digest = hashlib.sha256(json.dumps(identity, ensure_ascii=False).encode()).hexdigest()[:10]
    return f"{readable}-{digest}"


def inspect_input(lane, manifest_dir):
    path = (manifest_dir / lane["path"]).resolve()
    with Image.open(path) as stack:
        if stack.n_frames != len(lane["frames"]):
            raise ValueError(f"Frame count disagrees with the manifest for {path}")
        if stack.n_frames == 1:
            raise ValueError(f"ACORBA v1.3 Scanner expects a time stack; {path} contains only one TIFF page.")
        size = stack.size
        for index in range(stack.n_frames):
            stack.seek(index)
            if stack.mode != "L" or stack.tag_v2.get(259, 1) != 1 or stack.size != size:
                raise ValueError(f"Expected an uncompressed, single-channel 8-bit TIFF stack: {path}")
    return path


def evidence(folder, nframes):
    log = folder / "native.log"
    mask = folder / "Saved Segmentations" / "lane.tif_segmented.tif"
    records = final_lists(log)
    if set(records) != set(range(nframes)) or not mask.is_file():
        raise ValueError("Missing final angle lists or native segmentation")
    with Image.open(mask) as image:
        if image.n_frames != nframes:
            raise ValueError("Native mask frame count differs from the source")
    files = [log, mask]
    for index, record in records.items():
        suffix = "_recalculated" if record["stage"] == "recalculated" else ""
        files.append(folder / f"lane.tif_plateplot_{index + 1}{suffix}.png")
    if not all(p.is_file() for p in files):
        raise ValueError("Missing final native skeleton plots")
    return files


def valid_cache(record, signature, nframes):
    if not record.is_file():
        return False
    try:
        metadata = json.loads(record.read_text(encoding="utf-8"))
        if metadata["signature"] != signature or metadata["status"] not in ("complete", "native_export_failed"):
            return False
        attempt = record.parent / metadata["attempt"]
        return all(metadata["sha256"].get(p.relative_to(attempt).as_posix()) == file_hash(p)
                   for p in evidence(attempt, nframes))
    except (KeyError, ValueError, OSError):
        return False


def native_export(folder, exit_code):
    workbook = folder / "lane.tif__output.xlsx"
    log = (folder / "native.log").read_bytes().decode("utf-8", errors="replace")
    error = next((line.strip() for line in log.splitlines() if "Error on line" in line), "")
    if exit_code or error or not workbook.is_file():
        return {"status": "native_export_failed", "native_workbook": workbook.is_file(),
                "export_detail": error or f"Native exit {exit_code}; workbook present={workbook.is_file()}"}
    return {"status": "complete", "native_workbook": True, "export_detail": ""}


def run_lane(lane, manifest_dir, output, executable, options, runtime, env):
    source = inspect_input(lane, manifest_dir)
    key = lane_key(lane)
    root = output / key
    root.mkdir(parents=True, exist_ok=True)
    record = root / "run.json"
    nframes = len(lane["frames"])
    signature = dict(input_sha256=file_hash(source), options=options, runtime=runtime)
    if valid_cache(record, signature, nframes):
        prior = json.loads(record.read_text(encoding="utf-8"))
        folder = root / prior["attempt"]
        result = native_export(folder, prior["exit_code"])
        if any(prior.get(k) != v for k, v in result.items()):
            prior.update(result)
            write_json(folder / "metadata.json", prior)
            write_json(record, prior)
        print(f"{key}: verified previous run ({prior['status']})", flush=True)
        return True
    # Unique attempts protect concurrent/retried native processes from overwriting each other's TIFFs.
    attempt = Path(tempfile.mkdtemp(prefix="attempt-", dir=root))
    shutil.copyfile(source, attempt / "lane.tif")
    if file_hash(attempt / "lane.tif") != signature["input_sha256"]:
        raise ValueError(f"Prepared stack changed during copying: {source}")
    args = [str(executable), "--input_folder", str(attempt)]
    for name, value in options.items():
        args.extend(("--" + name, value))
    meta = dict(signature=signature, attempt=attempt.name, command=args, started=time.time(),
                source=str(source), accession=lane["accession"], plate_id=lane["plate_id"],
                lane=lane["lane"], status="running")
    write_json(attempt / "metadata.json", meta)
    print(f"{key}: started {attempt.name}", flush=True)
    with (attempt / "native.log").open("wb") as out:
        result = subprocess.run(args, cwd=executable.parent, env=env, stdout=out,
                                stderr=subprocess.STDOUT, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    meta.update(exit_code=result.returncode, finished=time.time())
    try:
        files = evidence(attempt, nframes)
        meta["sha256"] = {p.relative_to(attempt).as_posix(): file_hash(p) for p in files}
        # ACORBA can fail when forming its workbook after processing all frames (unequal counts).
        # Keep this state separate from a successful native exit; QC can still use complete native evidence.
        meta.update(native_export(attempt, result.returncode))
    except (ValueError, OSError) as err:
        meta.update(status="failed", detail=str(err))
    write_json(attempt / "metadata.json", meta)
    if not valid_cache(record, signature, nframes):
        write_json(record, meta)
    print(f"{key}: {meta['status']} in {meta['finished'] - meta['started']:.0f}s", flush=True)
    return meta["status"] in ("complete", "native_export_failed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Folder for native runs")
    parser.add_argument("--scanner", type=Path, required=True, help="Path to installed scanner.exe")
    parser.add_argument("--gpu-libs", type=Path, action="append", default=[],
                        help="Installed CUDA/cuDNN dependency directory. Repeat for additional directories.")
    parser.add_argument("--accession")
    parser.add_argument("--plate-id")
    parser.add_argument("--lane", type=int)
    parser.add_argument("--vector", type=int, default=5)
    parser.add_argument("--accuracy", choices=("super", "regular"), default="super")
    args = parser.parse_args()
    exe = args.scanner.resolve(strict=True)
    if args.vector < 1:
        parser.error("Vector length must be positive")
    opts = dict(OPTIONS, vector=str(args.vector), superaccuracy=str(args.accuracy == "super"))
    env = dict(os.environ, PYTHONIOENCODING="utf-8", TF_FORCE_GPU_ALLOW_GROWTH="true")
    if args.gpu_libs:
        env["PATH"] = os.pathsep.join([str(p.resolve(strict=True)) for p in args.gpu_libs]
                                      + [env.get("PATH", "")])
    runtime = {"scanner_sha256": file_hash(exe), "models": {
        p.name: file_hash(p) for p in sorted((exe.parent / "models").glob("unetnocontrast*.h5"))}}
    manifest = json.loads(args.manifest.read_text(encoding="utf-8-sig"))
    lanes = [lane for lane in manifest["lanes"] if
             (args.accession is None or lane["accession"] == args.accession) and
             (args.plate_id is None or str(lane["plate_id"]) == args.plate_id) and
             (args.lane is None or int(lane["lane"]) == args.lane)]
    if not lanes or len({lane_key(lane) for lane in lanes}) != len(lanes):
        parser.error("Selection is empty or the manifest contains duplicate lane identities")
    statuses = []
    for lane in lanes:
        try:
            ok = run_lane(lane, args.manifest.resolve().parent, args.output.resolve(), exe, opts, runtime, env)
        except (OSError, ValueError) as err:
            destination = args.output.resolve() / lane_key(lane)
            destination.mkdir(parents=True, exist_ok=True)
            write_json(destination / "run.json", {"status": "failed", "detail": str(err)})
            print(f"{lane_key(lane)}: {err}", file=sys.stderr, flush=True)
            ok = False
        statuses.append(ok)
    return 0 if all(statuses) else 1


if __name__ == "__main__":
    sys.exit(main())
