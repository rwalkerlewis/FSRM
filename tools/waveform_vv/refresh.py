#!/usr/bin/env python3
"""ObsPy-based IRIS waveform refresh tool for the FSRM pass-8
waveform V&V infrastructure.

Reads tools/waveform_vv/events.yaml, queries the IRIS DMC for each
station/channel triple per event, deconvolves the instrument
response to displacement (m), and writes one SAC file per
(event, station, channel) under tools/waveform_vv/cache/<EventName>/
along with a metadata.yaml recording the IRIS query parameters and
the download timestamp.

This is a runtime tool. It is NOT a build-time dependency. The
fsrm-ci:local Docker image does not need ObsPy. The C++
iris_validation tests GTEST_SKIP cleanly when the cache is empty.

Usage:
    python tools/waveform_vv/refresh.py --event Salmon1964 \\
        --cache-root tools/waveform_vv/cache

Optional args:
    --all         refresh every event in the manifest
    --manifest    path to alternative manifest (default events.yaml)
    --dry-run     print queries but skip download

Requirements (install in your local venv, NOT in the FSRM Docker
image):
    pip install obspy>=1.4 pyyaml

References:
    IRIS DMC FDSN web services:
        https://service.iris.edu/fdsnws/
    ObsPy:
        https://docs.obspy.org/
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
DEFAULT_MANIFEST = THIS_DIR / "events.yaml"
DEFAULT_CACHE_ROOT = THIS_DIR / "cache"


def _require(modname):
    try:
        return __import__(modname)
    except ImportError as e:
        print(f"ERROR: this script requires {modname!r}.", file=sys.stderr)
        print(f"       pip install {modname}", file=sys.stderr)
        raise SystemExit(2) from e


def load_manifest(path: Path):
    yaml = _require("yaml")
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def fetch_event(event_name: str, event_cfg: dict, cache_root: Path,
                dry_run: bool) -> None:
    obspy = _require("obspy")
    from obspy.clients.fdsn import Client  # type: ignore[import]

    event_dir = cache_root / event_name
    event_dir.mkdir(parents=True, exist_ok=True)

    origin = obspy.UTCDateTime(event_cfg["origin_time"])
    pre = float(event_cfg.get("pre_origin_s", 30.0))
    post = float(event_cfg.get("post_origin_s", 600.0))
    t_begin = origin - pre
    t_end = origin + post

    networks = event_cfg.get("networks", ["II", "IU"])
    networks_str = ",".join(networks)

    metadata = {
        "event_name": event_name,
        "origin_time": str(origin),
        "fetched_utc": datetime.utcnow().isoformat() + "Z",
        "iris_networks": networks,
        "stations": [],
    }

    client = Client("IRIS") if not dry_run else None

    for s in event_cfg.get("stations", []):
        sta = s["code"]
        for ch in s.get("channels", ["BHZ"]):
            print(f"[{event_name}] {sta}.{ch} "
                  f"{t_begin} -- {t_end}")
            if dry_run:
                continue
            try:
                st = client.get_waveforms(
                    network=networks_str, station=sta,
                    location="*", channel=ch,
                    starttime=t_begin, endtime=t_end,
                    attach_response=True)
                # Pre-process: detrend, taper, deconvolve to disp.
                for tr in st:
                    tr.detrend("demean")
                    tr.detrend("linear")
                    tr.taper(0.05, type="cosine")
                    try:
                        tr.remove_response(output="DISP",
                                           pre_filt=(0.05, 0.1, 20.0, 25.0))
                    except Exception as e:
                        print(f"  WARN: response deconv failed: {e}")
                # Write each trace separately.
                for tr in st:
                    out = event_dir / f"{sta}.{ch}.sac"
                    tr.write(str(out), format="SAC")
                    print(f"  wrote {out}")
                metadata["stations"].append({
                    "code": sta, "channel": ch,
                    "n_traces": len(st),
                    "iris_status": "ok"})
            except Exception as e:
                print(f"  ERROR: {e}")
                metadata["stations"].append({
                    "code": sta, "channel": ch,
                    "iris_status": f"fail: {e}"})

    yaml = _require("yaml")
    with open(event_dir / "metadata.yaml", "w", encoding="utf-8") as f:
        yaml.safe_dump(metadata, f)
    print(f"[{event_name}] metadata.yaml updated")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--event", help="Single event name to refresh")
    p.add_argument("--all", action="store_true",
                   help="Refresh every event in the manifest")
    p.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    p.add_argument("--cache-root", default=str(DEFAULT_CACHE_ROOT))
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args(argv)

    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        print(f"manifest not found: {manifest_path}", file=sys.stderr)
        return 2

    manifest = load_manifest(manifest_path)
    events = manifest.get("events", {}) or {}
    if not events:
        print("manifest has no 'events' section", file=sys.stderr)
        return 2

    cache_root = Path(args.cache_root)
    cache_root.mkdir(parents=True, exist_ok=True)

    if args.event:
        if args.event not in events:
            print(f"event {args.event!r} not in manifest", file=sys.stderr)
            return 2
        fetch_event(args.event, events[args.event], cache_root,
                    args.dry_run)
    elif args.all:
        for name, cfg in events.items():
            fetch_event(name, cfg, cache_root, args.dry_run)
    else:
        p.print_help()
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
