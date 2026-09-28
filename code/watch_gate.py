#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""watch_gate.py -- live progress monitor for 15_stratified_gate.py.

Read-only: polls the checkpoint JSONs and prints one updating status line.
Run this in any terminal while the gate is running; it never touches the run:

  python watch_gate.py
"""
import json
import pathlib
import sys
import time

OUT = (pathlib.Path(__file__).resolve().parent.parent
       / "results" / "stratified_gate")
CKPT = OUT / "checkpoints"

TOTAL = 400  # default --reps; only used to compute % of the arm


def snapshot():
    items = []
    for ck in sorted(CKPT.glob("*.json")):
        try:
            blob = json.loads(ck.read_text())
        except (json.JSONDecodeError, OSError):
            continue  # mid-write
        items.append((ck.stem, int(blob.get("done", 0))))
    return items


def main():
    print(f"watching {CKPT}  (Ctrl+C to stop; the run is unaffected)")
    last = None
    while True:
        items = snapshot()
        stamps = {p.name: p.stat().st_mtime for p in CKPT.glob("*.json")}
        now = time.time()
        if items != last:
            last = items
            line = " | ".join(f"{tag}: {done}/{TOTAL}" for tag, done in items)
            if line:
                print(f"\r\033[K[{time.strftime('%H:%M:%S')}] {line}",
                      end="", flush=True)
        if not (OUT / "strat_gate_mag_r400_summary.json").exists():
            time.sleep(5)
        else:
            print("\nsummary written -- run finished:")
            print((OUT / "strat_gate_mag_r400_summary.json").read_text())
            return 0


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nstopped watching (run keeps going)")
