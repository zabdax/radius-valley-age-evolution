"""make_data_bundle.py — Build colab_gate_data.zip for the Colab experiment.

Zips the five frozen input CSVs (from ../data) into colab_gate_data.zip,
which is the only file that needs to be uploaded to Google Colab to run
colab_highrep_gate.py. See colab/README.md.
"""
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = HERE / "colab_gate_data.zip"

FILES = [
    "planet_sample_FGK.csv",
    "hosts_kinematics_FGK.csv",
    "completeness_weights_FGK.csv",
    "planet_sample_M.csv",
    "hosts_kinematics_M.csv",
]

with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as zf:
    for name in FILES:
        zf.write(DATA / name, arcname=name)
        print(f"added {name}")

print(f"wrote {OUT} ({OUT.stat().st_size:,} bytes)")