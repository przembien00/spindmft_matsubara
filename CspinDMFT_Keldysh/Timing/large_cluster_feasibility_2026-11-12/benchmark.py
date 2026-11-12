"""Full-grid N=7/8 feasibility probes; no production sweep changes.

Four samples complete each stage. Production-block memory is accounted for
separately because four samples permit only two equal statistical blocks.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import time

OUT = Path(__file__).resolve().parent
SOLVER = OUT.parents[1]
EXE = SOLVER / "executable_DOUBLE.out"
fingerprint = hashlib.sha256(EXE.read_bytes()).hexdigest()
env = dict(os.environ)
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "BLAZE_NUM_THREADS"):
    env[key] = "1"
record = OUT / "measurements.json"
rows = json.loads(record.read_text()) if record.exists() else []
for sites, symmetry in ((7, "A"), (7, "C"), (8, "A"), (8, "C")):
    if any(r["sites"] == sites and r["symmetry"] == symmetry
           and r["executable_sha256"] == fingerprint and r["returncode"] == 0
           and r["peak_rss_bytes"] is not None for r in rows):
        continue
    destination = OUT / f"N{sites}_{symmetry}"
    destination.mkdir(exist_ok=True)
    command = [str(EXE),
        f"--config={SOLVER}/Configuration_Data/CubicLattice/Cubic_3D_N={sites}_NN_J=-0.40825.hdf5",
        "--project=CubicLattice", "--dstproject=Results", f"--cstype={symmetry}",
        "--beta=5", "--numImagTimeSteps=50", "--numRealTimeSteps=150", "--Tmax=10",
        "--realTimeSubsteps=1", "--timeTranslationInvariant=true",
        "--samplingStrategy=independent", "--numSamplesPerSet=1", "--numBlocks=32",
        "--numSamplesPerCore=4", "--preparationSamplesPerCore=4", "--noselfcons", "--seed=613"]
    if symmetry == "C":
        command += ["--unifield=0.5"]
    print(f"Starting N={sites}, {symmetry}", flush=True)
    start = time.monotonic()
    logfile = destination / "run.log"
    if logfile.exists():
        logfile.rename(destination / f"run_previous_{time.time_ns()}.log")
    with logfile.open("w") as stream:
        result = subprocess.run(["/usr/bin/time", "-l", *command], cwd=destination,
                                env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=3600)
    log = re.sub(r"\x1b\[[0-9;]*m", "", logfile.read_text())
    row = dict(sites=sites, symmetry=symmetry, samples=4, requested_blocks=32,
               actual_blocks=2, returncode=result.returncode,
               elapsed_seconds=time.monotonic() - start,
               executable_sha256=fingerprint, command=command)
    peak = re.search(r"(\d+)\s+maximum resident set size", log)
    row["peak_rss_bytes"] = int(peak[1]) if peak else None
    for label in ("initialization", "Matsubara field distribution", "Matsubara sampling",
                  "joint field distribution", "mixed real-time sampling"):
        values = re.findall(re.escape(label) + r" done \(0\).*?it took me ([\deE.+-]+) seconds", log)
        row[label] = max(map(float, values)) if values else None
    row["largest_factor_blocks"] = list(map(int,
        re.findall(r"largest Gaussian factor block\s+: (\d+)", log)))
    rows.append(row)
    record.write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps({k: v for k, v in row.items() if k != "command"}), flush=True)
    if result.returncode:
        raise SystemExit(result.returncode)
