#!/usr/bin/env python3
"""One-rank full-grid timing probes and sample budgets for the N=2,4 sweep.

Probes disable feedback to measure one iteration of each stage. Production
always uses self consistency. Only the Python standard library is required.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
BETAS = ("0.2", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5")
CONFIGURATION_SIZES = (2, 4)


def configurations():
    paths = sorted((ROOT / "Configuration_Data/CubicLattice").glob("Cubic_3D_N=*_NN_J=*0.40825.hdf5"))
    paths = [p for p in paths if int(re.search(r"N=(\d+)", p.stem)[1]) in CONFIGURATION_SIZES]
    if len(paths) != 4:
        raise RuntimeError("Expected both coupling signs for N=2 and N=4")
    return paths


def cases(config):
    return (("zero", "A", "unifield", "0"),
            ("hz0p5", "C", "unifield", "0.5"),
            ("ordering", "C", "unifield" if "J=-" in config else "stagfield", "0.005"))


def common(config, beta, case):
    _, symmetry, field, amplitude = case
    return [f"--config={config}", "--project=CubicLattice", "--spinspinmodel=ISO",
            f"--beta={beta}", "--numImagTimeSteps=50", "--numRealTimeSteps=150",
            "--Tmax=10", "--realTimeSubsteps=1", "--samplingStrategy=independent",
            f"--cstype={symmetry}", f"--{field}={amplitude}"]


def duration(log, label):
    values = re.findall(re.escape(label) + r" done \(0\).*?it took me ([\deE.+-]+) seconds", log)
    if not values:
        raise RuntimeError(f"Missing timing: {label}")
    return max(map(float, values))


def make_plan(rows, output):
    # Inflate measured times 2x and reserve another 20% for drift, MPI and I/O.
    safety, available = 2.0, 6 * 3600 * .8
    fields = ["config", "case", "symmetry", "field", "amplitude", "samples_per_core",
              "preparation_samples_per_core", "preparation_iterlimit", "estimated_seconds",
              "real_sample_seconds", "preparation_sample_seconds", "real_setup_seconds",
              "preparation_setup_seconds", "initialization_seconds", "timing_source", "status"]
    temporary = output.with_suffix(".tsv.tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for config in configurations():
            for case in cases(config.stem):
                sites = re.search(r"N=(\d+)", config.stem)[1]
                # Pool both coupling signs and C-field amplitudes conservatively.
                probes = [r for r in rows if re.search(r"N=(\d+)", r["config"])[1] == sites
                          and (r["case"] == "zero") == (case[0] == "zero")]
                if not probes:
                    raise RuntimeError(f"No timing for N={sites}, symmetry {case[1]}")
                rt = safety * max(r["real_sample_seconds"] for r in probes)
                prep = safety * max(r["preparation_sample_seconds"] for r in probes)
                rs = safety * max(r["real_setup_seconds"] for r in probes)
                ps = safety * max(r["preparation_setup_seconds"] for r in probes)
                initialization = safety * max(r["initialization_seconds"] for r in probes)
                remaining = available - initialization
                samples = math.floor((remaining * .8 / 10 - rs) / rt)
                quantum = 100 if samples >= 1000 else 10 if samples >= 100 else 1
                samples = samples // quantum * quantum
                prep_limit = 20
                prep_samples = min(samples // 10, math.floor((remaining * .2 / prep_limit - ps) / prep))
                status = "ready"
                if min(samples, prep_samples) < 4:
                    samples = prep_samples = 4
                    status = "over_budget"
                if max(r["peak_child_rss_bytes"] for r in probes) > 5_000_000_000:
                    status = "over_memory"
                estimate = initialization + 10 * (rs + rt * samples) + prep_limit * (ps + prep * prep_samples)
                writer.writerow(dict(zip(fields, [config.stem, *case, samples, prep_samples, prep_limit,
                    round(estimate, 1), rt, prep, rs, ps, initialization,
                    ";".join(sorted(set(r["config"] + "/" + r["case"] for r in probes))), status])))
    temporary.replace(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=4)
    parser.add_argument("--betas", nargs="+", default=["5"])
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--sites", nargs="+", type=int, choices=CONFIGURATION_SIZES)
    parser.add_argument("--cases", nargs="+", choices=["zero", "hz0p5", "ordering"], default=["zero", "hz0p5"])
    parser.add_argument("--sign", choices=["FM", "AFM"], default="FM")
    parser.add_argument("--output", type=Path, default=ROOT / "Timing/cubic_beta_sweep")
    parser.add_argument("--executable", type=Path, default=ROOT / "executable_DOUBLE.out")
    args = parser.parse_args()
    if args.samples < 4:
        parser.error("--samples must be at least 4")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    record = args.output / "measurements.json"
    rows = json.loads(record.read_text()) if record.exists() else []
    if not args.plan_only:
        executable = args.executable.resolve()
        fingerprint = hashlib.sha256(executable.read_bytes()).hexdigest()
        env = dict(os.environ)
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "BLAZE_NUM_THREADS"):
            env[key] = "1"
        for config in configurations():
            if ("FM" if "J=-" in config.stem else "AFM") != args.sign:
                continue
            if args.sites and int(re.search(r"N=(\d+)", config.stem)[1]) not in args.sites:
                continue
            for case in cases(config.stem):
                if case[0] not in args.cases:
                    continue
                for beta in args.betas:
                    tag = f"{config.stem}_{case[0]}_beta{beta}_ni50_nr150_n{args.samples}_{fingerprint[:8]}"
                    if any(r["tag"] == tag for r in rows):
                        continue
                    command = [str(executable), *common(str(config), beta, case),
                        "--dstproject=Timing_cubic_sweep", f"--numSamplesPerCore={args.samples}",
                        f"--preparationSamplesPerCore={args.samples}", "--numBlocks=2",
                        "--noselfcons", "--seed=613", f"--fileext={tag}"]
                    print(f"Timing {tag}", flush=True)
                    start = time.monotonic()
                    logfile = args.output / (tag + ".log")
                    with logfile.open("w") as log:
                        subprocess.run(command, cwd=ROOT, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, check=True, timeout=args.timeout)
                    log = re.sub(r"\x1b\[[0-9;]*m", "", logfile.read_text())
                    rows.append(dict(tag=tag, config=config.stem, case=case[0], beta=beta,
                        samples=args.samples, num_imag_time_steps=50, num_real_time_steps=150,
                        executable_sha256=fingerprint,
                        elapsed_seconds=time.monotonic() - start,
                        initialization_seconds=duration(log, "initialization"),
                        # RUSAGE_CHILDREN is the largest completed-child peak so far.
                        peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * (1 if sys.platform == "darwin" else 1024),
                        real_sample_seconds=duration(log, "mixed real-time sampling") / args.samples,
                        preparation_sample_seconds=duration(log, "Matsubara sampling") / args.samples,
                        real_setup_seconds=duration(log, "joint field distribution"),
                        preparation_setup_seconds=duration(log, "Matsubara field distribution")))
                    record.write_text(json.dumps(rows, indent=2) + "\n")
                    print(f"  real {rows[-1]['real_sample_seconds']:.4g} s/sample; prep {rows[-1]['preparation_sample_seconds']:.4g} s/sample", flush=True)
    if args.plan_only or all(any(re.search(r"N=(\d+)", r["config"])[1] == str(n)
        and (r["case"] == "zero") == zero for r in rows) for n in CONFIGURATION_SIZES for zero in (True, False)):
        make_plan(rows, args.output / "sample_plan.tsv")
        print(f"Wrote {args.output / 'sample_plan.tsv'}")


if __name__ == "__main__":
    main()
