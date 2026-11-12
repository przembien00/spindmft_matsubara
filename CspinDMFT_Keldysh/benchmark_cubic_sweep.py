#!/usr/bin/env python3
"""One-rank full-grid timing probes and sample budgets for the N=2,4,7,8 sweep.

Probes normally measure one iteration of each stage; --feedback-iterations
also measures later self-consistent field laws. Only the standard library is required.
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
CONFIGURATION_SIZES = (2, 4, 7, 8)


def configurations():
    paths = sorted((ROOT / "Configuration_Data/CubicLattice").glob("Cubic_3D_N=*_NN_J=*0.40825.hdf5"))
    paths = [p for p in paths if int(re.search(r"N=(\d+)", p.stem)[1]) in CONFIGURATION_SIZES]
    if len(paths) != 2 * len(CONFIGURATION_SIZES):
        raise RuntimeError("Expected both coupling signs for N=2,4,7,8")
    return paths


def cases(config):
    return (("zero", "A", "unifield", "0"),
            ("hz0p5", "C", "unifield", "0.5"),
            ("ordering", "C", "unifield" if "J=-" in config else "stagfield", "0.005"))


def common(config, beta, case):
    _, symmetry, field, amplitude = case
    return [f"--config={config}", "--project=CubicLattice", "--spinspinmodel=ISO",
            f"--beta={beta}", "--numImagTimeSteps=50", "--numRealTimeSteps=150",
            "--Tmax=10", "--realTimeSubsteps=1", "--timeTranslationInvariant=true",
            "--samplingStrategy=independent", "--numSamplesPerSet=1",
            f"--cstype={symmetry}", f"--{field}={amplitude}"]


def duration(log, label):
    values = re.findall(re.escape(label) + r" done \(0\).*?it took me ([\deE.+-]+) seconds", log)
    if not values:
        raise RuntimeError(f"Missing timing: {label}")
    return max(map(float, values))


def statistical_memory_allowance(log, samples, sites, symmetry):
    if samples >= 64 and samples % 32 == 0:
        return 0  # The probe already allocates all 32 blocks in both stages.
    pairs = re.search(r"numbering starts at 1\): ([^\n]+)", log)
    if not pairs:
        raise RuntimeError("Missing correlation categories for memory accounting")
    categories = len(re.findall(r"\d+-\d+", pairs[1]))
    directions, moments = (1, 0) if symmetry == "A" else (4, 1)
    width = 1 + 151 * categories * directions * 51 + 151 * sites * moments + 151
    preparation_width = 2 + categories * directions * 51 + sites * moments
    # Every valid probe has at least two blocks; allowing 30 extra blocks
    # conservatively covers production even if its sample count differs.
    return 30 * 16 * (width + preparation_width)


def make_plan(rows, output, hours):
    # Inflate measured times 2x and normally reserve another 20% for MPI/I/O.
    # If the minimum four samples need that reserve, retain them only when
    # their doubled-cost estimate still fits the full wall-time allocation.
    safety, available = 2.0, hours * 3600 * .8
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
                prep_limit = 20

                def preparation_count(count):
                    return max(4, count // 20 * 2)

                def total_seconds(count):
                    return (initialization + 10 * (rs + rt * count)
                            + prep_limit * (ps + prep * preparation_count(count)))

                # Search even counts so both stages admit equal statistical
                # blocks, including the low-count seven/eight-spin cases.
                low, high = 2, max(2, math.floor(available / (20 * rt)))
                while low < high:
                    middle = (low + high + 1) // 2
                    if total_seconds(2 * middle) <= available:
                        low = middle
                    else:
                        high = middle - 1
                samples = 2 * low
                quantum = 100 if samples >= 1000 else 10 if samples >= 100 else 2
                samples = samples // quantum * quantum
                prep_samples = preparation_count(samples)
                status = "ready"
                estimate = total_seconds(samples)
                if estimate > hours * 3600:
                    status = "over_budget"
                if max(r.get("production_peak_rss_bytes", r["peak_child_rss_bytes"])
                       for r in probes) > 4_000_000_000:
                    status = "over_memory"
                writer.writerow(dict(zip(fields, [config.stem, *case, samples, prep_samples, prep_limit,
                    round(estimate, 1), rt, prep, rs, ps, initialization,
                    ";".join(sorted(set(r["config"] + "/" + r["case"] for r in probes))), status])))
    temporary.replace(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int,
                        help="samples per stage; defaults to 256 for N=2/4 and 4 for N=7/8")
    parser.add_argument("--hours", type=float, default=8,
                        help="job wall-time budget; normally reserves 20%% for MPI/I/O")
    parser.add_argument("--betas", nargs="+", default=["5"])
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--feedback-iterations", type=int, default=0,
                        help="enable self consistency with this real-time iteration cap")
    parser.add_argument("--sites", nargs="+", type=int, choices=CONFIGURATION_SIZES)
    parser.add_argument("--cases", nargs="+", choices=["zero", "hz0p5", "ordering"], default=["zero", "hz0p5", "ordering"])
    parser.add_argument("--sign", choices=["FM", "AFM", "both"], default="both")
    parser.add_argument("--output", type=Path, default=ROOT / "Timing/cubic_beta_sweep")
    parser.add_argument("--executable", type=Path, default=ROOT / "executable_DOUBLE.out")
    args = parser.parse_args()
    if args.samples is not None and (args.samples < 4 or args.samples % 2):
        parser.error("--samples must be even and at least 4")
    if args.hours <= 0 or not math.isfinite(args.hours):
        parser.error("--hours must be positive and finite")
    if args.feedback_iterations < 0:
        parser.error("--feedback-iterations must be nonnegative")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    record = args.output / "measurements.json"
    rows = json.loads(record.read_text()) if record.exists() else []
    executable = args.executable.resolve()
    fingerprint = hashlib.sha256(executable.read_bytes()).hexdigest()
    # Never pool timings from older binaries, grids or sampling settings.
    def current(row):
        return (row.get("executable_sha256") == fingerprint
                and row.get("num_imag_time_steps") == 50
                and row.get("num_real_time_steps") == 150
                and row.get("num_blocks") == 32
                and row.get("time_translation_invariant") is True)
    if not args.plan_only:
        env = dict(os.environ)
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "BLAZE_NUM_THREADS"):
            env[key] = "1"
        for config in configurations():
            sites = int(re.search(r"N=(\d+)", config.stem)[1])
            probe_samples = args.samples if args.samples is not None else (4 if sites >= 7 else 256)
            if args.sign != "both" and ("FM" if "J=-" in config.stem else "AFM") != args.sign:
                continue
            if args.sites and int(re.search(r"N=(\d+)", config.stem)[1]) not in args.sites:
                continue
            for case in cases(config.stem):
                if case[0] not in args.cases:
                    continue
                for beta in args.betas:
                    tag = f"{config.stem}_{case[0]}_beta{beta}_ni50_nr150_b32_tti_n{probe_samples}_{fingerprint[:8]}"
                    if args.feedback_iterations:
                        tag += f"_feedback{args.feedback_iterations}"
                    if any(r["tag"] == tag and current(r) for r in rows):
                        continue
                    command = [str(executable), *common(str(config), beta, case),
                        "--dstproject=Timing_cubic_sweep", f"--numSamplesPerCore={probe_samples}",
                        f"--preparationSamplesPerCore={probe_samples}", "--numBlocks=32",
                        "--seed=613", f"--fileext={tag}"]
                    if args.feedback_iterations:
                        command += [f"--iterlimit={args.feedback_iterations}",
                                    "--preparationIterlimit=20", "--reliterror=5"]
                    else:
                        command += ["--noselfcons"]
                    print(f"Timing {tag}", flush=True)
                    start = time.monotonic()
                    logfile = args.output / (tag + ".log")
                    with logfile.open("w") as log:
                        result = subprocess.run(command, cwd=ROOT, env=env, stdout=log,
                                                stderr=subprocess.STDOUT, timeout=args.timeout)
                    if result.returncode != 0 and not (args.feedback_iterations and result.returncode == 3):
                        raise RuntimeError(f"Probe failed with status {result.returncode}: {logfile}")
                    log = re.sub(r"\x1b\[[0-9;]*m", "", logfile.read_text())
                    rows.append(dict(tag=tag, config=config.stem, case=case[0], beta=beta,
                        samples=probe_samples, num_imag_time_steps=50, num_real_time_steps=150,
                        num_blocks=32, time_translation_invariant=True,
                        command=command, returncode=result.returncode,
                        self_consistency=bool(args.feedback_iterations),
                        preparation_iterations=len(re.findall(r"--- Matsubara iteration", log)),
                        real_time_iterations=len(re.findall(r"--- Mixed/real-time iteration", log)),
                        executable_sha256=fingerprint,
                        elapsed_seconds=time.monotonic() - start,
                        initialization_seconds=duration(log, "initialization"),
                        # RUSAGE_CHILDREN is the largest completed-child peak so far.
                        peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * (1 if sys.platform == "darwin" else 1024),
                        real_sample_seconds=duration(log, "mixed real-time sampling") / probe_samples,
                        preparation_sample_seconds=duration(log, "Matsubara sampling") / probe_samples,
                        real_setup_seconds=duration(log, "joint field distribution"),
                        preparation_setup_seconds=duration(log, "Matsubara field distribution")))
                    rows[-1]["production_peak_rss_bytes"] = (
                        rows[-1]["peak_child_rss_bytes"]
                        + statistical_memory_allowance(log, probe_samples, sites, case[1]))
                    record.write_text(json.dumps(rows, indent=2) + "\n")
                    print(f"  real {rows[-1]['real_sample_seconds']:.4g} s/sample; prep {rows[-1]['preparation_sample_seconds']:.4g} s/sample", flush=True)
    current_rows = [r for r in rows if current(r)]
    if args.plan_only or all(any(re.search(r"N=(\d+)", r["config"])[1] == str(n)
        and (r["case"] == "zero") == zero for r in current_rows) for n in CONFIGURATION_SIZES for zero in (True, False)):
        make_plan(current_rows, args.output / "sample_plan.tsv", args.hours)
        print(f"Wrote {args.output / 'sample_plan.tsv'}")


if __name__ == "__main__":
    main()
