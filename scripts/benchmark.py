#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
from pathlib import Path
import random
import statistics
import time

from qiskit import QuantumCircuit, transpile

from laph.qiskit_backend import LAPHBackend, build_laph_from_circuit, sample_laph_counts


def parse_int_list(value: str) -> list[int]:
    return [int(x.strip()) for x in value.split(",") if x.strip()]


def random_clifford_t_circuit(
    nqubits: int,
    depth: int,
    *,
    seed: int,
    t_probability: float,
    entangle: str,
    measure: bool,
) -> QuantumCircuit:
    rng = random.Random(seed)
    circuit = QuantumCircuit(nqubits, nqubits if measure else 0)
    clifford_1q = ("h", "s", "z", "x")

    for layer in range(depth):
        for q in range(nqubits):
            if rng.random() < t_probability:
                circuit.t(q)
            else:
                gate = rng.choice(clifford_1q)
                getattr(circuit, gate)(q)

        if nqubits > 1:
            if entangle == "chain":
                start = layer & 1
                for q in range(start, nqubits - 1, 2):
                    circuit.cx(q, q + 1)
            elif entangle == "ring":
                start = layer & 1
                for q in range(start, nqubits - 1, 2):
                    circuit.cx(q, q + 1)
                if nqubits > 2 and layer % 2 == 0:
                    circuit.cx(nqubits - 1, 0)
            elif entangle == "random":
                qubits = list(range(nqubits))
                rng.shuffle(qubits)
                for a, b in zip(qubits[0::2], qubits[1::2]):
                    circuit.cx(a, b)
            elif entangle == "none":
                pass
            else:
                raise ValueError(f"unknown entangler pattern: {entangle}")

    if measure:
        circuit.measure(range(nqubits), range(nqubits))

    return circuit


def timed(callable_):
    start = time.perf_counter()
    value = callable_()
    return value, time.perf_counter() - start


def percentile(values: list[float], pct: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    idx = min(len(ordered) - 1, max(0, round((len(ordered) - 1) * pct)))
    return ordered[idx]


def measure_sample_latencies(sim, measure_map, num_clbits: int, shots: int, seed: int):
    rng = random.Random(seed)
    latencies = []
    counts = {}

    for _ in range(shots):
        shot_seed = rng.getrandbits(64)
        start = time.perf_counter()
        sample = sim.exact_sample(shot_seed)
        latencies.append(time.perf_counter() - start)
        bitstring = LAPHBackend.format_memory(sample, measure_map, num_clbits)
        counts[bitstring] = counts.get(bitstring, 0) + 1

    return latencies, counts


def run_case(args, nqubits: int, depth: int, trial: int):
    seed = args.seed + 1_000_003 * trial + 10_007 * nqubits + depth
    circuit, generate_s = timed(lambda: random_clifford_t_circuit(
        nqubits,
        depth,
        seed=seed,
        t_probability=args.t_probability,
        entangle=args.entangle,
        measure=True,
    ))

    backend = LAPHBackend(num_qubits=nqubits)
    compiled, transpile_s = timed(lambda: transpile(circuit, backend, optimization_level=0))
    (sim, measure_map), build_compress_s = timed(lambda: build_laph_from_circuit(compiled))

    stats = sim.stats()
    latencies, counts = measure_sample_latencies(
        sim,
        measure_map,
        compiled.num_clbits,
        args.shots,
        seed + 17,
    )

    _, backend_run_s = timed(lambda: backend.run(compiled, shots=args.shots, seed_simulator=seed + 29).result())

    row = {
        "nqubits": nqubits,
        "depth": depth,
        "trial": trial,
        "shots": args.shots,
        "t_probability": args.t_probability,
        "entangle": args.entangle,
        "circuit_ops": len(compiled.data),
        "generate_s": generate_s,
        "transpile_s": transpile_s,
        "build_compress_s": build_compress_s,
        "backend_run_s": backend_run_s,
        "sample_total_s": sum(latencies),
        "sample_mean_s": statistics.fmean(latencies) if latencies else 0.0,
        "sample_p50_s": percentile(latencies, 0.50),
        "sample_p95_s": percentile(latencies, 0.95),
        "distinct_counts": len(counts),
    }
    row.update(stats)
    return row


def print_row(row: dict):
    print(
        "n={nqubits:4d} depth={depth:4d} trial={trial:2d} "
        "tau={non_clifford_cut:3d} m={latent_variables:6d} "
        "constraints={constraints:6d} terms={phase_terms:6d} comps={components:4d} "
        "build={build_compress_s:8.4f}s sample_mean={sample_mean_s:8.4f}s "
        "backend={backend_run_s:8.4f}s".format(**row),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description="Benchmark the LAPH Qiskit backend.")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument(
        "--family",
        choices=["connected_clifford", "connected_sparse_t", "disconnected_sparse_t"],
    )
    parser.add_argument("--qubits", nargs="+", type=int)
    parser.add_argument("--widths", default="4,8,12,16,24,32")
    parser.add_argument("--depth", type=int)
    parser.add_argument("--depths", default="4,8,16")
    parser.add_argument("--trials", type=int, default=1)
    parser.add_argument("--shots", type=int, default=1)
    parser.add_argument("--t-probability", type=float, default=0.15)
    parser.add_argument("--entangle", choices=["none", "chain", "ring", "random"], default="chain")
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--max-case-seconds", type=float, default=30.0)
    parser.add_argument("--case-timeout", type=float, default=0.0)
    parser.add_argument("--stop-after-slow", type=int, default=2)
    parser.add_argument("--output", default="benchmark_results/laph_benchmark.csv")
    default_output = parser.get_default("output")
    args = parser.parse_args()

    if args.smoke:
        args.widths = "4,8"
        args.depths = "4"
        args.trials = 1
        args.shots = 1
        args.t_probability = 0.1
        args.entangle = "chain"
        if args.output == default_output:
            args.output = "benchmark_results/smoke.csv"

    if args.family == "connected_clifford":
        args.t_probability = 0.0
        args.entangle = "chain"
    elif args.family == "connected_sparse_t":
        args.t_probability = 0.05
        args.entangle = "chain"
    elif args.family == "disconnected_sparse_t":
        args.t_probability = 0.05
        args.entangle = "none"

    widths = args.qubits if args.qubits is not None else parse_int_list(args.widths)
    depths = [args.depth] if args.depth is not None else parse_int_list(args.depths)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    slow_cases = 0

    for nqubits in widths:
        for depth in depths:
            for trial in range(args.trials):
                start = time.perf_counter()
                row = run_case_with_timeout(args, nqubits, depth, trial)

                row["case_total_s"] = time.perf_counter() - start
                rows.append(row)
                write_csv(output, rows)

                if row.get("error"):
                    print(json.dumps(row, sort_keys=True), flush=True)
                    slow_cases += 1
                else:
                    print_row(row)
                    if row["case_total_s"] > args.max_case_seconds:
                        slow_cases += 1
                    else:
                        slow_cases = 0

                if slow_cases >= args.stop_after_slow:
                    print("stopping after consecutive slow/failing cases", flush=True)
                    return

    write_csv(output, rows)


def run_case_with_timeout(args, nqubits: int, depth: int, trial: int):
    if args.case_timeout <= 0:
        try:
            row = run_case(args, nqubits, depth, trial)
            row["error"] = ""
            return row
        except Exception as exc:  # noqa: BLE001 - benchmark records failures.
            return error_row(args, nqubits, depth, trial, repr(exc))

    queue: mp.Queue = mp.Queue(maxsize=1)
    process = mp.Process(target=run_case_child, args=(queue, args, nqubits, depth, trial))
    process.start()
    process.join(args.case_timeout)

    if process.is_alive():
        process.terminate()
        process.join(5)
        return error_row(args, nqubits, depth, trial, f"timeout after {args.case_timeout}s")

    if queue.empty():
        return error_row(args, nqubits, depth, trial, f"child exited with code {process.exitcode}")

    return queue.get()


def run_case_child(queue, args, nqubits: int, depth: int, trial: int):
    try:
        row = run_case(args, nqubits, depth, trial)
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001 - benchmark records failures.
        row = error_row(args, nqubits, depth, trial, repr(exc))
    queue.put(row)


def error_row(args, nqubits: int, depth: int, trial: int, error: str):
    return {
        "nqubits": nqubits,
        "depth": depth,
        "trial": trial,
        "shots": args.shots,
        "t_probability": args.t_probability,
        "entangle": args.entangle,
        "error": error,
    }


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
