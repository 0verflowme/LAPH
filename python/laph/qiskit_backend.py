from __future__ import annotations

from collections import Counter
import random
from typing import Iterable

from qiskit.circuit import QuantumCircuit
from qiskit.circuit import Measure
from qiskit.circuit.library import CCZGate, CXGate, CZGate, HGate, SGate, TGate, XGate, ZGate
from qiskit.providers import BackendV2, JobStatus, JobV1, Options
from qiskit.result import Result
from qiskit.result.models import ExperimentResult, ExperimentResultData
from qiskit.transpiler import Target

from ._laph import LAPH


class LAPHJob(JobV1):
    """Synchronous Qiskit job wrapper for an already-computed LAPH result."""

    def __init__(self, backend: BackendV2, job_id: str, result: Result):
        super().__init__(backend, job_id)
        self._result = result

    def submit(self):
        return None

    def result(self):
        return self._result

    def status(self):
        return JobStatus.DONE


class LAPHBackend(BackendV2):
    """Qiskit BackendV2 adapter for the LAPH exact sampler."""

    def __init__(self, num_qubits: int | None = None):
        super().__init__(
            name="laph_simulator",
            description="Lifted Affine-Phase Hypergraph exact simulator",
            backend_version="0.3.0",
        )
        target_kwargs = {}
        if num_qubits is not None:
            target_kwargs["num_qubits"] = num_qubits

        self._target = Target(**target_kwargs)
        self._target.add_instruction(XGate())
        self._target.add_instruction(HGate())
        self._target.add_instruction(TGate())
        self._target.add_instruction(SGate())
        self._target.add_instruction(ZGate())
        self._target.add_instruction(CXGate())
        self._target.add_instruction(CZGate())
        self._target.add_instruction(CCZGate())
        self._target.add_instruction(Measure())

    @property
    def target(self):
        return self._target

    @property
    def max_circuits(self):
        return None

    @classmethod
    def _default_options(cls):
        return Options(shots=1024, seed_simulator=None)

    def run(self, run_input, **options):
        run_options = dict(self.options)
        run_options.update(options)

        circuits = list(run_input) if isinstance(run_input, Iterable) and not isinstance(run_input, QuantumCircuit) else [run_input]
        seed = run_options.get("seed_simulator")
        rng = random.Random(seed)
        shots = int(run_options.get("shots", 1024))

        results = [
            self._run_circuit(circuit, shots=shots, rng=rng, seed=seed)
            for circuit in circuits
        ]

        result = Result(
            backend_name=self.name,
            backend_version=self.backend_version,
            job_id="laph-job-0",
            success=True,
            results=results,
        )
        return LAPHJob(self, "laph-job-0", result)

    def _run_circuit(self, circuit: QuantumCircuit, shots: int, rng: random.Random, seed):
        sim = LAPH(circuit.num_qubits)
        measure_map: list[tuple[int, int]] = []

        for instruction in circuit.data:
            operation = instruction.operation
            name = operation.name
            qubits = [circuit.find_bit(qubit).index for qubit in instruction.qubits]
            clbits = [circuit.find_bit(clbit).index for clbit in instruction.clbits]

            if name == "x":
                sim.x(qubits[0])
            elif name == "h":
                sim.h(qubits[0])
            elif name == "t":
                sim.t(qubits[0])
            elif name == "s":
                sim.s(qubits[0])
            elif name == "z":
                sim.z(qubits[0])
            elif name == "cx":
                sim.cnot(qubits[0], qubits[1])
            elif name == "cz":
                sim.cz(qubits[0], qubits[1])
            elif name == "ccz":
                sim.ccz(qubits[0], qubits[1], qubits[2])
            elif name == "measure":
                measure_map.append((qubits[0], clbits[0]))
            elif name in {"barrier", "delay"}:
                continue
            else:
                raise NotImplementedError(f"LAPHBackend does not support instruction {name!r}")

        sim.compress()
        counts = Counter()

        for _ in range(shots):
            shot_seed = rng.getrandbits(64)
            sample = sim.exact_sample(shot_seed)
            counts[self._format_memory(sample, measure_map, circuit.num_clbits)] += 1

        return ExperimentResult(
            shots=shots,
            success=True,
            data=ExperimentResultData(counts=dict(counts)),
            seed=seed,
            header={"name": circuit.name},
        )

    @staticmethod
    def _format_memory(sample: list[int], measure_map: list[tuple[int, int]], num_clbits: int):
        if not measure_map:
            return "".join(str(bit) for bit in reversed(sample))

        classical = [0] * num_clbits
        for qubit, clbit in measure_map:
            classical[clbit] = sample[qubit]

        return "".join(str(bit) for bit in reversed(classical))
