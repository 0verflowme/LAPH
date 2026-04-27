from qiskit import QuantumCircuit, transpile

from laph.qiskit_backend import LAPHBackend


def test_bell_counts():
    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    qc.measure([0, 1], [0, 1])

    backend = LAPHBackend(num_qubits=2)
    compiled = transpile(qc, backend)
    result = backend.run(compiled, shots=20, seed_simulator=123).result()
    counts = result.get_counts()

    assert set(counts) <= {"00", "11"}
    assert sum(counts.values()) == 20


def test_unmeasured_circuit_samples_all_qubits():
    qc = QuantumCircuit(1)
    qc.x(0)

    backend = LAPHBackend(num_qubits=1)
    result = backend.run(qc, shots=3, seed_simulator=123).result()

    assert result.get_counts() == {"1": 3}
