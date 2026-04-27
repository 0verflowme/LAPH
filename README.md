# LAPH

Lifted Affine-Phase Hypergraph (LAPH) is an exact, matrix-free quantum simulation
kernel. It represents a state as a sparse affine map, sparse GF(2) constraints,
and a sparse mod-8 phase hypergraph. Dense affine phases are lifted into explicit
latent variables so the exponential driver is the residual non-Clifford cut, not
the visible qubit count.

This repository promotes the C++ v0.2 prototype into a reusable C++17 library
with optional Python bindings.

## Layout

```text
include/laph/types.hpp       value types and sparse algebra containers
include/laph/phase.hpp       phase polynomial and Clifford summation API
include/laph/solver.hpp      GF(2) solving and partition sums
include/laph/optimizer.hpp   components, remapping, and backend selection
include/laph/laph.hpp        public simulator API
src/*.cpp                    compiled laph_core implementation
bindings/python_module.cpp   pybind11 extension module
python/laph/__init__.py      Python package facade
```

## Current library surface

- sparse GF(2) affine forms and constraints
- sparse mod-8 phase polynomial with modulo folding
- lifted affine variables and lift caching
- exact amplitude queries
- exact prefix-probability queries using the doubled-density construction
- exact sequential sampling
- component-local exact sampling
- greedy non-Clifford cut selection
- connected-component detection
- component-factorized partition queries
- scaled complex numbers with sqrt(2) exponents

## Build

```bash
make
make test
make run
```

Or with CMake:

```bash
cmake -S . -B build-cmake
cmake --build build-cmake
ctest --test-dir build-cmake --output-on-failure
```

Install the Python package with:

```bash
pip install .
python -c "import laph; print(laph.LAPH(1).h(0).stats())"
```

Install the optional Qiskit adapter with:

```bash
pip install ".[qiskit]"
```

## Minimal use

```cpp
#include "laph/laph.hpp"

int main() {
    laph::LAPH st(3);
    st.h(0).h(1).cnot(0, 2).t(2).cz(0, 1).h(2).ccz(0, 1, 2);
    st.compress();

    laph::ScaledComplex amp = st.amplitude_factorized(0b101);
    long double p0 = st.probability_prefix_factorized({{0, false}});
}
```

The default query backend is component-factorized. The monolithic backend remains
available for debugging:

```cpp
auto a = st.amplitude(0b101, {laph::PartitionBackend::Monolithic});
auto b = st.amplitude(0b101, {laph::PartitionBackend::Factorized});
```

For circuits that split into independent latent components, sample each component
locally and concatenate the visible bits:

```cpp
std::mt19937_64 rng(123);
auto sample = st.exact_sample_by_component(rng);
```

## Python Use

```python
import laph

st = laph.LAPH(3)
st.h(0).h(1).cnot(0, 2).t(2).cz(0, 1).h(2).ccz(0, 1, 2)
st.compress()

print(st.stats())
print(st.exact_sample(seed=123))
amp = st.amplitude(0b001)
print(amp.real(), amp.imag())
```

## Qiskit Use

```python
from qiskit import QuantumCircuit, transpile
from laph.qiskit_backend import LAPHBackend

qc = QuantumCircuit(2, 2)
qc.h(0)
qc.cx(0, 1)
qc.measure([0, 1], [0, 1])

backend = LAPHBackend(num_qubits=2)
compiled = transpile(qc, backend)
job = backend.run(compiled, shots=1000, seed_simulator=123)
print(job.result().get_counts())
```

## Notes

The current implementation is still intentionally sparse and clarity-oriented.
The next production steps are component-local variable renumbering at mutation
time, Gray-code cut enumeration, SIMD/dense GF(2) kernels for dense regions,
thread-parallel cut loops, and stronger compression passes for equivalent lifts
and singleton substitutions.
