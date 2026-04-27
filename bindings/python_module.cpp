#include <cstdint>
#include <random>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "laph/laph.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(_laph, m) {
    m.doc() = "LAPH C++ core backend";

    py::class_<laph::ScaledComplex>(m, "ScaledComplex")
        .def("real", &laph::ScaledComplex::real_ld)
        .def("imag", &laph::ScaledComplex::imag_ld)
        .def_property_readonly("sqrt2_power", [](const laph::ScaledComplex& z) {
            return z.k;
        });

    py::class_<laph::LAPH>(m, "LAPH")
        .def(py::init<int>(), py::arg("nqubits"))
        .def("x", &laph::LAPH::x, py::arg("q"), py::return_value_policy::reference_internal)
        .def("h", &laph::LAPH::h, py::arg("q"), py::return_value_policy::reference_internal)
        .def("t", &laph::LAPH::t, py::arg("q"), py::return_value_policy::reference_internal)
        .def("s", &laph::LAPH::s, py::arg("q"), py::return_value_policy::reference_internal)
        .def("z", &laph::LAPH::z, py::arg("q"), py::return_value_policy::reference_internal)
        .def(
            "cnot",
            &laph::LAPH::cnot,
            py::arg("control"),
            py::arg("target"),
            py::return_value_policy::reference_internal
        )
        .def(
            "cz",
            &laph::LAPH::cz,
            py::arg("a"),
            py::arg("b"),
            py::return_value_policy::reference_internal
        )
        .def(
            "ccz",
            &laph::LAPH::ccz,
            py::arg("a"),
            py::arg("b"),
            py::arg("c"),
            py::return_value_policy::reference_internal
        )
        .def("compress", &laph::LAPH::compress)
        .def("cut_set", &laph::LAPH::cut_set)
        .def("amplitude", &laph::LAPH::amplitude_factorized, py::arg("basis_state"))
        .def("probability_prefix", &laph::LAPH::probability_prefix_factorized, py::arg("prefix"))
        .def(
            "exact_sample",
            [](const laph::LAPH& self, std::uint64_t seed) {
                std::mt19937_64 rng(seed);
                return self.exact_sample(rng);
            },
            py::arg("seed")
        )
        .def(
            "exact_sample",
            [](const laph::LAPH& self) {
                std::random_device rd;
                std::mt19937_64 rng(
                    (static_cast<std::uint64_t>(rd()) << 32) ^
                    static_cast<std::uint64_t>(rd())
                );
                return self.exact_sample(rng);
            }
        )
        .def(
            "exact_sample_density_oracle",
            [](const laph::LAPH& self, std::uint64_t seed) {
                std::mt19937_64 rng(seed);
                return self.exact_sample_density_oracle(rng);
            },
            py::arg("seed")
        )
        .def("stats", [](const laph::LAPH& self) {
            laph::Stats s = self.stats();
            return py::dict(
                "qubits"_a = s.qubits,
                "latent_variables"_a = s.latent_variables,
                "constraints"_a = s.constraints,
                "phase_terms"_a = s.phase_terms,
                "hadamard_scale"_a = s.hadamard_scale,
                "non_clifford_cut"_a = s.non_clifford_cut,
                "hidden_interference_rank"_a = s.hidden_interference_rank,
                "fabqnf_rho"_a = s.fabqnf_rho,
                "fabqnf_out_dim"_a = s.fabqnf_out_dim,
                "fabqnf_hid_dim"_a = s.fabqnf_hid_dim,
                "fabqnf_table_size"_a = s.fabqnf_table_size,
                "fiber_active_terms"_a = s.fiber_active_terms,
                "fiber_active_vars"_a = s.fiber_active_vars,
                "max_component_rho"_a = s.max_component_rho,
                "sum_component_table_size"_a = s.sum_component_table_size,
                "max_component_kappa"_a = s.max_component_kappa,
                "component_table_failures"_a = s.component_table_failures,
                "components"_a = s.components
            );
        });
}
