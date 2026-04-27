#pragma once

#include "laph/solver.hpp"

#include <vector>

namespace laph {

struct DSU {
    std::vector<int> p;
    std::vector<int> r;

    explicit DSU(int n = 0);
    int find(int x);
    void unite(int a, int b);
};

std::vector<std::vector<int>> partition_components_for_problem(
    int total_vars,
    const PhasePoly& phase,
    const std::vector<Constraint>& constraints
);
Monomial remap_monomial(const Monomial& mon, const std::vector<int>& local_index);
XorSet remap_xorset(const XorSet& xs, const std::vector<int>& local_index);
AffineForm remap_affine_form(const AffineForm& f, const std::vector<int>& local_index);
ScaledComplex exact_partition_sum_factorized(
    int total_vars,
    const PhasePoly& phase,
    const std::vector<Constraint>& constraints
);

enum class PartitionBackend {
    Monolithic,
    Factorized
};

struct QueryOptions {
    PartitionBackend backend = PartitionBackend::Factorized;
};

ScaledComplex exact_partition_sum_with_backend(
    int total_vars,
    const PhasePoly& phase,
    const std::vector<Constraint>& constraints,
    PartitionBackend backend
);

struct Stats {
    int qubits = 0;
    int latent_variables = 0;
    size_t constraints = 0;
    size_t phase_terms = 0;
    int hadamard_scale = 0;
    size_t non_clifford_cut = 0;
    size_t hidden_interference_rank = 0;
    int fabqnf_rho = 0;
    int fabqnf_out_dim = 0;
    int fabqnf_hid_dim = 0;
    int fabqnf_table_size = 0;
    size_t fiber_active_terms = 0;
    size_t fiber_active_vars = 0;
    size_t components = 0;
};

struct StateComponent {
    std::vector<int> variables;
    std::vector<int> qubits;
};

} // namespace laph
