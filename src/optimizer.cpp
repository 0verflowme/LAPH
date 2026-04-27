#include "laph/optimizer.hpp"

#include <algorithm>
#include <map>
#include <numeric>

namespace laph {

DSU::DSU(int n) : p(n), r(n, 0) {
    std::iota(p.begin(), p.end(), 0);
}

int DSU::find(int x) {
    return p[x] == x ? x : p[x] = find(p[x]);
}

void DSU::unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a == b) return;
    if (r[a] < r[b]) std::swap(a, b);
    p[b] = a;
    if (r[a] == r[b]) ++r[a];
}

std::vector<std::vector<int>> partition_components_for_problem(
    int total_vars,
    const PhasePoly& phase,
    const std::vector<Constraint>& constraints
) {
    if (total_vars <= 0) return {};

    DSU dsu(total_vars);
    auto connect_all = [&](const std::vector<int>& xs) {
        if (xs.empty()) return;
        int root = xs[0];
        for (size_t i = 1; i < xs.size(); ++i) dsu.unite(root, xs[i]);
    };

    for (const auto& kv : phase) connect_all(kv.first.v);
    for (const Constraint& c : constraints) connect_all(c.lhs.v);

    std::map<int, std::vector<int>> grouped;
    for (int x = 0; x < total_vars; ++x) grouped[dsu.find(x)].push_back(x);

    std::vector<std::vector<int>> components;
    components.reserve(grouped.size());
    for (auto& kv : grouped) components.push_back(std::move(kv.second));
    return components;
}

Monomial remap_monomial(const Monomial& mon, const std::vector<int>& local_index) {
    Monomial out;
    out.v.reserve(mon.v.size());
    for (int x : mon.v) out.v.push_back(local_index[x]);
    out.normalize_or();
    return out;
}

XorSet remap_xorset(const XorSet& xs, const std::vector<int>& local_index) {
    XorSet out;
    out.v.reserve(xs.v.size());
    for (int x : xs.v) out.v.push_back(local_index[x]);
    out.normalize_xor();
    return out;
}

AffineForm remap_affine_form(const AffineForm& f, const std::vector<int>& local_index) {
    return AffineForm(remap_xorset(f.vars, local_index), f.c);
}

ScaledComplex exact_partition_sum_factorized(
    int total_vars,
    const PhasePoly& phase,
    const std::vector<Constraint>& constraints
) {
    std::vector<std::vector<int>> components =
        partition_components_for_problem(total_vars, phase, constraints);

    if (components.size() <= 1) {
        return exact_partition_sum(total_vars, phase, constraints);
    }

    std::vector<int> var_to_component(total_vars, -1);
    std::vector<int> local_index(total_vars, -1);
    for (size_t ci = 0; ci < components.size(); ++ci) {
        for (size_t li = 0; li < components[ci].size(); ++li) {
            int x = components[ci][li];
            var_to_component[x] = static_cast<int>(ci);
            local_index[x] = static_cast<int>(li);
        }
    }

    std::vector<PhasePoly> component_phase(components.size());
    std::vector<std::vector<Constraint>> component_constraints(components.size());
    int constant_phase = 0;

    for (const auto& kv : phase) {
        int coeff = mod8(kv.second);
        if (coeff == 0) continue;

        if (kv.first.empty()) {
            constant_phase = mod8(constant_phase + coeff);
            continue;
        }

        int ci = var_to_component[kv.first.v[0]];
        add_mod8(component_phase[ci], remap_monomial(kv.first, local_index), coeff);
    }

    for (const Constraint& c : constraints) {
        if (c.lhs.empty()) {
            if (c.rhs) return ScaledComplex::zero();
            continue;
        }
        int ci = var_to_component[c.lhs.v[0]];
        component_constraints[ci].push_back({remap_xorset(c.lhs, local_index), c.rhs});
    }

    ScaledComplex total = ScaledComplex::omega_phase(constant_phase);
    for (size_t ci = 0; ci < components.size(); ++ci) {
        ScaledComplex z = exact_partition_sum(
            static_cast<int>(components[ci].size()),
            component_phase[ci],
            component_constraints[ci]
        );
        total.mul(z);
        if (total.is_zero()) break;
    }
    return total;
}

ScaledComplex exact_partition_sum_with_backend(
    int total_vars,
    const PhasePoly& phase,
    const std::vector<Constraint>& constraints,
    PartitionBackend backend
) {
    if (backend == PartitionBackend::Factorized) {
        return exact_partition_sum_factorized(total_vars, phase, constraints);
    }
    return exact_partition_sum(total_vars, phase, constraints);
}

} // namespace laph
