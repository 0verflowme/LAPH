#include "laph/laph.hpp"

#include "laph/density_kernel.hpp"
#include "laph/math/born_clifford_sampler.hpp"
#include "laph/math/fabqnf.hpp"
#include "laph/math/hidden_rank.hpp"
#include "laph/phase.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>

namespace laph {

namespace {

struct ComponentStatsSummary {
    int max_rho = 0;
    long long sum_table_size = 0;
    int max_kappa = 0;
    int table_failures = 0;
};

ComponentStatsSummary component_fabqnf_summary(
    const LAPH& st,
    const std::vector<StateComponent>& components
) {
    ComponentStatsSummary out;

    for (const StateComponent& component : components) {
        if (component.qubits.empty()) continue;

        LAPH local = st.component_view(component);
        std::vector<int> cut = local.cut_set();
        int kappa = hidden_interference_rank(local, cut);
        out.max_kappa = std::max(out.max_kappa, kappa);

        FABQNFStats fs;
        if (is_clifford_poly(local.phase)) {
            FiberChart ch = build_fiber_chart(local);
            fs.out_dim = ch.out_dim;
            fs.hid_dim = ch.hid_dim;
            fs.rho = 0;
            fs.table_size = 1;
        } else {
            fs = fabqnf_stats(local);
        }

        out.max_rho = std::max(out.max_rho, fs.rho);
        if (fs.table_size < 0) ++out.table_failures;
        else out.sum_table_size += fs.table_size;
    }

    return out;
}

} // namespace

LAPH::LAPH(int nqubits)
    : n(nqubits),
      visible(nqubits, AffineForm::constant(false)),
      tableau(nqubits) {}

int LAPH::new_var() {
    return m++;
}

std::pair<bool, int> LAPH::lift(const AffineForm& f) {
    if (f.vars.empty()) return {false, f.c ? 1 : 0};
    if (!f.c && f.vars.size() == 1) return {true, f.vars.v[0]};

    auto it = lift_cache.find(f);
    if (it != lift_cache.end()) return {true, it->second};

    int u = new_var();
    XorSet lhs = f.vars;
    lhs.xor_with(XorSet::singleton(u));
    constraints.push_back({lhs, f.c});
    lift_cache[f] = u;
    return {true, u};
}

void LAPH::add_phase_monomial(int coeff, const Monomial& mon) {
    add_mod8(phase, mon, coeff);
}

void LAPH::add_phase_product_forms(int coeff, const std::vector<AffineForm>& forms) {
    Monomial mon;
    for (const AffineForm& f : forms) {
        auto [is_var, val] = lift(f);
        if (!is_var) {
            if (val == 0) return;
            continue;
        }
        mon.or_with_var(val);
    }
    add_phase_monomial(coeff, mon);
}

LAPH& LAPH::x(int q) {
    if (tableau_valid) tableau.x_gate(q);
    visible[q].c ^= true;
    return *this;
}

LAPH& LAPH::cnot(int c, int t) {
    if (tableau_valid) tableau.cnot(c, t);
    visible[t].xor_with(visible[c]);
    return *this;
}

LAPH& LAPH::h(int q) {
    if (tableau_valid) tableau.h(q);
    AffineForm old = visible[q];
    int s = new_var();
    add_phase_product_forms(4, {AffineForm::variable(s), old});
    visible[q] = AffineForm::variable(s);
    scale += 1;
    return *this;
}

LAPH& LAPH::t(int q) {
    tableau_valid = false;
    add_phase_product_forms(1, {visible[q]});
    return *this;
}

LAPH& LAPH::s(int q) {
    if (tableau_valid) tableau.s(q);
    add_phase_product_forms(2, {visible[q]});
    return *this;
}

LAPH& LAPH::z(int q) {
    if (tableau_valid) tableau.z_gate(q);
    add_phase_product_forms(4, {visible[q]});
    return *this;
}

LAPH& LAPH::cz(int a, int b) {
    if (tableau_valid) tableau.cz(a, b);
    add_phase_product_forms(4, {visible[a], visible[b]});
    return *this;
}

LAPH& LAPH::ccz(int a, int b, int c) {
    tableau_valid = false;
    add_phase_product_forms(4, {visible[a], visible[b], visible[c]});
    return *this;
}

std::vector<Constraint> LAPH::output_constraints(uint64_t xbits) const {
    std::vector<Constraint> out = constraints;
    out.reserve(out.size() + n);
    for (int q = 0; q < n; ++q) {
        bool bit = ((xbits >> q) & 1ull) != 0;
        out.push_back({visible[q].vars, static_cast<bool>(bit ^ visible[q].c)});
    }
    return out;
}

std::vector<int> LAPH::cut_set() const {
    return greedy_cut_set_for_phase(phase);
}

ScaledComplex LAPH::amplitude(uint64_t xbits, QueryOptions options) const {
    std::vector<Constraint> cs = output_constraints(xbits);
    ScaledComplex z = exact_partition_sum_with_backend(m, phase, cs, options.backend);
    z.mul_sqrt2_power(-scale);
    return z;
}

ScaledComplex LAPH::amplitude_monolithic(uint64_t xbits) const {
    return amplitude(xbits, QueryOptions{PartitionBackend::Monolithic});
}

ScaledComplex LAPH::amplitude_factorized(uint64_t xbits) const {
    return amplitude(xbits, QueryOptions{PartitionBackend::Factorized});
}

PhasePoly LAPH::density_phase() const {
    PhasePoly ph;
    for (const auto& kv : phase) {
        add_mod8(ph, kv.first, kv.second);
        add_mod8(ph, kv.first.shifted(m), -kv.second);
    }
    return ph;
}

std::vector<Constraint> LAPH::density_constraints(
    const std::vector<std::pair<int, bool>>& prefix
) const {
    std::vector<Constraint> cs;
    cs.reserve(2 * constraints.size() + n + prefix.size());

    for (const auto& c : constraints) {
        cs.push_back(c);
        cs.push_back({c.lhs.shifted(m), c.rhs});
    }

    for (int q = 0; q < n; ++q) {
        XorSet lhs = visible[q].vars;
        lhs.xor_with(visible[q].vars.shifted(m));
        cs.push_back({lhs, false});
    }

    for (auto [q, bit] : prefix) {
        cs.push_back({visible[q].vars, static_cast<bool>(bit ^ visible[q].c)});
    }

    return cs;
}

ScaledComplex LAPH::density_prefix_sum(
    const std::vector<std::pair<int, bool>>& prefix,
    QueryOptions options
) const {
    PhasePoly ph = density_phase();
    std::vector<Constraint> cs = density_constraints(prefix);
    return exact_partition_sum_with_backend(2 * m, ph, cs, options.backend);
}

long double LAPH::density_prefix_weight(
    const std::vector<std::pair<int, bool>>& prefix,
    QueryOptions options
) const {
    ScaledComplex z = density_prefix_sum(prefix, options);
    long double p = z.real_ld();
    if (p < 0 && p > -1e-12L) p = 0;
    return p;
}

long double LAPH::probability_prefix(
    const std::vector<std::pair<int, bool>>& prefix,
    QueryOptions options
) const {
    ScaledComplex z = density_prefix_sum(prefix, options);
    z.mul_sqrt2_power(-2LL * scale);
    long double p = z.real_ld();
    if (p < 0 && p > -1e-12L) p = 0;
    if (p > 1 && p < 1 + 1e-12L) p = 1;
    return p;
}

long double LAPH::probability_prefix_monolithic(
    const std::vector<std::pair<int, bool>>& prefix
) const {
    return probability_prefix(prefix, QueryOptions{PartitionBackend::Monolithic});
}

long double LAPH::probability_prefix_factorized(
    const std::vector<std::pair<int, bool>>& prefix
) const {
    return probability_prefix(prefix, QueryOptions{PartitionBackend::Factorized});
}

std::vector<int> LAPH::exact_sample(
    std::mt19937_64& rng,
    QueryOptions options
) const {
    if (is_clifford_poly(phase)) return sample_born_clifford(*this, rng);
    if (options.backend == PartitionBackend::Factorized &&
        connected_components().size() > 1) {
        return exact_sample_by_component(rng, options);
    }
    return build_fabqnf(*this).sample(rng);
}

std::vector<int> LAPH::exact_sample_density_oracle(std::mt19937_64& rng) const {
    if (is_clifford_poly(phase)) return sample_born_clifford(*this, rng);
    return exact_sample_cached_density(*this, rng);
}

std::vector<int> LAPH::exact_sample_monolithic(std::mt19937_64& rng) const {
    return exact_sample(rng, QueryOptions{PartitionBackend::Monolithic});
}

std::vector<int> LAPH::exact_sample_factorized(std::mt19937_64& rng) const {
    return exact_sample(rng, QueryOptions{PartitionBackend::Factorized});
}

std::vector<StateComponent> LAPH::state_components() const {
    std::vector<std::vector<int>> vars = connected_components();
    std::vector<StateComponent> out(vars.size());
    std::vector<int> var_to_component(std::max(0, m), -1);

    for (size_t ci = 0; ci < vars.size(); ++ci) {
        out[ci].variables = std::move(vars[ci]);
        for (int x : out[ci].variables) var_to_component[x] = static_cast<int>(ci);
    }

    for (int q = 0; q < n; ++q) {
        if (visible[q].vars.empty()) continue;
        int ci = var_to_component[visible[q].vars.v[0]];
        out[ci].qubits.push_back(q);
    }

    return out;
}

LAPH LAPH::component_view(const StateComponent& component) const {
    LAPH local(static_cast<int>(component.qubits.size()));
    local.m = static_cast<int>(component.variables.size());
    local.tableau_valid = false;

    std::vector<int> local_index(std::max(0, m), -1);
    for (size_t i = 0; i < component.variables.size(); ++i) {
        local_index[component.variables[i]] = static_cast<int>(i);
    }

    for (size_t i = 0; i < component.qubits.size(); ++i) {
        local.visible[i] = remap_affine_form(visible[component.qubits[i]], local_index);
    }

    for (const Constraint& c : constraints) {
        if (c.lhs.empty()) {
            if (c.rhs) throw std::runtime_error("inconsistent LAPH constraints");
            continue;
        }
        if (local_index[c.lhs.v[0]] >= 0) {
            local.constraints.push_back({remap_xorset(c.lhs, local_index), c.rhs});
        }
    }

    for (const auto& kv : phase) {
        if (kv.first.empty()) continue;
        if (local_index[kv.first.v[0]] >= 0) {
            add_mod8(local.phase, remap_monomial(kv.first, local_index), kv.second);
        }
    }

    return local;
}

std::vector<int> LAPH::exact_sample_by_component(
    std::mt19937_64& rng,
    QueryOptions
) const {
    if (tableau_valid) return tableau.sample_all(rng);

    for (const Constraint& c : constraints) {
        if (c.lhs.empty() && c.rhs) throw std::runtime_error("state has zero norm");
    }

    std::vector<int> sample(n, 0);
    for (int q = 0; q < n; ++q) {
        if (visible[q].vars.empty()) sample[q] = visible[q].c ? 1 : 0;
    }

    for (const StateComponent& component : state_components()) {
        if (component.qubits.empty()) continue;

        LAPH local = component_view(component);
        std::vector<int> local_sample =
            local.exact_sample(rng, QueryOptions{PartitionBackend::Monolithic});

        for (int local_q = 0; local_q < local.n; ++local_q) {
            sample[component.qubits[local_q]] = local_sample[local_q];
        }
    }

    return sample;
}

void LAPH::fold_phase_terms() {
    std::vector<Monomial> kill;
    for (auto& kv : phase) {
        kv.second = mod8(kv.second);
        if (kv.second == 0) kill.push_back(kv.first);
    }
    for (const Monomial& mon : kill) phase.erase(mon);
}

void LAPH::canonicalize_constraints() {
    std::map<int, Constraint> basis;

    for (Constraint row : constraints) {
        row.lhs.normalize_xor();
        while (!row.lhs.empty()) {
            int p = row.lhs.first();
            auto it = basis.find(p);
            if (it == basis.end()) {
                basis[p] = std::move(row);
                goto inserted;
            }
            row.lhs.xor_with(it->second.lhs);
            row.rhs ^= it->second.rhs;
        }
        if (row.rhs) throw std::runtime_error("inconsistent LAPH constraints");
        inserted:;
    }

    constraints.clear();
    constraints.reserve(basis.size());
    for (auto& kv : basis) constraints.push_back(std::move(kv.second));
}

void LAPH::compress() {
    canonicalize_constraints();
    fold_phase_terms();
}

std::vector<std::vector<int>> LAPH::connected_components() const {
    DSU dsu(std::max(1, m));
    auto connect_all = [&](const std::vector<int>& xs) {
        if (xs.empty()) return;
        int root = xs[0];
        for (size_t i = 1; i < xs.size(); ++i) dsu.unite(root, xs[i]);
    };

    for (const auto& c : constraints) connect_all(c.lhs.v);
    for (const auto& kv : phase) connect_all(kv.first.v);
    for (const auto& f : visible) connect_all(f.vars.v);

    std::map<int, std::vector<int>> grouped;
    for (int x = 0; x < m; ++x) grouped[dsu.find(x)].push_back(x);

    std::vector<std::vector<int>> components;
    components.reserve(grouped.size());
    for (auto& kv : grouped) components.push_back(std::move(kv.second));
    return components;
}

Stats LAPH::stats() const {
    std::vector<int> cut = cut_set();
    std::vector<std::vector<int>> components = connected_components();
    bool disconnected = components.size() > 1;

    std::vector<StateComponent> state_comps;
    ComponentStatsSummary component_summary;
    if (disconnected) {
        state_comps = state_components();
        component_summary = component_fabqnf_summary(*this, state_comps);
    }

    FABQNFStats fs;
    int kappa = -1;
    if (disconnected) {
        fs.rho = -1;
        fs.out_dim = -1;
        fs.hid_dim = -1;
        fs.table_size = -1;
    } else if (is_clifford_poly(phase)) {
        FiberChart ch = build_fiber_chart(*this);
        fs.out_dim = ch.out_dim;
        fs.hid_dim = ch.hid_dim;
        fs.rho = 0;
        fs.table_size = 1;
        kappa = 0;
    } else if (m <= 5000) {
        fs = fabqnf_stats(*this);
        kappa = hidden_interference_rank(*this, cut);
    } else {
        fs.rho = -1;
        fs.out_dim = -1;
        fs.hid_dim = -1;
        fs.table_size = -1;
    }

    return {
        n,
        m,
        constraints.size(),
        phase.size(),
        scale,
        cut.size(),
        kappa,
        fs.rho,
        fs.out_dim,
        fs.hid_dim,
        fs.table_size,
        static_cast<size_t>(fs.fiber_active_terms),
        static_cast<size_t>(fs.fiber_active_vars),
        component_summary.max_rho,
        component_summary.sum_table_size,
        component_summary.max_kappa,
        component_summary.table_failures,
        components.size()
    };
}

void LAPH::print_stats() const {
    Stats s = stats();
    std::cout << "qubits=" << n
              << " latent_vars=" << s.latent_variables
              << " constraints=" << s.constraints
              << " phase_terms=" << s.phase_terms
              << " hadamard_scale=" << s.hadamard_scale
              << " non_clifford_cut=" << s.non_clifford_cut
              << " hidden_interference_rank=" << s.hidden_interference_rank
              << " fabqnf_rho=" << s.fabqnf_rho
              << " fabqnf_out_dim=" << s.fabqnf_out_dim
              << " fabqnf_hid_dim=" << s.fabqnf_hid_dim
              << " fabqnf_table_size=" << s.fabqnf_table_size
              << " fiber_active_terms=" << s.fiber_active_terms
              << " max_component_rho=" << s.max_component_rho
              << " sum_component_table_size=" << s.sum_component_table_size
              << " max_component_kappa=" << s.max_component_kappa
              << " component_table_failures=" << s.component_table_failures
              << " components=" << s.components
              << "\n";
}

} // namespace laph
