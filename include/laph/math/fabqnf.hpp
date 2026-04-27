#pragma once

#include "laph/math/fiber_chart.hpp"
#include "laph/solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace laph {

static constexpr int FABQNF_MAX_EXACT_RHO = 20;

struct FABQNFStats {
    int out_dim = 0;
    int hid_dim = 0;
    int rho = 0;
    int table_size = 0;
    int fiber_active_terms = 0;
    int fiber_active_vars = 0;
};

struct FABQNF {
    int out_dim = 0;
    int rho = 0;
    BitRow x0;
    std::vector<BitRow> out_x;
    std::vector<BitRow> quotient_rows;
    std::vector<long double> cell_weight;
    int hidden_partition_evaluations = 0;
    int hidden_partition_cache_hits = 0;
    int hidden_partition_cache_entries = 0;

    BitRow output_from_t(const BitRow& t) const {
        BitRow x = x0;
        for (int i = 0; i < out_dim; ++i) {
            if (t.test(i)) x.xor_with(out_x[i]);
        }
        return x;
    }

    uint64_t quotient_value(const BitRow& t) const {
        uint64_t y = 0;
        for (int i = 0; i < rho; ++i) {
            if (quotient_rows[i].dot(t)) y |= (1ull << i);
        }
        return y;
    }

    BitRow sample_t_given_quotient(uint64_t y, std::mt19937_64& rng) const {
        std::vector<BitRow> rows;
        std::vector<uint8_t> rhs;
        rows.reserve(rho);
        rhs.reserve(rho);

        for (int i = 0; i < rho; ++i) {
            rows.push_back(quotient_rows[i]);
            rhs.push_back(static_cast<uint8_t>((y >> i) & 1ull));
        }

        AffineSystemSolution sol = solve_affine(rows, rhs, out_dim);
        if (!sol.ok) throw std::runtime_error("empty FABQNF quotient cell");
        return random_affine_point(sol, rng);
    }

    std::vector<int> sample(std::mt19937_64& rng) const {
        long double total = 0.0L;
        for (long double w : cell_weight) total += w;
        if (!(total > 0.0L)) throw std::runtime_error("FABQNF has zero total weight");

        std::uniform_real_distribution<long double> U(0.0L, total);
        long double r = U(rng);
        uint64_t y = 0;

        for (; y < cell_weight.size(); ++y) {
            r -= cell_weight[y];
            if (r <= 0.0L) break;
        }
        if (y >= cell_weight.size()) y = static_cast<uint64_t>(cell_weight.size() - 1);

        BitRow t = sample_t_given_quotient(y, rng);
        BitRow x = output_from_t(t);

        std::vector<int> out(x.n, 0);
        for (int i = 0; i < x.n; ++i) out[i] = x.test(i) ? 1 : 0;
        return out;
    }
};

inline bool is_fiber_active(const Monomial& monomial, const FiberChart& ch) {
    for (int v : monomial.v) {
        if (!ch.var_hid[v].empty()) return true;
    }
    return false;
}

template <class State>
std::vector<int> fabqnf_active_factor_variables(const State& st, const FiberChart& ch) {
    std::vector<uint8_t> seen(ch.m, 0);
    std::vector<int> vars;

    for (const auto& kv : st.phase) {
        const Monomial& monomial = kv.first;
        int coeff = mod8(kv.second);
        if (!coeff || monomial.empty()) continue;
        if (!is_fiber_active(monomial, ch)) continue;

        for (int v : monomial.v) {
            if (!seen[v]) {
                seen[v] = 1;
                vars.push_back(v);
            }
        }
    }

    return vars;
}

inline std::vector<BitRow> quotient_rows_from_active_factors(
    const std::vector<int>& active_vars,
    const FiberChart& ch
) {
    if (active_vars.empty()) return {};

    std::vector<BitRow> hidden_equations;
    hidden_equations.reserve(ch.hid_dim);

    for (int h = 0; h < ch.hid_dim; ++h) {
        BitRow eq(static_cast<int>(active_vars.size()));
        for (size_t i = 0; i < active_vars.size(); ++i) {
            if (ch.var_hid[active_vars[i]].test(h)) eq.set(static_cast<int>(i));
        }
        if (!eq.empty()) hidden_equations.push_back(std::move(eq));
    }

    std::vector<BitRow> hidden_left_nullspace =
        nullspace(hidden_equations, static_cast<int>(active_vars.size()));

    std::vector<BitRow> rows;
    rows.reserve(hidden_left_nullspace.size());

    for (const BitRow& alpha : hidden_left_nullspace) {
        BitRow q(ch.out_dim);
        for (size_t i = 0; i < active_vars.size(); ++i) {
            if (alpha.test(static_cast<int>(i))) q.xor_with(ch.var_out[active_vars[i]]);
        }
        if (!q.empty()) rows.push_back(std::move(q));
    }

    std::vector<uint8_t> rhs(rows.size(), 0);
    RrefResult rr = bit_rref(rows, rhs, ch.out_dim);
    if (!rr.ok) throw std::runtime_error("quotient-row RREF failed");

    std::vector<BitRow> independent;
    independent.reserve(rr.rows.size());
    for (const BitRow& r : rr.rows) independent.push_back(r);
    return independent;
}

template <class State>
std::vector<BitRow> fabqnf_quotient_rows(const State& st, const FiberChart& ch) {
    return quotient_rows_from_active_factors(fabqnf_active_factor_variables(st, ch), ch);
}

struct HiddenPhaseModel {
    int vars = 0;
    std::vector<Constraint> constraints;
    PhasePoly phase;
    std::unordered_map<AffineForm, int, AffineHash> lift_cache;

    explicit HiddenPhaseModel(int hidden_vars) : vars(hidden_vars) {}

    int new_var() { return vars++; }

    int lift_affine(const BitRow& row, bool c) {
        XorSet xs = row.to_xorset();

        if (xs.empty()) return c ? -2 : -1;
        if (!c && xs.size() == 1) return xs.v[0];

        AffineForm key(xs, c);
        auto it = lift_cache.find(key);
        if (it != lift_cache.end()) return it->second;

        int u = new_var();
        XorSet lhs = xs;
        lhs.xor_with(XorSet::singleton(u));
        constraints.push_back({lhs, c});
        lift_cache.emplace(std::move(key), u);
        return u;
    }

    void add_phase_monomial(int coeff, const std::vector<int>& factors) {
        add_mod8(phase, Monomial(factors), coeff);
    }

    ScaledComplex partition_sum() const {
        return exact_partition_sum(vars, phase, constraints);
    }

    std::string canonical_key() const {
        std::vector<std::string> constraint_keys;
        constraint_keys.reserve(constraints.size());
        for (const Constraint& c : constraints) {
            std::string key = c.rhs ? "1:" : "0:";
            for (int v : c.lhs.v) {
                key += std::to_string(v);
                key.push_back(',');
            }
            constraint_keys.push_back(std::move(key));
        }
        std::sort(constraint_keys.begin(), constraint_keys.end());

        std::vector<std::string> phase_keys;
        phase_keys.reserve(phase.size());
        for (const auto& kv : phase) {
            int coeff = mod8(kv.second);
            if (!coeff) continue;

            std::string key = std::to_string(coeff);
            key.push_back(':');
            for (int v : kv.first.v) {
                key += std::to_string(v);
                key.push_back(',');
            }
            phase_keys.push_back(std::move(key));
        }
        std::sort(phase_keys.begin(), phase_keys.end());

        std::string key = "v=" + std::to_string(vars) + ";c=";
        for (const std::string& c : constraint_keys) {
            key += c;
            key.push_back(';');
        }
        key += "p=";
        for (const std::string& p : phase_keys) {
            key += p;
            key.push_back(';');
        }
        return key;
    }
};

template <class State>
HiddenPhaseModel build_hidden_phase_model_for_t(
    const State& st,
    const FiberChart& ch,
    const BitRow& t
) {
    HiddenPhaseModel model(ch.hid_dim);

    for (const auto& kv : st.phase) {
        const Monomial& monomial = kv.first;
        int coeff = mod8(kv.second);
        if (!coeff || monomial.empty()) continue;
        if (!is_fiber_active(monomial, ch)) continue;

        bool killed = false;
        std::vector<int> hidden_factors;
        hidden_factors.reserve(monomial.v.size());

        for (int v : monomial.v) {
            bool c = static_cast<bool>(ch.var_const[v]) ^ ch.var_out[v].dot(t);
            int lifted = model.lift_affine(ch.var_hid[v], c);

            if (lifted == -1) {
                killed = true;
                break;
            }
            if (lifted == -2) continue;
            hidden_factors.push_back(lifted);
        }

        if (killed || hidden_factors.empty()) continue;
        model.add_phase_monomial(coeff, hidden_factors);
    }

    return model;
}

template <class State>
ScaledComplex hidden_partition_for_t(
    const State& st,
    const FiberChart& ch,
    const BitRow& t
) {
    return build_hidden_phase_model_for_t(st, ch, t).partition_sum();
}

template <class State>
FABQNF build_fabqnf(const State& st) {
    FiberChart ch = build_fiber_chart(st);

    FABQNF law;
    law.out_dim = ch.out_dim;
    law.x0 = ch.x0;
    law.out_x = ch.out_x;
    law.quotient_rows = fabqnf_quotient_rows(st, ch);
    law.rho = static_cast<int>(law.quotient_rows.size());

    if (law.rho >= FABQNF_MAX_EXACT_RHO) {
        throw std::runtime_error("FABQNF quotient too large");
    }

    uint64_t table_size = 1ull << law.rho;
    law.cell_weight.assign(table_size, 0.0L);

    int cell_free_dim = law.out_dim - law.rho;
    long double cell_size = std::ldexp(1.0L, cell_free_dim);

    std::unordered_map<std::string, ScaledComplex> hidden_cache;

    for (uint64_t y = 0; y < table_size; ++y) {
        std::vector<BitRow> rows;
        std::vector<uint8_t> rhs;
        rows.reserve(law.rho);
        rhs.reserve(law.rho);

        for (int i = 0; i < law.rho; ++i) {
            rows.push_back(law.quotient_rows[i]);
            rhs.push_back(static_cast<uint8_t>((y >> i) & 1ull));
        }

        AffineSystemSolution sol = solve_affine(rows, rhs, law.out_dim);
        if (!sol.ok) {
            law.cell_weight[y] = 0.0L;
            continue;
        }

        HiddenPhaseModel model = build_hidden_phase_model_for_t(st, ch, sol.particular);
        std::string key = model.canonical_key();
        ScaledComplex raw;

        auto cached = hidden_cache.find(key);
        if (cached == hidden_cache.end()) {
            raw = model.partition_sum();
            hidden_cache.emplace(std::move(key), raw);
            ++law.hidden_partition_evaluations;
        } else {
            raw = cached->second;
            ++law.hidden_partition_cache_hits;
        }

        long double re = raw.real_ld();
        long double im = raw.imag_ld();
        law.cell_weight[y] = (re * re + im * im) * cell_size;
    }

    law.hidden_partition_cache_entries = static_cast<int>(hidden_cache.size());
    return law;
}

template <class State>
FABQNFStats fabqnf_stats(const State& st) {
    FiberChart ch = build_fiber_chart(st);
    std::vector<BitRow> qrows = fabqnf_quotient_rows(st, ch);

    FABQNFStats s;
    s.out_dim = ch.out_dim;
    s.hid_dim = ch.hid_dim;
    s.rho = static_cast<int>(qrows.size());
    s.table_size = s.rho < FABQNF_MAX_EXACT_RHO ? (1 << s.rho) : -1;

    std::vector<int> active_vars = fabqnf_active_factor_variables(st, ch);
    s.fiber_active_vars = static_cast<int>(active_vars.size());

    for (const auto& kv : st.phase) {
        const Monomial& monomial = kv.first;
        int coeff = mod8(kv.second);
        if (coeff && !monomial.empty() && is_fiber_active(monomial, ch)) {
            ++s.fiber_active_terms;
        }
    }

    return s;
}

} // namespace laph
