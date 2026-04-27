#include "laph/solver.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace laph {

SolveResult solve_affine_system(const std::vector<Constraint>& rows_in, const std::vector<int>& variables) {
    std::map<int, Constraint> basis;

    for (Constraint row : rows_in) {
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
        if (row.rhs) return {};
        inserted:;
    }

    std::unordered_set<int> varset(variables.begin(), variables.end());
    std::unordered_set<int> pivots;
    for (const auto& kv : basis) pivots.insert(kv.first);

    std::vector<int> free_vars;
    for (int x : variables) if (!pivots.count(x)) free_vars.push_back(x);
    std::sort(free_vars.begin(), free_vars.end());

    SolveResult res;
    res.ok = true;
    res.nfree = static_cast<int>(free_vars.size());

    for (int i = 0; i < static_cast<int>(free_vars.size()); ++i) {
        res.expr[free_vars[i]] = AffineForm::variable(i);
    }

    std::vector<int> pivot_list;
    for (const auto& kv : basis) pivot_list.push_back(kv.first);
    std::sort(pivot_list.begin(), pivot_list.end(), std::greater<int>());

    for (int p : pivot_list) {
        const Constraint& row = basis[p];
        AffineForm e = AffineForm::constant(row.rhs);
        for (int x : row.lhs.v) {
            if (x == p) continue;
            auto it = res.expr.find(x);
            if (it == res.expr.end()) {
                // A variable outside the requested variable list appeared. Treat it as a new free var.
                int idx = res.nfree++;
                res.expr[x] = AffineForm::variable(idx);
                e.xor_with(res.expr[x]);
            } else {
                e.xor_with(it->second);
            }
        }
        res.expr[p] = e;
    }

    return res;
}

PhasePoly restrict_poly(const PhasePoly& poly, const std::unordered_map<int, bool>& assign) {
    PhasePoly out;
    for (const auto& kv : poly) {
        bool killed = false;
        Monomial nm;
        for (int x : kv.first.v) {
            auto it = assign.find(x);
            if (it == assign.end()) nm.v.push_back(x);
            else if (!it->second) { killed = true; break; }
            // if assigned one, variable is removed from product
        }
        if (!killed) add_mod8(out, nm, kv.second);
    }
    return out;
}

std::vector<Constraint> restrict_constraints(const std::vector<Constraint>& cs, const std::unordered_map<int, bool>& assign) {
    std::vector<Constraint> out;
    out.reserve(cs.size());
    for (const Constraint& c : cs) {
        Constraint nc;
        nc.rhs = c.rhs;
        for (int x : c.lhs.v) {
            auto it = assign.find(x);
            if (it == assign.end()) nc.lhs.v.push_back(x);
            else nc.rhs ^= it->second;
        }
        nc.lhs.normalize_xor();
        out.push_back(std::move(nc));
    }
    return out;
}

std::vector<int> greedy_cut_set_for_phase(const PhasePoly& phase) {
    struct Demand { Monomial m; int quota; };
    std::vector<Demand> demands;
    for (const auto& kv : phase) {
        const Monomial& m = kv.first;
        int c = mod8(kv.second);
        if (c == 0 || m.degree() == 0) continue;
        int allowed;
        if (c % 2 == 1) allowed = 0;
        else if (c % 4 == 2) allowed = 1;
        else allowed = 2;
        int quota = std::max<int>(0, static_cast<int>(m.degree()) - allowed);
        if (quota > 0) demands.push_back({m, quota});
    }

    std::unordered_set<int> cut;
    while (true) {
        bool done = true;
        std::unordered_map<int, int64_t> score;
        for (const Demand& d : demands) {
            int have = 0;
            for (int x : d.m.v) if (cut.count(x)) ++have;
            if (have >= d.quota) continue;
            done = false;
            int need = d.quota - have;
            int avail = 0;
            for (int x : d.m.v) if (!cut.count(x)) ++avail;
            int64_t weight = (need >= avail) ? 1000 : 1;
            for (int x : d.m.v) if (!cut.count(x)) score[x] += weight;
        }
        if (done) break;
        int best = -1;
        int64_t best_score = -1;
        for (const auto& kv : score) {
            if (kv.second > best_score) { best = kv.first; best_score = kv.second; }
        }
        if (best < 0) throw std::runtime_error("cut-set failed");
        cut.insert(best);
    }
    std::vector<int> out(cut.begin(), cut.end());
    std::sort(out.begin(), out.end());
    return out;
}

ScaledComplex exact_partition_sum(
    int total_vars,
    const PhasePoly& phase,
    const std::vector<Constraint>& constraints,
    const std::vector<int>* provided_cut
) {
    std::vector<int> cut = provided_cut ? *provided_cut : greedy_cut_set_for_phase(phase);
    if (cut.size() >= static_cast<size_t>(8 * sizeof(uint64_t))) {
        throw std::runtime_error("cut too large for this debug enumerator; use Gray-code / parallel backend");
    }

    ScaledComplex total = ScaledComplex::zero();
    uint64_t lim = 1ull << cut.size();

    for (uint64_t mask = 0; mask < lim; ++mask) {
        std::unordered_map<int, bool> assignment;
        assignment.reserve(cut.size() * 2 + 1);
        for (size_t i = 0; i < cut.size(); ++i) assignment[cut[i]] = ((mask >> i) & 1ull) != 0;

        PhasePoly ph = restrict_poly(phase, assignment);
        std::vector<Constraint> cs = restrict_constraints(constraints, assignment);

        std::vector<int> remaining;
        remaining.reserve(total_vars - static_cast<int>(cut.size()));
        std::unordered_set<int> cutset(cut.begin(), cut.end());
        for (int x = 0; x < total_vars; ++x) if (!cutset.count(x)) remaining.push_back(x);

        SolveResult sol = solve_affine_system(cs, remaining);
        if (!sol.ok) continue;

        PhasePoly residual = substitute_clifford_poly(ph, sol.expr);
        if (!is_clifford_poly(residual)) throw std::runtime_error("cut failed to Cliffordize residual");
        ScaledComplex term = clifford_sum(residual, sol.nfree);
        total.add(term);
    }
    return total;
}

} // namespace laph
