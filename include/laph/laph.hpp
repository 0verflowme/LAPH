#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <set>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace laph {

using std::int64_t;

static inline int mod8(int x) {
    x %= 8;
    if (x < 0) x += 8;
    return x;
}

struct XorSet {
    std::vector<int> v;

    XorSet() = default;
    explicit XorSet(std::vector<int> xs) : v(std::move(xs)) { normalize_xor(); }

    static XorSet singleton(int x) { return XorSet({x}); }

    void normalize_xor() {
        std::sort(v.begin(), v.end());
        std::vector<int> out;
        for (size_t i = 0; i < v.size();) {
            size_t j = i + 1;
            while (j < v.size() && v[j] == v[i]) ++j;
            if (((j - i) & 1u) != 0) out.push_back(v[i]);
            i = j;
        }
        v.swap(out);
    }

    bool empty() const { return v.empty(); }
    size_t size() const { return v.size(); }
    int first() const { return v.front(); }

    bool contains(int x) const {
        return std::binary_search(v.begin(), v.end(), x);
    }

    void xor_with(const XorSet& b) {
        std::vector<int> out;
        out.reserve(v.size() + b.v.size());
        size_t i = 0, j = 0;
        while (i < v.size() || j < b.v.size()) {
            if (j == b.v.size() || (i < v.size() && v[i] < b.v[j])) {
                out.push_back(v[i++]);
            } else if (i == v.size() || b.v[j] < v[i]) {
                out.push_back(b.v[j++]);
            } else {
                ++i;
                ++j;
            }
        }
        v.swap(out);
    }

    XorSet xored(const XorSet& b) const {
        XorSet r = *this;
        r.xor_with(b);
        return r;
    }

    XorSet shifted(int delta) const {
        XorSet r;
        r.v.reserve(v.size());
        for (int x : v) r.v.push_back(x + delta);
        return r;
    }

    XorSet without(int x) const {
        XorSet r;
        r.v.reserve(v.size());
        for (int y : v) if (y != x) r.v.push_back(y);
        return r;
    }

    bool operator==(const XorSet& o) const { return v == o.v; }
};

struct XorHash {
    size_t operator()(const XorSet& s) const noexcept {
        uint64_t h = 1469598103934665603ull;
        for (int x : s.v) {
            uint64_t y = static_cast<uint64_t>(x) + 0x9e3779b97f4a7c15ull;
            h ^= y;
            h *= 1099511628211ull;
        }
        return static_cast<size_t>(h);
    }
};

struct AffineForm {
    XorSet vars;
    bool c = false;

    AffineForm() = default;
    AffineForm(XorSet xs, bool cc=false) : vars(std::move(xs)), c(cc) {}
    static AffineForm constant(bool b) { return AffineForm({}, b); }
    static AffineForm variable(int x) { return AffineForm(XorSet::singleton(x), false); }

    void xor_with(const AffineForm& b) {
        vars.xor_with(b.vars);
        c ^= b.c;
    }

    AffineForm xored(const AffineForm& b) const {
        AffineForm r = *this;
        r.xor_with(b);
        return r;
    }

    AffineForm shifted(int delta) const {
        return AffineForm(vars.shifted(delta), c);
    }

    bool operator==(const AffineForm& o) const { return c == o.c && vars == o.vars; }
};

struct AffineHash {
    size_t operator()(const AffineForm& a) const noexcept {
        return XorHash{}(a.vars) ^ (a.c ? 0x517cc1b727220a95ull : 0ull);
    }
};

struct Monomial {
    std::vector<int> v;

    Monomial() = default;
    explicit Monomial(std::vector<int> xs) : v(std::move(xs)) { normalize_or(); }

    static Monomial constant() { return Monomial(); }
    static Monomial singleton(int x) { return Monomial({x}); }

    void normalize_or() {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    size_t degree() const { return v.size(); }
    bool empty() const { return v.empty(); }

    bool contains(int x) const {
        return std::binary_search(v.begin(), v.end(), x);
    }

    Monomial without(int x) const {
        Monomial r;
        r.v.reserve(v.size());
        for (int y : v) if (y != x) r.v.push_back(y);
        return r;
    }

    Monomial shifted(int delta) const {
        Monomial r;
        r.v.reserve(v.size());
        for (int x : v) r.v.push_back(x + delta);
        return r;
    }

    void or_with_var(int x) {
        if (!contains(x)) {
            v.insert(std::lower_bound(v.begin(), v.end(), x), x);
        }
    }

    static Monomial boolean_mul(const Monomial& a, const Monomial& b) {
        Monomial r;
        r.v.reserve(a.v.size() + b.v.size());
        std::set_union(a.v.begin(), a.v.end(), b.v.begin(), b.v.end(), std::back_inserter(r.v));
        return r;
    }

    bool operator==(const Monomial& o) const { return v == o.v; }
};

struct MonomialHash {
    size_t operator()(const Monomial& m) const noexcept {
        uint64_t h = 1469598103934665603ull;
        for (int x : m.v) {
            uint64_t y = static_cast<uint64_t>(x) + 0x9e3779b97f4a7c15ull;
            h ^= y;
            h *= 1099511628211ull;
        }
        return static_cast<size_t>(h);
    }
};

using PhasePoly = std::unordered_map<Monomial, int, MonomialHash>;

struct Constraint {
    XorSet lhs;
    bool rhs = false;
};

static void add_mod8(PhasePoly& p, const Monomial& m, int coeff) {
    coeff = mod8(coeff);
    if (coeff == 0) return;
    auto it = p.find(m);
    if (it == p.end()) {
        p.emplace(m, coeff);
    } else {
        int nc = mod8(it->second + coeff);
        if (nc == 0) p.erase(it);
        else it->second = nc;
    }
}

static Monomial mon_from_vars(const std::vector<int>& xs) {
    return Monomial(xs);
}

static void add_affine_phase(PhasePoly& poly, int coeff, const XorSet& mask, bool constant) {
    coeff = mod8(coeff);
    if (coeff == 0) return;

    if (constant) {
        add_mod8(poly, Monomial::constant(), coeff);
        coeff = mod8(-coeff);
    }

    const auto& vs = mask.v;
    const int n = static_cast<int>(vs.size());

    int c1 = mod8(coeff);
    if (c1) {
        for (int i = 0; i < n; ++i) {
            add_mod8(poly, Monomial::singleton(vs[i]), c1);
        }
    }

    int c2 = mod8(coeff * -2);
    if (c2) {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                add_mod8(poly, mon_from_vars({vs[i], vs[j]}), c2);
            }
        }
    }

    int c3 = mod8(coeff * 4);
    if (c3) {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = j + 1; k < n; ++k) {
                    add_mod8(poly, mon_from_vars({vs[i], vs[j], vs[k]}), c3);
                }
            }
        }
    }
}

using ANF = std::unordered_set<Monomial, MonomialHash>;

static ANF anf_affine(const AffineForm& a) {
    ANF out;
    if (a.c) out.insert(Monomial::constant());
    for (int x : a.vars.v) out.insert(Monomial::singleton(x));
    return out;
}

static void anf_toggle(ANF& s, const Monomial& m) {
    auto it = s.find(m);
    if (it == s.end()) s.insert(m);
    else s.erase(it);
}

static ANF anf_mul(const ANF& a, const ANF& b) {
    ANF out;
    for (const auto& x : a) {
        for (const auto& y : b) {
            anf_toggle(out, Monomial::boolean_mul(x, y));
        }
    }
    return out;
}

struct ScaledComplex {
    // value = (re + i*im) * (sqrt(2)^k)
    long double re = 0.0L;
    long double im = 0.0L;
    int64_t k = 0;

    static constexpr long double SQRT2 = 1.41421356237309504880168872420969807857L;
    static constexpr long double INV_SQRT2 = 0.70710678118654752440084436210484903928L;

    static ScaledComplex zero() { return {}; }

    static ScaledComplex omega_phase(int phase, int64_t sqrt2_power = 0) {
        phase = mod8(phase);
        switch (phase) {
            case 0: return {1.0L, 0.0L, sqrt2_power};
            case 1: return {INV_SQRT2, INV_SQRT2, sqrt2_power};
            case 2: return {0.0L, 1.0L, sqrt2_power};
            case 3: return {-INV_SQRT2, INV_SQRT2, sqrt2_power};
            case 4: return {-1.0L, 0.0L, sqrt2_power};
            case 5: return {-INV_SQRT2, -INV_SQRT2, sqrt2_power};
            case 6: return {0.0L, -1.0L, sqrt2_power};
            case 7: return {INV_SQRT2, -INV_SQRT2, sqrt2_power};
        }
        return zero();
    }

    bool is_zero() const { return re == 0.0L && im == 0.0L; }

    static long double pow_sqrt2(int64_t d) {
        if (d > 4096) return std::numeric_limits<long double>::infinity();
        if (d < -4096) return 0.0L;
        return std::powl(SQRT2, static_cast<long double>(d));
    }

    void normalize() {
        if (is_zero()) { k = 0; return; }
        long double mag2 = re * re + im * im;
        if (mag2 == 0.0L) { re = im = 0.0L; k = 0; return; }
        long double log2mag = 0.5L * std::log2(mag2);
        int64_t shift = static_cast<int64_t>(std::floor(2.0L * log2mag));
        if (shift != 0) {
            long double f = pow_sqrt2(-shift);
            re *= f;
            im *= f;
            k += shift;
        }
    }

    void add(const ScaledComplex& b) {
        if (b.is_zero()) return;
        if (is_zero()) { *this = b; return; }
        int64_t K = std::max(k, b.k);
        long double fa = pow_sqrt2(k - K);
        long double fb = pow_sqrt2(b.k - K);
        re = re * fa + b.re * fb;
        im = im * fa + b.im * fb;
        k = K;
        normalize();
    }

    void mul(const ScaledComplex& b) {
        if (is_zero() || b.is_zero()) {
            *this = zero();
            return;
        }
        long double nr = re * b.re - im * b.im;
        long double ni = re * b.im + im * b.re;
        re = nr;
        im = ni;
        k += b.k;
        normalize();
    }

    void mul_sqrt2_power(int64_t d) {
        if (!is_zero()) k += d;
    }

    long double real_ld() const {
        return re * pow_sqrt2(k);
    }

    long double imag_ld() const {
        return im * pow_sqrt2(k);
    }
};

static bool is_clifford_poly(const PhasePoly& poly) {
    for (const auto& kv : poly) {
        const Monomial& m = kv.first;
        int c = mod8(kv.second);
        if (c == 0) continue;
        size_t deg = m.degree();
        if (deg == 0) continue;
        if (deg == 1) {
            if (c % 2 != 0) return false;
        } else if (deg == 2) {
            if (c % 4 != 0) return false;
        } else {
            return false;
        }
    }
    return true;
}

static PhasePoly substitute_pivot_in_clifford(const PhasePoly& poly, int pivot, bool constant, const XorSet& mask) {
    PhasePoly out;
    for (const auto& kv : poly) {
        const Monomial& m = kv.first;
        int c = mod8(kv.second);
        if (!m.contains(pivot)) {
            add_mod8(out, m, c);
            continue;
        }
        Monomial rest = m.without(pivot);
        if (rest.empty()) {
            add_affine_phase(out, c, mask, constant);
        } else {
            if (rest.degree() != 1 || c % 4 != 0) {
                throw std::runtime_error("bad Clifford pivot substitution");
            }
            int j = rest.v[0];
            AffineForm a(mask, constant);
            ANF prod = anf_mul(anf_affine(a), ANF{Monomial::singleton(j)});
            for (const auto& mon : prod) add_mod8(out, mon, c);
        }
    }
    return out;
}

static ScaledComplex clifford_sum(PhasePoly poly, int nvars) {
    if (!is_clifford_poly(poly)) throw std::runtime_error("not Clifford");

    int phase = 0;
    auto it0 = poly.find(Monomial::constant());
    if (it0 != poly.end()) {
        phase = mod8(it0->second);
        poly.erase(it0);
    }

    std::set<int> active;
    for (int i = 0; i < nvars; ++i) active.insert(i);

    int64_t log_sqrt2 = 0;

    while (!active.empty()) {
        int chosen = -1;

        for (int v : active) {
            auto it = poly.find(Monomial::singleton(v));
            int c = (it == poly.end()) ? 0 : mod8(it->second);
            if (c % 2 == 0 && ((c / 2) % 2 == 1)) {
                chosen = v;
                break;
            }
        }

        if (chosen == -1) {
            for (int v : active) {
                bool has_quad = false;
                for (const auto& kv : poly) {
                    if (kv.first.degree() == 2 && kv.first.contains(v) && mod8(kv.second) != 0) {
                        has_quad = true;
                        break;
                    }
                }
                if (has_quad) { chosen = v; break; }
            }
        }

        if (chosen == -1) chosen = *active.begin();

        int v = chosen;
        auto lin_it = poly.find(Monomial::singleton(v));
        int lin = (lin_it == poly.end()) ? 0 : mod8(lin_it->second);
        int l = (lin / 2) % 4;

        XorSet neigh;
        for (const auto& kv : poly) {
            const Monomial& m = kv.first;
            int c = mod8(kv.second);
            if (m.degree() == 2 && m.contains(v) && c != 0) {
                if (c % 4 != 0) throw std::runtime_error("non-Clifford quadratic");
                if (c % 8 != 0) {
                    int other = (m.v[0] == v) ? m.v[1] : m.v[0];
                    neigh.xor_with(XorSet::singleton(other));
                }
            }
        }

        PhasePoly p0;
        for (const auto& kv : poly) {
            if (!kv.first.contains(v)) add_mod8(p0, kv.first, kv.second);
        }

        active.erase(v);

        if (l == 1 || l == 3) {
            log_sqrt2 += 1;
            phase = mod8(phase + (l == 1 ? 1 : -1));
            add_affine_phase(p0, (l == 1 ? 6 : 2), neigh, false);
            poly.swap(p0);
        } else {
            bool required = (l == 2);
            if (neigh.empty()) {
                if (required) return ScaledComplex::zero();
                log_sqrt2 += 2;
                poly.swap(p0);
            } else {
                int pivot = neigh.first();
                XorSet rest = neigh.without(pivot);
                poly = substitute_pivot_in_clifford(p0, pivot, required, rest);
                active.erase(pivot);
                log_sqrt2 += 2;
            }
        }

        std::vector<Monomial> kill;
        for (const auto& kv : poly) if (mod8(kv.second) == 0) kill.push_back(kv.first);
        for (const auto& m : kill) poly.erase(m);
    }

    auto it = poly.find(Monomial::constant());
    if (it != poly.end()) phase = mod8(phase + it->second);
    return ScaledComplex::omega_phase(phase, log_sqrt2);
}

struct SolveResult {
    bool ok = false;
    int nfree = 0;
    std::unordered_map<int, AffineForm> expr; // original var -> affine over compact free vars 0..nfree-1
};

static SolveResult solve_affine_system(const std::vector<Constraint>& rows_in, const std::vector<int>& variables) {
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

static PhasePoly substitute_clifford_poly(const PhasePoly& poly, const std::unordered_map<int, AffineForm>& exprs) {
    PhasePoly out;
    for (const auto& kv : poly) {
        const Monomial& m = kv.first;
        int c = mod8(kv.second);
        if (c == 0) continue;
        size_t deg = m.degree();
        if (deg == 0) {
            add_mod8(out, Monomial::constant(), c);
        } else if (deg == 1) {
            auto it = exprs.find(m.v[0]);
            if (it == exprs.end()) throw std::runtime_error("missing expression");
            add_affine_phase(out, c, it->second.vars, it->second.c);
        } else if (deg == 2 && c % 4 == 0) {
            auto it0 = exprs.find(m.v[0]);
            auto it1 = exprs.find(m.v[1]);
            if (it0 == exprs.end() || it1 == exprs.end()) throw std::runtime_error("missing expression");
            ANF prod = anf_mul(anf_affine(it0->second), anf_affine(it1->second));
            for (const auto& mon : prod) add_mod8(out, mon, c);
        } else {
            throw std::runtime_error("non-Clifford residual term after cut");
        }
    }
    return out;
}

static PhasePoly restrict_poly(const PhasePoly& poly, const std::unordered_map<int, bool>& assign) {
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

static std::vector<Constraint> restrict_constraints(const std::vector<Constraint>& cs, const std::unordered_map<int, bool>& assign) {
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

static std::vector<int> greedy_cut_set_for_phase(const PhasePoly& phase) {
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

static ScaledComplex exact_partition_sum(
    int total_vars,
    const PhasePoly& phase,
    const std::vector<Constraint>& constraints,
    const std::vector<int>* provided_cut = nullptr
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

struct DSU {
    std::vector<int> p, r;
    explicit DSU(int n=0) : p(n), r(n,0) { std::iota(p.begin(), p.end(), 0); }
    int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return;
        if (r[a] < r[b]) std::swap(a, b);
        p[b] = a;
        if (r[a] == r[b]) ++r[a];
    }
};

static std::vector<std::vector<int>> partition_components_for_problem(
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

static Monomial remap_monomial(const Monomial& mon, const std::vector<int>& local_index) {
    Monomial out;
    out.v.reserve(mon.v.size());
    for (int x : mon.v) out.v.push_back(local_index[x]);
    out.normalize_or();
    return out;
}

static XorSet remap_xorset(const XorSet& xs, const std::vector<int>& local_index) {
    XorSet out;
    out.v.reserve(xs.v.size());
    for (int x : xs.v) out.v.push_back(local_index[x]);
    out.normalize_xor();
    return out;
}

static AffineForm remap_affine_form(const AffineForm& f, const std::vector<int>& local_index) {
    return AffineForm(remap_xorset(f.vars, local_index), f.c);
}

static ScaledComplex exact_partition_sum_factorized(
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

enum class PartitionBackend {
    Monolithic,
    Factorized
};

struct QueryOptions {
    PartitionBackend backend = PartitionBackend::Factorized;
};

static ScaledComplex exact_partition_sum_with_backend(
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

struct Stats {
    int qubits = 0;
    int latent_variables = 0;
    size_t constraints = 0;
    size_t phase_terms = 0;
    int hadamard_scale = 0;
    size_t non_clifford_cut = 0;
    size_t components = 0;
};

struct StateComponent {
    std::vector<int> variables;
    std::vector<int> qubits;
};

class LAPH {
public:
    int n = 0;
    int m = 0;
    std::vector<AffineForm> visible;
    std::vector<Constraint> constraints;
    PhasePoly phase;
    int scale = 0;
    std::unordered_map<AffineForm, int, AffineHash> lift_cache;

    explicit LAPH(int nqubits) : n(nqubits), visible(nqubits, AffineForm::constant(false)) {}

    int new_var() { return m++; }

    std::pair<bool,int> lift(const AffineForm& f) {
        // returns {is_var, value}; if !is_var, value is constant 0/1
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

    void add_phase_monomial(int coeff, const Monomial& mon) {
        add_mod8(phase, mon, coeff);
    }

    void add_phase_product_forms(int coeff, const std::vector<AffineForm>& forms) {
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

    LAPH& x(int q) {
        visible[q].c ^= true;
        return *this;
    }

    LAPH& cnot(int c, int t) {
        visible[t].xor_with(visible[c]);
        return *this;
    }

    LAPH& h(int q) {
        AffineForm old = visible[q];
        int s = new_var();
        add_phase_product_forms(4, {AffineForm::variable(s), old});
        visible[q] = AffineForm::variable(s);
        scale += 1;
        return *this;
    }

    LAPH& t(int q) { add_phase_product_forms(1, {visible[q]}); return *this; }
    LAPH& s(int q) { add_phase_product_forms(2, {visible[q]}); return *this; }
    LAPH& z(int q) { add_phase_product_forms(4, {visible[q]}); return *this; }

    LAPH& cz(int a, int b) {
        add_phase_product_forms(4, {visible[a], visible[b]});
        return *this;
    }

    LAPH& ccz(int a, int b, int c) {
        add_phase_product_forms(4, {visible[a], visible[b], visible[c]});
        return *this;
    }

    std::vector<Constraint> output_constraints(uint64_t xbits) const {
        std::vector<Constraint> out = constraints;
        out.reserve(out.size() + n);
        for (int q = 0; q < n; ++q) {
            bool bit = ((xbits >> q) & 1ull) != 0;
            out.push_back({visible[q].vars, static_cast<bool>(bit ^ visible[q].c)});
        }
        return out;
    }

    std::vector<int> cut_set() const { return greedy_cut_set_for_phase(phase); }

    ScaledComplex amplitude(
        uint64_t xbits,
        QueryOptions options = {}
    ) const {
        std::vector<Constraint> cs = output_constraints(xbits);
        ScaledComplex z = exact_partition_sum_with_backend(m, phase, cs, options.backend);
        z.mul_sqrt2_power(-scale);
        return z;
    }

    ScaledComplex amplitude_monolithic(uint64_t xbits) const {
        return amplitude(xbits, QueryOptions{PartitionBackend::Monolithic});
    }

    ScaledComplex amplitude_factorized(uint64_t xbits) const {
        return amplitude(xbits, QueryOptions{PartitionBackend::Factorized});
    }

    PhasePoly density_phase() const {
        PhasePoly ph;
        for (const auto& kv : phase) {
            add_mod8(ph, kv.first, kv.second);
            add_mod8(ph, kv.first.shifted(m), -kv.second);
        }
        return ph;
    }

    std::vector<Constraint> density_constraints(const std::vector<std::pair<int,bool>>& prefix) const {
        std::vector<Constraint> cs;
        cs.reserve(2 * constraints.size() + n + prefix.size());

        for (const auto& c : constraints) {
            cs.push_back(c);
            cs.push_back({c.lhs.shifted(m), c.rhs});
        }

        // visible equality: A_q z = A_q z'
        for (int q = 0; q < n; ++q) {
            XorSet lhs = visible[q].vars;
            lhs.xor_with(visible[q].vars.shifted(m));
            cs.push_back({lhs, false});
        }

        // prefix constraints on unprimed visible bits.
        for (auto [q, bit] : prefix) {
            cs.push_back({visible[q].vars, static_cast<bool>(bit ^ visible[q].c)});
        }

        return cs;
    }

    ScaledComplex density_prefix_sum(
        const std::vector<std::pair<int,bool>>& prefix,
        QueryOptions options = {}
    ) const {
        PhasePoly ph = density_phase();
        std::vector<Constraint> cs = density_constraints(prefix);
        return exact_partition_sum_with_backend(2 * m, ph, cs, options.backend);
    }

    long double density_prefix_weight(
        const std::vector<std::pair<int,bool>>& prefix,
        QueryOptions options = {}
    ) const {
        ScaledComplex z = density_prefix_sum(prefix, options);
        long double p = z.real_ld();
        if (p < 0 && p > -1e-12L) p = 0;
        return p;
    }

    long double probability_prefix(
        const std::vector<std::pair<int,bool>>& prefix,
        QueryOptions options = {}
    ) const {
        ScaledComplex z = density_prefix_sum(prefix, options);
        z.mul_sqrt2_power(-2LL * scale); // probability factor = 2^{-scale}
        long double p = z.real_ld();
        if (p < 0 && p > -1e-12L) p = 0;
        if (p > 1 && p < 1 + 1e-12L) p = 1;
        return p;
    }

    long double probability_prefix_monolithic(
        const std::vector<std::pair<int,bool>>& prefix
    ) const {
        return probability_prefix(prefix, QueryOptions{PartitionBackend::Monolithic});
    }

    long double probability_prefix_factorized(
        const std::vector<std::pair<int,bool>>& prefix
    ) const {
        return probability_prefix(prefix, QueryOptions{PartitionBackend::Factorized});
    }

    std::vector<int> exact_sample(
        std::mt19937_64& rng,
        QueryOptions options = {}
    ) const {
        std::vector<std::pair<int,bool>> prefix;
        prefix.reserve(n);
        std::vector<int> sample(n, 0);
        long double prefix_prob = probability_prefix(prefix, options);
        if (prefix_prob <= 0) throw std::runtime_error("state has zero norm");

        std::uniform_real_distribution<long double> U(0.0L, 1.0L);
        for (int q = 0; q < n; ++q) {
            auto pref0 = prefix;
            pref0.push_back({q, false});
            long double p0_abs = probability_prefix(pref0, options);
            long double p0_cond = p0_abs / prefix_prob;
            if (p0_cond < 0) p0_cond = 0;
            if (p0_cond > 1) p0_cond = 1;
            bool bit = !(U(rng) < p0_cond);
            prefix.push_back({q, bit});
            sample[q] = bit ? 1 : 0;
            prefix_prob = bit ? (prefix_prob - p0_abs) : p0_abs;
            if (prefix_prob < 0 && prefix_prob > -1e-12L) prefix_prob = 0;
        }
        return sample;
    }

    std::vector<int> exact_sample_monolithic(std::mt19937_64& rng) const {
        return exact_sample(rng, QueryOptions{PartitionBackend::Monolithic});
    }

    std::vector<int> exact_sample_factorized(std::mt19937_64& rng) const {
        return exact_sample(rng, QueryOptions{PartitionBackend::Factorized});
    }

    std::vector<StateComponent> state_components() const {
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

    LAPH component_view(const StateComponent& component) const {
        LAPH local(static_cast<int>(component.qubits.size()));
        local.m = static_cast<int>(component.variables.size());

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

    std::vector<int> exact_sample_by_component(
        std::mt19937_64& rng,
        QueryOptions options = {}
    ) const {
        for (const Constraint& c : constraints) {
            if (c.lhs.empty() && c.rhs) throw std::runtime_error("state has zero norm");
        }

        std::vector<int> sample(n, 0);
        for (int q = 0; q < n; ++q) {
            if (visible[q].vars.empty()) sample[q] = visible[q].c ? 1 : 0;
        }

        std::uniform_real_distribution<long double> U(0.0L, 1.0L);

        for (const StateComponent& component : state_components()) {
            if (component.qubits.empty()) continue;

            LAPH local = component_view(component);
            std::vector<std::pair<int, bool>> prefix;
            prefix.reserve(local.n);

            long double prefix_weight = local.density_prefix_weight(prefix, options);
            if (prefix_weight <= 0) throw std::runtime_error("component has zero norm");

            for (int local_q = 0; local_q < local.n; ++local_q) {
                auto pref0 = prefix;
                pref0.push_back({local_q, false});

                long double p0_weight = local.density_prefix_weight(pref0, options);
                long double p0_cond = p0_weight / prefix_weight;
                if (p0_cond < 0) p0_cond = 0;
                if (p0_cond > 1) p0_cond = 1;

                bool bit = !(U(rng) < p0_cond);
                prefix.push_back({local_q, bit});
                sample[component.qubits[local_q]] = bit ? 1 : 0;
                prefix_weight = bit ? (prefix_weight - p0_weight) : p0_weight;
                if (prefix_weight < 0 && prefix_weight > -1e-12L) prefix_weight = 0;
            }
        }

        return sample;
    }

    void fold_phase_terms() {
        std::vector<Monomial> kill;
        for (auto& kv : phase) {
            kv.second = mod8(kv.second);
            if (kv.second == 0) kill.push_back(kv.first);
        }
        for (const Monomial& m : kill) phase.erase(m);
    }

    void canonicalize_constraints() {
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

    void compress() {
        canonicalize_constraints();
        fold_phase_terms();
    }

    std::vector<std::vector<int>> connected_components() const {
        DSU dsu(std::max(1, m));
        auto connect_all = [&](const std::vector<int>& xs) {
            if (xs.empty()) return;
            int root = xs[0];
            for (size_t i = 1; i < xs.size(); ++i) dsu.unite(root, xs[i]);
        };
        for (const auto& c : constraints) connect_all(c.lhs.v);
        for (const auto& kv : phase) connect_all(kv.first.v);
        for (const auto& f : visible) connect_all(f.vars.v);

        std::map<int, std::vector<int>> mp;
        for (int x = 0; x < m; ++x) mp[dsu.find(x)].push_back(x);
        std::vector<std::vector<int>> comps;
        for (auto& kv : mp) comps.push_back(std::move(kv.second));
        return comps;
    }

    Stats stats() const {
        return {
            n,
            m,
            constraints.size(),
            phase.size(),
            scale,
            cut_set().size(),
            connected_components().size()
        };
    }

    void print_stats() const {
        Stats s = stats();
        std::cout << "qubits=" << n
                  << " latent_vars=" << s.latent_variables
                  << " constraints=" << s.constraints
                  << " phase_terms=" << s.phase_terms
                  << " hadamard_scale=" << s.hadamard_scale
                  << " non_clifford_cut=" << s.non_clifford_cut
                  << " components=" << s.components
                  << "\n";
    }
};

} // namespace laph
