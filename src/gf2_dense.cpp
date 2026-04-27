#include "laph/gf2_dense.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

namespace laph {

DenseRow::DenseRow(int n) : nbits(n), w((n + 63) >> 6, 0) {}

DenseRow DenseRow::from_xorset(const XorSet& xs, int nbits) {
    DenseRow r(nbits);
    for (int x : xs.v) {
        if (x < 0 || x >= nbits) {
            throw std::runtime_error("DenseRow variable index out of range");
        }
        r.w[static_cast<size_t>(x) >> 6] ^= (1ull << (x & 63));
    }
    return r;
}

bool DenseRow::test(int i) const {
    return (w[static_cast<size_t>(i) >> 6] >> (i & 63)) & 1ull;
}

void DenseRow::xor_with(const DenseRow& b) {
    assert(nbits == b.nbits);
    for (size_t i = 0; i < w.size(); ++i) w[i] ^= b.w[i];
}

bool DenseRow::empty() const {
    for (uint64_t x : w) {
        if (x) return false;
    }
    return true;
}

int DenseRow::first_set() const {
    for (size_t i = 0; i < w.size(); ++i) {
        uint64_t x = w[i];
        if (!x) continue;
        int b = __builtin_ctzll(x);
        int idx = static_cast<int>(i * 64 + b);
        if (idx < nbits) return idx;
    }
    return -1;
}

XorSet DenseRow::to_xorset() const {
    XorSet xs;
    for (int i = 0; i < nbits; ++i) {
        if (test(i)) xs.v.push_back(i);
    }
    return xs;
}

DenseRrefState::DenseRrefState(int n)
    : nvars(n), row(n), rhs(n, 0), has(n, 0) {
    for (int i = 0; i < n; ++i) row[i] = DenseRow(n);
}

bool DenseRrefState::insert(DenseRow r, bool b) {
    if (r.nbits != nvars) throw std::runtime_error("RREF row width mismatch");

    for (int p : pivots) {
        if (r.test(p)) {
            r.xor_with(row[p]);
            b ^= static_cast<bool>(rhs[p]);
        }
    }

    if (r.empty()) return !b;

    int p = r.first_set();
    if (p < 0) return !b;

    for (int q : pivots) {
        if (row[q].test(p)) {
            row[q].xor_with(r);
            rhs[q] ^= static_cast<unsigned char>(b);
        }
    }

    has[p] = 1;
    row[p] = std::move(r);
    rhs[p] = static_cast<unsigned char>(b);
    pivots.insert(std::lower_bound(pivots.begin(), pivots.end(), p), p);
    return true;
}

bool DenseRrefState::insert_constraint(const Constraint& c) {
    return insert(DenseRow::from_xorset(c.lhs, nvars), c.rhs);
}

DenseSolveResult DenseRrefState::solve() const {
    DenseSolveResult res;
    res.ok = true;
    res.expr.assign(nvars, AffineForm::constant(false));

    std::vector<int> free_index(nvars, -1);
    for (int c = 0; c < nvars; ++c) {
        if (!has[c]) {
            free_index[c] = res.nfree++;
            res.expr[c] = AffineForm::variable(free_index[c]);
        }
    }

    for (int p : pivots) {
        AffineForm e = AffineForm::constant(static_cast<bool>(rhs[p]));
        for (int c = 0; c < nvars; ++c) {
            if (!has[c] && row[p].test(c)) {
                e.vars.v.push_back(free_index[c]);
            }
        }
        e.vars.normalize_xor();
        res.expr[p] = std::move(e);
    }

    return res;
}

DenseSolveResult solve_affine_system_dense(
    const std::vector<Constraint>& rows,
    int nvars
) {
    DenseRrefState st(nvars);
    for (const Constraint& c : rows) {
        if (!st.insert_constraint(c)) return {};
    }
    return st.solve();
}

AffineForm compose_affine(
    const AffineForm& outer,
    const std::vector<AffineForm>& inner
) {
    AffineForm out = AffineForm::constant(outer.c);
    for (int x : outer.vars.v) {
        if (x < 0 || x >= static_cast<int>(inner.size())) {
            throw std::runtime_error("compose variable out of range");
        }
        out.xor_with(inner[x]);
    }
    return out;
}

AffineForm project_xorset(
    const XorSet& xs,
    const std::vector<AffineForm>& expr
) {
    AffineForm out = AffineForm::constant(false);
    for (int x : xs.v) {
        if (x < 0 || x >= static_cast<int>(expr.size())) {
            throw std::runtime_error("project variable out of range");
        }
        out.xor_with(expr[x]);
    }
    return out;
}

Constraint affine_equals_bit(const AffineForm& f, bool bit) {
    return Constraint{f.vars, static_cast<bool>(bit ^ f.c)};
}

} // namespace laph
