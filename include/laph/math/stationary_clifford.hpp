#pragma once

#include "laph/types.hpp"

#include <algorithm>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

namespace laph {

struct BitRow {
    int n = 0;
    std::vector<uint64_t> w;

    BitRow() = default;
    explicit BitRow(int nbits) : n(nbits), w((nbits + 63) >> 6, 0) {}

    static BitRow from_xorset(const XorSet& xs, int nbits) {
        BitRow r(nbits);
        for (int x : xs.v) {
            if (x < 0 || x >= nbits) {
                throw std::runtime_error("BitRow variable index out of range");
            }
            r.flip(x);
        }
        return r;
    }

    void flip(int i) {
        w[static_cast<size_t>(i) >> 6] ^= 1ull << (i & 63);
    }

    void set(int i) {
        w[static_cast<size_t>(i) >> 6] |= 1ull << (i & 63);
    }

    bool test(int i) const {
        return (w[static_cast<size_t>(i) >> 6] >> (i & 63)) & 1ull;
    }

    void xor_with(const BitRow& b) {
        if (n != b.n) throw std::runtime_error("BitRow width mismatch");
        for (size_t i = 0; i < w.size(); ++i) w[i] ^= b.w[i];
    }

    bool empty() const {
        for (uint64_t x : w) {
            if (x) return false;
        }
        return true;
    }

    int first_set() const {
        for (size_t i = 0; i < w.size(); ++i) {
            uint64_t x = w[i];
            if (!x) continue;

            int bit = __builtin_ctzll(x);
            int idx = static_cast<int>(64 * i + bit);
            if (idx < n) return idx;
        }
        return -1;
    }

    bool dot(const BitRow& b) const {
        if (n != b.n) throw std::runtime_error("BitRow dot width mismatch");
        uint64_t acc = 0;
        for (size_t i = 0; i < w.size(); ++i) acc ^= (w[i] & b.w[i]);
        return __builtin_parityll(acc) != 0;
    }
};

struct AffineSystemSolution {
    bool ok = false;
    BitRow particular;
    std::vector<BitRow> kernel;
};

struct RrefResult {
    bool ok = false;
    std::vector<BitRow> rows;
    std::vector<uint8_t> rhs;
    std::vector<int> pivots;
};

inline RrefResult bit_rref(
    std::vector<BitRow> rows,
    std::vector<uint8_t> rhs,
    int nvars
) {
    if (rows.size() != rhs.size()) {
        throw std::runtime_error("RREF row/rhs size mismatch");
    }

    RrefResult out;
    int r = 0;

    for (int c = 0; c < nvars; ++c) {
        int pivot = -1;
        for (int i = r; i < static_cast<int>(rows.size()); ++i) {
            if (rows[i].test(c)) {
                pivot = i;
                break;
            }
        }
        if (pivot < 0) continue;

        std::swap(rows[r], rows[pivot]);
        std::swap(rhs[r], rhs[pivot]);

        for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
            if (i != r && rows[i].test(c)) {
                rows[i].xor_with(rows[r]);
                rhs[i] ^= rhs[r];
            }
        }

        out.pivots.push_back(c);
        ++r;
    }

    for (int i = r; i < static_cast<int>(rows.size()); ++i) {
        if (rows[i].empty() && rhs[i]) {
            out.ok = false;
            return out;
        }
    }

    out.rows.assign(rows.begin(), rows.begin() + r);
    out.rhs.assign(rhs.begin(), rhs.begin() + r);
    out.ok = true;
    return out;
}

inline std::vector<BitRow> nullspace(
    const std::vector<BitRow>& input_rows,
    int nvars
) {
    std::vector<uint8_t> rhs(input_rows.size(), 0);
    RrefResult rr = bit_rref(input_rows, rhs, nvars);
    if (!rr.ok) {
        throw std::runtime_error("homogeneous system unexpectedly inconsistent");
    }

    std::vector<uint8_t> is_pivot(nvars, 0);
    for (int p : rr.pivots) is_pivot[p] = 1;

    std::vector<BitRow> basis;
    for (int f = 0; f < nvars; ++f) {
        if (is_pivot[f]) continue;

        BitRow v(nvars);
        v.set(f);
        for (size_t row = 0; row < rr.rows.size(); ++row) {
            int p = rr.pivots[row];
            if (rr.rows[row].test(f)) v.flip(p);
        }
        basis.push_back(std::move(v));
    }
    return basis;
}

inline AffineSystemSolution solve_affine(
    const std::vector<BitRow>& input_rows,
    const std::vector<uint8_t>& input_rhs,
    int nvars
) {
    RrefResult rr = bit_rref(input_rows, input_rhs, nvars);

    AffineSystemSolution sol;
    sol.particular = BitRow(nvars);
    if (!rr.ok) return sol;

    sol.ok = true;

    std::vector<uint8_t> is_pivot(nvars, 0);
    for (int p : rr.pivots) is_pivot[p] = 1;

    for (size_t i = 0; i < rr.pivots.size(); ++i) {
        if (rr.rhs[i]) sol.particular.set(rr.pivots[i]);
    }

    for (int f = 0; f < nvars; ++f) {
        if (is_pivot[f]) continue;

        BitRow v(nvars);
        v.set(f);
        for (size_t row = 0; row < rr.rows.size(); ++row) {
            int p = rr.pivots[row];
            if (rr.rows[row].test(f)) v.flip(p);
        }
        sol.kernel.push_back(std::move(v));
    }

    return sol;
}

inline BitRow random_affine_point(
    const AffineSystemSolution& sol,
    std::mt19937_64& rng
) {
    if (!sol.ok) throw std::runtime_error("cannot sample inconsistent affine space");

    BitRow z = sol.particular;
    for (const BitRow& k : sol.kernel) {
        if (rng() & 1ull) z.xor_with(k);
    }
    return z;
}

} // namespace laph
