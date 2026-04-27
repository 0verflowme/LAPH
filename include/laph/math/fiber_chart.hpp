#pragma once

#include "laph/math/stationary_clifford.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace laph {

struct FiberChart {
    int m = 0;
    int n = 0;
    int out_dim = 0;
    int hid_dim = 0;

    BitRow z0;
    BitRow x0;

    std::vector<BitRow> out_lift;
    std::vector<BitRow> out_x;
    std::vector<BitRow> hid_lift;

    std::vector<uint8_t> var_const;
    std::vector<BitRow> var_out;
    std::vector<BitRow> var_hid;
};

template <class State>
BitRow visible_delta(const State& st, const BitRow& dz) {
    BitRow out(st.n);
    for (int q = 0; q < st.n; ++q) {
        if (BitRow::from_xorset(st.visible[q].vars, st.m).dot(dz)) out.set(q);
    }
    return out;
}

template <class State>
BitRow visible_point(const State& st, const BitRow& z) {
    BitRow out(st.n);
    for (int q = 0; q < st.n; ++q) {
        bool bit = BitRow::from_xorset(st.visible[q].vars, st.m).dot(z) ^ st.visible[q].c;
        if (bit) out.set(q);
    }
    return out;
}

struct WitnessOutputBasis {
    int n = 0;
    int m = 0;
    std::vector<BitRow> row;
    std::vector<BitRow> witness;
    std::vector<uint8_t> has;
    std::vector<int> pivots;
    std::vector<BitRow> hidden;

    WitnessOutputBasis(int nbits, int mbits)
        : n(nbits), m(mbits), row(nbits), witness(nbits), has(nbits, 0) {
        for (int i = 0; i < n; ++i) {
            row[i] = BitRow(n);
            witness[i] = BitRow(m);
        }
    }

    void insert(BitRow x, BitRow z) {
        for (int p : pivots) {
            if (x.test(p)) {
                x.xor_with(row[p]);
                z.xor_with(witness[p]);
            }
        }

        if (x.empty()) {
            if (!z.empty()) hidden.push_back(std::move(z));
            return;
        }

        int p = x.first_set();
        for (int q : pivots) {
            if (row[q].test(p)) {
                row[q].xor_with(x);
                witness[q].xor_with(z);
            }
        }

        has[p] = 1;
        row[p] = std::move(x);
        witness[p] = std::move(z);
        pivots.insert(std::lower_bound(pivots.begin(), pivots.end(), p), p);
    }

    std::vector<BitRow> output_rows() const {
        std::vector<BitRow> out;
        out.reserve(pivots.size());
        for (int p : pivots) out.push_back(row[p]);
        return out;
    }

    std::vector<BitRow> output_lifts() const {
        std::vector<BitRow> out;
        out.reserve(pivots.size());
        for (int p : pivots) out.push_back(witness[p]);
        return out;
    }
};

template <class State>
FiberChart build_fiber_chart(const State& st) {
    std::vector<BitRow> rows;
    std::vector<uint8_t> rhs;
    rows.reserve(st.constraints.size());
    rhs.reserve(st.constraints.size());

    for (const Constraint& c : st.constraints) {
        rows.push_back(BitRow::from_xorset(c.lhs, st.m));
        rhs.push_back(static_cast<uint8_t>(c.rhs));
    }

    AffineSystemSolution sol = solve_affine(rows, rhs, st.m);
    if (!sol.ok) throw std::runtime_error("LAPH constraint system is inconsistent");

    WitnessOutputBasis basis(st.n, st.m);
    for (const BitRow& k : sol.kernel) basis.insert(visible_delta(st, k), k);

    FiberChart ch;
    ch.m = st.m;
    ch.n = st.n;
    ch.z0 = sol.particular;
    ch.x0 = visible_point(st, sol.particular);
    ch.out_x = basis.output_rows();
    ch.out_lift = basis.output_lifts();
    ch.hid_lift = basis.hidden;
    ch.out_dim = static_cast<int>(ch.out_lift.size());
    ch.hid_dim = static_cast<int>(ch.hid_lift.size());

    ch.var_const.assign(ch.m, 0);
    ch.var_out.reserve(ch.m);
    ch.var_hid.reserve(ch.m);

    for (int v = 0; v < ch.m; ++v) {
        ch.var_const[v] = static_cast<uint8_t>(ch.z0.test(v));

        BitRow out_row(ch.out_dim);
        for (int j = 0; j < ch.out_dim; ++j) {
            if (ch.out_lift[j].test(v)) out_row.set(j);
        }
        ch.var_out.push_back(std::move(out_row));

        BitRow hid_row(ch.hid_dim);
        for (int j = 0; j < ch.hid_dim; ++j) {
            if (ch.hid_lift[j].test(v)) hid_row.set(j);
        }
        ch.var_hid.push_back(std::move(hid_row));
    }

    return ch;
}

} // namespace laph
