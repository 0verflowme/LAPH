#pragma once

#include "laph/math/stationary_clifford.hpp"

#include <stdexcept>
#include <vector>

namespace laph {

template <class State>
int hidden_interference_rank(
    const State& st,
    const std::vector<int>& cut_variables
) {
    const int m = st.m;
    const int n = st.n;

    std::vector<BitRow> homogeneous;
    homogeneous.reserve(st.constraints.size() + st.visible.size());

    for (const Constraint& c : st.constraints) {
        homogeneous.push_back(BitRow::from_xorset(c.lhs, m));
    }

    for (int q = 0; q < n; ++q) {
        homogeneous.push_back(BitRow::from_xorset(st.visible[q].vars, m));
    }

    std::vector<BitRow> K = nullspace(homogeneous, m);

    std::vector<BitRow> projection_rows;
    projection_rows.reserve(cut_variables.size());

    for (int v : cut_variables) {
        if (v < 0 || v >= m) {
            throw std::runtime_error("cut variable out of range");
        }

        BitRow row(static_cast<int>(K.size()));
        for (size_t j = 0; j < K.size(); ++j) {
            if (K[j].test(v)) row.set(static_cast<int>(j));
        }
        projection_rows.push_back(std::move(row));
    }

    std::vector<uint8_t> zero_rhs(projection_rows.size(), 0);
    RrefResult rr = bit_rref(
        projection_rows,
        zero_rhs,
        static_cast<int>(K.size())
    );

    if (!rr.ok) throw std::runtime_error("rank computation failed");
    return static_cast<int>(rr.pivots.size());
}

} // namespace laph
