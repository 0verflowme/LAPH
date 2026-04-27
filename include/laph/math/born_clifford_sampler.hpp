#pragma once

#include "laph/math/clifford_phase.hpp"

#include <random>
#include <stdexcept>
#include <vector>

namespace laph {

template <class State>
std::vector<int> sample_born_clifford(
    const State& st,
    std::mt19937_64& rng
) {
    const int m = st.m;
    const int n = st.n;

    CliffordPhase ph = extract_clifford_phase(st.phase, m);

    std::vector<BitRow> homogeneous;
    homogeneous.reserve(st.constraints.size() + st.visible.size());

    for (const Constraint& c : st.constraints) {
        homogeneous.push_back(BitRow::from_xorset(c.lhs, m));
    }

    for (int q = 0; q < n; ++q) {
        homogeneous.push_back(BitRow::from_xorset(st.visible[q].vars, m));
    }

    std::vector<BitRow> K = nullspace(homogeneous, m);
    std::vector<BitRow> B = polar_rows(ph);

    std::vector<BitRow> BK;
    BK.reserve(K.size());
    for (const BitRow& k : K) BK.push_back(apply_polar(B, k));

    std::vector<BitRow> restricted;
    restricted.reserve(K.size());
    for (size_t row = 0; row < K.size(); ++row) {
        BitRow rr(static_cast<int>(K.size()));
        for (size_t col = 0; col < K.size(); ++col) {
            if (K[row].dot(BK[col])) rr.set(static_cast<int>(col));
        }
        restricted.push_back(std::move(rr));
    }

    std::vector<BitRow> radical_coordinates =
        nullspace(restricted, static_cast<int>(K.size()));

    std::vector<BitRow> rows;
    std::vector<uint8_t> rhs;
    rows.reserve(st.constraints.size() + radical_coordinates.size());
    rhs.reserve(st.constraints.size() + radical_coordinates.size());

    for (const Constraint& c : st.constraints) {
        rows.push_back(BitRow::from_xorset(c.lhs, m));
        rhs.push_back(static_cast<uint8_t>(c.rhs));
    }

    for (const BitRow& alpha : radical_coordinates) {
        BitRow r(m);
        for (size_t i = 0; i < K.size(); ++i) {
            if (alpha.test(static_cast<int>(i))) r.xor_with(K[i]);
        }

        DerivativeEquation eq = derivative_zero_equation(ph, r);
        if (!eq.possible) {
            throw std::runtime_error("zero-norm or inconsistent Clifford branch");
        }

        rows.push_back(std::move(eq.row));
        rhs.push_back(eq.rhs);
    }

    AffineSystemSolution sol = solve_affine(rows, rhs, m);
    if (!sol.ok) throw std::runtime_error("empty Clifford support");

    BitRow z = random_affine_point(sol, rng);
    std::vector<int> x(n, 0);

    for (int q = 0; q < n; ++q) {
        bool bit = BitRow::from_xorset(st.visible[q].vars, m).dot(z) ^ st.visible[q].c;
        x[q] = bit ? 1 : 0;
    }

    return x;
}

} // namespace laph
