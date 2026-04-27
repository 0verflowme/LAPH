#pragma once

#include "laph/math/stationary_clifford.hpp"
#include "laph/phase.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace laph {

struct CliffordPhase {
    int m = 0;
    uint8_t constant = 0;
    std::vector<uint8_t> linear;
    std::vector<std::pair<int, int>> quad4;
};

inline CliffordPhase extract_clifford_phase(const PhasePoly& phase, int m) {
    CliffordPhase ph;
    ph.m = m;
    ph.linear.assign(m, 0);

    for (const auto& kv : phase) {
        const Monomial& mon = kv.first;
        uint8_t coeff = static_cast<uint8_t>(mod8(kv.second));
        if (!coeff) continue;

        size_t deg = mon.degree();
        if (deg == 0) {
            ph.constant = static_cast<uint8_t>((ph.constant + coeff) & 7);
        } else if (deg == 1) {
            if (coeff & 1u) {
                throw std::runtime_error("not Clifford: odd linear phase");
            }
            int i = mon.v[0];
            if (i < 0 || i >= m) {
                throw std::runtime_error("linear phase variable out of range");
            }
            ph.linear[i] = static_cast<uint8_t>((ph.linear[i] + coeff) & 7);
        } else if (deg == 2) {
            if (coeff != 4) {
                throw std::runtime_error("not Clifford: non-4 quadratic phase");
            }
            int i = mon.v[0];
            int j = mon.v[1];
            if (i < 0 || j < 0 || i >= m || j >= m || i == j) {
                throw std::runtime_error("bad quadratic monomial");
            }
            if (i > j) std::swap(i, j);
            ph.quad4.push_back({i, j});
        } else {
            throw std::runtime_error("not Clifford: degree >= 3 phase");
        }
    }

    return ph;
}

inline std::vector<BitRow> polar_rows(const CliffordPhase& ph) {
    std::vector<BitRow> B;
    B.reserve(ph.m);
    for (int i = 0; i < ph.m; ++i) B.emplace_back(ph.m);

    for (int i = 0; i < ph.m; ++i) {
        uint8_t a = ph.linear[i] & 7u;
        if ((a & 3u) == 2u) B[i].flip(i);
    }

    for (auto [i, j] : ph.quad4) {
        B[i].flip(j);
        B[j].flip(i);
    }

    return B;
}

inline BitRow apply_polar(
    const std::vector<BitRow>& B,
    const BitRow& v
) {
    int m = v.n;
    BitRow out(m);
    for (int i = 0; i < m; ++i) {
        if (B[i].dot(v)) out.set(i);
    }
    return out;
}

struct DerivativeEquation {
    bool possible = true;
    BitRow row;
    uint8_t rhs = 0;
};

inline DerivativeEquation derivative_zero_equation(
    const CliffordPhase& ph,
    const BitRow& r
) {
    DerivativeEquation eq;
    eq.row = BitRow(ph.m);

    uint8_t c0 = 0;

    for (int i = 0; i < ph.m; ++i) {
        if (!r.test(i)) continue;

        uint8_t a = ph.linear[i] & 7u;
        c0 = static_cast<uint8_t>((c0 + a) & 7u);

        if ((a & 3u) == 2u) eq.row.flip(i);
    }

    for (auto [i, j] : ph.quad4) {
        bool ri = r.test(i);
        bool rj = r.test(j);

        if (ri && rj) c0 = static_cast<uint8_t>((c0 + 4u) & 7u);
        if (ri) eq.row.flip(j);
        if (rj) eq.row.flip(i);
    }

    if (c0 & 3u) {
        eq.possible = false;
        return eq;
    }

    eq.rhs = static_cast<uint8_t>((c0 >> 2) & 1u);
    return eq;
}

} // namespace laph
