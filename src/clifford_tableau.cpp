#include "laph/clifford_tableau.hpp"

#include <algorithm>

namespace laph {

CliffordTableau::CliffordTableau(int nqubits) {
    reset(nqubits);
}

void CliffordTableau::reset(int nqubits) {
    n = nqubits;
    words = (n + 63) >> 6;
    row.assign(2 * n, Row{});

    for (Row& r : row) {
        r.x.assign(words, 0);
        r.z.assign(words, 0);
        r.phase = 0;
    }

    for (int q = 0; q < n; ++q) {
        set_x_row(q, q);
        set_z_row(n + q, q);
    }
}

bool CliffordTableau::test_x(int r, int q) const {
    return (row[r].x[static_cast<size_t>(q) >> 6] >> (q & 63)) & 1ull;
}

bool CliffordTableau::test_z(int r, int q) const {
    return (row[r].z[static_cast<size_t>(q) >> 6] >> (q & 63)) & 1ull;
}

void CliffordTableau::flip_x(int r, int q) {
    row[r].x[static_cast<size_t>(q) >> 6] ^= 1ull << (q & 63);
}

void CliffordTableau::flip_z(int r, int q) {
    row[r].z[static_cast<size_t>(q) >> 6] ^= 1ull << (q & 63);
}

void CliffordTableau::set_identity(int r) {
    std::fill(row[r].x.begin(), row[r].x.end(), 0);
    std::fill(row[r].z.begin(), row[r].z.end(), 0);
    row[r].phase = 0;
}

void CliffordTableau::set_x_row(int r, int q, unsigned char phase) {
    set_identity(r);
    flip_x(r, q);
    row[r].phase = phase & 3u;
}

void CliffordTableau::set_z_row(int r, int q, unsigned char phase) {
    set_identity(r);
    flip_z(r, q);
    row[r].phase = phase & 3u;
}

void CliffordTableau::x_gate(int q) {
    for (int r = 0; r < 2 * n; ++r) {
        if (test_z(r, q)) row[r].phase = (row[r].phase + 2) & 3u;
    }
}

void CliffordTableau::z_gate(int q) {
    for (int r = 0; r < 2 * n; ++r) {
        if (test_x(r, q)) row[r].phase = (row[r].phase + 2) & 3u;
    }
}

void CliffordTableau::h(int q) {
    for (int r = 0; r < 2 * n; ++r) {
        bool x = test_x(r, q);
        bool z = test_z(r, q);
        if (x && z) row[r].phase = (row[r].phase + 2) & 3u;
        if (x != z) {
            flip_x(r, q);
            flip_z(r, q);
        }
    }
}

void CliffordTableau::s(int q) {
    for (int r = 0; r < 2 * n; ++r) {
        bool x = test_x(r, q);
        if (x) row[r].phase = (row[r].phase + 1) & 3u;
        if (x) flip_z(r, q);
    }
}

void CliffordTableau::cnot(int control, int target) {
    for (int r = 0; r < 2 * n; ++r) {
        bool xc = test_x(r, control);
        bool zt = test_z(r, target);
        if (xc) flip_x(r, target);
        if (zt) flip_z(r, control);
    }
}

void CliffordTableau::cz(int a, int b) {
    h(b);
    cnot(a, b);
    h(b);
}

void CliffordTableau::row_multiply(int target, int source) {
    unsigned phase = row[target].phase + row[source].phase;

    for (int w = 0; w < words; ++w) {
        phase += 2u * static_cast<unsigned>(
            __builtin_popcountll(row[target].z[w] & row[source].x[w]) & 1u
        );
    }

    for (int w = 0; w < words; ++w) {
        row[target].x[w] ^= row[source].x[w];
        row[target].z[w] ^= row[source].z[w];
    }
    row[target].phase = static_cast<unsigned char>(phase & 3u);
}

bool CliffordTableau::measure_z(int q, std::mt19937_64& rng) {
    int pivot = -1;
    for (int r = n; r < 2 * n; ++r) {
        if (test_x(r, q)) {
            pivot = r;
            break;
        }
    }

    if (pivot >= 0) {
        for (int r = 0; r < 2 * n; ++r) {
            if (r != pivot && test_x(r, q)) row_multiply(r, pivot);
        }

        row[pivot - n] = row[pivot];
        bool bit = (rng() & 1ull) != 0;
        set_z_row(pivot, q, bit ? 2 : 0);
        return bit;
    }

    Row scratch;
    scratch.x.assign(words, 0);
    scratch.z.assign(words, 0);
    scratch.phase = 0;
    row.push_back(std::move(scratch));
    int scratch_index = static_cast<int>(row.size()) - 1;

    for (int r = 0; r < n; ++r) {
        if (test_x(r, q)) row_multiply(scratch_index, n + r);
    }

    bool bit = ((row[scratch_index].phase >> 1) & 1u) != 0;
    row.pop_back();
    return bit;
}

std::vector<int> CliffordTableau::sample_all(std::mt19937_64& rng) const {
    CliffordTableau tmp = *this;
    std::vector<int> sample(n, 0);
    for (int q = 0; q < n; ++q) {
        sample[q] = tmp.measure_z(q, rng) ? 1 : 0;
    }
    return sample;
}

} // namespace laph
