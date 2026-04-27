#pragma once

#include <cstdint>
#include <random>
#include <vector>

namespace laph {

class CliffordTableau {
public:
    struct Row {
        std::vector<uint64_t> x;
        std::vector<uint64_t> z;
        unsigned char phase = 0; // i^phase * X^x Z^z
    };

    int n = 0;
    int words = 0;
    std::vector<Row> row;

    CliffordTableau() = default;
    explicit CliffordTableau(int nqubits);

    void reset(int nqubits);

    void x_gate(int q);
    void z_gate(int q);
    void h(int q);
    void s(int q);
    void cnot(int control, int target);
    void cz(int a, int b);

    std::vector<int> sample_all(std::mt19937_64& rng) const;

private:
    bool test_x(int r, int q) const;
    bool test_z(int r, int q) const;
    void flip_x(int r, int q);
    void flip_z(int r, int q);
    void set_identity(int r);
    void set_x_row(int r, int q, unsigned char phase = 0);
    void set_z_row(int r, int q, unsigned char phase = 0);
    void row_multiply(int target, int source);
    bool measure_z(int q, std::mt19937_64& rng);
};

} // namespace laph
