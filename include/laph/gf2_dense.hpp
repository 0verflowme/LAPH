#pragma once

#include "laph/types.hpp"

#include <cstdint>
#include <vector>

namespace laph {

struct DenseRow {
    int nbits = 0;
    std::vector<uint64_t> w;

    DenseRow() = default;
    explicit DenseRow(int n);

    static DenseRow from_xorset(const XorSet& xs, int nbits);

    bool test(int i) const;
    void xor_with(const DenseRow& b);
    bool empty() const;
    int first_set() const;
    XorSet to_xorset() const;
};

struct DenseSolveResult {
    bool ok = false;
    int nfree = 0;
    std::vector<AffineForm> expr;
};

struct DenseRrefState {
    int nvars = 0;
    std::vector<DenseRow> row;
    std::vector<unsigned char> rhs;
    std::vector<unsigned char> has;
    std::vector<int> pivots;

    DenseRrefState() = default;
    explicit DenseRrefState(int n);

    bool insert(DenseRow r, bool b);
    bool insert_constraint(const Constraint& c);
    DenseSolveResult solve() const;
};

DenseSolveResult solve_affine_system_dense(
    const std::vector<Constraint>& rows,
    int nvars
);
AffineForm compose_affine(
    const AffineForm& outer,
    const std::vector<AffineForm>& inner
);
AffineForm project_xorset(
    const XorSet& xs,
    const std::vector<AffineForm>& expr
);
Constraint affine_equals_bit(const AffineForm& f, bool bit);

} // namespace laph
