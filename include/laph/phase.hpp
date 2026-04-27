#pragma once

#include "laph/types.hpp"

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace laph {

using ANF = std::unordered_set<Monomial, MonomialHash>;

void add_mod8(PhasePoly& p, const Monomial& m, int coeff);
void add_affine_phase(PhasePoly& poly, int coeff, const XorSet& mask, bool constant);
bool is_clifford_poly(const PhasePoly& poly);
PhasePoly substitute_clifford_poly(
    const PhasePoly& poly,
    const std::unordered_map<int, AffineForm>& exprs
);
PhasePoly substitute_clifford_poly_vec(
    const PhasePoly& poly,
    const std::vector<AffineForm>& exprs
);
ScaledComplex clifford_sum(PhasePoly poly, int nvars);

} // namespace laph
