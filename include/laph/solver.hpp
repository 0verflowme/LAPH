#pragma once

#include "laph/phase.hpp"

#include <unordered_map>
#include <vector>

namespace laph {

struct SolveResult {
    bool ok = false;
    int nfree = 0;
    std::unordered_map<int, AffineForm> expr;
};

SolveResult solve_affine_system(
    const std::vector<Constraint>& rows_in,
    const std::vector<int>& variables
);
PhasePoly restrict_poly(
    const PhasePoly& poly,
    const std::unordered_map<int, bool>& assign
);
std::vector<Constraint> restrict_constraints(
    const std::vector<Constraint>& cs,
    const std::unordered_map<int, bool>& assign
);
std::vector<int> greedy_cut_set_for_phase(const PhasePoly& phase);
ScaledComplex exact_partition_sum(
    int total_vars,
    const PhasePoly& phase,
    const std::vector<Constraint>& constraints,
    const std::vector<int>* provided_cut = nullptr
);

} // namespace laph
