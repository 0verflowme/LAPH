#pragma once

#include "laph/gf2_dense.hpp"
#include "laph/optimizer.hpp"

#include <random>
#include <vector>

namespace laph {

class LAPH;

class DensityKernel {
public:
    const LAPH& st;
    const int total_vars;
    PhasePoly phase;
    DenseSolveResult base;
    std::vector<int> cut;
    std::vector<AffineForm> cut_forms_base;
    std::vector<AffineForm> prefix_forms_base;

    explicit DensityKernel(const LAPH& s);

    ScaledComplex partition_from_prefix_state(
        const DenseRrefState& prefix_state
    ) const;
    long double probability_from_prefix_state(
        const DenseRrefState& prefix_state
    ) const;
    std::vector<int> sample(std::mt19937_64& rng) const;
};

std::vector<int> exact_sample_cached_density(
    const LAPH& st,
    std::mt19937_64& rng
);

} // namespace laph
