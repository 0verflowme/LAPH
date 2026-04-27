#include "laph/density_kernel.hpp"

#include "laph/laph.hpp"

#include <stdexcept>
#include <unordered_map>

namespace laph {

DensityKernel::DensityKernel(const LAPH& s)
    : st(s), total_vars(2 * s.m), phase(s.density_phase()) {
    std::vector<Constraint> base_constraints = st.density_constraints({});
    base = solve_affine_system_dense(base_constraints, total_vars);
    if (!base.ok) {
        throw std::runtime_error("inconsistent base density system");
    }

    cut = greedy_cut_set_for_phase(phase);
    if (cut.size() >= 63) {
        throw std::runtime_error(
            "density cut too large for single-kernel sampler; split components first"
        );
    }

    cut_forms_base.reserve(cut.size());
    for (int x : cut) cut_forms_base.push_back(base.expr[x]);

    prefix_forms_base.reserve(st.n);
    for (int q = 0; q < st.n; ++q) {
        AffineForm f = AffineForm::constant(st.visible[q].c);
        f.xor_with(project_xorset(st.visible[q].vars, base.expr));
        prefix_forms_base.push_back(std::move(f));
    }
}

ScaledComplex DensityKernel::partition_from_prefix_state(
    const DenseRrefState& prefix_state
) const {
    ScaledComplex total = ScaledComplex::zero();
    const uint64_t lim = 1ull << cut.size();

    for (uint64_t mask = 0; mask < lim; ++mask) {
        DenseRrefState qs = prefix_state;
        std::unordered_map<int, bool> assignment;
        assignment.reserve(cut.size() * 2 + 1);

        bool ok = true;
        for (size_t i = 0; i < cut.size(); ++i) {
            bool bit = ((mask >> i) & 1ull) != 0;
            assignment[cut[i]] = bit;
            if (!qs.insert_constraint(affine_equals_bit(cut_forms_base[i], bit))) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;

        DenseSolveResult qsol = qs.solve();
        if (!qsol.ok) continue;

        std::vector<AffineForm> composed(total_vars);
        for (int x = 0; x < total_vars; ++x) {
            composed[x] = compose_affine(base.expr[x], qsol.expr);
        }

        PhasePoly ph = restrict_poly(phase, assignment);
        PhasePoly residual = substitute_clifford_poly_vec(ph, composed);
        if (!is_clifford_poly(residual)) {
            throw std::runtime_error("density cut failed");
        }

        ScaledComplex term = clifford_sum(residual, qsol.nfree);
        total.add(term);
    }

    total.mul_sqrt2_power(-2LL * st.scale);
    return total;
}

long double DensityKernel::probability_from_prefix_state(
    const DenseRrefState& prefix_state
) const {
    ScaledComplex z = partition_from_prefix_state(prefix_state);
    long double p = z.real_ld();
    if (p < 0 && p > -1e-12L) p = 0;
    if (p > 1 && p < 1 + 1e-12L) p = 1;
    return p;
}

std::vector<int> DensityKernel::sample(std::mt19937_64& rng) const {
    DenseRrefState prefix_state(base.nfree);
    std::vector<int> sample(st.n, 0);

    long double prefix_prob = probability_from_prefix_state(prefix_state);
    if (prefix_prob <= 0) throw std::runtime_error("state has zero norm");

    std::uniform_real_distribution<long double> U(0.0L, 1.0L);

    for (int q = 0; q < st.n; ++q) {
        const AffineForm& fq = prefix_forms_base[q];

        DenseRrefState cand0 = prefix_state;
        bool ok0 = cand0.insert_constraint(affine_equals_bit(fq, false));
        long double p0_abs = ok0 ? probability_from_prefix_state(cand0) : 0.0L;

        long double p0_cond = p0_abs / prefix_prob;
        if (p0_cond < 0) p0_cond = 0;
        if (p0_cond > 1) p0_cond = 1;

        bool bit = !(U(rng) < p0_cond);
        sample[q] = bit ? 1 : 0;

        bool ok = prefix_state.insert_constraint(affine_equals_bit(fq, bit));
        if (!ok) throw std::runtime_error("sampled zero-probability branch");

        prefix_prob = bit ? (prefix_prob - p0_abs) : p0_abs;
        if (prefix_prob < 0 && prefix_prob > -1e-12L) prefix_prob = 0;
        if (prefix_prob > 1 && prefix_prob < 1 + 1e-12L) prefix_prob = 1;
    }

    return sample;
}

std::vector<int> exact_sample_cached_density(
    const LAPH& st,
    std::mt19937_64& rng
) {
    DensityKernel kernel(st);
    return kernel.sample(rng);
}

} // namespace laph
