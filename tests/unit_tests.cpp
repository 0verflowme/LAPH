#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "laph/laph.hpp"
#include "laph/math/fabqnf.hpp"

namespace {

bool near(long double a, long double b, long double eps = 1e-9L) {
    return std::fabsl(a - b) <= eps;
}

void assert_complex_near(
    const laph::ScaledComplex& z,
    long double re,
    long double im,
    long double eps = 1e-9L
) {
    assert(near(z.real_ld(), re, eps));
    assert(near(z.imag_ld(), im, eps));
}

laph::LAPH demo_state() {
    laph::LAPH st(3);
    st.h(0).h(1).cnot(0, 2).t(2).cz(0, 1).h(2).ccz(0, 1, 2);
    st.compress();
    return st;
}

uint64_t bitrow_value(const laph::BitRow& row) {
    uint64_t value = 0;
    for (int i = 0; i < row.n; ++i) {
        if (row.test(i)) value |= (1ull << i);
    }
    return value;
}

laph::BitRow bitrow_from_value(int nbits, uint64_t value) {
    laph::BitRow row(nbits);
    for (int i = 0; i < nbits; ++i) {
        if ((value >> i) & 1ull) row.set(i);
    }
    return row;
}

std::vector<long double> exact_distribution(laph::LAPH& st) {
    std::vector<long double> dist(1ull << st.n, 0.0L);
    long double total = 0.0L;

    for (uint64_t x = 0; x < dist.size(); ++x) {
        laph::ScaledComplex amp = st.amplitude_factorized(x);
        long double p = amp.real_ld() * amp.real_ld() + amp.imag_ld() * amp.imag_ld();
        dist[x] = p;
        total += p;
    }

    for (long double& p : dist) p /= total;
    return dist;
}

std::vector<long double> fabqnf_distribution(const laph::LAPH& st) {
    laph::FABQNF law = laph::build_fabqnf(st);
    std::vector<long double> dist(1ull << st.n, 0.0L);
    long double total = 0.0L;
    long double cell_size = std::ldexp(1.0L, law.out_dim - law.rho);

    for (uint64_t tv = 0; tv < (1ull << law.out_dim); ++tv) {
        laph::BitRow t = bitrow_from_value(law.out_dim, tv);
        uint64_t y = law.quotient_value(t);
        long double w = law.cell_weight[y] / cell_size;
        uint64_t x = bitrow_value(law.output_from_t(t));
        dist[x] += w;
        total += w;
    }

    for (long double& p : dist) p /= total;
    return dist;
}

void assert_distributions_near(
    const std::vector<long double>& a,
    const std::vector<long double>& b,
    long double eps = 1e-8L
) {
    assert(a.size() == b.size());
    long double l1 = 0.0L;
    for (size_t i = 0; i < a.size(); ++i) l1 += std::fabsl(a[i] - b[i]);
    assert(l1 <= eps);
}

int conservative_fabqnf_rho(const laph::LAPH& st) {
    laph::FiberChart ch = laph::build_fiber_chart(st);
    std::vector<laph::BitRow> rows;

    for (const auto& kv : st.phase) {
        int coeff = laph::mod8(kv.second);
        if (!coeff || kv.first.empty()) continue;
        if (!laph::is_fiber_active(kv.first, ch)) continue;

        for (int v : kv.first.v) {
            if (!ch.var_out[v].empty()) rows.push_back(ch.var_out[v]);
        }
    }

    std::vector<uint8_t> rhs(rows.size(), 0);
    laph::RrefResult rr = laph::bit_rref(rows, rhs, ch.out_dim);
    assert(rr.ok);
    return static_cast<int>(rr.rows.size());
}

laph::LAPH hidden_translation_reduction_state() {
    laph::LAPH st(3);
    st.tableau_valid = false;
    st.m = 5;

    st.visible[0] = laph::AffineForm(laph::XorSet({0}), false);
    st.visible[1] = laph::AffineForm(laph::XorSet({0, 1, 2}), false);
    st.visible[2] = laph::AffineForm(laph::XorSet({1, 2, 4}), false);

    st.add_phase_monomial(1, laph::Monomial::singleton(4));
    st.add_phase_monomial(1, laph::Monomial::singleton(1));
    st.add_phase_monomial(2, laph::Monomial::singleton(3));

    st.compress();
    return st;
}

laph::LAPH repeated_hidden_packet_state() {
    laph::LAPH st(3);
    st.tableau_valid = false;

    int hidden = st.new_var();
    std::vector<int> out_vars;
    for (int q = 0; q < 3; ++q) {
        int v = st.new_var();
        out_vars.push_back(v);
        st.visible[q] = laph::AffineForm::variable(v);
    }

    for (int v : out_vars) {
        laph::XorSet xs({hidden, v});
        st.add_phase_product_forms(1, {laph::AffineForm(xs, false)});
    }

    st.compress();
    return st;
}

void test_demo_amplitudes() {
    laph::LAPH st = demo_state();
    long double s = 0.35355339059327376220L;

    assert_complex_near(st.amplitude_factorized(0), s, 0.0L);
    assert_complex_near(st.amplitude_factorized(1), 0.25L, 0.25L);
    assert_complex_near(st.amplitude_factorized(2), s, 0.0L);
    assert_complex_near(st.amplitude_factorized(3), -0.25L, -0.25L);
    assert_complex_near(st.amplitude_factorized(4), s, 0.0L);
    assert_complex_near(st.amplitude_factorized(5), -0.25L, -0.25L);
    assert_complex_near(st.amplitude_factorized(6), s, 0.0L);
    assert_complex_near(st.amplitude_factorized(7), -0.25L, -0.25L);

    assert(st.stats().non_clifford_cut == 1);
}

void test_factorized_matches_monolithic() {
    laph::LAPH st(6);
    st.h(0).t(0);
    st.h(1).s(1);
    st.h(2).z(2);
    st.h(3).t(3);
    st.h(4);
    st.h(5).cz(4, 5);
    st.compress();

    for (uint64_t x = 0; x < 64; ++x) {
        laph::ScaledComplex a = st.amplitude_monolithic(x);
        laph::ScaledComplex b = st.amplitude_factorized(x);
        assert(near(a.real_ld(), b.real_ld()));
        assert(near(a.imag_ld(), b.imag_ld()));
    }

    std::vector<std::pair<int, bool>> none;
    assert(near(st.probability_prefix_monolithic(none), 1.0L));
    assert(near(st.probability_prefix_factorized(none), 1.0L));
    assert(near(
        st.probability_prefix_monolithic({{0, false}, {3, true}}),
        st.probability_prefix_factorized({{0, false}, {3, true}})
    ));
    assert(st.connected_components().size() > 1);
}

void test_clifford_cut_zero() {
    laph::LAPH st(5);
    st.h(0).h(1).h(2).cnot(0, 3).cz(1, 2).s(4).z(3);
    st.compress();

    assert(st.cut_set().empty());
    std::vector<std::pair<int, bool>> none;
    assert(near(st.probability_prefix_factorized(none), 1.0L));
}

void test_exact_sample_shape() {
    laph::LAPH st = demo_state();
    std::mt19937_64 rng(123);
    auto sample = st.exact_sample_factorized(rng);

    assert(sample.size() == 3);
    for (int bit : sample) assert(bit == 0 || bit == 1);

    std::mt19937_64 rng2(123);
    auto component_sample = st.exact_sample_by_component(rng2);
    assert(component_sample == sample);
}

void test_component_sampler_keeps_constants() {
    laph::LAPH st(4);
    st.h(0).t(0);
    st.h(1).t(1);
    st.x(3);
    st.compress();

    assert(st.state_components().size() == 2);

    std::mt19937_64 rng(7);
    auto sample = st.exact_sample_by_component(rng);
    assert(sample.size() == 4);
    assert(sample[2] == 0);
    assert(sample[3] == 1);
}

void test_tableau_bell_sampler() {
    for (uint64_t seed = 0; seed < 16; ++seed) {
        laph::LAPH st(2);
        st.h(0).cnot(0, 1);

        assert(st.tableau_valid);
        std::mt19937_64 rng(seed);
        auto sample = st.exact_sample(rng);

        assert(sample.size() == 2);
        assert(sample[0] == sample[1]);
    }
}

void test_tableau_deterministic_paulis() {
    laph::LAPH st(3);
    st.x(0).h(1).cnot(1, 2);

    for (uint64_t seed = 0; seed < 16; ++seed) {
        std::mt19937_64 rng(seed);
        auto sample = st.exact_sample(rng);

        assert(sample[0] == 1);
        assert(sample[1] == sample[2]);
    }
}

void test_non_clifford_invalidates_tableau() {
    laph::LAPH st(1);
    st.h(0).t(0);
    assert(!st.tableau_valid);

    std::mt19937_64 rng(4);
    auto sample = st.exact_sample(rng);
    assert(sample.size() == 1);
    assert(sample[0] == 0 || sample[0] == 1);
}

void test_hidden_interference_rank_distinguishes_output_phase() {
    laph::LAPH output_phase(1);
    output_phase.h(0).t(0);
    output_phase.compress();

    laph::Stats output_stats = output_phase.stats();
    assert(output_stats.non_clifford_cut == 1);
    assert(output_stats.hidden_interference_rank == 0);

    laph::LAPH hidden_phase(1);
    hidden_phase.h(0).t(0).h(0);
    hidden_phase.compress();

    laph::Stats hidden_stats = hidden_phase.stats();
    assert(hidden_stats.non_clifford_cut == 1);
    assert(hidden_stats.hidden_interference_rank == 1);
}

void test_fabqnf_output_only_t_has_rho_zero() {
    laph::LAPH st(1);
    st.h(0).t(0);
    st.compress();

    laph::FABQNFStats stats = laph::fabqnf_stats(st);
    assert(stats.rho == 0);
    assert(stats.fiber_active_terms == 0);

    assert_distributions_near(exact_distribution(st), fabqnf_distribution(st));
}

void test_fabqnf_hidden_t_matches_amplitudes() {
    laph::LAPH st(2);
    st.h(0).h(1).cnot(0, 1).t(1).h(1);
    st.compress();

    laph::FABQNFStats stats = laph::fabqnf_stats(st);
    assert(stats.rho < 4);

    assert_distributions_near(exact_distribution(st), fabqnf_distribution(st));
}

void test_fabqnf_hidden_translation_reduces_rho() {
    laph::LAPH st = hidden_translation_reduction_state();

    int old_rho = conservative_fabqnf_rho(st);
    laph::FABQNF law = laph::build_fabqnf(st);
    laph::FABQNFStats stats = laph::fabqnf_stats(st);

    assert(old_rho == 1);
    assert(stats.rho == 0);
    assert(law.rho == 0);
    assert(law.rho < old_rho);

    assert_distributions_near(exact_distribution(st), fabqnf_distribution(st));
}

void test_fabqnf_hidden_packet_cache_hits() {
    laph::LAPH st = repeated_hidden_packet_state();
    laph::FABQNF law = laph::build_fabqnf(st);

    assert(law.rho == 2);
    assert(law.hidden_partition_evaluations < static_cast<int>(law.cell_weight.size()));
    assert(law.hidden_partition_cache_hits > 0);
    assert(law.hidden_partition_cache_entries == law.hidden_partition_evaluations);

    assert_distributions_near(exact_distribution(st), fabqnf_distribution(st));
}

void test_disconnected_stats_use_component_local_fabqnf() {
    laph::LAPH st(4);
    st.h(0).t(0);
    st.h(1).t(1);
    st.compress();

    laph::Stats stats = st.stats();
    assert(stats.components == 2);
    assert(stats.hidden_interference_rank == -1);
    assert(stats.fabqnf_rho == -1);
    assert(stats.max_component_rho == 0);
    assert(stats.max_component_kappa == 0);
    assert(stats.sum_component_table_size == 2);
    assert(stats.component_table_failures == 0);
}

void test_fabqnf_random_small_circuits_match_amplitudes() {
    for (uint64_t seed = 0; seed < 25; ++seed) {
        std::mt19937_64 rng(seed);
        laph::LAPH st(4);

        for (int layer = 0; layer < 6; ++layer) {
            for (int q = 0; q < 4; ++q) {
                uint64_t g = rng() % 5;
                if (g == 0) st.h(q);
                else if (g == 1) st.s(q);
                else if (g == 2) st.z(q);
                else if (g == 3) st.x(q);
                else st.t(q);
            }

            int c = static_cast<int>(rng() % 4);
            int t = static_cast<int>((c + 1 + (rng() % 3)) % 4);
            st.cnot(c, t);
        }

        st.compress();
        assert_distributions_near(exact_distribution(st), fabqnf_distribution(st), 1e-7L);
    }
}

} // namespace

int main() {
    test_demo_amplitudes();
    test_factorized_matches_monolithic();
    test_clifford_cut_zero();
    test_exact_sample_shape();
    test_component_sampler_keeps_constants();
    test_tableau_bell_sampler();
    test_tableau_deterministic_paulis();
    test_non_clifford_invalidates_tableau();
    test_hidden_interference_rank_distinguishes_output_phase();
    test_fabqnf_output_only_t_has_rho_zero();
    test_fabqnf_hidden_t_matches_amplitudes();
    test_fabqnf_hidden_translation_reduces_rho();
    test_fabqnf_hidden_packet_cache_hits();
    test_disconnected_stats_use_component_local_fabqnf();
    test_fabqnf_random_small_circuits_match_amplitudes();

    std::cout << "all LAPH tests passed\n";
}
