#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "laph/laph.hpp"

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

    std::cout << "all LAPH tests passed\n";
}
