#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "laph/laph.hpp"

int main() {
    laph::LAPH st(3);
    st.h(0).h(1).cnot(0, 2).t(2).cz(0, 1).h(2).ccz(0, 1, 2);
    st.compress();

    st.print_stats();

    std::cout << std::setprecision(12);
    for (uint64_t x = 0; x < 8; ++x) {
        laph::ScaledComplex a = st.amplitude_factorized(x);
        long double re = a.real_ld();
        long double im = a.imag_ld();
        if (std::fabs(static_cast<double>(re)) > 1e-10 ||
            std::fabs(static_cast<double>(im)) > 1e-10) {
            std::cout << "|";
            for (int q = 2; q >= 0; --q) std::cout << ((x >> q) & 1ull);
            std::cout << "> -> " << static_cast<double>(re)
                      << " + " << static_cast<double>(im) << "i\n";
        }
    }

    std::vector<std::pair<int, bool>> none;
    std::cout << "norm_prob=" << static_cast<double>(st.probability_prefix_factorized(none)) << "\n";
    for (int q = 0; q < 3; ++q) {
        std::cout << "P(q" << q << "=0)="
                  << static_cast<double>(st.probability_prefix_factorized({{q, false}})) << "\n";
    }

    std::mt19937_64 rng(123);
    auto sample = st.exact_sample_by_component(rng);
    std::cout << "sample=";
    for (int q = 2; q >= 0; --q) std::cout << sample[q];
    std::cout << "\n";
}
