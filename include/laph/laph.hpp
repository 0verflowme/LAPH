#pragma once

#include "laph/clifford_tableau.hpp"
#include "laph/optimizer.hpp"

#include <cstdint>
#include <iosfwd>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

namespace laph {

class LAPH {
public:
    int n = 0;
    int m = 0;
    std::vector<AffineForm> visible;
    std::vector<Constraint> constraints;
    PhasePoly phase;
    int scale = 0;
    CliffordTableau tableau;
    bool tableau_valid = true;
    std::unordered_map<AffineForm, int, AffineHash> lift_cache;

    explicit LAPH(int nqubits);

    int new_var();
    std::pair<bool, int> lift(const AffineForm& f);
    void add_phase_monomial(int coeff, const Monomial& mon);
    void add_phase_product_forms(int coeff, const std::vector<AffineForm>& forms);

    LAPH& x(int q);
    LAPH& cnot(int c, int t);
    LAPH& h(int q);
    LAPH& t(int q);
    LAPH& s(int q);
    LAPH& z(int q);
    LAPH& cz(int a, int b);
    LAPH& ccz(int a, int b, int c);

    std::vector<Constraint> output_constraints(uint64_t xbits) const;
    std::vector<int> cut_set() const;

    ScaledComplex amplitude(uint64_t xbits, QueryOptions options = {}) const;
    ScaledComplex amplitude_monolithic(uint64_t xbits) const;
    ScaledComplex amplitude_factorized(uint64_t xbits) const;

    PhasePoly density_phase() const;
    std::vector<Constraint> density_constraints(
        const std::vector<std::pair<int, bool>>& prefix
    ) const;
    ScaledComplex density_prefix_sum(
        const std::vector<std::pair<int, bool>>& prefix,
        QueryOptions options = {}
    ) const;
    long double density_prefix_weight(
        const std::vector<std::pair<int, bool>>& prefix,
        QueryOptions options = {}
    ) const;
    long double probability_prefix(
        const std::vector<std::pair<int, bool>>& prefix,
        QueryOptions options = {}
    ) const;
    long double probability_prefix_monolithic(
        const std::vector<std::pair<int, bool>>& prefix
    ) const;
    long double probability_prefix_factorized(
        const std::vector<std::pair<int, bool>>& prefix
    ) const;

    std::vector<int> exact_sample(std::mt19937_64& rng, QueryOptions options = {}) const;
    std::vector<int> exact_sample_monolithic(std::mt19937_64& rng) const;
    std::vector<int> exact_sample_factorized(std::mt19937_64& rng) const;

    std::vector<StateComponent> state_components() const;
    LAPH component_view(const StateComponent& component) const;
    std::vector<int> exact_sample_by_component(
        std::mt19937_64& rng,
        QueryOptions options = {}
    ) const;

    void fold_phase_terms();
    void canonicalize_constraints();
    void compress();

    std::vector<std::vector<int>> connected_components() const;
    Stats stats() const;
    void print_stats() const;
};

} // namespace laph
