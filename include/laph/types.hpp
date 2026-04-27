#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace laph {

using std::int64_t;
static inline int mod8(int x) {
    x %= 8;
    if (x < 0) x += 8;
    return x;
}

struct XorSet {
    std::vector<int> v;

    XorSet() = default;
    explicit XorSet(std::vector<int> xs) : v(std::move(xs)) { normalize_xor(); }

    static XorSet singleton(int x) { return XorSet({x}); }

    void normalize_xor() {
        std::sort(v.begin(), v.end());
        std::vector<int> out;
        for (size_t i = 0; i < v.size();) {
            size_t j = i + 1;
            while (j < v.size() && v[j] == v[i]) ++j;
            if (((j - i) & 1u) != 0) out.push_back(v[i]);
            i = j;
        }
        v.swap(out);
    }

    bool empty() const { return v.empty(); }
    size_t size() const { return v.size(); }
    int first() const { return v.front(); }

    bool contains(int x) const {
        return std::binary_search(v.begin(), v.end(), x);
    }

    void xor_with(const XorSet& b) {
        std::vector<int> out;
        out.reserve(v.size() + b.v.size());
        size_t i = 0, j = 0;
        while (i < v.size() || j < b.v.size()) {
            if (j == b.v.size() || (i < v.size() && v[i] < b.v[j])) {
                out.push_back(v[i++]);
            } else if (i == v.size() || b.v[j] < v[i]) {
                out.push_back(b.v[j++]);
            } else {
                ++i;
                ++j;
            }
        }
        v.swap(out);
    }

    XorSet xored(const XorSet& b) const {
        XorSet r = *this;
        r.xor_with(b);
        return r;
    }

    XorSet shifted(int delta) const {
        XorSet r;
        r.v.reserve(v.size());
        for (int x : v) r.v.push_back(x + delta);
        return r;
    }

    XorSet without(int x) const {
        XorSet r;
        r.v.reserve(v.size());
        for (int y : v) if (y != x) r.v.push_back(y);
        return r;
    }

    bool operator==(const XorSet& o) const { return v == o.v; }
};

struct XorHash {
    size_t operator()(const XorSet& s) const noexcept {
        uint64_t h = 1469598103934665603ull;
        for (int x : s.v) {
            uint64_t y = static_cast<uint64_t>(x) + 0x9e3779b97f4a7c15ull;
            h ^= y;
            h *= 1099511628211ull;
        }
        return static_cast<size_t>(h);
    }
};

struct AffineForm {
    XorSet vars;
    bool c = false;

    AffineForm() = default;
    AffineForm(XorSet xs, bool cc=false) : vars(std::move(xs)), c(cc) {}
    static AffineForm constant(bool b) { return AffineForm({}, b); }
    static AffineForm variable(int x) { return AffineForm(XorSet::singleton(x), false); }

    void xor_with(const AffineForm& b) {
        vars.xor_with(b.vars);
        c ^= b.c;
    }

    AffineForm xored(const AffineForm& b) const {
        AffineForm r = *this;
        r.xor_with(b);
        return r;
    }

    AffineForm shifted(int delta) const {
        return AffineForm(vars.shifted(delta), c);
    }

    bool operator==(const AffineForm& o) const { return c == o.c && vars == o.vars; }
};

struct AffineHash {
    size_t operator()(const AffineForm& a) const noexcept {
        return XorHash{}(a.vars) ^ (a.c ? 0x517cc1b727220a95ull : 0ull);
    }
};

struct Monomial {
    std::vector<int> v;

    Monomial() = default;
    explicit Monomial(std::vector<int> xs) : v(std::move(xs)) { normalize_or(); }

    static Monomial constant() { return Monomial(); }
    static Monomial singleton(int x) { return Monomial({x}); }

    void normalize_or() {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    size_t degree() const { return v.size(); }
    bool empty() const { return v.empty(); }

    bool contains(int x) const {
        return std::binary_search(v.begin(), v.end(), x);
    }

    Monomial without(int x) const {
        Monomial r;
        r.v.reserve(v.size());
        for (int y : v) if (y != x) r.v.push_back(y);
        return r;
    }

    Monomial shifted(int delta) const {
        Monomial r;
        r.v.reserve(v.size());
        for (int x : v) r.v.push_back(x + delta);
        return r;
    }

    void or_with_var(int x) {
        if (!contains(x)) {
            v.insert(std::lower_bound(v.begin(), v.end(), x), x);
        }
    }

    static Monomial boolean_mul(const Monomial& a, const Monomial& b) {
        Monomial r;
        r.v.reserve(a.v.size() + b.v.size());
        std::set_union(a.v.begin(), a.v.end(), b.v.begin(), b.v.end(), std::back_inserter(r.v));
        return r;
    }

    bool operator==(const Monomial& o) const { return v == o.v; }
};

struct MonomialHash {
    size_t operator()(const Monomial& m) const noexcept {
        uint64_t h = 1469598103934665603ull;
        for (int x : m.v) {
            uint64_t y = static_cast<uint64_t>(x) + 0x9e3779b97f4a7c15ull;
            h ^= y;
            h *= 1099511628211ull;
        }
        return static_cast<size_t>(h);
    }
};

using PhasePoly = std::unordered_map<Monomial, int, MonomialHash>;

struct Constraint {
    XorSet lhs;
    bool rhs = false;
};

struct ScaledComplex {
    // value = (re + i*im) * (sqrt(2)^k)
    long double re = 0.0L;
    long double im = 0.0L;
    int64_t k = 0;

    static constexpr long double SQRT2 = 1.41421356237309504880168872420969807857L;
    static constexpr long double INV_SQRT2 = 0.70710678118654752440084436210484903928L;

    static ScaledComplex zero() { return {}; }

    static ScaledComplex omega_phase(int phase, int64_t sqrt2_power = 0) {
        phase = mod8(phase);
        switch (phase) {
            case 0: return {1.0L, 0.0L, sqrt2_power};
            case 1: return {INV_SQRT2, INV_SQRT2, sqrt2_power};
            case 2: return {0.0L, 1.0L, sqrt2_power};
            case 3: return {-INV_SQRT2, INV_SQRT2, sqrt2_power};
            case 4: return {-1.0L, 0.0L, sqrt2_power};
            case 5: return {-INV_SQRT2, -INV_SQRT2, sqrt2_power};
            case 6: return {0.0L, -1.0L, sqrt2_power};
            case 7: return {INV_SQRT2, -INV_SQRT2, sqrt2_power};
        }
        return zero();
    }

    bool is_zero() const { return re == 0.0L && im == 0.0L; }

    static long double pow_sqrt2(int64_t d) {
        if (d > 4096) return std::numeric_limits<long double>::infinity();
        if (d < -4096) return 0.0L;
        return std::powl(SQRT2, static_cast<long double>(d));
    }

    void normalize() {
        if (is_zero()) { k = 0; return; }
        long double mag2 = re * re + im * im;
        if (mag2 == 0.0L) { re = im = 0.0L; k = 0; return; }
        long double log2mag = 0.5L * std::log2(mag2);
        int64_t shift = static_cast<int64_t>(std::floor(2.0L * log2mag));
        if (shift != 0) {
            long double f = pow_sqrt2(-shift);
            re *= f;
            im *= f;
            k += shift;
        }
    }

    void add(const ScaledComplex& b) {
        if (b.is_zero()) return;
        if (is_zero()) { *this = b; return; }
        int64_t K = std::max(k, b.k);
        long double fa = pow_sqrt2(k - K);
        long double fb = pow_sqrt2(b.k - K);
        re = re * fa + b.re * fb;
        im = im * fa + b.im * fb;
        k = K;
        normalize();
    }

    void mul(const ScaledComplex& b) {
        if (is_zero() || b.is_zero()) {
            *this = zero();
            return;
        }
        long double nr = re * b.re - im * b.im;
        long double ni = re * b.im + im * b.re;
        re = nr;
        im = ni;
        k += b.k;
        normalize();
    }

    void mul_sqrt2_power(int64_t d) {
        if (!is_zero()) k += d;
    }

    long double real_ld() const {
        return re * pow_sqrt2(k);
    }

    long double imag_ld() const {
        return im * pow_sqrt2(k);
    }
};

} // namespace laph
