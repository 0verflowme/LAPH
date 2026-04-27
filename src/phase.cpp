#include "laph/phase.hpp"

#include <set>
#include <stdexcept>

namespace laph {

void add_mod8(PhasePoly& p, const Monomial& m, int coeff) {
    coeff = mod8(coeff);
    if (coeff == 0) return;
    auto it = p.find(m);
    if (it == p.end()) {
        p.emplace(m, coeff);
    } else {
        int nc = mod8(it->second + coeff);
        if (nc == 0) p.erase(it);
        else it->second = nc;
    }
}

static Monomial mon_from_vars(const std::vector<int>& xs) {
    return Monomial(xs);
}

void add_affine_phase(PhasePoly& poly, int coeff, const XorSet& mask, bool constant) {
    coeff = mod8(coeff);
    if (coeff == 0) return;

    if (constant) {
        add_mod8(poly, Monomial::constant(), coeff);
        coeff = mod8(-coeff);
    }

    const auto& vs = mask.v;
    const int n = static_cast<int>(vs.size());

    int c1 = mod8(coeff);
    if (c1) {
        for (int i = 0; i < n; ++i) {
            add_mod8(poly, Monomial::singleton(vs[i]), c1);
        }
    }

    int c2 = mod8(coeff * -2);
    if (c2) {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                add_mod8(poly, mon_from_vars({vs[i], vs[j]}), c2);
            }
        }
    }

    int c3 = mod8(coeff * 4);
    if (c3) {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = j + 1; k < n; ++k) {
                    add_mod8(poly, mon_from_vars({vs[i], vs[j], vs[k]}), c3);
                }
            }
        }
    }
}


static ANF anf_affine(const AffineForm& a) {
    ANF out;
    if (a.c) out.insert(Monomial::constant());
    for (int x : a.vars.v) out.insert(Monomial::singleton(x));
    return out;
}

static void anf_toggle(ANF& s, const Monomial& m) {
    auto it = s.find(m);
    if (it == s.end()) s.insert(m);
    else s.erase(it);
}

static ANF anf_mul(const ANF& a, const ANF& b) {
    ANF out;
    for (const auto& x : a) {
        for (const auto& y : b) {
            anf_toggle(out, Monomial::boolean_mul(x, y));
        }
    }
    return out;
}

bool is_clifford_poly(const PhasePoly& poly) {
    for (const auto& kv : poly) {
        const Monomial& m = kv.first;
        int c = mod8(kv.second);
        if (c == 0) continue;
        size_t deg = m.degree();
        if (deg == 0) continue;
        if (deg == 1) {
            if (c % 2 != 0) return false;
        } else if (deg == 2) {
            if (c % 4 != 0) return false;
        } else {
            return false;
        }
    }
    return true;
}

static PhasePoly substitute_pivot_in_clifford(const PhasePoly& poly, int pivot, bool constant, const XorSet& mask) {
    PhasePoly out;
    for (const auto& kv : poly) {
        const Monomial& m = kv.first;
        int c = mod8(kv.second);
        if (!m.contains(pivot)) {
            add_mod8(out, m, c);
            continue;
        }
        Monomial rest = m.without(pivot);
        if (rest.empty()) {
            add_affine_phase(out, c, mask, constant);
        } else {
            if (rest.degree() != 1 || c % 4 != 0) {
                throw std::runtime_error("bad Clifford pivot substitution");
            }
            int j = rest.v[0];
            AffineForm a(mask, constant);
            ANF prod = anf_mul(anf_affine(a), ANF{Monomial::singleton(j)});
            for (const auto& mon : prod) add_mod8(out, mon, c);
        }
    }
    return out;
}

ScaledComplex clifford_sum(PhasePoly poly, int nvars) {
    if (!is_clifford_poly(poly)) throw std::runtime_error("not Clifford");

    int phase = 0;
    auto it0 = poly.find(Monomial::constant());
    if (it0 != poly.end()) {
        phase = mod8(it0->second);
        poly.erase(it0);
    }

    std::set<int> active;
    for (int i = 0; i < nvars; ++i) active.insert(i);

    int64_t log_sqrt2 = 0;

    while (!active.empty()) {
        int chosen = -1;

        for (int v : active) {
            auto it = poly.find(Monomial::singleton(v));
            int c = (it == poly.end()) ? 0 : mod8(it->second);
            if (c % 2 == 0 && ((c / 2) % 2 == 1)) {
                chosen = v;
                break;
            }
        }

        if (chosen == -1) {
            for (int v : active) {
                bool has_quad = false;
                for (const auto& kv : poly) {
                    if (kv.first.degree() == 2 && kv.first.contains(v) && mod8(kv.second) != 0) {
                        has_quad = true;
                        break;
                    }
                }
                if (has_quad) { chosen = v; break; }
            }
        }

        if (chosen == -1) chosen = *active.begin();

        int v = chosen;
        auto lin_it = poly.find(Monomial::singleton(v));
        int lin = (lin_it == poly.end()) ? 0 : mod8(lin_it->second);
        int l = (lin / 2) % 4;

        XorSet neigh;
        for (const auto& kv : poly) {
            const Monomial& m = kv.first;
            int c = mod8(kv.second);
            if (m.degree() == 2 && m.contains(v) && c != 0) {
                if (c % 4 != 0) throw std::runtime_error("non-Clifford quadratic");
                if (c % 8 != 0) {
                    int other = (m.v[0] == v) ? m.v[1] : m.v[0];
                    neigh.xor_with(XorSet::singleton(other));
                }
            }
        }

        PhasePoly p0;
        for (const auto& kv : poly) {
            if (!kv.first.contains(v)) add_mod8(p0, kv.first, kv.second);
        }

        active.erase(v);

        if (l == 1 || l == 3) {
            log_sqrt2 += 1;
            phase = mod8(phase + (l == 1 ? 1 : -1));
            add_affine_phase(p0, (l == 1 ? 6 : 2), neigh, false);
            poly.swap(p0);
        } else {
            bool required = (l == 2);
            if (neigh.empty()) {
                if (required) return ScaledComplex::zero();
                log_sqrt2 += 2;
                poly.swap(p0);
            } else {
                int pivot = neigh.first();
                XorSet rest = neigh.without(pivot);
                poly = substitute_pivot_in_clifford(p0, pivot, required, rest);
                active.erase(pivot);
                log_sqrt2 += 2;
            }
        }

        std::vector<Monomial> kill;
        for (const auto& kv : poly) if (mod8(kv.second) == 0) kill.push_back(kv.first);
        for (const auto& m : kill) poly.erase(m);
    }

    auto it = poly.find(Monomial::constant());
    if (it != poly.end()) phase = mod8(phase + it->second);
    return ScaledComplex::omega_phase(phase, log_sqrt2);
}

PhasePoly substitute_clifford_poly(const PhasePoly& poly, const std::unordered_map<int, AffineForm>& exprs) {
    PhasePoly out;
    for (const auto& kv : poly) {
        const Monomial& m = kv.first;
        int c = mod8(kv.second);
        if (c == 0) continue;
        size_t deg = m.degree();
        if (deg == 0) {
            add_mod8(out, Monomial::constant(), c);
        } else if (deg == 1) {
            auto it = exprs.find(m.v[0]);
            if (it == exprs.end()) throw std::runtime_error("missing expression");
            add_affine_phase(out, c, it->second.vars, it->second.c);
        } else if (deg == 2 && c % 4 == 0) {
            auto it0 = exprs.find(m.v[0]);
            auto it1 = exprs.find(m.v[1]);
            if (it0 == exprs.end() || it1 == exprs.end()) throw std::runtime_error("missing expression");
            ANF prod = anf_mul(anf_affine(it0->second), anf_affine(it1->second));
            for (const auto& mon : prod) add_mod8(out, mon, c);
        } else {
            throw std::runtime_error("non-Clifford residual term after cut");
        }
    }
    return out;
}

PhasePoly substitute_clifford_poly_vec(
    const PhasePoly& poly,
    const std::vector<AffineForm>& exprs
) {
    PhasePoly out;
    for (const auto& kv : poly) {
        const Monomial& m = kv.first;
        int c = mod8(kv.second);
        if (c == 0) continue;
        size_t deg = m.degree();
        if (deg == 0) {
            add_mod8(out, Monomial::constant(), c);
        } else if (deg == 1) {
            int x = m.v[0];
            if (x < 0 || x >= static_cast<int>(exprs.size())) {
                throw std::runtime_error("missing expression");
            }
            add_affine_phase(out, c, exprs[x].vars, exprs[x].c);
        } else if (deg == 2 && c % 4 == 0) {
            int x = m.v[0];
            int y = m.v[1];
            if (x < 0 || y < 0 ||
                x >= static_cast<int>(exprs.size()) ||
                y >= static_cast<int>(exprs.size())) {
                throw std::runtime_error("missing expression");
            }
            ANF prod = anf_mul(anf_affine(exprs[x]), anf_affine(exprs[y]));
            for (const auto& mon : prod) add_mod8(out, mon, c);
        } else {
            throw std::runtime_error("non-Clifford residual term after cut");
        }
    }
    return out;
}

} // namespace laph
