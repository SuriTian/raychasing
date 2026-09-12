#ifndef BRAID_H
#define BRAID_H

#include <vector>
#include <optional>
#include <map>
#include <string>
#include <sstream>
#include <numeric>
#include <cmath>
#include <utility>
#include <stdexcept>
#include "polynomial.h"

using namespace std;

// A braid word: signed integers. +i = strand i crosses OVER i+1 (sigma_i), -i = strand i crosses UNDER i+1 (sigma_i^-1)
// Strand count is implicit: max(|entry|) + 1
using BraidWord = vector<int>;

// Number of strands = highest strand index touched by the word, plus one.
inline int strand_count(const BraidWord& w) {
    int m = 0;
    for (int x : w) m = max(m, abs(x));
    return m + 1;
}

// Generates the standard torus knot braid word for T(p,q):
// (sigma_1 sigma_2 ... sigma_(p-1)) repeated q times, on p strands.
inline BraidWord generate_torus_braid(int p, int q) {
    BraidWord w;
    for (int rep = 0; rep < q; rep++)
        for (int i = 1; i <= p - 1; i++)
            w.push_back(i);
    return w;
}

// Braid closures are invariant under cyclic rotation of the word.
inline bool is_cyclic_equivalent(const BraidWord& a, const BraidWord& b) {
    if (a.size() != b.size()) return false;
    size_t n = a.size();
    if (n == 0) return true;
    for (size_t shift = 0; shift < n; shift++) {
        bool match = true;
        for (size_t i = 0; i < n; i++) {
            if (a[i] != b[(i + shift) % n]) { match = false; break; }
        }
        if (match) return true;
    }
    return false;
}

// NOTE: this only recognizes trefoil-*shaped* words (3 crossings, uniform sign),
// it does not verify the closure is actually a knot (1 component). E.g. [1,2,1]
// has this shape but its closure is a 2-component Hopf link after Markov
// destabilization, not the trefoil. Used only for renderability/torus-knot
// hints elsewhere -- conway_polynomial() below does not trust it.
inline bool is_trefoil_like(const BraidWord& word) {
    if (word.empty()) return false;

    if (word.size() == 3 && strand_count(word) == 2) {
        bool all_positive = true;
        bool all_negative = true;
        for (int x : word) {
            if (x <= 0) all_positive = false;
            if (x >= 0) all_negative = false;
        }
        if (all_positive || all_negative) {
            return true;
        }
    }

    if (word.size() == 3 && strand_count(word) == 3) {
        const vector<BraidWord> candidates = {
            {1, 2, 1},
            {1, -2, 1},
            {-1, -2, -1},
            {-1, 2, -1}
        };
        for (const auto& candidate : candidates) {
            if (is_cyclic_equivalent(word, candidate)) {
                return true;
            }
        }
    }

    return false;
}

// Checks whether `word` matches a (p, q) torus knot braid pattern
// (up to cyclic rotation). Returns {p, q} if found.
inline optional<pair<int,int>> match_torus_knot(const BraidWord& word, int max_p = 10, int max_q = 10) {
    if (word.empty()) return nullopt;
    int p = strand_count(word);
    if (p < 2 || p > max_p) return nullopt;

    for (int q = 1; q <= max_q; q++) {
        if (gcd(p, q) != 1) continue;
        BraidWord candidate = generate_torus_braid(p, q);
        if (is_cyclic_equivalent(word, candidate)) return make_pair(p, q);
    }

    if (is_trefoil_like(word)) {
        return make_pair(2, 3);
    }

    return nullopt;
}

// On 2 strands there's only one generator (sigma_1), so a word is fully
// described by its net exponent e = (#positive crossings) - (#negative
// crossings) -- order doesn't matter. V_e follows a simple recursion:
//   V_0 = 0 (unlink), V_1 = 1 (unknot), V_e = z*V_{e-1} + V_{e-2}
// Negative e (more under-crossings than over) just runs that recursion
// backwards: V_k = V_{k+2} - z*V_{k+1}.
inline Polynomial conway_polynomial_two_strand(int e) {
    vector<Polynomial> v;
    v.push_back(Polynomial::zero()); // V_0
    v.push_back(Polynomial::one());  // V_1

    if (e >= 0) {
        for (int k = 2; k <= e; k++) {
            v.push_back(Polynomial::z() * v[k - 1] + v[k - 2]);
        }
        return v[e];
    }

    // same recursion, run backwards from k=-1 down to e
    Polynomial v_next2 = v[1]; // V_1
    Polynomial v_next1 = v[0]; // V_0
    Polynomial result = Polynomial::zero();
    for (int k = -1; k >= e; k--) {
        result = v_next2 - Polynomial::z() * v_next1;
        v_next2 = v_next1;
        v_next1 = result;
    }
    return result;
}

// Only 2-strand braids are supported right now -- see the note above
// conway_polynomial_two_strand(). Words on 3+ strands need the
// Burau/Alexander-polynomial machinery, which isn't implemented yet
// (a naive attempt at this lived here before and was silently wrong for
// every multi-strand input, so it was removed rather than patched).
inline Polynomial conway_polynomial(const BraidWord& word) {
    if (word.empty()) {
        return (strand_count(word) == 1) ? Polynomial::one() : Polynomial::zero();
    }

    if (strand_count(word) != 2) {
        throw std::logic_error("conway_polynomial: only 2-strand braids are currently supported");
    }

    int e = 0;
    for (int x : word) e += (x > 0) ? 1 : -1;
    return conway_polynomial_two_strand(e);
}

#endif
