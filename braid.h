#ifndef BRAID_H
#define BRAID_H

#include <vector>
#include <optional>
#include <string>
#include <sstream>
#include <numeric>
#include <cmath>
#include <utility>
#include <stdexcept>
#include "polynomial.h"
#include "burau.h"

using namespace std;

// A braid word: signed integers. +i = strand i crosses OVER i+1 (sigma_i), -i = strand i crosses UNDER i+1 (sigma_i^-1)
// Strand count is implicit: max(|entry|) + 1
using BraidWord = vector<int>;

// Parses a comma-separated list like "1,-2,1" into a BraidWord. Throws
// invalid_argument on anything that isn't a list of nonzero whole numbers --
// 0 is not a generator, and letting it through indexes off the end of the
// strand array further down.
inline BraidWord parse_braid_word(const string& text) {
    BraidWord word;
    stringstream ss(text);
    string token;
    while (getline(ss, token, ',')) {
        const size_t first = token.find_first_not_of(" \t");
        if (first == string::npos) continue;
        token = token.substr(first, token.find_last_not_of(" \t") - first + 1);

        size_t consumed = 0;
        int value = 0;
        try {
            value = stoi(token, &consumed);
        } catch (const exception&) {
            throw invalid_argument("'" + token + "' is not a whole number");
        }
        if (consumed != token.size()) {
            throw invalid_argument("'" + token + "' is not a whole number");
        }
        if (value == 0) {
            throw invalid_argument("0 is not a braid generator -- use 1, -1, 2, -2, ...");
        }
        word.push_back(value);
    }
    return word;
}

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

// Same as the burau.h version, with the strand count filled in from the word.
inline int closure_component_count(const BraidWord& w) {
    return closure_component_count(w, strand_count(w));
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

// Two strands keep the closed form above -- it is cheap, and it doubles as an
// independent check on the matrices. Everything wider goes through the Burau
// representation in burau.h, knots and links alike.
inline Polynomial conway_polynomial(const BraidWord& word) {
    if (word.empty()) {
        return (strand_count(word) == 1) ? Polynomial::one() : Polynomial::zero();
    }

    const int strands = strand_count(word);

    if (strands == 2) {
        int e = 0;
        for (int x : word) e += (x > 0) ? 1 : -1;
        return conway_polynomial_two_strand(e);
    }

    return conway_polynomial_via_burau(word, strands);
}

// True when the closure of `word` is, as far as its Conway polynomial can
// tell, the trefoil. This asks an invariant of the closed-up knot rather than
// pattern-matching the word, so [1,2,1] no longer counts: it has trefoil shape
// but closes to a 2-component link. Two caveats worth knowing:
//   - the Conway polynomial does not separate every knot, so an exotic knot
//     sharing z^2 + 1 would slip through;
//   - it is blind to mirror images, so both handednesses pass -- which suits
//     the renderer, since it draws one fixed chirality either way.
inline bool is_trefoil_like(const BraidWord& word) {
    if (word.empty()) return false;
    try {
        if (closure_component_count(word) != 1) return false;
        return conway_polynomial(word).coeffs == vector<long long>{1, 0, 1}; // z^2 + 1
    } catch (const exception&) {
        // A word we cannot evaluate is not one we can claim is a trefoil, and a
        // predicate the renderer calls should answer rather than throw.
        return false;
    }
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

#endif
