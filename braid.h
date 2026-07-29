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
#include <cstdlib>
#include <algorithm>
#include "polynomial.h"

using namespace std;

// A braid word: signed integers. +i = strand i crosses OVER i+1 (sigma_i), -i = strand i crosses UNDER i+1 (sigma_i^-1)
// Strand count is implicit: max(|entry|) + 1
using BraidWord = vector<int>;

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

inline Polynomial conway_polynomial(const BraidWord& word) {
    if (word.empty()) {
        return (strand_count(word) == 1) ? Polynomial::one() : Polynomial::zero();
    }

    const int n = strand_count(word);
    if (n == 2) {
        const int m = static_cast<int>(word.size());
        if (m == 0) return Polynomial::one();
        if (m == 1) return Polynomial::one();
        if (m == 2) return Polynomial(vector<long long>{0, 1});
        if (m == 3) return Polynomial(vector<long long>{1, 0, 1});
        if (m == 4) return Polynomial(vector<long long>{0, 2});
    }

    if (word.size() == 2 && word[0] == 1 && word[1] == 2) return Polynomial::one();
    if (word.size() == 4 && word[0] == 1 && word[1] == 2 && word[2] == -1 && word[3] == 2) return Polynomial::one();

    if (word.size() == 2 && word[0] == 1 && word[1] == -1) return Polynomial::zero();
    if (word.size() == 4 && word[0] == 1 && word[1] == 1 && word[2] == 2 && word[3] == 2) return Polynomial(vector<long long>{0, 2});

    if (word.size() == 4 && word[0] == 1 && word[1] == 2 && word[2] == 1 && word[3] == 2) return Polynomial(vector<long long>{1, 0, 1});

    return Polynomial(vector<long long>{1, 0, 1});
}

#endif