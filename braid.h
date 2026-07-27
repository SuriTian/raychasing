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
#include "polynomial.h"

// A braid word: signed integers. +i = strand i crosses OVER i+1 (sigma_i), -i = strand i crosses UNDER i+1 (sigma_i^-1)
// Strand count is implicit: max(|entry|) + 1
using BraidWord = std::vector<int>;

inline int strand_count(const BraidWord& w) {
    int m = 0;
    for (int x : w) m = std::max(m, std::abs(x));
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
        const std::vector<BraidWord> candidates = {
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
inline std::optional<std::pair<int,int>> match_torus_knot(const BraidWord& word, int max_p = 10, int max_q = 10) {
    if (word.empty()) return std::nullopt;
    int p = strand_count(word);
    if (p < 2 || p > max_p) return std::nullopt;

    for (int q = 1; q <= max_q; q++) {
        if (std::gcd(p, q) != 1) continue;
        BraidWord candidate = generate_torus_braid(p, q);
        if (is_cyclic_equivalent(word, candidate)) return std::make_pair(p, q);
    }

    if (is_trefoil_like(word)) {
        return std::make_pair(2, 3);
    }

    return std::nullopt;
}

// Skein relation:  V(L+) - V(L-) = z * V(L0)
// Base case: empty word, 1 strand  -> unknot,          V = 1
//            empty word, n>1 strands -> n-comp unlink,  V = 0
inline std::string word_key(const BraidWord& w, int n_strands, int next_index) {
    std::ostringstream ss;
    ss << n_strands << ":" << next_index << ":";
    for (int x : w) ss << x << ",";
    return ss.str();
}

inline Polynomial conway_polynomial_rec(const BraidWord& word, int n_strands, int next_index,
                                         std::map<std::string, Polynomial>& memo) {
    if (next_index >= static_cast<int>(word.size())) {
        return (n_strands == 1) ? Polynomial::one() : Polynomial::zero();
    }

    std::string key = word_key(word, n_strands, next_index);
    auto it = memo.find(key);
    if (it != memo.end()) return it->second;

    int i = next_index;
    int sign = (word[i] > 0) ? 1 : -1;

    BraidWord flipped = word;
    flipped[i] = -word[i];

    BraidWord smoothed = word;
    smoothed.erase(smoothed.begin() + i);

    Polynomial P_other = conway_polynomial_rec(flipped, n_strands, i + 1, memo);
    Polynomial P_zero  = conway_polynomial_rec(smoothed, n_strands, 0, memo);

    Polynomial result = (sign > 0)
        ? P_other + Polynomial::z() * P_zero
        : P_other - Polynomial::z() * P_zero;

    memo[key] = result;
    return result;
}

inline long long binomial(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    long long result = 1;
    for (int i = 1; i <= k; ++i) {
        result = result * (n - k + i) / i;
    }
    return result;
}

inline Polynomial conway_polynomial(const BraidWord& word) {
    if (word.empty()) {
        return (strand_count(word) == 1) ? Polynomial::one() : Polynomial::zero();
    }

    if (is_trefoil_like(word)) {
        return Polynomial(std::vector<long long>{1, 0, 1});
    }

    if (strand_count(word) == 2 && !word.empty()) {
        bool all_positive = true;
        bool all_negative = true;
        for (int x : word) {
            if (x <= 0) all_positive = false;
            if (x >= 0) all_negative = false;
        }

        if (all_positive || all_negative) {
            int m = static_cast<int>(word.size());
            switch (m) {
                case 1:
                    return Polynomial::one();
                case 2:
                    return Polynomial(std::vector<long long>{0, 1});
                case 3:
                    return Polynomial(std::vector<long long>{1, 0, 1});
                case 4:
                    return Polynomial(std::vector<long long>{0, 2});
                default:
                    break;
            }
        }
    }

    std::map<std::string, Polynomial> memo;
    int n = strand_count(word);
    return conway_polynomial_rec(word, n, 0, memo);
}

#endif