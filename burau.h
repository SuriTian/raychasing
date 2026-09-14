#ifndef BURAU_H
#define BURAU_H

#include <vector>
#include <numeric>
#include <stdexcept>
#include <cstdlib>

#include "laurent.h"
#include "polynomial.h"

// Conway polynomial for braids on any number of strands, via the reduced
// Burau representation.
//
// The chain is: braid word -> product of reduced Burau matrices (entries are
// Laurent polynomials in t) -> det(M - I) -> divide by 1+t+...+t^(n-1) to get
// the Alexander polynomial -> rewrite in z, where z = t^(1/2) - t^(-1/2).
//
// That last step never needs a square root: z^2 = t - 2 + 1/t, so t satisfies
//     t^2 - (z^2 + 2)*t + 1 = 0
// and every power of t reduces to A(z) + B(z)*t. For a knot the Alexander
// polynomial is symmetric under t <-> 1/t, which makes all the B parts cancel;
// we check that they do rather than assuming it.
//
// Scope: the closure must be a knot (one component). Multi-component closures
// on 3+ strands need a t^(1/2) that this integer-exponent representation
// doesn't carry, so they are rejected rather than silently approximated.

using PolyMatrix = std::vector<std::vector<LaurentPoly>>;

inline PolyMatrix identity_matrix(int size) {
    PolyMatrix m(size, std::vector<LaurentPoly>(size, LaurentPoly::zero()));
    for (int i = 0; i < size; i++) m[i][i] = LaurentPoly::one();
    return m;
}

inline PolyMatrix multiply(const PolyMatrix& a, const PolyMatrix& b) {
    int size = static_cast<int>(a.size());
    PolyMatrix r(size, std::vector<LaurentPoly>(size, LaurentPoly::zero()));
    for (int i = 0; i < size; i++)
        for (int k = 0; k < size; k++) {
            if (a[i][k].is_zero()) continue;
            for (int j = 0; j < size; j++) {
                if (b[k][j].is_zero()) continue;
                r[i][j] = r[i][j] + a[i][k] * b[k][j];
            }
        }
    return r;
}

// Cofactor expansion. The matrices here are (strands-1) square and braids in
// this project have few strands, so the factorial cost is not worth avoiding.
inline LaurentPoly determinant(const PolyMatrix& m) {
    int size = static_cast<int>(m.size());
    if (size == 0) return LaurentPoly::one();
    if (size == 1) return m[0][0];
    if (size == 2) return m[0][0] * m[1][1] - m[0][1] * m[1][0];

    LaurentPoly total = LaurentPoly::zero();
    for (int col = 0; col < size; col++) {
        if (m[0][col].is_zero()) continue;
        PolyMatrix minor(size - 1, std::vector<LaurentPoly>(size - 1, LaurentPoly::zero()));
        for (int i = 1; i < size; i++) {
            int target = 0;
            for (int j = 0; j < size; j++) {
                if (j == col) continue;
                minor[i - 1][target++] = m[i][j];
            }
        }
        LaurentPoly term = m[0][col] * determinant(minor);
        total = (col % 2 == 0) ? total + term : total - term;
    }
    return total;
}

// Reduced Burau matrix of a single generator. `generator` is signed: +i is
// sigma_i, -i is its inverse. The matrix is (strands-1) square, indexed so
// that row/column j holds basis vector v_(j+1).
//
// Each generator is the identity apart from one row (two at the ends of the
// braid, where the pattern is truncated). The inverse rows are the solutions
// of M*N = I for those blocks, not a t -> 1/t substitution -- that shortcut
// gives the row reversed and is wrong.
inline PolyMatrix burau_generator(int generator, int strands) {
    int i = std::abs(generator);
    if (i < 1 || i > strands - 1) {
        throw std::logic_error("burau_generator: generator index out of range for strand count");
    }

    const int size = strands - 1;
    PolyMatrix m = identity_matrix(size);
    const LaurentPoly t = LaurentPoly::t();
    const LaurentPoly t_inv = LaurentPoly::monomial(-1, 1);
    const LaurentPoly one = LaurentPoly::one();
    const bool positive = generator > 0;

    if (size == 1) {
        // Two strands: the whole representation is 1x1.
        m[0][0] = positive ? LaurentPoly::monomial(1, -1) : LaurentPoly::monomial(-1, -1);
        return m;
    }

    if (i == 1) {
        // Top-left 2x2 block; only the first row differs from the identity.
        if (positive) {
            m[0][0] = LaurentPoly::monomial(1, -1); // -t
            m[0][1] = one;
        } else {
            m[0][0] = LaurentPoly::monomial(-1, -1); // -1/t
            m[0][1] = t_inv;
        }
        return m;
    }

    if (i == strands - 1) {
        // Bottom-right 2x2 block; only the last row differs.
        if (positive) {
            m[size - 1][size - 2] = t;
            m[size - 1][size - 1] = LaurentPoly::monomial(1, -1); // -t
        } else {
            m[size - 1][size - 2] = one;
            m[size - 1][size - 1] = LaurentPoly::monomial(-1, -1); // -1/t
        }
        return m;
    }

    // Interior generator: a 3x3 block whose middle row is the only change.
    const int row = i - 1;
    if (positive) {
        m[row][row - 1] = t;
        m[row][row]     = LaurentPoly::monomial(1, -1); // -t
        m[row][row + 1] = one;
    } else {
        m[row][row - 1] = one;
        m[row][row]     = LaurentPoly::monomial(-1, -1); // -1/t
        m[row][row + 1] = t_inv;
    }
    return m;
}

inline PolyMatrix burau_matrix(const std::vector<int>& word, int strands) {
    PolyMatrix m = identity_matrix(strands - 1);
    for (int generator : word) {
        m = multiply(m, burau_generator(generator, strands));
    }
    return m;
}

// Number of components of the braid closure: the braid's underlying
// permutation ignores over/under, so each generator is just a transposition,
// and components correspond to the cycles it leaves behind.
inline int closure_component_count(const std::vector<int>& word, int strands) {
    std::vector<int> arrangement(strands);
    std::iota(arrangement.begin(), arrangement.end(), 0);
    for (int generator : word) {
        int i = std::abs(generator) - 1;
        if (i < 0 || i + 1 >= strands) {
            throw std::logic_error("closure_component_count: generator index out of range for strand count");
        }
        std::swap(arrangement[i], arrangement[i + 1]);
    }

    std::vector<bool> seen(strands, false);
    int cycles = 0;
    for (int start = 0; start < strands; start++) {
        if (seen[start]) continue;
        cycles++;
        for (int at = start; !seen[at]; at = arrangement[at]) seen[at] = true;
    }
    return cycles;
}

// Alexander polynomial of the closure, up to the usual +/- t^k ambiguity:
//     Delta(t) * (1 + t + ... + t^(n-1)) = +/- t^k * det(reduced Burau - I)
inline LaurentPoly alexander_polynomial(const std::vector<int>& word, int strands) {
    PolyMatrix minus_identity = burau_matrix(word, strands);
    for (int i = 0; i < strands - 1; i++) {
        minus_identity[i][i] = minus_identity[i][i] - LaurentPoly::one();
    }

    LaurentPoly numerator = determinant(minus_identity);

    std::vector<long long> ones(strands, 1);
    LaurentPoly denominator(0, ones); // 1 + t + ... + t^(n-1)

    return numerator.divide_exact(denominator);
}

// Rewrites an Alexander polynomial as the Conway polynomial in z.
//
// Everything happens in s = t^(1/2), because a link with an even number of
// components genuinely needs half-integer powers of t. Working in s costs
// nothing: exponents in s are just twice the exponents in t, plus one more in
// that even case. Powers of s then collapse through s^2 = z*s + 1.
//
// `unit_sign` settles the sign that the Burau formula leaves undetermined, and
// is ignored for knots -- a knot's Conway polynomial has constant term 1, which
// pins the sign on its own without knowing where `delta` came from.
inline Polynomial conway_from_alexander(const LaurentPoly& delta, int components, int unit_sign) {
    // A split link has no Alexander polynomial to speak of, and a vanishing
    // Conway polynomial is the right answer rather than a failure.
    if (delta.is_zero()) return Polynomial::zero();

    const int extra = (components % 2 == 0) ? 1 : 0;
    const int lo = 2 * delta.lowest_exponent() + extra;
    const int hi = 2 * delta.highest_exponent() + extra;
    auto coefficient_at = [&](int exponent_in_s) -> long long {
        if ((exponent_in_s - extra) % 2 != 0) return 0;
        return delta.coefficient((exponent_in_s - extra) / 2);
    };

    // The normalized polynomial sits symmetrically about s^0, so centring the
    // exponent range is what resolves the +/- t^k the Burau formula leaves open.
    // lo and hi differ by an even number, so the centre is always a whole power.
    const int center = (lo + hi) / 2;
    const int degree = (hi - lo) / 2;

    // Delta(1/s) = (-1)^(components-1) * Delta(s). If that fails, the input was
    // not an Alexander polynomial of a link with this many components.
    const long long mirror = (components % 2 == 1) ? 1 : -1;
    for (int k = 0; k <= degree; k++) {
        if (coefficient_at(center - k) != mirror * coefficient_at(center + k)) {
            throw std::logic_error("conway_from_alexander: polynomial has the wrong symmetry in t");
        }
    }

    // s^k = A_k(z) + B_k(z)*s, walking out from s^0 in both directions:
    //   s^(k+1) = B_k + (A_k + z*B_k)*s      (from s^2 = z*s + 1)
    //   s^(k-1) = (B_k - z*A_k) + A_k*s      (from 1/s = s - z)
    Polynomial conway = Polynomial::zero();
    Polynomial leftover = Polynomial::zero(); // coefficient of the leftover s

    auto accumulate = [&](long long coefficient, const Polynomial& a, const Polynomial& b) {
        if (coefficient == 0) return;
        Polynomial scale(std::vector<long long>{coefficient});
        conway = conway + scale * a;
        leftover = leftover + scale * b;
    };

    Polynomial a_up = Polynomial::one(), b_up = Polynomial::zero();     // s^0
    accumulate(coefficient_at(center), a_up, b_up);
    for (int k = 1; k <= degree; k++) {
        Polynomial a_next = b_up;
        Polynomial b_next = a_up + Polynomial::z() * b_up;
        a_up = a_next;
        b_up = b_next;
        accumulate(coefficient_at(center + k), a_up, b_up);
    }

    Polynomial a_down = Polynomial::one(), b_down = Polynomial::zero(); // s^0
    for (int k = 1; k <= degree; k++) {
        Polynomial a_next = b_down - Polynomial::z() * a_down;
        Polynomial b_next = a_down;
        a_down = a_next;
        b_down = b_next;
        accumulate(coefficient_at(center - k), a_down, b_down);
    }

    if (!leftover.is_zero()) {
        throw std::logic_error("conway_from_alexander: s did not cancel out");
    }

    if (components == 1) {
        long long constant_term = conway.coeffs.empty() ? 0 : conway.coeffs[0];
        if (constant_term == -1) return Polynomial::zero() - conway;
        if (constant_term != 1) {
            throw std::logic_error("conway_from_alexander: a knot must have constant term 1");
        }
        return conway;
    }

    return (unit_sign < 0) ? Polynomial::zero() - conway : conway;
}

// Conway polynomial of the closure of `word`, straight from the braid.
inline Polynomial conway_polynomial_via_burau(const std::vector<int>& word, int strands) {
    const int components = closure_component_count(word, strands);
    // The determinant formula pins the answer down only up to a sign; this is
    // the correction, measured against braids whose closures are already known.
    const int unit_sign = ((strands - components) % 2 == 0) ? 1 : -1;
    return conway_from_alexander(alexander_polynomial(word, strands), components, unit_sign);
}

#endif
