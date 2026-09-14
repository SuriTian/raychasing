#ifndef CHECKED_INT_H
#define CHECKED_INT_H

#include <limits>
#include <stdexcept>

// Polynomial coefficients here grow fast -- a torus knot around T(10,13)
// already runs past 64 bits. Wrapping around would hand back a polynomial that
// looks plausible and is wrong, so these throw instead. The comparisons happen
// before the arithmetic, since signed overflow itself is undefined behaviour
// and cannot be detected after the fact.

inline long long checked_add(long long a, long long b) {
    if (b > 0 && a > std::numeric_limits<long long>::max() - b) {
        throw std::overflow_error("polynomial coefficient is too large for 64-bit integers");
    }
    if (b < 0 && a < std::numeric_limits<long long>::min() - b) {
        throw std::overflow_error("polynomial coefficient is too large for 64-bit integers");
    }
    return a + b;
}

inline long long checked_mul(long long a, long long b) {
    if (a == 0 || b == 0) return 0;

    const long long max = std::numeric_limits<long long>::max();
    const long long min = std::numeric_limits<long long>::min();
    bool overflows = false;
    if (a > 0) {
        overflows = (b > 0) ? (a > max / b) : (b < min / a);
    } else {
        overflows = (b > 0) ? (a < min / b) : (a < max / b);
    }
    if (overflows) {
        throw std::overflow_error("polynomial coefficient is too large for 64-bit integers");
    }
    return a * b;
}

#endif
