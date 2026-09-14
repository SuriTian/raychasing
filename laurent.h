#ifndef LAURENT_H
#define LAURENT_H

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include "checked_int.h"

// Laurent polynomial in one variable t with integer coefficients.
// Burau matrix entries need negative powers of t (an inverse generator
// contributes 1/t), so a plain coefficient vector isn't enough.
// coeffs[i] is the coefficient of t^(min_exp + i).
class LaurentPoly {
public:
    LaurentPoly() : min_exp(0), coeffs{0} {}
    explicit LaurentPoly(long long constant) : min_exp(0), coeffs{constant} { trim(); }
    LaurentPoly(int exponent, long long coefficient) : min_exp(exponent), coeffs{coefficient} { trim(); }
    LaurentPoly(int lowest_exponent, std::vector<long long> c)
        : min_exp(lowest_exponent), coeffs(std::move(c)) { trim(); }

    static LaurentPoly zero() { return LaurentPoly(0LL); }
    static LaurentPoly one()  { return LaurentPoly(1LL); }
    static LaurentPoly t()    { return LaurentPoly(1, 1); }
    static LaurentPoly monomial(int exponent, long long coefficient) {
        return LaurentPoly(exponent, coefficient);
    }

    bool is_zero() const { return coeffs.size() == 1 && coeffs[0] == 0; }
    int lowest_exponent()  const { return min_exp; }
    int highest_exponent() const { return min_exp + static_cast<int>(coeffs.size()) - 1; }

    long long coefficient(int exponent) const {
        int i = exponent - min_exp;
        if (i < 0 || i >= static_cast<int>(coeffs.size())) return 0;
        return coeffs[i];
    }

    LaurentPoly operator+(const LaurentPoly& o) const { return combine(o, 1); }
    LaurentPoly operator-(const LaurentPoly& o) const { return combine(o, -1); }

    LaurentPoly operator*(const LaurentPoly& o) const {
        if (is_zero() || o.is_zero()) return zero();
        std::vector<long long> r(coeffs.size() + o.coeffs.size() - 1, 0);
        for (std::size_t i = 0; i < coeffs.size(); i++)
            for (std::size_t j = 0; j < o.coeffs.size(); j++)
                r[i + j] = checked_add(r[i + j], checked_mul(coeffs[i], o.coeffs[j]));
        return LaurentPoly(min_exp + o.min_exp, std::move(r));
    }

    bool operator==(const LaurentPoly& o) const {
        return min_exp == o.min_exp && coeffs == o.coeffs;
    }

    // Exact division. Throws if `divisor` does not divide this evenly --
    // every caller here divides by a factor that is mathematically
    // guaranteed to come out even, so a remainder means a bug upstream.
    LaurentPoly divide_exact(const LaurentPoly& divisor) const {
        if (divisor.is_zero()) throw std::logic_error("LaurentPoly: division by zero");
        if (is_zero()) return zero();

        std::vector<long long> num = coeffs;
        const std::vector<long long>& den = divisor.coeffs;
        int quotient_degree = static_cast<int>(num.size()) - static_cast<int>(den.size());
        if (quotient_degree < 0) throw std::logic_error("LaurentPoly: inexact division");

        std::vector<long long> q(quotient_degree + 1, 0);
        for (int i = quotient_degree; i >= 0; i--) {
            long long lead = num[i + den.size() - 1];
            if (lead % den.back() != 0) throw std::logic_error("LaurentPoly: inexact division");
            long long f = lead / den.back();
            q[i] = f;
            if (f == 0) continue;
            for (std::size_t j = 0; j < den.size(); j++) {
                num[i + j] = checked_add(num[i + j], -checked_mul(f, den[j]));
            }
        }
        for (long long r : num) {
            if (r != 0) throw std::logic_error("LaurentPoly: inexact division");
        }
        return LaurentPoly(min_exp - divisor.min_exp, std::move(q));
    }

    std::string to_string(const std::string& var = "t") const {
        std::ostringstream out;
        bool first = true;
        for (int e = highest_exponent(); e >= min_exp; e--) {
            long long c = coefficient(e);
            if (c == 0) continue;
            if (!first) out << (c > 0 ? " + " : " - ");
            else if (c < 0) out << "-";
            long long mag = std::llabs(c);
            if (e == 0) out << mag;
            else {
                if (mag != 1) out << mag << "*";
                out << var;
                if (e != 1) out << "^" << e;
            }
            first = false;
        }
        if (first) out << "0";
        return out.str();
    }

private:
    int min_exp;
    std::vector<long long> coeffs;

    void trim() {
        while (coeffs.size() > 1 && coeffs.back() == 0) coeffs.pop_back();
        while (coeffs.size() > 1 && coeffs.front() == 0) {
            coeffs.erase(coeffs.begin());
            min_exp++;
        }
        if (coeffs.empty()) { coeffs.push_back(0); min_exp = 0; }
        if (coeffs.size() == 1 && coeffs[0] == 0) min_exp = 0;
    }

    LaurentPoly combine(const LaurentPoly& o, long long sign) const {
        int lo = std::min(min_exp, o.min_exp);
        int hi = std::max(highest_exponent(), o.highest_exponent());
        std::vector<long long> r(hi - lo + 1, 0);
        for (int e = min_exp; e <= highest_exponent(); e++) {
            r[e - lo] = checked_add(r[e - lo], coefficient(e));
        }
        for (int e = o.min_exp; e <= o.highest_exponent(); e++) {
            r[e - lo] = checked_add(r[e - lo], checked_mul(sign, o.coefficient(e)));
        }
        return LaurentPoly(lo, std::move(r));
    }
};

#endif
