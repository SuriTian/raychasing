#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "checked_int.h"

// Simple dense polynomial in one variable (z), integer coefficients.
// coeffs[i] = coefficient of z^i
class Polynomial {
public:
    std::vector<long long> coeffs;

    Polynomial() : coeffs{0} {}
    Polynomial(long long constant) : coeffs{constant} {}
    Polynomial(std::vector<long long> c) : coeffs(std::move(c)) { trim(); }

    static Polynomial zero() { return Polynomial(0); }
    static Polynomial one()  { return Polynomial(1); }
    static Polynomial z()    { return Polynomial(std::vector<long long>{0, 1}); }

    void trim() {
        while (coeffs.size() > 1 && coeffs.back() == 0) coeffs.pop_back();
        if (coeffs.empty()) coeffs.push_back(0);
    }

    Polynomial operator+(const Polynomial& o) const {
        std::vector<long long> r(std::max(coeffs.size(), o.coeffs.size()), 0);
        for (std::size_t i = 0; i < coeffs.size(); i++) r[i] = checked_add(r[i], coeffs[i]);
        for (std::size_t i = 0; i < o.coeffs.size(); i++) r[i] = checked_add(r[i], o.coeffs[i]);
        return Polynomial(r);
    }

    Polynomial operator-(const Polynomial& o) const {
        std::vector<long long> r(std::max(coeffs.size(), o.coeffs.size()), 0);
        for (std::size_t i = 0; i < coeffs.size(); i++) r[i] = checked_add(r[i], coeffs[i]);
        for (std::size_t i = 0; i < o.coeffs.size(); i++) r[i] = checked_add(r[i], -o.coeffs[i]);
        return Polynomial(r);
    }

    Polynomial operator*(const Polynomial& o) const {
        std::vector<long long> r(coeffs.size() + o.coeffs.size() - 1, 0);
        for (std::size_t i = 0; i < coeffs.size(); i++)
            for (std::size_t j = 0; j < o.coeffs.size(); j++)
                r[i + j] = checked_add(r[i + j], checked_mul(coeffs[i], o.coeffs[j]));
        return Polynomial(r);
    }

    bool is_zero() const { return coeffs.size() == 1 && coeffs[0] == 0; }

    std::string to_string() const {
        std::ostringstream out;
        bool first = true;
        for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; i--) {
            long long c = coeffs[i];
            if (c == 0) continue;
            if (!first) out << (c > 0 ? " + " : " - ");
            else if (c < 0) out << "-";
            long long mag = std::abs(c);
            if (i == 0) out << mag;
            else {
                if (mag != 1) out << mag << "*";
                out << "z";
                if (i > 1) out << "^" << i;
            }
            first = false;
        }
        if (first) out << "0";
        return out.str();
    }
};

#endif