#include <iostream>
#include <string>
#include <vector>
#include <random>
#include "braid.h"

using namespace std;

static int failures = 0;

static void fail(const string& label, const string& detail) {
    cerr << "FAIL: " << label << " -- " << detail << "\n";
    failures++;
}

static void check(bool ok, const string& label) {
    if (!ok) fail(label, "check returned false");
}

static string word_to_string(const BraidWord& word) {
    string text;
    for (size_t i = 0; i < word.size(); i++) {
        if (i) text += ",";
        text += to_string(word[i]);
    }
    return text.empty() ? "(empty)" : text;
}

// The Conway polynomial cannot tell a knot from its mirror image, so these
// cases pin down the size and shape of the answer, not its handedness.
static void check_conway(const BraidWord& word, const string& expected, const string& label) {
    try {
        string actual = conway_polynomial(word).to_string();
        cout << "  " << word_to_string(word) << " => " << actual << "   (" << label << ")\n";
        if (actual != expected) fail(label, "expected " + expected + ", got " + actual);
    } catch (const exception& e) {
        fail(label, string("threw: ") + e.what());
    }
}

static void check_rejected(const BraidWord& word, const string& label) {
    try {
        conway_polynomial(word);
        fail(label, "expected a rejection but got a polynomial");
    } catch (const logic_error&) {
        cout << "  " << word_to_string(word) << " => rejected   (" << label << ")\n";
    }
}

static bool matrices_equal(const PolyMatrix& a, const PolyMatrix& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        for (size_t j = 0; j < a.size(); j++)
            if (!(a[i][j] == b[i][j])) return false;
    return true;
}

// If these fail, the Burau matrices are not a braid group representation at
// all and every polynomial computed from them is meaningless.
static void check_burau_is_a_representation() {
    cout << "Burau representation:\n";
    for (int n = 2; n <= 6; n++) {
        for (int i = 1; i <= n - 1; i++) {
            PolyMatrix id = identity_matrix(n - 1);
            PolyMatrix positive = burau_generator(i, n);
            PolyMatrix negative = burau_generator(-i, n);
            check(matrices_equal(multiply(positive, negative), id) &&
                      matrices_equal(multiply(negative, positive), id),
                  "sigma_" + to_string(i) + " inverse on " + to_string(n) + " strands");
        }
    }

    for (int n = 3; n <= 6; n++) {
        for (int i = 1; i <= n - 2; i++) {
            PolyMatrix a = burau_generator(i, n);
            PolyMatrix b = burau_generator(i + 1, n);
            check(matrices_equal(multiply(multiply(a, b), a), multiply(multiply(b, a), b)),
                  "braid relation at sigma_" + to_string(i) + " on " + to_string(n) + " strands");
        }
    }

    for (int n = 4; n <= 6; n++) {
        for (int i = 1; i <= n - 1; i++) {
            for (int j = i + 2; j <= n - 1; j++) {
                PolyMatrix a = burau_generator(i, n);
                PolyMatrix b = burau_generator(j, n);
                check(matrices_equal(multiply(a, b), multiply(b, a)),
                      "sigma_" + to_string(i) + " commutes with sigma_" + to_string(j));
            }
        }
    }
    cout << "  generator inverses, braid relations and far commutativity hold\n";
}

// Two strands go through a closed-form recursion instead of the matrices, so
// the two routes must agree -- for links as well as knots.
static void check_two_strand_routes_agree() {
    cout << "Two-strand fast path vs Burau:\n";
    const vector<BraidWord> words = {
        {1}, {1, 1, 1}, {1, 1, 1, 1, 1}, {-1, -1, -1}, {1, 1, -1, 1, 1}, {1, 1, 1, 1, 1, 1, 1},
        {1, 1}, {-1, -1}, {1, 1, 1, 1}, {1, -1, 1, 1},   // links too
    };
    for (const BraidWord& word : words) {
        string closed_form = conway_polynomial(word).to_string();
        string via_burau = conway_polynomial_via_burau(word, 2).to_string();
        cout << "  " << word_to_string(word) << " => " << closed_form << "\n";
        if (closed_form != via_burau) {
            fail("two-strand agreement for " + word_to_string(word),
                 "closed form " + closed_form + " vs Burau " + via_burau);
        }
    }
}

// Torus knots have a closed-form Alexander polynomial that owes nothing to
// Burau matrices:
//     Delta(t) = (t^(pq) - 1)(t - 1) / ((t^p - 1)(t^q - 1))
// Running both routes over a spread of (p, q) checks the whole pipeline --
// matrices, determinant, division, rewrite in z -- against outside arithmetic.
static void check_against_closed_form_torus_knots() {
    cout << "Torus knots vs their closed form:\n";
    auto t_power_minus_one = [](int k) {
        vector<long long> c(k + 1, 0);
        c[0] = -1;
        c[k] = 1;
        return LaurentPoly(0, c);
    };

    for (int p = 2; p <= 5; p++) {
        for (int q = 2; q <= 7; q++) {
            if (gcd(p, q) != 1) continue;
            string label = "T(" + to_string(p) + "," + to_string(q) + ")";

            LaurentPoly closed_form =
                (t_power_minus_one(p * q) * t_power_minus_one(1))
                    .divide_exact(t_power_minus_one(p) * t_power_minus_one(q));
            // A torus knot is a knot, so the sign comes from the polynomial
            // itself and the unit argument is not consulted.
            string expected = conway_from_alexander(closed_form, 1, 1).to_string();

            BraidWord word = generate_torus_braid(p, q);
            string actual = conway_polynomial(word).to_string();

            cout << "  " << label << " => " << actual << "\n";
            if (actual != expected) fail(label, "closed form says " + expected + ", braid says " + actual);
        }
    }
}

// The Conway polynomial belongs to the closure, not to the word, so it has to
// survive every move that leaves the closure alone: rotating the word,
// conjugating it, and stabilizing it onto an extra strand. Random words here,
// so nothing is hand-picked to pass -- and stabilizing a 2-strand word pushes
// it through the matrices, checking the closed-form path against them.
static void check_closure_moves_do_not_change_the_answer() {
    cout << "Invariance under Markov moves:\n";
    mt19937 rng(20260913);
    int knots = 0;
    int comparisons = 0;

    for (int trial = 0; trial < 2000 && knots < 60; trial++) {
        const int strands = 2 + static_cast<int>(rng() % 4);
        const int length = 2 + static_cast<int>(rng() % 7);

        BraidWord word;
        for (int i = 0; i < length; i++) {
            int generator = 1 + static_cast<int>(rng() % (strands - 1));
            if (rng() % 2) generator = -generator;
            word.push_back(generator);
        }
        if (strand_count(word) != strands) continue;
        knots++;

        const string expected = conway_polynomial(word).to_string();

        auto expect_same = [&](const BraidWord& moved, const string& move) {
            comparisons++;
            string actual = conway_polynomial(moved).to_string();
            if (actual != expected) {
                fail(move + " changed the polynomial",
                     word_to_string(word) + " => " + expected + " but " +
                         word_to_string(moved) + " => " + actual);
            }
        };

        BraidWord rotated(word.begin() + 1, word.end());
        rotated.push_back(word.front());
        if (strand_count(rotated) == strands) expect_same(rotated, "cyclic rotation");

        for (int generator = 1; generator <= strands - 1; generator++) {
            BraidWord conjugated;
            conjugated.push_back(generator);
            conjugated.insert(conjugated.end(), word.begin(), word.end());
            conjugated.push_back(-generator);
            if (strand_count(conjugated) == strands) expect_same(conjugated, "conjugation");
        }

        for (int sign : {1, -1}) {
            BraidWord stabilized = word;
            stabilized.push_back(sign * strands);
            expect_same(stabilized, "stabilization");
        }
    }

    cout << "  " << knots << " random braids survived " << comparisons << " closure moves\n";
}

// The Conway polynomial multiplies under connected sum, and the closure of
// sigma_1^a sigma_2^b sigma_3^c is exactly the connected sum of the three
// two-strand closures. That gives an answer worked out without any matrices,
// for links of up to four components.
static void check_connected_sums_multiply() {
    cout << "Connected sums:\n";
    int cases = 0;
    for (int a = -3; a <= 3; a++) {
        for (int b = -3; b <= 3; b++) {
            for (int c = -3; c <= 3; c++) {
                if (a == 0 || b == 0) continue;

                BraidWord word;
                for (int i = 0; i < abs(a); i++) word.push_back(a > 0 ? 1 : -1);
                for (int i = 0; i < abs(b); i++) word.push_back(b > 0 ? 2 : -2);
                for (int i = 0; i < abs(c); i++) word.push_back(c > 0 ? 3 : -3);

                Polynomial expected = conway_polynomial_two_strand(a) * conway_polynomial_two_strand(b);
                if (c != 0) expected = expected * conway_polynomial_two_strand(c);

                string actual = conway_polynomial(word).to_string();
                cases++;
                if (actual != expected.to_string()) {
                    fail("connected sum " + word_to_string(word),
                         "expected " + expected.to_string() + ", got " + actual);
                }
            }
        }
    }
    cout << "  " << cases << " connected sums matched the product of their pieces\n";
}

// Coefficients outgrow 64 bits well before the knots get interesting, and a
// wrapped coefficient looks like a perfectly ordinary answer. Big inputs have
// to refuse rather than quietly lie.
static void check_oversized_input_is_refused() {
    cout << "Coefficient overflow:\n";

    // The largest torus knots that still fit, and still correct.
    for (int p : {9, 10}) {
        BraidWord word = generate_torus_braid(p, p + 1);
        try {
            Polynomial conway = conway_polynomial(word);
            bool shaped_right = conway.coeffs.front() == 1 && conway.coeffs.back() == 1;
            for (size_t i = 1; i < conway.coeffs.size(); i += 2) {
                if (conway.coeffs[i] != 0) shaped_right = false;   // knots have no odd powers
            }
            for (size_t i = 0; i < conway.coeffs.size(); i += 2) {
                if (conway.coeffs[i] <= 0) shaped_right = false;   // torus knots stay positive
            }
            cout << "  T(" << p << "," << p + 1 << ") degree " << (conway.coeffs.size() - 1)
                 << " computed cleanly\n";
            check(shaped_right, "T(" + to_string(p) + "," + to_string(p + 1) + ") coefficient shape");
        } catch (const exception& e) {
            fail("T(" + to_string(p) + "," + to_string(p + 1) + ")", string("threw: ") + e.what());
        }
    }

    // Past that, refusing is the only honest answer.
    for (int p : {11, 12}) {
        BraidWord word = generate_torus_braid(p, p + 1);
        try {
            Polynomial conway = conway_polynomial(word);
            fail("T(" + to_string(p) + "," + to_string(p + 1) + ")",
                 "expected an overflow refusal, got degree " + to_string(conway.coeffs.size() - 1));
        } catch (const overflow_error&) {
            cout << "  T(" << p << "," << p + 1 << ") refused, as it should be\n";
        }
    }

    // A predicate must not throw, whatever it is handed.
    try {
        check(!is_trefoil_like(generate_torus_braid(12, 13)), "huge braid is not a trefoil");
    } catch (const exception& e) {
        fail("is_trefoil_like on a huge braid", string("threw: ") + e.what());
    }
}

// The renderer draws a trefoil, so its gate has to agree with the knot the
// braid actually closes to -- not with the shape of the word.
static void check_trefoil_recognition() {
    cout << "Trefoil recognition:\n";
    struct Expectation { BraidWord word; bool trefoil; const char* note; };
    const vector<Expectation> expectations = {
        {{1, 1, 1}, true, "T(2,3)"},
        {{-1, -1, -1}, true, "mirror image, which Conway cannot see"},
        {{1, 2, 1, 2}, true, "T(3,2), the same knot on 3 strands"},
        {{1, 2, 1}, false, "trefoil-shaped but closes to a 2-component link"},
        {{1, 2}, false, "unknot"},
        {{1, 1}, false, "Hopf link"},
        {{1, -2, 1, -2}, false, "figure-eight knot"},
        {{1, 1, 1, 1, 1}, false, "cinquefoil"},
    };
    for (const Expectation& e : expectations) {
        bool actual = is_trefoil_like(e.word);
        cout << "  " << word_to_string(e.word) << " => " << (actual ? "trefoil" : "not a trefoil")
             << "   (" << e.note << ")\n";
        if (actual != e.trefoil) {
            fail("trefoil recognition for " + word_to_string(e.word),
                 string("expected ") + (e.trefoil ? "trefoil" : "not a trefoil"));
        }
    }

    if (match_torus_knot({1, 2, 1})) {
        fail("match_torus_knot for 1,2,1", "a 2-component closure was reported as a torus knot");
    }
}

static void check_component_counts() {
    cout << "Closure component counts:\n";
    struct Expectation { BraidWord word; int strands; int components; };
    const vector<Expectation> expectations = {
        {{1, 1, 1}, 2, 1},        // trefoil
        {{1, 1}, 2, 2},           // Hopf link
        {{1, 2, 1}, 3, 2},        // looks trefoil-shaped, closes to a 2-component link
        {{1, 2, 1, 2}, 3, 1},     // T(3,2), a knot
        {{1, 2}, 3, 1},           // unknot
        {{1, 1, 2, 2}, 3, 3},
    };
    for (const Expectation& e : expectations) {
        int actual = closure_component_count(e.word, e.strands);
        cout << "  " << word_to_string(e.word) << " => " << actual << " component(s)\n";
        if (actual != e.components) {
            fail("component count for " + word_to_string(e.word),
                 "expected " + to_string(e.components) + ", got " + to_string(actual));
        }
    }
}

int main() {
    cout << "Two-strand words:\n";
    check_conway({1, 1, 1}, "z^2 + 1", "trefoil, T(2,3)");
    check_conway({1, 1, 1, 1}, "z^3 + 2*z", "T(2,4) torus link");
    check_conway({1}, "1", "unknot");
    check_conway({1, -1}, "0", "2-component unlink");
    check_conway({-1, -1, -1}, "z^2 + 1", "mirror trefoil");
    check_conway({-1, -1}, "-z", "negative Hopf link");
    check_conway({1, 1, 1, -1}, "z", "net exponent 2, order does not matter");

    cout << "Wider braids closing to knots:\n";
    check_conway({1, 2}, "1", "unknot as a 3-strand braid");
    check_conway({1, 2, 1, 2}, "z^2 + 1", "trefoil as T(3,2)");
    check_conway({1, -2, 1, -2}, "-z^2 + 1", "figure-eight knot");
    check_conway({1, 1, 1, 1, 1}, "z^4 + 3*z^2 + 1", "cinquefoil T(2,5)");
    check_conway({1, 2, 1, 2, 1, 2, 1, 2}, "z^6 + 5*z^4 + 5*z^2 + 1", "T(3,4) on 3 strands");
    check_conway({1, 2, 3, 1, 2, 3, 1, 2, 3}, "z^6 + 5*z^4 + 5*z^2 + 1", "T(4,3) on 4 strands");
    check_conway({1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4},
                 "z^12 + 11*z^10 + 44*z^8 + 77*z^6 + 56*z^4 + 15*z^2 + 1", "T(5,4) on 5 strands");

    cout << "Wider braids closing to links:\n";
    check_conway({1, 2, 1}, "z", "Hopf link on 3 strands, same answer as 1,1 on two");
    check_conway({1, 1, 2, 2}, "z^2", "3-chain: two Hopf links joined");
    check_conway({1, 1, 1, 2, 2}, "z^3 + z", "trefoil joined to a Hopf link");
    check_conway({1, -2, 1, -2, 1, -2}, "z^4", "Borromean rings");
    check_conway({1, 2, 1, 2, 1, 2}, "z^4 + 3*z^2", "T(3,3), a 3-component torus link");

    cout << "Split links vanish:\n";
    check_conway({1, -1, 3}, "0", "a link that falls apart into separate pieces");

    cout << "Bad input:\n";
    check_rejected({0}, "0 is not a generator");
    check_rejected({1, 0, 2}, "0 in the middle of a word");

    check_burau_is_a_representation();
    check_two_strand_routes_agree();
    check_against_closed_form_torus_knots();
    check_component_counts();
    check_trefoil_recognition();
    check_closure_moves_do_not_change_the_answer();
    check_connected_sums_multiply();
    check_oversized_input_is_refused();

    if (failures) {
        cerr << failures << " check(s) failed\n";
        return 1;
    }
    cout << "All braid polynomial checks passed.\n";
    return 0;
}
