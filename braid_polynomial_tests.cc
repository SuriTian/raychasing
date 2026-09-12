#include <iostream>
#include <string>
#include "braid.h"

using namespace std;

int main() {
    // Only 2-strand braids are currently supported (see braid.h TODO on
    // general n-strand support via Burau matrices). The 3-strand cases that
    // used to live here ({1,2,1}, {1,1,2,2}, {1,2,1,2}, {1,2}, {1,2,-1,2})
    // are removed until that lands -- the old code "passed" them via a
    // hardcoded lookup table, not by actually computing them.
    const vector<pair<BraidWord, string>> cases = {
        {{1, 1, 1}, "z^2 + 1"},        // (2,3) torus knot = trefoil
        {{1, 1, 1, 1}, "z^3 + 2*z"},   // (2,4) torus link
        {{1}, "1"},                    // unknot
        {{1, -1}, "0"},                // 2-component unlink
        {{-1, -1, -1}, "z^2 + 1"},     // mirror trefoil (Conway poly is even, so unchanged)
        {{-1, -1}, "-z"},              // negative Hopf link
        {{1, 1, 1, -1}, "z"},          // net exponent 2 -> Hopf link, regardless of order
    };

    for (const auto& [word, expected] : cases) {
        Polynomial poly = conway_polynomial(word);
        cout << "word=";
        for (int x : word) {
            cout << x << (x == word.back() ? "" : ",");
        }
        cout << " => " << poly.to_string() << "\n";
        if (poly.to_string() != expected) {
            cerr << "Expected " << expected << " but got " << poly.to_string() << "\n";
            return 1;
        }
    }

    return 0;
}
