#include <iostream>
#include <string>
#include "braid.h"

using namespace std;

int main() {
    const vector<pair<BraidWord, string>> cases = {
        {{1, 1, 1}, "z^2 + 1"},
        {{1, 1, 1, 1}, "2*z"},
        {{1, 2, 1}, "z^2 + 1"},
        {{1}, "1"},
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
