#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "braid.h"

using namespace std; 

// Parses a comma-separated list like "1,-2,1" into a BraidWord.
static BraidWord parse_braid_word(const std::string& text) {
    BraidWord word;
    stringstream ss(text);
    string token;
    while (getline(ss, token, ',')) {
        if (token.empty()) continue;
        int value = stoi(token);
        word.push_back(value);
    }
    return word;
}

int main(int argc, char** argv) {
    string input = argc > 1 ? argv[1] : "1,1,1";
    BraidWord word = parse_braid_word(input);

    cout << "Input braid: " << input << "\n";
    cout << "Strands: " << strand_count(word) << "\n";

    try {
        Polynomial poly = conway_polynomial(word);
        cout << "Conway polynomial: " << poly.to_string() << "\n";
    } catch (const std::logic_error& e) {
        cout << "Conway polynomial: unavailable (" << e.what() << ")\n";
    }

    auto torus = match_torus_knot(word);
    if (torus) {
        cout << "Matches torus knot pattern: T(" << torus->first << ", " << torus->second << ")\n";
    } else {
        cout << "Does not match a simple torus knot pattern.\n";
    }

    if (word == BraidWord{1, 1, 1}) {
        cout << "Trefoil-like input detected; current renderer supports trefoil geometry.\n";
    }

    return 0;
}
