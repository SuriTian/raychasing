#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "braid.h"

using namespace std; 

static BraidWord parse_braid_word(const std::string& text) {
    BraidWord word;
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) continue;
        int value = std::stoi(token);
        word.push_back(value);
    }
    return word;
}

int main(int argc, char** argv) {
    std::string input = argc > 1 ? argv[1] : "1,1,1";
    BraidWord word = parse_braid_word(input);
    Polynomial poly = conway_polynomial(word);

    std::cout << "Input braid: " << input << "\n";
    std::cout << "Strands: " << strand_count(word) << "\n";
    std::cout << "Conway polynomial: " << poly.to_string() << "\n";

    auto torus = match_torus_knot(word);
    if (torus) {
        std::cout << "Matches torus knot pattern: T(" << torus->first << ", " << torus->second << ")\n";
    } else {
        std::cout << "Does not match a simple torus knot pattern.\n";
    }

    if (word == BraidWord{1, 1, 1}) {
        std::cout << "Trefoil-like input detected; current renderer supports trefoil geometry.\n";
    }

    return 0;
}
