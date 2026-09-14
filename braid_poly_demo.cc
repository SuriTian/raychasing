#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "braid.h"

using namespace std; 

int main(int argc, char** argv) {
    string input = argc > 1 ? argv[1] : "1,1,1";

    BraidWord word;
    try {
        word = parse_braid_word(input);
    } catch (const invalid_argument& e) {
        cerr << "Could not read that braid word: " << e.what() << "\n"
             << "Expected a comma-separated list of nonzero integers, e.g. 1,-2,1\n";
        return 1;
    }

    cout << "Input braid: " << input << "\n";
    cout << "Strands: " << strand_count(word) << "\n";
    cout << "Closes to: " << closure_component_count(word) << " component(s)\n";

    try {
        Polynomial poly = conway_polynomial(word);
        cout << "Conway polynomial: " << poly.to_string() << "\n";
    } catch (const std::exception& e) {
        cout << "Conway polynomial: unavailable (" << e.what() << ")\n";
    }

    auto torus = match_torus_knot(word);
    if (torus) {
        cout << "Matches torus knot pattern: T(" << torus->first << ", " << torus->second << ")\n";
    } else {
        cout << "Does not match a simple torus knot pattern.\n";
    }

    if (is_trefoil_like(word)) {
        cout << "Closure is a trefoil; the renderer can draw this one.\n";
    }

    return 0;
}
