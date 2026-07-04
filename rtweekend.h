#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>

using std::make_shared;
using std::shared_ptr;
using namespace std;

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// utility functions
inline double degrees_to_radians(double degree) {
    return degree * pi / 180.0;
}

inline double random_double() {
    // real in [0, 1)
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
    // real in [min, max)
    return min + (max - min) * random_double(); 
}

#include "color.h"
#include "interval.h"
#include "ray.h"
#include "vec3.h"

#endif