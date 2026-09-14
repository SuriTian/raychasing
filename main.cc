#include <iostream>
#include <sstream>
#include <string>

#include "rtweekend.h"

#include "camera.h"
#include "hittable.h"
#include "hittable_list.h"
#include "material.h"
#include "sphere.h"
#include "trefoil.h"
#include "braid.h"

// The renderer below only knows how to draw a trefoil, so we only attempt a
// render when the braid actually closes to one.
static bool is_supported_braid(const BraidWord& word) {
    if (word.empty()) return false;
    return is_trefoil_like(word);
}

int main(int argc, char** argv) {
    std::string input = (argc > 1) ? argv[1] : "1,1,1";

    BraidWord braid;
    try {
        braid = parse_braid_word(input);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Could not read that braid word: " << e.what() << "\n"
                  << "Expected a comma-separated list of nonzero integers, e.g. 1,-2,1\n";
        return 1;
    }

    // stdout carries the PPM and nothing else -- the usual way to run this is
    // `RayChasing.exe > image.ppm`, and anything else printed there lands in
    // the middle of the image file.
    std::cerr << "Input braid: " << input << "\n";
    try {
        Polynomial polynomial = conway_polynomial(braid);
        std::cerr << "Conway polynomial: " << polynomial.to_string() << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Conway polynomial: unavailable (" << e.what() << ")\n";
    }

    if (!is_supported_braid(braid)) {
        std::cerr << "That braid does not close to a trefoil, which is the only knot "
                     "this renderer can draw so far.\n";
        return 0;
    }

    // World
    hittable_list world;

    auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    // auto material_center = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    // auto material_left   = make_shared<dielectric>(1.50);
    // auto material_bubble = make_shared<dielectric>(1.00 / 1.50);
    // auto material_right  = make_shared<metal>(color(0.8, 0.6, 0.2), 1.0);

    // world.add(make_shared<sphere>(point3( 0.0, -100.5, -1.0), 100.0, material_ground));
    // world.add(make_shared<sphere>(point3( 0.0, 0.0, -1.2), 0.5, material_center));
    // world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), 0.5, material_left));
    // world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), 0.4, material_bubble));
    // world.add(make_shared<sphere>(point3( 1.0, 0.0, -1.0), 0.5, material_right));

    // TREFOIL
    auto material_trefoil = make_shared<dielectric>(1.5);
    world.add(make_shared<trefoil>(point3(0, 0, -1), 0.15, 0.04, 100, material_trefoil));

    // auto R = cos(pi/4);

    // auto material_left  = make_shared<lambertian>(color(0,0,1));
    // auto material_right = make_shared<lambertian>(color(1,0,0));

    // world.add(make_shared<sphere>(point3(-R, 0, -1), R, material_left));
    // world.add(make_shared<sphere>(point3( R, 0, -1), R, material_right));

    // camera
    camera cam;
    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 400;
    cam.samples_per_pixel = 100;
    cam.max_depth = 50;

    cam.vfov = 20; // large is farther, small is closer
    cam.lookfrom = point3(-2, 2, 1);
    cam.lookat = point3(0, 0, -1);
    cam.vup = vec3(0, 1, 0);
    cam.defocus_angle = 0;
    // cam.focus_distance = 3.4;

    cam.render(world);
    return 0;
}