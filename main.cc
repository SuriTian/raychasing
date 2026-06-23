#include "color.h"
#include "ray.h"
#include "vec3.h"

#include <iostream>

using namespace std; 

color ray_color(const ray& r) {
    vec3 unit_direction = unit_vector(r.direction());
    auto a = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - a) * color(1.0, 1.0, 1.0) + a * color(0.5, 0.7, 1.0);
}

int main() {
    // Image

    auto aspect_ratio = 16.0 / 9.0;

    int image_width = 400;
    int image_height = int(image_width / aspect_ratio); // respect the aspect ratio
    image_height = (image_height < 1) ? 1 : image_height; // for some reason cannot have less than 1, possibly because it would not render

    // Camera
    auto focal_length = 1.0;
    auto viewport_height = 2.0;
    auto viewport_width = viewport_height * (double(image_width)/image_height);
    auto camera_center = point3(0, 0, 0);

    // Calculate the vectors across the horizontal and down the vertical viewpoint edges
    auto viewport_u = vec3(viewport_width, 0, 0);
    auto viewport_v = vec3(0, -viewport_height, 0);

    // Calculate the horizontal and vertical delta vectors form pixel to pixel 
    auto pixel_delta_u = viewport_u / image_width;
    auto pixel_delta_v = viewport_v / image_height;

    // Calculate location of upper left pixel
    auto viewport_upper_left = camera_center - vec3(0, 0, focal_length) - viewport_u/2 - viewport_v/2;
    auto pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v); 

    // Render
    cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = 0; j < image_height; j++) { // stacks the rows top to bottom
        clog << "\rScanlines remaining: " << (image_height - j) << ' ' << flush; // progress indicter
        for (int i = 0; i < image_width; i++) { // row
            // auto pixel_color = color(double(i) / (image_width - 1), double(j) / (image_height - 1), 0); 
            
            auto pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
            auto ray_direction = pixel_center - camera_center; 
            ray r(camera_center, ray_direction);
            
            color pixel_color = ray_color(r);
            write_color(cout, pixel_color); 

            // auto r = double(i) / (image_width - 1);
            // auto g = double(j) / (image_height - 1); 
            // auto b = 0.0; 

            // int ir = int(255.999 * r);
            // int ig = int(255.999 * g);
            // int ib = int(255.999 * b);

            // cout << ir << ' ' << ig << ' ' << ib << '\n'; 
        }
    }

    clog << "\rDone.            \n";
}