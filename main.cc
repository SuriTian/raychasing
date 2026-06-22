#include "color.h"
#include "vec3.h"

#include <iostream>

using namespace std; 

int main() {
    // Image
    int image_width = 256;
    int image_height = 256; 

    // Render
    cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = 0; j < image_height; j++) { // stacks the rows top to bottom
        clog << "\rScanlines remaining: " << (image_height - j) << ' ' << flush; // progress indicter
        for (int i = 0; i < image_width; i++) { // row
            auto pixel_color = color(double(i) / (image_width - 1), double(j) / (image_height - 1), 0); 
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