#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"

class sphere : public hittable {
    public:
        sphere(const point3& center, double radius) : center(center), radius(fmax(0, radius)) {}

        bool hit(const ray& r, double ray_tmin, double ray_tmax, hit_record& rec) const override {
            vec3 oc = center - r.origin();
            auto a = r.direction().length_squared();
            auto h = dot(r.direction(), oc);
            auto c = oc.length_squared() - radius * radius;

            auto discriminant = h * h - a * c;
            if (discriminant < 0) return false;

            auto sqrtdis = sqrt(discriminant);

            // find nearest root that lies in acceptable range 
            auto root = (h - sqrtdis) / a; // simplified form of [-b - sqrt(b^2 - 4ac)] / 2a
            if (root <= ray_tmin || root >= ray_tmax) { // check - root
                root = (h + sqrtdis) / a; // check + root
                if (root <= ray_tmin || root >= ray_tmax) return false; 
            }

            rec.t = root; // took the root that was in bounds 
            rec.p = r.at(rec.t);
            // rec.normal = (rec.p - center) / radius;
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);

            return true; 
        }

    private:
        point3 center;
        double radius;
};

#endif