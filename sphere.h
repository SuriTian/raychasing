#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"

class sphere : public hittable {
    public:
        sphere(const point3& center, double radius, shared_ptr<material> mat) : center(center), radius(fmax(0, radius)), mat(mat) {}

        bool hit(const ray& r, interval ray_t, hit_record& rec) const override {
            vec3 oc = center - r.origin();
            auto a = r.direction().length_squared();
            auto h = dot(r.direction(), oc);
            auto c = oc.length_squared() - radius * radius;

            auto discriminant = h * h - a * c;
            if (discriminant < 0) return false;

            auto sqrtdis = sqrt(discriminant);

            // find nearest root that lies in acceptable range 
            auto root = (h - sqrtdis) / a; // simplified form of [-b - sqrt(b^2 - 4ac)] / 2a
            if (!ray_t.surrounds(root)) { // check - root
                root = (h + sqrtdis) / a; // check + root
                if (!ray_t.surrounds(root)) return false; 
            }

            rec.t = root; // took the root that was in bounds 
            rec.p = r.at(rec.t);
            // rec.normal = (rec.p - center) / radius;
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);
            rec.mat = mat;

            return true; 
        }

    private:
        point3 center;
        double radius;
        shared_ptr<material> mat;
};

#endif