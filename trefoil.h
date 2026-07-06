#ifndef TREFOIL_H
#define TREFOIL_H

#include "hittable.h"
#include <vector>
#include <cmath>

class trefoil : public hittable {
public:
    trefoil(point3 center, double scale, double tube_radius,
            int segments, shared_ptr<material> mat)
        : center(center), tube_radius(tube_radius), mat(mat)
    {
        // Oversample the raw curve, then resample evenly by arc length
        const int raw_n = segments * 8;
        std::vector<point3> raw(raw_n);
        for (int i = 0; i < raw_n; i++) {
            double t = 2.0 * pi * i / raw_n;
            raw[i] = center + scale * curve(t);
        }

        std::vector<double> arclen(raw_n + 1, 0.0);
        for (int i = 0; i < raw_n; i++) {
            point3 a = raw[i];
            point3 b = raw[(i + 1) % raw_n];
            arclen[i + 1] = arclen[i] + (b - a).length();
        }
        double total = arclen[raw_n];

        points.resize(segments);
        int j = 0;
        for (int i = 0; i < segments; i++) {
            double target = total * i / segments;
            while (j < raw_n && arclen[j + 1] < target) j++;
            double segStart = arclen[j], segEnd = arclen[j + 1];
            double frac = (segEnd > segStart) ? (target - segStart) / (segEnd - segStart) : 0.0;
            point3 a = raw[j], b = raw[(j + 1) % raw_n];
            points[i] = a + frac * (b - a);
        }

        // Precompute segment direction vectors and squared lengths ONCE,
        // instead of recomputing them on every sdf() call.
        int n = points.size();
        seg_dir.resize(n);
        seg_len2.resize(n);
        for (int i = 0; i < n; i++) {
            vec3 ab = points[(i + 1) % n] - points[i];
            seg_dir[i] = ab;
            seg_len2[i] = fmax(dot(ab, ab), 1e-12); // avoid div-by-zero
        }

        bbox_radius = scale * 3.5 + tube_radius; // rough bound for early-out
    }

    bool hit(const ray& r, interval ray_t, hit_record& rec) const override {
        // The ray's direction vector is NOT guaranteed to be unit length in
        // this codebase, but our SDF marching works in real-world distances.
        // So we normalize direction for marching, then convert t back at the end.
        double dir_len = r.direction().length();
        if (dir_len < 1e-12) return false;
        vec3 unit_dir = r.direction() / dir_len;

        double dist_min = ray_t.min * dir_len;
        double dist_max = ray_t.max * dir_len;

        // bounding sphere reject (uses normalized direction)
        vec3 oc = center - r.origin();
        double tca = dot(oc, unit_dir);
        double d2 = dot(oc, oc) - tca * tca;
        if (d2 > bbox_radius * bbox_radius) return false;

        double dist = dist_min;
        for (int i = 0; i < max_steps; i++) {
            point3 p = r.origin() + dist * unit_dir;
            double d = sdf(p);
            if (d < epsilon) {
                double t = dist / dir_len; // back to original ray parameterization
                if (t < ray_t.min || t > ray_t.max) return false;
                rec.t = t;
                rec.p = p;
                vec3 outward_normal = calc_normal(p);
                rec.set_face_normal(r, outward_normal);
                rec.mat = mat;
                return true;
            }
            dist += fmax(d, 1e-3); // clamp step so we never stall/move backwards
            if (dist > dist_max) return false;
        }
        return false;
    }

private:
    point3 center;
    double tube_radius;
    shared_ptr<material> mat;
    std::vector<point3> points;
    std::vector<vec3> seg_dir;    // precomputed b - a per segment
    std::vector<double> seg_len2; // precomputed |b - a|^2 per segment
    double bbox_radius;
    static constexpr int max_steps = 96;
    static constexpr double epsilon = 5e-4;

    static point3 curve(double t) {
        double x = sin(t) + 2 * sin(2 * t);
        double y = cos(t) - 2 * cos(2 * t);
        double z = -sin(3 * t);
        return point3(x, y, z);
    }

    // Returns squared distance to segment i (skip sqrt until we know the winner)
    double dist2_to_segment(const point3& p, int i) const {
        const point3& a = points[i];
        vec3 ap = p - a;
        double t = dot(ap, seg_dir[i]) / seg_len2[i];
        t = fmax(0.0, fmin(1.0, t));
        point3 closest = a + t * seg_dir[i];
        vec3 diff = p - closest;
        return dot(diff, diff);
    }

    double sdf(const point3& p) const {
        double min_d2 = infinity;
        int n = points.size();
        for (int i = 0; i < n; i++) {
            double d2 = dist2_to_segment(p, i);
            if (d2 < min_d2) min_d2 = d2;
        }
        return sqrt(min_d2) - tube_radius;
    }

    vec3 calc_normal(const point3& p) const {
        const double h = 1e-4;
        return unit_vector(vec3(
            sdf(p + vec3(h,0,0)) - sdf(p - vec3(h,0,0)),
            sdf(p + vec3(0,h,0)) - sdf(p - vec3(0,h,0)),
            sdf(p + vec3(0,0,h)) - sdf(p - vec3(0,0,h))
        ));
    }
};

#endif