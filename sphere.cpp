#include "sphere.h"
#include "shapes.h"
#include <cmath>
#include <math.h>
double refraction_index = 1.0 / 1.52;

Sphere::Sphere(Vector _C, double _R, Vector _albedo, bool _mirror,
               bool _transparent, bool _is_interior) {
    C = _C;
    R = _R;
    albedo = _albedo;
    is_mirror_b = _mirror;
    is_transparent_b = _transparent;
    is_interior_b = _is_interior;
}
bool Sphere::intersect(const Ray &r, double &t, Vector &N) {
    Vector dif = C - r.get_origin();
    double dot_prod = dot(r.get_dir(), dif);
    double delta = dot_prod * dot_prod - dif.norm2() + R * R;

    bool intersected = false;
    if (delta >= 0) {
        double sqrt_delta = sqrt(delta);
        if (dot_prod + sqrt_delta > 0) {
            intersected = true;
            t = dot_prod > sqrt_delta ? dot_prod - sqrt_delta
                                      : dot_prod + sqrt_delta;
        } else {
            t = dot_prod + sqrt_delta;
        }
        Vector P = r.get_origin() + t * r.get_dir();
        N = (is_interior_b ? -1 : 1) * (P - C) / R;
    }
    return intersected;
}
