#ifndef SPHERE_H
#define SPHERE_H
#include "shapes.h"

class Sphere : public Shape {
  public:
    Sphere(Vector _C, double _R, Vector _albedo, bool _mirror = false,
           bool _transparent = false, bool _is_interior = false);
    Vector get_center() { return C; }
    double get_radius() { return R; }
    bool intersect(const Ray &r, double &t, Vector &N) override;
    // Vector get_normal(const Vector &P) override;

  private:
    Vector C;
    double R;
};
#endif
