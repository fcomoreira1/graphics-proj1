#ifndef SHAPES_H
#define SHAPES_H
#include "vector.h"

class Ray {
  public:
    Ray(Vector o, Vector u) {
        O = o;
        U = u;
        U.normalize();
    }
    Vector get_origin() const { return O; }
    Vector get_dir() const { return U; }

  private:
    Vector O;
    Vector U;
};

class Shape {
  public:
    virtual Vector get_normal(const Vector &P) = 0;
    virtual bool intersect(const Ray &r, double &t) = 0;
    bool is_mirror() { return is_mirror_b; }
    bool is_transparent() { return is_transparent_b; }
    bool is_interior() { return is_interior_b; }
    Vector get_color() { return albedo; }

  protected:
    Vector albedo;
    bool is_mirror_b;
    bool is_transparent_b;
    bool is_interior_b;
};
class Sphere : public Shape {
  public:
    Sphere(Vector _C, double _R, Vector _albedo, bool _mirror = false,
           bool _transparent = false, bool _is_interior = false);
    Vector get_center() { return C; }
    double get_radius() { return R; }
    bool intersect(const Ray &r, double &t) override;
    Vector get_normal(const Vector &P) override;

  private:
    Vector C;
    double R;
};
#endif
