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
    // virtual Vector get_normal(const Vector &P) = 0;
    virtual bool intersect(const Ray &r, double &t, Vector &N) = 0;
    bool is_mirror() { return is_mirror_b; }
    bool is_transparent() { return is_transparent_b; }
    bool is_interior() { return is_interior_b; }
    Vector get_color() { return albedo; }
    void set_color(Vector c) { albedo = c; }
    void set_mirror(bool m) { is_mirror_b = m; }
    void set_transparent(bool t) { is_transparent_b = t; }

  protected:
    Vector albedo;
    bool is_mirror_b;
    bool is_transparent_b;
    bool is_interior_b;
};
#endif
