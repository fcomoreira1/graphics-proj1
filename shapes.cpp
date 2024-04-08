#include <cmath>
#include <math.h>
#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "./stb/stb_image.h"

class Vector {
  public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const { return sqrt(norm2()); }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double &operator[](int i) { return data[i]; };
    double data[3];
    void print() {
        std::cout << data[0] << " " << data[1] << " " << data[2] << std::endl;
    }
};

Vector operator+(const Vector &a, const Vector &b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector &b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &a, const double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator/(const Vector &a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector &a, const Vector &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector &a, const Vector &b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                  a[0] * b[1] - a[1] * b[0]);
}

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

class Sphere {
  public:
    Sphere(Vector _C, double _R, Vector _albedo = Vector(1.0, 1.0, 1.0),
           bool _mirror = false) {
        C = _C;
        R = _R;
        albedo = _albedo;
        mirror = _mirror;
    }
    bool intersect(const Ray &r, double &t) {
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
        }
        return intersected;
    }
    Vector get_color() { return albedo; }
    Vector get_center() { return C; }
    double get_radius() { return R; }
    bool is_mirror() { return mirror; }

  private:
    Vector C;
    double R;
    Vector albedo;
    bool mirror;
};

class Scene {
  public:
    Scene(Vector _center, int _H, int _W, double _alpha,
          Vector _light_source = Vector(-10, 20, 40)) {
        center = _center;
        H = _H;
        W = _W;
        alpha = _alpha;
        // f = (double)_W / (2.0 * tan(_alpha / 2.0));
        light_source = _light_source;
    }
    void add_sphere(Sphere S) { spheres.emplace_back(S); }
    void render() {
        std::vector<unsigned char> image(W * H * 3, 0);
        for (size_t i = 0; i < H; i++) {
            for (size_t j = 0; j < W; j++) {
                double z = -(double)W / (2.0 * tan(alpha / 2.0));
                Vector u(j + 0.5 - W / 2.0, H / 2.0 - i - 0.5, z);
                Ray r(center, u);
                Vector color = raytrace(r, 10);
                image[(i * W + j) * 3 + 0] =
                    std::min(255.0, pow(color.data[0], 1 / 2.2));
                image[(i * W + j) * 3 + 1] =
                    std::min(255.0, pow(color.data[1], 1 / 2.2));
                image[(i * W + j) * 3 + 2] =
                    std::min(255.0, pow(color.data[2], 1 / 2.2));
            }
        }
        stbi_write_png("image.png", W, H, 3, &image[0], 0);
    }

  private:
    Vector raytrace(const Ray &r, int ray_depth = 0) {
        if (ray_depth < 0) {
            return Vector(0.0, 0.0, 0.0);
        }
        double min_t = MAXFLOAT;
        double t;
        Sphere *S;
        for (auto &s : spheres) {
            if (s.intersect(r, t) && t < min_t) {
                S = &s;
                min_t = t;
            }
        }
        Vector rho = Vector();
        if (min_t < MAXFLOAT) {
            Vector P = r.get_origin() + (min_t * r.get_dir());
            Vector N = (P - S->get_center()) / S->get_radius();
            Vector V = light_source - P;
            P = P + 1e-8 * N;
            if (S->is_mirror()) {
                std::cout << "This is a fucking mirror" << std::endl;
                Vector u = r.get_dir() - 2 * dot(r.get_dir(), N) * N;
                return raytrace(Ray(P, u), ray_depth - 1);
            }
            double t0;
            bool non_direct_path = false;
            for (auto s : spheres) {
                if (s.intersect(Ray(P, light_source - P), t0)) {
                    if (t0 < (light_source - P).norm()) {
                        non_direct_path = true;
                        break;
                    }
                }
            }
            if (!non_direct_path) {
                rho = (intensity / (4 * M_PI * V.norm2())) *
                      dot(N, V / V.norm()) * (S->get_color() / M_PI);
            }
        }
        return rho;
    }
    std::vector<Sphere> spheres;
    size_t H, W;
    double alpha, f;
    Vector center;
    Vector light_source;
    double intensity = 1e10;
};

int main() {
    int W = 512;
    int H = 512;
    float alpha = 60.0 * M_PI / 180.0;

    Vector Q(0, 0, 55);
    Scene scene(Q, H, W, alpha);

    Sphere S(Vector(0, 0, 0), 10.0, Vector(1, 1, 1), true);
    scene.add_sphere(S);
    scene.add_sphere(Sphere(Vector(0.0, -1000, 0.0), 990.0, Vector(0, 0, 1)));
    scene.add_sphere(Sphere(Vector(0.0, 1000, 0.0), 940.0, Vector(1, 0, 0)));
    scene.add_sphere(Sphere(Vector(0.0, 0.0, -1000.0), 940.0, Vector(0, 1, 0)));
    scene.add_sphere(Sphere(Vector(0.0, 0.0, 1000.0), 940.0, Vector(1, 0, 1)));
    scene.add_sphere(Sphere(Vector(1000.0, 0.0, 0.0), 940.0, Vector(1, 1, 0)));
    scene.add_sphere(Sphere(Vector(-1000.0, 0.0, 0.0), 940.0, Vector(0, 1, 1)));
    scene.render();
}
