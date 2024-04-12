#include "scene.h"
#include "stb/stb_image_write.h"
#include <cmath>
#include <iostream>

const double refractive_index = 1.0 / 1.52;

Scene::Scene(Vector _center, int _H, int _W, double _alpha,
             Vector _light_source) {
    center = _center;
    H = _H;
    W = _W;
    alpha = _alpha;
    light_source = _light_source;
}
void Scene::add_shape(Shape *S) { shapes.emplace_back(S); }
void Scene::render() {
    std::vector<unsigned char> image(W * H * 3, 0);
    for (size_t i = 0; i < H; i++) {
        for (size_t j = 0; j < W; j++) {
            double z = -(double)W / (2.0 * std::tan(alpha / 2.0));
            Vector u(j + 0.5 - W / 2.0, H / 2.0 - i - 0.5, z);
            Ray r(center, u);
            Vector color = raytrace(r);
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
Vector Scene::raytrace(const Ray &r) {
    Vector color;
    double sigma = 0.001;
    for (size_t i = 0; i < number_of_rays; i++) {
        // Antialising
        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;
        double x = sigma * sqrt(-2.0 * log(r1)) * cos(2 * M_PI * r2);
        double y = sigma * sqrt(-2.0 * log(r1)) * sin(2 * M_PI * r2);
        Vector dir = r.get_dir() + Vector(x, y, 0);
        color = color + pathtrace(Ray(r.get_origin(), dir), max_ray_bounces);
    }
    return color / number_of_rays;
}

Vector Scene::random_direction(const Vector &N) {
    double r1 = (double)rand() / RAND_MAX;
    double r2 = (double)rand() / RAND_MAX;
    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);
    if (std::max(abs(N.data[1]), abs(N.data[0])) > 1e-4) {
        Vector T(-N.data[1], N.data[0], 0);
        T.normalize();
        Vector T2 = cross(N, T);
        return x * T + y * T2 + z * N;
    } else {
        Vector T(0, -N.data[2], N.data[1]);
        T.normalize();
        Vector T2 = cross(N, T);
        return x * T + y * T2 + z * N;
    }
}

Vector Scene::pathtrace(const Ray &r, int ray_depth = 0) {
    if (ray_depth < 0) {
        return Vector(0.0, 0.0, 0.0);
    }
    double min_t = MAXFLOAT;
    double t;
    Shape *S;
    for (auto &s : shapes) {
        if (s->intersect(r, t) && t < min_t) {
            S = s;
            min_t = t;
        }
    }
    if (min_t == MAXFLOAT) {
        return Vector(0.0, 0.0, 0.0);
    }
    // Compute the intersection point P and the normal N
    Vector P = r.get_origin() + (min_t * r.get_dir());
    Vector N = S->get_normal(P);
    double dot_prod = dot(r.get_dir(), N);
    Vector V = light_source - P;

    // Check for special surfaces
    if (S->is_mirror()) {
        Vector u = r.get_dir() - 2 * dot_prod * N;
        P = P + 1e-8 * N;
        return pathtrace(Ray(P, u), ray_depth - 1);
    } else if (S->is_transparent()) {
        double r_index = refractive_index;
        if (dot_prod > 0) {
            r_index = 1.0 / r_index;
            N = -1.0 * N;
            dot_prod = -dot_prod;
        }
        double sq = 1 - r_index * r_index * (1 - dot_prod * dot_prod);
        if (sq >= 0) {
            P = P - 1e-8 * N;
            Vector t = r_index * (r.get_dir() - dot_prod * N) - sqrt(sq) * N;
            return pathtrace(Ray(P, t), ray_depth - 1);
        } else {
            P = P + 1e-8 * N;
            Vector t = r.get_dir() - 2 * dot_prod * N;
            return pathtrace(Ray(P, t), ray_depth - 1);
        }
    }

    // Get direct light
    // Avoid self-intersection
    P = P + 1e-8 * N;
    double t0;
    bool non_direct_path = false;
    Vector direct_light;
    for (auto s : shapes) {
        if (s->intersect(Ray(P, light_source - P), t0)) {
            if (t0 < (light_source - P).norm()) {
                non_direct_path = true;
                break;
            }
        }
    }
    // Check for non-direct path
    if (!non_direct_path) {
        direct_light = (intensity / (4 * M_PI * V.norm2())) *
                       dot(N, V / V.norm()) * (S->get_color() / M_PI);
    }

    // Get indirect light
    Vector omega_i = random_direction(N);
    Vector indirect_light =
        S->get_color() * pathtrace(Ray(P, omega_i), ray_depth - 1);
    return direct_light + indirect_light;
}
