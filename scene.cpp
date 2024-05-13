#include "scene.h"
#include "random.h"
#define _CRT_SECURE_NO_WARNINGS 1

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "./stb/stb_image.h"
#include <float.h>
#include <iostream>
#include <omp.h>

const double refractive_index = 1.0 / 1.52;

Scene::Scene(Vector _center) : center(_center) {}
void Scene::add_shape(Shape *S) { shapes.emplace_back(S); }
void Scene::render() {
    std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for
    for (size_t i = 0; i < H; i++) {
        for (size_t j = 0; j < W; j++) {
            double z = -(double)W / (2.0 * std::tan(alpha / 2.0));
            Vector u(j + 0.5 - W / 2.0, H / 2.0 - i - 0.5, z);
            Ray r(center, u);
            Vector color = get_color(r);
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

Vector Scene::get_color(const Ray &r) {
    Vector color;
    double sigma = 0.001;
    Ray og_ray = r;
    if (DoF) {
        // Depth of Field
        Vector P = og_ray.get_origin() +
                   (D / abs(og_ray.get_dir().data[2])) * og_ray.get_dir();

        double r = 5e-2 * uniform_distribution(),
               theta = 2.0 * M_PI * uniform_distribution();
        Vector O = og_ray.get_origin() +
                   1 * Vector(r * sin(theta), r * cos(theta), 0.0);
        og_ray = Ray(O, P - O);
    }

    for (size_t i = 0; i < rays_per_pixel; i++) {
        Ray fi_ray = og_ray;
        if (ANTIALIASING) {
            // Antialising
            std::pair<double, double> noise = box_muller(0, sigma);
            fi_ray =
                Ray(og_ray.get_origin(),
                    og_ray.get_dir() + Vector(noise.first, noise.second, 0));
        }
        color = color + raytrace(fi_ray, max_ray_bounces);
    }
    return color / rays_per_pixel;
}

Vector Scene::random_cos(const Vector &N) {
    double r1 = uniform_distribution();
    double r2 = uniform_distribution();
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

Vector Scene::raytrace(const Ray &r, int ray_depth = 0) {
    if (ray_depth < 0) {
        return Vector(0.0, 0.0, 0.0);
    }
    double min_t = DBL_MAX;
    double t;
    Vector N0, N;
    Shape *S;
    for (auto &s : shapes) {
        if (s->intersect(r, t, N0) && t < min_t) {
            S = s;
            min_t = t;
            N = N0;
        }
    }
    if (min_t == DBL_MAX) {
        return Vector(0.0, 0.0, 0.0);
    }
    // Compute the intersection point P and the normal N
    Vector P = r.get_origin() + (min_t * r.get_dir());
    double dot_prod = dot(r.get_dir(), N);
    Vector V = light_source - P;

    // Check for special surfaces
    if (S->is_mirror()) {
        Vector u = r.get_dir() - 2 * dot_prod * N;
        P = P + 1e-8 * N;
        return raytrace(Ray(P, u), ray_depth - 1);
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
            return raytrace(Ray(P, t), ray_depth - 1);
        } else {
            // Total Refraction
            P = P + 1e-8 * N;
            Vector t = r.get_dir() - 2 * dot_prod * N;
            return raytrace(Ray(P, t), ray_depth - 1);
        }
    }

    // Get direct light
    // Avoid self-intersection
    P = P + 1e-8 * N;
    double t0;
    bool non_direct_path = false;
    Vector direct_light;
    for (auto s : shapes) {
        if (s->intersect(Ray(P, light_source - P), t0, N0)) {
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
    Vector indirect_light;
    if (INDIRECT_LIGHTING) {
        Vector omega_i = random_cos(N);
        indirect_light =
            S->get_color() * raytrace(Ray(P, omega_i), ray_depth - 1);
    }
    return direct_light + indirect_light;
}
