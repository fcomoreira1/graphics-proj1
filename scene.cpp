#include "scene.h"
#include "stb/stb_image_write.h"
#include <cmath>

double refractive_index = 1.0 / 1.52;

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

Vector Scene::raytrace(const Ray &r, int ray_depth = 0) {
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
    Vector rho = Vector();
    if (min_t < MAXFLOAT) {
        Vector P = r.get_origin() + (min_t * r.get_dir());
        Vector N = S->get_normal(P);
        double dot_prod = dot(r.get_dir(), N);
        Vector V = light_source - P;
        P = P + 1e-8 * N;
        if (S->is_mirror()) {
            Vector u = r.get_dir() - 2 * dot(r.get_dir(), N) * N;
            return raytrace(Ray(P, u), ray_depth - 1);
        } else if (S->is_transparent()) {
            double r_index =
                (dot_prod < 0) ? refractive_index : 1 / refractive_index;
            double sq = 1 - r_index * r_index * (1 - dot_prod * dot_prod);
            if (sq >= 0) {
                Vector t =
                    r_index * (r.get_dir() - dot_prod * N) - sqrt(sq) * N;
                return raytrace(Ray(P, t), ray_depth - 1);
            } else {
                Vector u = r.get_dir() - 2 * dot_prod * N;
                return raytrace(Ray(P, u), ray_depth - 1);
            }
        }
        double t0;
        bool non_direct_path = false;
        for (auto s : shapes) {
            if (s->intersect(Ray(P, light_source - P), t0)) {
                if (t0 < (light_source - P).norm()) {
                    non_direct_path = true;
                    break;
                }
            }
        }
        if (!non_direct_path) {
            rho = (intensity / (4 * M_PI * V.norm2())) * dot(N, V / V.norm()) *
                  (S->get_color() / M_PI);
        }
    }
    return rho;
}
