#ifndef SCENE_H
#define SCENE_H
#include "shapes.h"
#include <cmath>
#include <vector>

class Scene {
    Vector raytrace(const Ray &ray, int ray_depth);
    Vector get_color(const Ray &ray);
    Vector random_cos(const Vector &N);
    std::vector<Shape *> shapes;
    Vector center;
    double f;
    // Features used in the scene
    const bool ANTIALIASING = true;
    const bool INDIRECT_LIGHTING = true;
    const bool DoF = false;
    const bool Fresnel = true;
    // Parameters of the scene
    const double intensity = 3e10;
    const double rays_per_pixel = 100;
    const double max_ray_bounces = 10;
    const double focal_distance = 45;
    const double apperture_radius = 1;
    const size_t W = 512;
    const size_t H = 512;
    const float alpha = 60.0 * M_PI / 180.0;
    const Vector light_source = Vector(-10, 20, 40);

  public:
    Scene(Vector _center);
    void add_shape(Shape *S);
    void render();
};
#endif
