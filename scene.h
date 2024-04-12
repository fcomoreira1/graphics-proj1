#ifndef SCENE_H
#define SCENE_H
#include "shapes.h"
#include <vector>

class Scene {
  public:
    Scene(Vector _center, int _H, int _W, double _alpha,
          Vector _light_source = Vector(-10, 20, 40));
    void add_shape(Shape *S);
    void render();

  private:
    Vector pathtrace(const Ray &ray, int ray_depth);
    Vector raytrace(const Ray &ray);
    Vector random_direction(const Vector &N);
    std::vector<Shape *> shapes;
    Vector center, light_source;
    size_t H, W;
    double alpha, f, intensity = 1e10;
    double number_of_rays = 32;
    double max_ray_bounces = 3;
};
#endif
