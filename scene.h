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
    Vector raytrace(const Ray &ray, int ray_depth);
    std::vector<Shape *> shapes;
    Vector center, light_source;
    size_t H, W;
    double alpha, f, intensity = 1e10;
};
#endif
