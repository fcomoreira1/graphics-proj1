#include "scene.h"
#include "shapes.h"
#include <cmath>

int main() {
    int W = 512;
    int H = 512;
    float alpha = 60.0 * M_PI / 180.0;

    Vector Q(0, 0, 55);
    Scene scene(Q, H, W, alpha);

    // scene.add_shape(new Sphere(Vector(-20, 0, 0), 10.0, Vector(1, 1, 1),
    // true));
    scene.add_shape(new Sphere(Vector(0, 0, 0), 10.0, Vector(1, 1, 1), false,
                               false, false));
    /* scene.add_shape(new Sphere(Vector(20, 0, 0), 10.0, Vector(1, 1, 1),
    false, true, false)); scene.add_shape( new Sphere(Vector(20, 0, 0), 9.0,
    Vector(1, 1, 1), false, true, true)); */
    scene.add_shape(
        new Sphere(Vector(0.0, -1000, 0.0), 990.0, Vector(0, 0, 1)));
    scene.add_shape(new Sphere(Vector(0.0, 1000, 0.0), 940.0, Vector(1, 0, 0)));
    scene.add_shape(
        new Sphere(Vector(0.0, 0.0, -1000.0), 940.0, Vector(0, 1, 0)));
    scene.add_shape(
        new Sphere(Vector(0.0, 0.0, 1000.0), 940.0, Vector(1, 0, 1)));
    scene.add_shape(
        new Sphere(Vector(1000.0, 0.0, 0.0), 940.0, Vector(1, 1, 0)));
    scene.add_shape(
        new Sphere(Vector(-1000.0, 0.0, 0.0), 940.0, Vector(0, 1, 1)));
    scene.render();
}
