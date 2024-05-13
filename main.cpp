#include "mesh.h"
#include "scene.h"
#include "sphere.h"

int main() {
    Vector Q(0, 0, 55);
    Scene scene(Q);

    // Cat!
    TriangleMesh *mesh = new TriangleMesh();
    mesh->readOBJ("cat.obj");
    mesh->transform(0.6, Vector(0, -10, 0), -M_PI / 2.0);
    mesh->set_color(Vector(0.3, 0.2, 0.25));
    // mesh->set_color(Vector(1.0, 1.0, 1.0));
    scene.add_shape(mesh);

    // scene.add_shape(new Sphere(Vector(-20, 0, 0), 10.0, Vector(1, 1, 1),
    // true)); scene.add_shape(new Sphere(Vector(20, 0, 0), 10.0, Vector(1, 1,
    // 1), false,
    //                            true, false));
    // scene.add_shape(
    //     new Sphere(Vector(20, 0, 0), 9.5, Vector(1, 1, 1), false, true,
    //     true));

    // Nice Sphere rendering
    //
    // Center Sphere
    // scene.add_shape(new Sphere(Vector(0, 0, 0), 10.0, Vector(1, 1, 1), false,
    //                            false, false));
    // scene.add_shape(new Sphere(Vector(-10, 0, 20), 10.0, Vector(1, 1, 1),
    // false,
    //                            true, false));
    // scene.add_shape(new Sphere(Vector(-10, 0, 20), 9.5, Vector(1, 1, 1),
    // false,
    //                            true, true));
    // scene.add_shape(new Sphere(Vector(10, 0, -20), 10.0, Vector(1, 1, 1),
    // true,
    //                            false, false));

    // WALLS
    scene.add_shape(
        new Sphere(Vector(0.0, -1000, 0.0), 990.0, Vector(0, 0, 0.8)));
    scene.add_shape(
        new Sphere(Vector(0.0, 1000, 0.0), 940.0, Vector(0.8, 0, 0)));
    scene.add_shape(
        new Sphere(Vector(0.0, 0.0, -1000.0), 940.0, Vector(0, 0.8, 0)));
    scene.add_shape(
        new Sphere(Vector(0.0, 0.0, 1000.0), 940.0, Vector(0.8, 0, 0.8)));
    scene.add_shape(
        new Sphere(Vector(1000.0, 0.0, 0.0), 940.0, Vector(0.8, 0.8, 0)));
    scene.add_shape(
        new Sphere(Vector(-1000.0, 0.0, 0.0), 940.0, Vector(0, 0.8, 0.8)));
    scene.render();
}
