#include "mesh.h"
#include "scene.h"
#include "sphere.h"
#include <chrono>

int main() {
    Vector Q(0, 0, 55);
    Scene scene(Q);

    // Cat!
    TriangleMesh *mesh = new TriangleMesh();
    mesh->readOBJ("cat.obj");
    mesh->transform(0.6, Vector(0, -10, 0), -6.0 * M_PI / 16.0);
    mesh->transform(1, Vector(-8, 0, 6));
    mesh->set_color(Vector(0.6, 0.4, 0.5));
    // mesh->set_transparent(true);
    // mesh->set_mirror(true);
    scene.add_shape(mesh);

    // Center Sphere
    scene.add_shape(new Sphere(Vector(10, -2.0, 10), 8.0, Vector(0.2, 0.4, 0.1),
                               false, true, false));
    scene.add_shape(new Sphere(Vector(10, -2.0, 10), 7.7, Vector(0.2, 0.4, 0.1),
                               false, true, true));
    // scene.add_shape(
    //     new Sphere(Vector(0, -5.0, 25), 5.0, Vector(0.4, 0.1, 0.3)));

    // scene.add_shape(new Sphere(Vector(-20, 0, 0), 10.0, Vector(1, 1, 1),
    // true)); scene.add_shape(new Sphere(Vector(20, 0, 0), 10.0, Vector(1, 1,
    // 1), false,
    //                            true, false));
    // scene.add_shape(
    //     new Sphere(Vector(20, 0, 0), 9.5, Vector(1, 1, 1), false, true,
    //     true));

    // Nice Sphere rendering
    //
    // scene.add_shape(new Sphere(Vector(-10, 0, 20), 10.0, Vector(1, 1, 1),
    // false,
    //                            true, false));
    // scene.add_shape(new Sphere(Vector(-10, 0, 20), 9.5, Vector(1, 1, 1),
    // false,
    //                            true, true));
    // scene.add_shape(new Sphere(Vector(10, 0, -20), 10.0, Vector(1, 1, 1),
    // true,
    //                            false, false));

    // Nice DoF
    // scene.add_shape(
    //     new Sphere(Vector(0, 0, 0), 7.0, Vector(1, 1, 1), false, false,
    //     false));
    // scene.add_shape(new Sphere(Vector(-14, 0, 15), 7.0, Vector(0, 0.5, 0.5),
    //                            false, false, false));
    // scene.add_shape(new Sphere(Vector(-10, 0, 20), 9.5, Vector(1, 1, 1),
    // false,
    //                            true, true));
    // scene.add_shape(new Sphere(Vector(20, 0, -15), 7.0, Vector(0.5, 0, 0.2),
    //                            false, false, false));

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

    auto start = std::chrono::high_resolution_clock::now();
    scene.render();
    auto end = std::chrono::high_resolution_clock::now();
    auto time_elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << "Time elapsed: " << time_elapsed << "s" << std::endl;
}
