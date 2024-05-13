#include "mesh.h"
#include <cinttypes>
#include <limits>
#include <list>
#include <utility>

/**
 * @brief Computes the coefficient t in which the ray r intersects the plane
 * defined by normal vector N and point A.
 */
double BoundingBox::intersect_plane(const Ray &r, const Vector &N,
                                    const Vector &A) {
    return dot(A - r.get_origin(), N) / dot(r.get_dir(), N);
}

bool BoundingBox::intersect(const Ray &r) {
    double txmin = intersect_plane(r, Vector(1, 0, 0), Bmin);
    double tymin = intersect_plane(r, Vector(0, 1, 0), Bmin);
    double tzmin = intersect_plane(r, Vector(0, 0, 1), Bmin);
    double txmax = intersect_plane(r, Vector(1, 0, 0), Bmax);
    double tymax = intersect_plane(r, Vector(0, 1, 0), Bmax);
    double tzmax = intersect_plane(r, Vector(0, 0, 1), Bmax);

    if (txmin > txmax)
        std::swap(txmin, txmax);
    if (tymin > tymax)
        std::swap(tymin, tymax);
    if (tzmin > tzmax)
        std::swap(tzmin, tzmax);

    double tmin = std::max(txmin, std::max(tymin, tzmin));
    double tmax = std::min(txmax, std::min(tymax, tzmax));
    return tmin < tmax && 0 < tmax;
}

BVHNode::BVHNode(const std::vector<Vector> &vertices,
                 std::vector<TriangleIndices> &indices, int start, int end) {
    left = NULL;
    right = NULL;
    starting_triangle = start;
    ending_triangle = end;
    bbox = BoundingBox(vertices, indices, start, end);
    Vector box_center = 0.5 * (bbox.Bmax + bbox.Bmin);
    int longest_axis = (bbox.Bmax - bbox.Bmin).max_arg();
    int pivot_index = start;
    for (int i = start; i < end; ++i) {
        TriangleIndices t = indices[i];
        Vector barycenter =
            (vertices[t.vtxi] + vertices[t.vtxj] + vertices[t.vtxk]) / 3.0;

        if (barycenter[longest_axis] < box_center[longest_axis]) {
            std::swap(indices[i], indices[pivot_index]);
            pivot_index++;
        }
    }
    // std::cout << start << " - " << pivot_index << " - " << end << "; "
    //           << longest_axis << " -> " << bbox.Bmax.data[0] << " "
    //           << bbox.Bmin.data[0] << " " << bbox.Bmax.data[1] << " "
    //           << bbox.Bmin.data[1] << " " << bbox.Bmax.data[2] << " "
    //           << bbox.Bmin.data[2] << std::endl;
    if (pivot_index <= start || pivot_index >= end - 1 || end - start < 5) {
        return;
    }
    left = new BVHNode(vertices, indices, start, pivot_index);
    right = new BVHNode(vertices, indices, pivot_index, end);
}

void BVHNode::intersect(const Ray &r, std::vector<int> &reduced_indices) {
    if (!bbox.intersect(r))
        return;
    if (left == NULL && right == NULL) {
        for (size_t i = starting_triangle; i < ending_triangle; ++i) {
            reduced_indices.emplace_back(i);
        }
        return;
    }
    left->intersect(r, reduced_indices);
    right->intersect(r, reduced_indices);
}

bool TriangleMesh::intersect(const Ray &r, double &t, Vector &N) {
    Vector O = r.get_origin(), u = r.get_dir();
    double _t, alpha, beta, gamma, denom;
    Vector A, B, C, _N;
    t = std::numeric_limits<double>::max();
    std::vector<int> reduced_indices;
    if (do_BVH) {
        BVHroot->intersect(r, reduced_indices);
    } else {
        for (size_t i = 0; i < indices.size(); i++)
            reduced_indices.emplace_back(i);
    }
    for (int i : reduced_indices) {
        auto triang = indices[i];
        A = vertices[triang.vtxi];
        B = vertices[triang.vtxj];
        C = vertices[triang.vtxk];
        _N = cross(B - A, C - A);
        denom = dot(u, _N);
        beta = dot(C - A, cross(A - O, u)) / denom;
        gamma = -dot(B - A, cross(A - O, u)) / denom;
        alpha = 1 - beta - gamma;
        _t = dot(A - O, _N) / denom;
        if (0 < _t && _t < t && 0 <= alpha && 0 <= beta && 0 <= gamma) {
            t = _t;
            if (Phong_interpolation) {
                N = alpha * normals[triang.ni] + beta * normals[triang.nj] +
                    gamma * normals[triang.nk];
                N.normalize();
            } else {
                _N.normalize();
                N = _N;
            }
        }
    }
    return (t < std::numeric_limits<double>::max());
};
void TriangleMesh::transform(double scale_factor, Vector translation,
                             double rotation_angle) {
    // Vectors for rotation matrix (rotation on y-axis)
    Vector row1(cos(rotation_angle), 0, sin(rotation_angle));
    Vector row3(-sin(rotation_angle), 0, cos(rotation_angle));
    for (auto &v : vertices) {
        v = scale_factor * v + translation;
        v = Vector(dot(row1, v), v.data[1], dot(row3, v));
    }
    createBVH();
}
