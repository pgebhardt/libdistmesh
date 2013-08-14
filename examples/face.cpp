#include <distmesh/distmesh.h>
#include <fstream>

distmesh::Functional elliptical(
    distmesh::dtype::real f1, distmesh::dtype::real f2, distmesh::dtype::real r,
    distmesh::dtype::real x0, distmesh::dtype::real y0) {
    return DISTMESH_FUNCTIONAL({
        return (f1 * (points.col(0) - x0).square() + f2 * (points.col(1) - y0).square()).sqrt() - r;
    });
}

int main() {
    // bounding box in which the algorithm tries to create points
    distmesh::dtype::array<distmesh::dtype::real> bounding_box(2, 2);
    bounding_box << -1.0, 1.0, -1.0, 1.0;

    // fixed points at the corners of domain to guarantee convergence
    distmesh::dtype::array<distmesh::dtype::real> fixed_points(4, 2);
    fixed_points << -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;

    distmesh::dtype::array<distmesh::dtype::real> midpoints(2, 2);
    midpoints << -0.5, 0.5, 0.5, 0.5;

    auto distance_function =
        distmesh::distance_function::rectangular(bounding_box)
            .max(-distmesh::distance_function::circular(0.25, midpoints.row(0)))
            .max(-distmesh::distance_function::circular(0.25, midpoints.row(1)))
            .max(-elliptical(1.0, 10.0, 0.75, 0.0, -0.5))
            .max(-elliptical(3.0, 1.0, 0.15, 0.0, 0.1));

    // create mesh
    auto mesh = distmesh::distmesh(distance_function,
        0.2 - distance_function,
        0.02, bounding_box, fixed_points);

    // plot mesh
    std::ofstream points_file;
    std::ofstream triangulation_file;
    points_file.open("points.txt");
    triangulation_file.open("triangulation.txt");

    points_file << std::get<0>(mesh);
    triangulation_file << std::get<1>(mesh);

    points_file.close();
    triangulation_file.close();

    system("python plot_mesh.py");

    return 0;
}
