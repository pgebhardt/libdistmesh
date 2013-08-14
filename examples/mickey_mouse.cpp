#include <distmesh/distmesh.h>
#include <fstream>

int main() {
    // bounding box in which the algorithm tries to create points
    distmesh::dtype::array<distmesh::dtype::real> bounding_box(2, 2);
    bounding_box << -1.0, 1.0, -1.0, 1.0;

    // fixed points at the corners of domain to guarantee convergence
    distmesh::dtype::array<distmesh::dtype::real> fixed_points(4, 2);
    fixed_points << -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;

    distmesh::dtype::array<distmesh::dtype::real> midpoints(2, 2);
    midpoints << -0.8, 0.8, 0.8, 0.8;

    auto distance_function = distmesh::distance_function::circular()
        .min(distmesh::distance_function::circular(0.5, midpoints.row(0)))
        .min(distmesh::distance_function::circular(0.5, midpoints.row(1)));

    // create mesh
    auto mesh = distmesh::distmesh(distance_function,
        0.5 - distance_function, 0.05, bounding_box);

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
