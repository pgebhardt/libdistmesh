#include <distmesh/distmesh.h>
#include <fstream>

int main() {
    // midpoints of ears
    distmesh::dtype::array<distmesh::dtype::real> midpoints(2, 2);
    midpoints << -0.8, 0.8, 0.8, 0.8;

    // distance function as union of 3 circles
    auto distance_function = distmesh::distance_function::circular()
        .min(distmesh::distance_function::circular(0.5, midpoints.row(0)))
        .min(distmesh::distance_function::circular(0.5, midpoints.row(1)));

    // create mesh
    auto mesh = distmesh::distmesh(distance_function, 0.05,
        0.5 - distance_function, distmesh::bounding_box(2) * 2.0);

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
