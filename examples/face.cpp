#include <distmesh/distmesh.h>
#include <fstream>

int main() {
    // bounding box in which the algorithm tries to create points
    distmesh::dtype::array<distmesh::dtype::real> bounding_box(2, 2);
    bounding_box << -1.0, 1.0, -1.0, 1.0;

    // fixed points at the corners of domain to guarantee convergence
    distmesh::dtype::array<distmesh::dtype::real> fixed_points(4, 2);
    fixed_points << -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;

    auto distance_function = DISTMESH_FUNCTIONAL({
        auto out1 = points.col(0) + 1.0;
        auto out2 = -(points.col(0) - 1.0);
        auto out3 =  (points.col(1) + 1.0);
        auto out4 = -(points.col(1) - 1.0);
        auto eye1 = ((points.col(0) + 0.5).square() + (points.col(1) - 0.5).square()).sqrt() - 0.25;
        auto eye2 = ((points.col(0) - 0.5).square() + (points.col(1) - 0.5).square()).sqrt() - 0.25;
        auto mouth = (points.col(0).square() + 10.0 * (points.col(1) + 0.5).square()).sqrt() - 0.75;
        auto nose = (3.0 * points.col(0).square() + (points.col(1) - 0.1).square()).sqrt() - 0.15;

        return -out1.min(out2).min(out3).min(out4).min(eye1).min(eye2).min(mouth).min(nose);
    });

    // create mesh
    auto mesh = distmesh::distmesh(distance_function,
        DISTMESH_FUNCTIONAL({
            return 0.2 - distance_function(points);
        }), 0.02, bounding_box, fixed_points);

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
