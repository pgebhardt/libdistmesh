// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"
#include <vector>
#include <set>
#include <algorithm>
#include <sys/time.h>

class Time {
public:
    Time() {
        this->restart();
    }
    double elapsed() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return ((double)tv.tv_sec + (double)tv.tv_usec * 1e-6) - this->time;
    }
    void restart() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        this->time = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
    }

private:
    double time;
};

// apply the distmesh algorithm
std::tuple<std::shared_ptr<distmesh::dtype::array<distmesh::dtype::real>>,
    std::shared_ptr<distmesh::dtype::array<distmesh::dtype::index>>>
    distmesh::distmesh(
    std::function<dtype::real(dtype::array<dtype::real>)> distance_function,
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> edge_length_function,
    dtype::real initial_edge_length, dtype::array<dtype::real> bounding_box) {
    Time time;
    // create initial distribution in bounding_box
    auto points = utils::create_point_list(distance_function,
        edge_length_function, initial_edge_length, bounding_box);

    // create initial triangulation
    auto triangulation = triangulation::delaunay(points);

    // create points buffer for retriangulation and stop criterion
    dtype::array<dtype::real> buffer_retriangulation_criterion(points->rows(), points->cols());
    dtype::array<dtype::real> buffer_stop_criterion(points->rows(), points->cols());
    buffer_retriangulation_criterion.fill(INFINITY);

    // main distmesh loop
    std::shared_ptr<dtype::array<dtype::index>> bar_indices = nullptr;
//    for (dtype::index i = 0; i < 1; ++i) {
    while (true) {
        time.restart();
        // retriangulate if point movement is above tolerance
        auto retriangulation_criterion =
            (((*points) - buffer_retriangulation_criterion).square().rowwise().sum().sqrt() /
            initial_edge_length).maxCoeff();
        if (retriangulation_criterion > settings::retriangulation_tolerance) {
            buffer_retriangulation_criterion = *points;

            // reject triangles with circumcenter outside of the region
            triangulation = triangulation::delaunay(points);
            dtype::array<dtype::index> keep_triangulation(triangulation->rows(), triangulation->cols());
            dtype::index triangle_count = 0;
            for (dtype::index triangle = 0; triangle < triangulation->rows(); ++triangle) {
                // calculate circumcenter
                dtype::array<dtype::real> circumcenter(1, points->cols());
                circumcenter.fill(0.0);
                for (dtype::index point = 0; point < triangulation->cols(); ++point) {
                    circumcenter += points->row((*triangulation)(triangle, point)) / triangulation->cols();
                }
                if (distance_function(circumcenter) < -settings::general_precision * initial_edge_length) {
                    keep_triangulation.row(triangle_count) = triangulation->row(triangle);
                    triangle_count++;
                }
            }
            *triangulation = keep_triangulation;
            triangulation->conservativeResize(triangle_count, triangulation->cols());

            // find unique bar indices
            bar_indices = utils::find_unique_bars(points, triangulation);

            std::cout << "retriangulation: " << time.elapsed() * 1e3 << " ms" << std::endl;
            time.restart();
        }

        // calculate bar vector
        dtype::array<dtype::real> bar_vector(bar_indices->rows(), points->cols());
        bar_vector.fill(0.0);
        for (dtype::index bar = 0; bar < bar_indices->rows(); ++bar) {
            bar_vector.row(bar) += points->row((*bar_indices)(bar, 0)) -
                points->row((*bar_indices)(bar, 1));
        }

        std::cout << "calculate bar vector: " << time.elapsed() * 1e3 << " ms" << std::endl;
        time.restart();

        // calculate length of each bar
        dtype::array<dtype::real> bar_length(bar_indices->rows(), 1);
        bar_length = bar_vector.square().rowwise().sum().sqrt();

        std::cout << "calculate length of each bar: " << time.elapsed() * 1e3 << " ms" << std::endl;
        time.restart();

        // evaluate edge_length_function
        dtype::array<dtype::real> hbars(bar_indices->rows(), 1);
        dtype::array<dtype::real> bar_midpoints(bar_indices->rows(), points->cols());
        for (dtype::index bar = 0; bar < bar_indices->rows(); ++bar) {
            bar_midpoints.row(bar) = 0.5 * (points->row((*bar_indices)(bar, 0)) +
                points->row((*bar_indices)(bar, 1)));
        }
        hbars = edge_length_function(bar_midpoints);

        std::cout << "evaluate edge_length_function: " << time.elapsed() * 1e3 << " ms" << std::endl;
        time.restart();

        // calculate desired bar length
        dtype::array<dtype::real> desired_bar_length(bar_length.rows(), bar_length.cols());
        desired_bar_length = hbars * (1.0 + 0.4 / std::pow(2.0, points->cols() - 1)) *
            std::pow((bar_length.pow(points->cols()).sum() / hbars.pow(points->cols()).sum()), 1.0 / points->cols());

        std::cout << "calculate desired bar length: " << time.elapsed() * 1e3 << " ms" << std::endl;
        time.restart();

        // calculate force vector for each bar
        dtype::array<dtype::real> force(bar_indices->rows(), 1);
        dtype::array<dtype::real> force_vector(bar_indices->rows(), points->cols());
        force = ((desired_bar_length - bar_length) / bar_length).max(0.0);
        for (dtype::index dim = 0; dim < points->cols(); ++dim) {
            force_vector.col(dim) = force * bar_vector.col(dim);
        }

        std::cout << "calculate force vector for each bar: " << time.elapsed() * 1e3 << " ms" << std::endl;
        time.restart();

        // move points
        buffer_stop_criterion = *points;
        for (dtype::index bar = 0; bar < bar_indices->rows(); ++bar) {
            points->row((*bar_indices)(bar, 0)) += settings::deltaT *
                force_vector.row(bar);
            points->row((*bar_indices)(bar, 1)) -= settings::deltaT *
                force_vector.row(bar);
        }

        std::cout << "move points: " << time.elapsed() * 1e3 << " ms" << std::endl;
        time.restart();

        // project points outside of domain to boundary
        utils::project_points_to_function(distance_function,
            initial_edge_length, points);

        std::cout << "project points outside of domain to boundary: " << time.elapsed() * 1e3 << " ms" << std::endl;
        time.restart();

        // stop criterion
        auto stop_criterion = ((*points - buffer_stop_criterion).square().rowwise().sum().sqrt() / initial_edge_length).maxCoeff();
        if (stop_criterion < settings::point_movement_tolerance) {
            break;
        }
        std::cout << "stop criterion: " << time.elapsed() * 1e3 << " ms" << std::endl;
        time.restart();
    }

    return std::make_tuple(points, triangulation);
}

// determine boundary edges of given triangulation
std::shared_ptr<distmesh::dtype::array<distmesh::dtype::index>>
    distmesh::boundedges(
    std::shared_ptr<dtype::array<dtype::index>> triangulation) {
    std::set<std::vector<dtype::index>> edge_set;
    std::vector<std::vector<dtype::index>> boundary_edges;
    std::vector<dtype::index> facet(triangulation->cols());
    std::vector<dtype::index> edge(triangulation->cols() - 1);

    // find edges, which only appear once in triangulation
    for (dtype::index triangle = 0; triangle < triangulation->rows(); ++triangle) {
        // get current facet
        for (dtype::index vertex = 0; vertex < triangulation->cols(); ++vertex) {
            facet[vertex] = (*triangulation)(triangle, vertex);
        }
        std::sort(facet.begin(), facet.end());

        // check appearance for all permutation and take care of same permutation
        std::set<std::vector<dtype::index>> permutation_set;
        do {
            // use the first vertices of facet as edge
            for (dtype::index vertex = 0; vertex < triangulation->cols() - 1; ++vertex) {
                edge[vertex] = facet[vertex];
            }
            std::sort(edge.begin(), edge.end());
            if (!std::get<1>(permutation_set.insert(edge))) {
                continue;
            }

            // insert edge in set to get info about multiple appearance
            bool appearance = std::get<1>(edge_set.insert(edge));
            if (!appearance) {
                // find edge in vector and delete it
                auto it = std::find(boundary_edges.begin(), boundary_edges.end(), edge);
                if (it != boundary_edges.end()) {
                    boundary_edges.erase(it);
                }
            } else {
                boundary_edges.push_back(edge);
            }
        } while (std::next_permutation(facet.begin(), facet.end()));
    }

    // convert stl vector to eigen array
    auto boundary_array = std::make_shared<dtype::array<dtype::index>>(
        boundary_edges.size(), edge.size());
    for (dtype::index edge = 0; edge < boundary_array->rows(); ++edge)
    for (dtype::index vertex = 0; vertex < boundary_array->cols(); ++vertex) {
        (*boundary_array)(edge, vertex) = boundary_edges[edge][vertex];
    }

    return boundary_array;
}
