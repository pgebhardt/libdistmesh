// distmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"
#include <stdio.h>

std::shared_ptr<distmesh::dtype::array<distmesh::dtype::index, 2>>
    distmesh::triangulation::delaunay(std::shared_ptr<dtype::array<dtype::real, 2>> points) {
    // set flags for qhull
    std::string flags = "qhull d Qbb Qc Qz";

    // calculate delaunay triangulation
    qh_new_qhull(points->shape()[1], points->shape()[0], points->data(), False,
        (char*)flags.c_str(), nullptr, stderr);
    qh_triangulate();

    // count all upper delaunay facets
    dtype::index facet_count = 0;
    facetT* facet;
    FORALLfacets {
        if (!facet->upperdelaunay) {
            facet_count++;
        }
    }

    // extract point ids from delaunay triangulation
    auto triangulation = std::make_shared<dtype::array<dtype::index, 2>>(
        boost::extents[facet_count][3]);
    dtype::index facet_id = 0;
    dtype::index vertex_id = 0;
    vertexT* vertex, **vertexp;

    FORALLfacets {
        vertex_id = 0;
        if (!facet->upperdelaunay) {
            FOREACHvertex_(facet->vertices) {
                std::cout << qh_pointid(vertex->point) << " ";
                (*triangulation)[facet_id][vertex_id] = qh_pointid(vertex->point);
                vertex_id++;
            }
            std::cout << std::endl;
            facet_id++;
        }
    }

    return triangulation;
}
