// distmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include <stdio.h>
#include "distmesh/distmesh.h"

// qhull library used to calculate delaunay triangulation
extern "C" {
#define qh_QHimport
#include <qhull/qhull_a.h>
}

std::shared_ptr<distmesh::dtype::array<distmesh::dtype::index>>
    distmesh::triangulation::delaunay(
    std::shared_ptr<dtype::array<dtype::real>> points) {
    // set flags for qhull
    std::string flags = "qhull d Qbb Qc Qz";

    // calculate delaunay triangulation
    qh_new_qhull(points->cols(), points->rows(), points->data(), False,
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
    auto triangulation = std::make_shared<dtype::array<dtype::index>>(
        facet_count, points->cols() + 1);
    dtype::index facet_id = 0;
    dtype::index vertex_id = 0;
    vertexT* vertex, **vertexp;

    FORALLfacets {
        vertex_id = 0;
        if (!facet->upperdelaunay) {
            FOREACHvertex_(facet->vertices) {
                (*triangulation)(facet_id, vertex_id) = qh_pointid(vertex->point);
                vertex_id++;
            }
            facet_id++;
        }
    }

    return triangulation;
}
