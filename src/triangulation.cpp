// --------------------------------------------------------------------
// This file is part of libDistMesh.
//
// libDistMesh is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// libDistMesh is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libDistMesh. If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2015 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#include <stdio.h>

// qhull library used to calculate delaunay triangulation
extern "C" {
    #define qh_QHimport
    #include <qhull/qhull_a.h>
}

#include "distmesh/distmesh.h"
#include "distmesh/triangulation.h"

Eigen::ArrayXXi distmesh::triangulation::delaunay(
    Eigen::Ref<Eigen::ArrayXXd const> const points) {
    // reset qhull
    if (qh_qh) {
        qh_save_qhull();
    }

    // convert points array to row major format
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic,
        Eigen::RowMajor> pointsRowMajor = points;

    // calculate delaunay triangulation
    std::string flags = "qhull d Qt Qbb Qc Qz";
    qh_new_qhull(points.cols(), points.rows(), pointsRowMajor.data(), False,
        (char*)flags.c_str(), nullptr, stderr);
    qh_triangulate();

    // count all upper delaunay facets
    unsigned facetCount = 0;
    facetT* facet;
    FORALLfacets {
        if (!facet->upperdelaunay) {
            facetCount++;
        }
    }

    // extract point ids from delaunay triangulation
    Eigen::ArrayXXi triangulation(facetCount, points.cols() + 1);
    unsigned facetId = 0;
    unsigned vertexId = 0;
    vertexT* vertex, **vertexp;

    FORALLfacets {
        vertexId = 0;
        if (!facet->upperdelaunay) {
            qh_setsize(facet->vertices);
            FOREACHvertex_(facet->vertices) {
                triangulation(facetId, vertexId) = qh_pointid(vertex->point);
                vertexId++;
            }
            facetId++;
        }
    }

    return triangulation;
}
