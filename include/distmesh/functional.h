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
// along with libDistMesh.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2013 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#ifndef LIBDISTMESH_INCLUDE_FUNCTIONAL_H
#define LIBDISTMESH_INCLUDE_FUNCTIONAL_H

// macro for easies creation of distmesh lambda functions
#define DISTMESH_FUNCTIONAL(function_body) \
    ([=](const Eigen::Ref<distmesh::dtype::array<distmesh::dtype::real>>& points) -> \
    distmesh::dtype::array<distmesh::dtype::real> \
    function_body)

// namespace distmesh::functional
namespace distmesh {
namespace functional {
    // function type for edge length functions
    typedef std::function<dtype::array<dtype::real>(
        const Eigen::Ref<dtype::array<dtype::real>>&)> function_t;

    // generate new functional with difference of two
    function_t diff(function_t function1, function_t function2);

    // generates new functional with intersect of two
    function_t intersect(function_t function1, function_t function2);

    // generates new functional with union of two
    function_t union_(function_t function1, function_t function2);
}
}

#endif
