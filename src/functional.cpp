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

#include "distmesh/distmesh.h"

// generate new function with difference of two
distmesh::functional::function_t distmesh::functional::diff(
    function_t function1, function_t function2) {
    return DISTMESH_FUNCTIONAL({
        return function1(points).max(-function2(points));
    });
}

// generates new functional with intersect of two
distmesh::functional::function_t distmesh::functional::intersect(
    function_t function1, function_t function2) {
    return DISTMESH_FUNCTIONAL({
        return function1(points).max(function2(points));
    });
}

// generates new functional with union of two
distmesh::functional::function_t distmesh::functional::union_(
    function_t function1, function_t function2) {
    return DISTMESH_FUNCTIONAL({
        return function1(points).min(function2(points));
    });
}
