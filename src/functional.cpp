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

distmesh::functional::Function::Function(const function_t& func)
    : function_(func) { }

distmesh::functional::Function::Function(const Function& rhs)
    : function_(rhs.function()) { }

distmesh::dtype::array<distmesh::dtype::real> distmesh::functional::Function::operator()(
    const Eigen::Ref<dtype::array<dtype::real>>& points) const {
    return this->function()(points);
}

distmesh::functional::Function distmesh::functional::Function::operator-() {
    Function func(*this);
    return DISTMESH_FUNCTIONAL({
        return -func(points);
    });
}

distmesh::functional::Function& distmesh::functional::Function::operator+=(
    const Function& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points).max(rhs(points));
    });
    return *this;
}

distmesh::functional::Function& distmesh::functional::Function::operator-=(
    const Function& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points).max(-rhs(points));
    });
    return *this;
}

distmesh::functional::Function& distmesh::functional::Function::operator*=(
    const Function& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points).min(rhs(points));
    });
    return *this;
}

distmesh::functional::Function distmesh::functional::operator+(
    const Function& lhs, const Function& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points).max(rhs(points));
    });
}

distmesh::functional::Function distmesh::functional::operator-(
    const Function& lhs, const Function& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points).max(-rhs(points));
    });
}

distmesh::functional::Function distmesh::functional::operator*(
    const Function& lhs, const Function& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points).min(rhs(points));
    });
}
