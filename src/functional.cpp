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

// assignment operator
distmesh::Functional& distmesh::Functional::operator=(
    const Functional& rhs) {
    this->function() = rhs.function();
    return *this;
}
distmesh::Functional& distmesh::Functional::operator=(
    Functional&& rhs) {
    this->function() = std::move(rhs.function());
    return *this;
}

distmesh::dtype::array<distmesh::dtype::real> distmesh::Functional::operator()(
    const Eigen::Ref<dtype::array<dtype::real>>& points) const {
    return this->function()(points);
}

distmesh::Functional distmesh::Functional::operator-() const {
    function_t func = this->function();;
    return DISTMESH_FUNCTIONAL({
        return -func(points);
    });
}

distmesh::Functional& distmesh::Functional::operator+=(
    const Functional& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) + rhs(points);
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator+=(
    const dtype::real& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) + rhs;
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator-=(
    const Functional& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) - rhs(points);
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator-=(
    const dtype::real& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) - rhs;
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator*=(
    const Functional& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) * rhs(points);
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator*=(
    const dtype::real& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) * rhs;
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator/=(
    const Functional& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) / rhs(points);
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator/=(
    const dtype::real& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) / rhs;
    });
    return *this;
}

distmesh::Functional distmesh::operator+(
    const Functional& lhs, const Functional& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) + rhs(points);
    });
}

distmesh::Functional distmesh::operator+(
    const Functional& lhs, const dtype::real& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) + rhs;
    });
}

distmesh::Functional distmesh::operator+(
    const dtype::real& lhs, const Functional& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs + rhs(points);
    });
}

distmesh::Functional distmesh::operator-(
    const Functional& lhs, const Functional& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) - rhs(points);
    });
}

distmesh::Functional distmesh::operator-(
    const Functional& lhs, const dtype::real& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) - rhs;
    });
}

distmesh::Functional distmesh::operator-(
    const dtype::real& lhs, const Functional& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs - rhs(points);
    });
}

distmesh::Functional distmesh::operator*(
    const Functional& lhs, const Functional& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) * rhs(points);
    });
}

distmesh::Functional distmesh::operator*(
    const Functional& lhs, const dtype::real& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) * rhs;
    });
}

distmesh::Functional distmesh::operator*(
    const dtype::real& lhs, const Functional& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs * rhs(points);
    });
}

distmesh::Functional distmesh::operator/(
    const Functional& lhs, const Functional& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) / rhs(points);
    });
}

distmesh::Functional distmesh::operator/(
    const Functional& lhs, const dtype::real& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) / rhs;
    });
}

distmesh::Functional distmesh::operator/(
    const dtype::real& lhs, const Functional& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs / rhs(points);
    });
}

distmesh::Functional distmesh::Functional::min(
    const Functional& rhs) const {
    Functional res(this->function());
    return DISTMESH_FUNCTIONAL({
        return res(points).min(rhs(points));
    });
}

distmesh::Functional distmesh::Functional::max(
    const Functional& rhs) const {
    Functional res(this->function());
    return DISTMESH_FUNCTIONAL({
        return res(points).max(rhs(points));
    });
}

distmesh::Functional distmesh::Functional::abs() const {
    Functional res(this->function());
    return DISTMESH_FUNCTIONAL({
        return res(points).abs();
    });
}
