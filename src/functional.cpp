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
// Copyright (C) 2015 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#include "distmesh/distmesh.h"

// assignment operator
distmesh::Functional& distmesh::Functional::operator=(
    Functional const& rhs) {
    this->function() = rhs.function();
    return *this;
}
distmesh::Functional& distmesh::Functional::operator=(
    Functional&& rhs) {
    this->function() = std::move(rhs.function());
    return *this;
}

Eigen::ArrayXd distmesh::Functional::operator()(
    Eigen::Ref<Eigen::ArrayXXd const> const points) const {
    return this->function()(points);
}

distmesh::Functional distmesh::Functional::operator-() const {
    function_t func = this->function();;
    return DISTMESH_FUNCTIONAL({
        return -func(points);
    });
}

distmesh::Functional& distmesh::Functional::operator+=(
    Functional const& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) + rhs(points);
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator+=(
    double const rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) + rhs;
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator-=(
    Functional const& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) - rhs(points);
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator-=(
    double const rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) - rhs;
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator*=(
    Functional const& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) * rhs(points);
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator*=(
    double const rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) * rhs;
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator/=(
    Functional const& rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) / rhs(points);
    });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator/=(
    double const rhs) {
    function_t func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) / rhs;
    });
    return *this;
}

distmesh::Functional distmesh::operator+(
    Functional const& lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) + rhs(points);
    });
}

distmesh::Functional distmesh::operator+(
    Functional const& lhs, double const rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) + rhs;
    });
}

distmesh::Functional distmesh::operator+(
    double const lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs + rhs(points);
    });
}

distmesh::Functional distmesh::operator-(
    Functional const& lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) - rhs(points);
    });
}

distmesh::Functional distmesh::operator-(
    Functional const& lhs, double const rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) - rhs;
    });
}

distmesh::Functional distmesh::operator-(
    double const lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs - rhs(points);
    });
}

distmesh::Functional distmesh::operator*(
    Functional const& lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) * rhs(points);
    });
}

distmesh::Functional distmesh::operator*(
    Functional const& lhs, double const rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) * rhs;
    });
}

distmesh::Functional distmesh::operator*(
    double const lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs * rhs(points);
    });
}

distmesh::Functional distmesh::operator/(
    Functional const& lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) / rhs(points);
    });
}

distmesh::Functional distmesh::operator/(
    Functional const& lhs, double const rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) / rhs;
    });
}

distmesh::Functional distmesh::operator/(
    double const lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs / rhs(points);
    });
}

distmesh::Functional distmesh::Functional::min(
    Functional const& rhs) const {
    Functional res(this->function());
    return DISTMESH_FUNCTIONAL({
        return res(points).min(rhs(points));
    });
}

distmesh::Functional distmesh::Functional::max(
    Functional const& rhs) const {
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

// geometric transform
distmesh::Functional distmesh::Functional::shift(Eigen::Ref<Eigen::ArrayXd const> const offset) const {
    Functional res(this->function());
    return DISTMESH_FUNCTIONAL({
        return res(points.rowwise() - offset.transpose());
    });
}

distmesh::Functional distmesh::Functional::rotate2D(double const angle) const {
    Functional res(this->function());
    return DISTMESH_FUNCTIONAL({
        Eigen::ArrayXXd transformedPoints(points.rows(), points.cols());
        transformedPoints.col(0) = points.col(0) * std::cos(angle) + points.col(1) * std::sin(angle);
        transformedPoints.col(1) = -points.col(0) * std::sin(angle) + points.col(1) * std::cos(angle);

        return res(transformedPoints);
    });
}

