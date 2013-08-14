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
    (distmesh::Functional( \
        [=](const Eigen::Ref<distmesh::dtype::array<distmesh::dtype::real>>& points) -> \
        distmesh::dtype::array<distmesh::dtype::real> \
        function_body))

// namespace distmesh
namespace distmesh {
    // base class of all function expression for allowing easy function arithmetic
    class Functional {
    public:
        // function type of Functional callable
        typedef std::function<dtype::array<dtype::real>(
            const Eigen::Ref<dtype::array<dtype::real>>&)> function_t;

        // create class from function type
        Functional(const function_t& func);

        // copy constructor
        Functional(const Functional& rhs);

        // evaluate function by call
        dtype::array<dtype::real> operator()(
            const Eigen::Ref<dtype::array<dtype::real>>& points) const;

        // basic arithmetic operations
        Functional& operator+() { return *this; }
        Functional operator-();
        Functional& operator+=(const Functional& rhs);
        Functional& operator+=(const dtype::real& rhs);
        Functional& operator-=(const Functional& rhs);
        Functional& operator-=(const dtype::real& rhs);
        Functional& operator*=(const Functional& rhs);
        Functional& operator*=(const dtype::real& rhs);
        Functional& operator/=(const Functional& rhs);
        Functional& operator/=(const dtype::real& rhs);
        friend Functional operator+(const Functional& lhs, const Functional& rhs);
        friend Functional operator+(const Functional& lhs, const dtype::real& rhs);
        friend Functional operator+(const dtype::real& lhs, const Functional& rhs);
        friend Functional operator-(const Functional& lhs, const Functional& rhs);
        friend Functional operator-(const Functional& lhs, const dtype::real& rhs);
        friend Functional operator-(const dtype::real& lhs, const Functional& rhs);
        friend Functional operator*(const Functional& lhs, const Functional& rhs);
        friend Functional operator*(const Functional& lhs, const dtype::real& rhs);
        friend Functional operator*(const dtype::real& lhs, const Functional& rhs);
        friend Functional operator/(const Functional& lhs, const Functional& rhs);
        friend Functional operator/(const Functional& lhs, const dtype::real& rhs);
        friend Functional operator/(const dtype::real& lhs, const Functional& rhs);

        // comparison operations
        Functional min(const Functional& rhs);
        Functional max(const Functional& rhs);

        // enable easier compatibility with rest of distmesh
        operator function_t() { return this->function_; }

        // accessors
        function_t& function() { return this->function_; }
        const function_t& function() const { return this->function_; }

    private:
        // stores std function
        function_t function_;
    };
}

#endif
