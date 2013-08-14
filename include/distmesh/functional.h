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
    (distmesh::functional::Function( \
        [=](const Eigen::Ref<distmesh::dtype::array<distmesh::dtype::real>>& points) -> \
        distmesh::dtype::array<distmesh::dtype::real> \
        function_body))

// namespace distmesh::functional
namespace distmesh {
namespace functional {
    // function type for edge length functions
    typedef std::function<dtype::array<dtype::real>(
        const Eigen::Ref<dtype::array<dtype::real>>&)> function_t;

    // base class of all function expression for allowing easy function arithmetic
    class Function {
    public:
        // create class from function type
        Function(const function_t& func);

        // copy constructor
        Function(const Function& rhs);

        // evaluate function by call
        dtype::array<dtype::real> operator()(
            const Eigen::Ref<dtype::array<dtype::real>>& points) const;

        // basic arithmetic operations
        Function& operator+() { return *this; }
        Function operator-();
        Function& operator+=(const Function& rhs);
        Function& operator-=(const Function& rhs);
        Function& operator*=(const Function& rhs);
        friend Function operator+(const Function& lhs, const Function& rhs);
        friend Function operator-(const Function& lhs, const Function& rhs);
        friend Function operator*(const Function& lhs, const Function& rhs);

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
}

#endif
