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
    (distmesh::Functional([=](Eigen::Ref<const Eigen::ArrayXXd> points) -> Eigen::ArrayXXd \
        function_body))

// namespace distmesh
namespace distmesh {
    // base class of all function expression for allowing easy function arithmetic
    class Functional {
    public:
        // function type of Functional callable
        typedef std::function<Eigen::ArrayXXd(Eigen::Ref<const Eigen::ArrayXXd>)> function_t;

        // create class from function type
        Functional(const function_t& func) : function_(func) {}
        Functional(double constant) :
            Functional(DISTMESH_FUNCTIONAL({
                return Eigen::ArrayXXd::Constant(points.rows(), 1, constant);
            })) {}

        // copy constructor
        Functional(const Functional& rhs) : function_(rhs.function()) {}
        Functional(Functional&& rhs) : function_(std::move(rhs.function())) {}

        // assignment operator
        Functional& operator=(const Functional& rhs);
        Functional& operator=(Functional&& rhs);

        // evaluate function by call
        Eigen::ArrayXXd operator() (Eigen::Ref<const Eigen::ArrayXXd> points) const;

        // basic arithmetic operations
        Functional& operator+() { return *this; }
        Functional operator-() const;
        Functional& operator+=(const Functional& rhs);
        Functional& operator+=(const double& rhs);
        Functional& operator-=(const Functional& rhs);
        Functional& operator-=(const double& rhs);
        Functional& operator*=(const Functional& rhs);
        Functional& operator*=(const double& rhs);
        Functional& operator/=(const Functional& rhs);
        Functional& operator/=(const double& rhs);
        friend Functional operator+(const Functional& lhs, const Functional& rhs);
        friend Functional operator+(const Functional& lhs, const double& rhs);
        friend Functional operator+(const double& lhs, const Functional& rhs);
        friend Functional operator-(const Functional& lhs, const Functional& rhs);
        friend Functional operator-(const Functional& lhs, const double& rhs);
        friend Functional operator-(const double& lhs, const Functional& rhs);
        friend Functional operator*(const Functional& lhs, const Functional& rhs);
        friend Functional operator*(const Functional& lhs, const double& rhs);
        friend Functional operator*(const double& lhs, const Functional& rhs);
        friend Functional operator/(const Functional& lhs, const Functional& rhs);
        friend Functional operator/(const Functional& lhs, const double& rhs);
        friend Functional operator/(const double& lhs, const Functional& rhs);

        // mathematical methods
        Functional min(const Functional& rhs) const;
        Functional max(const Functional& rhs) const;
        Functional abs() const;

        // accessors
        function_t& function() { return this->function_; }
        const function_t& function() const { return this->function_; }

    private:
        // stores std function
        function_t function_;
    };
}

#endif
