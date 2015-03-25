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

#ifndef LIBDISTMESH_INCLUDE_UTILS_H
#define LIBDISTMESH_INCLUDE_UTILS_H

namespace distmesh {
namespace utils {
    // select array elements based on mask
    template <
        class type
    >
    Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> selectMaskedArrayElements(
        Eigen::Ref<Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> const> const array,
        Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> const> const mask) {
        Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> result(mask.count(), array.cols());

        int resultCount = 0;
        for (int row = 0; row < array.rows(); ++row) {
            if (mask(row, 0)) {
                result.row(resultCount) = array.row(row);
                resultCount++;
            }
        }

        return result;
    }

    // select array elements based on indices
    template <
        class type
    >
    Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> selectIndexedArrayElements(
        Eigen::Ref<Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> const> const array,
        Eigen::Ref<Eigen::ArrayXXi const> const indices) {
        Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> result(indices.rows(), array.cols());

        for (int row = 0; row < indices.rows(); ++row) {
            result.row(row) = array.row(indices(row, 0));
        }

        return result;
    }

    // calculate factorial recursively
    unsigned factorial(unsigned const n);

    // create initial points distribution
    Eigen::ArrayXXd createInitialPoints(Functional const& distanceFunction,
        double const baseEdgeLength, Functional const& edgeLengthFunction,
        Eigen::Ref<Eigen::ArrayXXd const> const boundingBox,
        Eigen::Ref<Eigen::ArrayXXd const> const fixedPoints);

    // create array with all unique combinations n over k
    Eigen::ArrayXXi nOverK(unsigned const n, unsigned const k);

    // find unique bars
    Eigen::ArrayXXi findUniqueBars(Eigen::Ref<Eigen::ArrayXXi const> const triangulation);

    // project points outside of boundary back to it
    void projectPointsToFunction(
        Functional const& distanceFunction, double const baseEdgeLength,
        Eigen::Ref<Eigen::ArrayXXd> points);

    // check whether points lies inside or outside of polygon
    Eigen::ArrayXXd pointsInsidePoly(
        Eigen::Ref<Eigen::ArrayXXd const> const points,
        Eigen::Ref<Eigen::ArrayXXd const> const polygon);
}
}

#endif
