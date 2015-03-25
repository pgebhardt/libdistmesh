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

#include <distmesh/distmesh.h>
#include <fstream>

// save eigen array to text file
template <typename type>
void savetxt(Eigen::Ref<Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> const> const array,
    std::string const& filename) {
    // open file
    std::ofstream file(filename);

    // save array to file with high precision
    Eigen::IOFormat const format(Eigen::FullPrecision, Eigen::DontAlignCols);
    file << array.format(format);

    file.close();
}

int main() {
    Eigen::ArrayXXd boundingBox(2, 2);
    boundingBox << 0.0, 1.0, 0.0, 1.0;

    // fixed points at corners of domain to guarantee convergence
    Eigen::ArrayXXd fixedPoints(4, 2);
    fixedPoints << 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0;

    // distance function for rectangular domain
    auto distanceFunction =
        distmesh::distanceFunction::rectangular(boundingBox);

    // corner points of polygon
    Eigen::ArrayXXd poly(2, 2);
    poly << 0.3, 0.7, 0.7, 0.5;

    // edge size function
    auto sizeFunction =
        (0.01 + 0.3 * distmesh::distanceFunction::circular(0.0).abs())
        .min(0.025 + 0.3 * distmesh::distanceFunction::polygon(poly).abs())
        .min(0.15);

    // create mesh
    auto mesh = distmesh::distmesh(distanceFunction, 0.01,
        sizeFunction, boundingBox, fixedPoints);

    // save mesh to file
    savetxt<double>(std::get<0>(mesh), "points.txt");
    savetxt<int>(std::get<1>(mesh), "triangulation.txt");

    // plot mesh using python
    system("python plot_mesh.py");

    return 0;
}
