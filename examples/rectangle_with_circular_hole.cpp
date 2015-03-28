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

#include <fstream>
#include <distmesh/distmesh.h>

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
    // fixed points at the corners of domain to guarantee convergence
    Eigen::ArrayXXd fixedPoints(4, 2);
    fixedPoints << -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;

    // create mesh
    auto mesh = distmesh::distmesh(
        distmesh::distanceFunction::rectangular(distmesh::boundingBox(2))
            .max(-distmesh::distanceFunction::circular(0.5)),
        0.05, 0.05 + 0.3 * distmesh::distanceFunction::circular(0.5),
        distmesh::boundingBox(2), fixedPoints);

    // save mesh to file
    savetxt<double>(std::get<0>(mesh), "points.txt");
    savetxt<int>(std::get<1>(mesh), "triangulation.txt");

    // plot mesh using python
    return system("python plot_mesh.py");
}
