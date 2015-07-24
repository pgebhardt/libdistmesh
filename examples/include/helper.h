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

#ifndef _eeb6ef75_b27b_4918_b448_3e996841ab4c
#define _eeb6ef75_b27b_4918_b448_3e996841ab4c

#include <fstream>
#include <chrono>

namespace distmesh {
namespace helper {
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

    class HighPrecisionTime {
    private:
        std::chrono::high_resolution_clock::time_point time;

    public:
        HighPrecisionTime() {
            this->restart();
        }
        
        void restart() {
            this->time = std::chrono::high_resolution_clock::now();
        }
        
        double elapsed() const {
            return std::chrono::duration_cast<std::chrono::duration<double>>(
                std::chrono::high_resolution_clock::now() - this->time).count();
        }
    };
}
}

#endif
