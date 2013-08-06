// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_DTYPE_H
#define LIBDISTMESH_INCLUDE_DTYPE_H

// namespace libdistmesh::dtype
namespace distmesh {
namespace dtype {
    // basic scalar types
    typedef double real;
    typedef unsigned int index;

    // basic linear algebra types
    template <
        class type
    >
    using matrix = Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>;
}
}

#endif
