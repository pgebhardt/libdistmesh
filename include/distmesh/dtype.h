// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_DTYPE_H
#define LIBDISTMESH_INCLUDE_DTYPE_H

// namespace libdistmesh::dtype
namespace mpFlow {
namespace dtype {
    // basic scalar types
    typedef double real;
    typedef boost::multi_array_types::index index;

    // basic array type
    template <
        int dims
    >
    using array = boost::multi_array<real, dims>;
}
}

#endif
