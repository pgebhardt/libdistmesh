// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_COMMON_H
#define LIBDISTMESH_INCLUDE_COMMON_H

// standard c++ includes
#include <iostream>
#include <memory>
#include <functional>

// boost multi array library for array handling
#include <boost/multi_array.hpp>

// qhull library used to calculate delaunay triangulation
extern "C" {
#define qh_QHimport
#include <qhull/qhull_a.h>
}

#endif
