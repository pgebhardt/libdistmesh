# --------------------------------------------------------------------
# This file is part of libdistmesh.
#
# libdistmesh is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# libdistmesh is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with libdistmesh. If not, see <http:#www.gnu.org/licenses/>.
#
# Copyright (C) 2015 Patrik Gebhardt
# Contact: patrik.gebhardt@rub.de
# --------------------------------------------------------------------

# The Makefile for libdistmesh
PROJECT := distmesh

##############################
# Directories
##############################
BUILD_DIR := build
prefix ?= /usr/local

##############################
# Compiler
##############################
AR := ar rcs
CXX := clang++

# Cross compile for arm architecture
ARM ?= 0
ifeq ($(ARM), 1)
	CXX := arm-linux-gnueabihf-g++
	TARGET_ARCH := armv7-linux-gnueabihf
endif

# Target build architecture
TARGET_ARCH ?= $(shell uname -m)-$(shell uname)
BUILD_DIR := $(BUILD_DIR)/$(TARGET_ARCH)

##############################
# The target shared library and static library names
##############################
LIB_BUILD_DIR := $(BUILD_DIR)/lib
NAME := $(LIB_BUILD_DIR)/lib$(PROJECT).so
STATIC_NAME := $(LIB_BUILD_DIR)/lib$(PROJECT)_static.a

##############################
# Includes and libraries
##############################
LIBRARIES := qhull
LIBRARY_DIRS :=
INCLUDE_DIRS := ./include

##############################
# Compiler Flags
##############################
GIT_VERSION := $(shell git describe --tags --long)
COMMON_FLAGS := $(addprefix -I, $(INCLUDE_DIRS)) -DGIT_VERSION=\"$(GIT_VERSION)\" -O3
CFLAGS := -std=c++11 -fPIC
LINKFLAGS := -O3 -fPIC -static-libgcc -static-libstdc++
LDFLAGS := $(addprefix -l, $(LIBRARIES)) $(addprefix -L, $(LIBRARY_DIRS))

##############################
# Source Files
##############################
CXX_SRCS := $(shell find src -name "*.cpp")
HXX_SRCS := $(shell find include -name "*.h")
EXAMPLE_SRCS := $(shell find examples -name "*.cpp")

# Object files
CXX_OBJS := $(addprefix $(BUILD_DIR)/objs/, $(CXX_SRCS:.cpp=.o))
EXAMPLE_OBJS := $(addprefix $(BUILD_DIR)/objs/, $(EXAMPLE_SRCS:.cpp=.o))
EXAMPLE_BINS := $(patsubst examples%.cpp, $(BUILD_DIR)/examples%, $(EXAMPLE_SRCS))

##############################
# Build targets
##############################
.PHONY: all install clean examples

all: $(NAME) $(STATIC_NAME) examples

examples: $(EXAMPLE_BINS)

$(EXAMPLE_BINS): $(BUILD_DIR)/examples/% : $(BUILD_DIR)/objs/examples/%.o $(STATIC_NAME)
	@mkdir -p $(BUILD_DIR)/examples
	@cp -u examples/plot_mesh.py $(BUILD_DIR)/examples
	$(CXX) $< $(STATIC_NAME) -o $@ $(LDFLAGS) $(LINKFLAGS)

$(NAME): $(CXX_OBJS) $(CU_OBJS)
	@mkdir -p $(LIB_BUILD_DIR)
	$(CXX) -shared -o $@ $(CXX_OBJS) $(CU_OBJS) $(LDFLAGS) $(LINKFLAGS)

$(STATIC_NAME): $(CXX_OBJS) $(CU_OBJS)
	@mkdir -p $(LIB_BUILD_DIR)
	$(AR) $@ $(CXX_OBJS) $(CU_OBJS)

$(BUILD_DIR)/objs/%.o: %.cpp $(HXX_SRCS)
	@$(foreach d, $(subst /, ,${@D}), mkdir -p $d && cd $d && ):
	$(CXX) $(CFLAGS) $(COMMON_FLAGS) -c -o $@ $<

install: $(NAME) $(STATIC_NAME) $(HXX_SRCS)
	@install -m 0644 $(NAME) $(prefix)/lib
	@install -m 0644 $(STATIC_NAME) $(prefix)/lib
	@$(foreach f, $(HXX_SRCS), install -D -m 0644 $f $(prefix)/$f && ):

clean:
	@rm -rf $(BUILD_DIR)
