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
# Load build configuration
##############################
CONFIG_FILE ?= Makefile.config
include $(CONFIG_FILE)

##############################
# Main output directories
##############################
prefix ?= /usr/local
ROOT_BUILD_DIR := build

# adjust build dir for debug configuration
DEBUG ?= 0
ifeq ($(DEBUG), 1)
	BUILD_DIR := $(ROOT_BUILD_DIR)/debug
else
	BUILD_DIR := $(ROOT_BUILD_DIR)/release
endif

##############################
# Compiler
##############################
AR := ar rcs
CXX ?= g++

# Target build architecture
TARGET_ARCH_NAME ?= $(shell $(CXX) -dumpmachine)
BUILD_DIR := $(BUILD_DIR)/$(TARGET_ARCH_NAME)

##############################
# The shared library and static library names
##############################
NAME := $(BUILD_DIR)/lib/lib$(PROJECT).so
STATIC_NAME := $(BUILD_DIR)/lib/lib$(PROJECT)_static.a

##############################
# Includes and libraries
##############################
LIBRARIES := qhull
LIBRARY_DIRS +=
INCLUDE_DIRS += ./include

##############################
# Compiler Flags
##############################
GIT_VERSION := $(shell git describe --tags --long)
COMMON_FLAGS := $(addprefix -I, $(INCLUDE_DIRS)) -DGIT_VERSION=\"$(GIT_VERSION)\"
CXXFLAGS := -std=c++11 -fPIC
LINKFLAGS := -fPIC -static-libstdc++
LDFLAGS := $(addprefix -l, $(LIBRARIES)) $(addprefix -L, $(LIBRARY_DIRS)) $(addprefix -Xlinker -rpath , $(LIBRARY_DIRS))

# Set compiler flags for debug configuration
ifeq ($(DEBUG), 1)
	COMMON_FLAGS += -g -O0 -DDEBUG
else
	COMMON_FLAGS += -O3 -DNDEBUG
endif

##############################
# Source Files
##############################
CXX_SRCS := $(shell find src -name "*.cpp")
HXX_SRCS := $(shell find include -name "*.h")
EXAMPLES_SRCS := $(shell find examples -name "*.cpp")
EXAMPLES_SRCS := $(filter-out $(UTILS_SRCS), $(EXAMPLES_SRCS))

# Object files
CXX_OBJS := $(addprefix $(BUILD_DIR)/objs/, $(CXX_SRCS:.cpp=.o))
EXAMPLES_OBJS := $(addprefix $(BUILD_DIR)/objs/, $(EXAMPLES_SRCS:.cpp=.o))
EXAMPLES_BINS := $(patsubst examples%.cpp, $(BUILD_DIR)/examples%, $(EXAMPLES_SRCS))

##############################
# Build targets
##############################
.PHONY: all install clean examples

all: $(NAME) $(STATIC_NAME) examples

examples: $(EXAMPLES_BINS)

$(EXAMPLES_BINS): $(BUILD_DIR)/examples/% : $(BUILD_DIR)/objs/examples/%.o $(UTILS_OBJS) $(STATIC_NAME)
	@echo [ Linking ] $@
	@mkdir -p $(BUILD_DIR)/examples
	@cp examples/plot_mesh.py $(BUILD_DIR)/examples
	@$(CXX) -o $@ $< $(UTILS_OBJS) $(STATIC_NAME) $(COMMON_FLAGS) $(LDFLAGS) $(LINKFLAGS)

$(NAME): $(CXX_OBJS)
	@echo [ Linking ] $@
	@mkdir -p $(BUILD_DIR)/lib
	@$(CXX) -shared -o $@ $(CXX_OBJS) $(COMMON_FLAGS) $(LDFLAGS) $(LINKFLAGS)

$(STATIC_NAME): $(CXX_OBJS)
	@echo [ Linking ] $@
	@mkdir -p $(BUILD_DIR)/lib
	@$(AR) $@ $(CXX_OBJS)

$(BUILD_DIR)/objs/%.o: %.cpp $(HXX_SRCS)
	@echo [ CXX ] $<
	@$(foreach d, $(subst /, ,${@D}), mkdir -p $d && cd $d && ):
	@$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) -c -o $@ $<

install: $(NAME) $(STATIC_NAME) $(HXX_SRCS) $(EXAMPLES_BINS)
	@install -m 0644 $(NAME) $(prefix)/lib
	@install -m 0644 $(STATIC_NAME) $(prefix)/lib
	@$(foreach f, $(HXX_SRCS), install -D -m 0644 $f $(prefix)/$f && ):
	@$(foreach f, $(EXAMPLES_BINS), install -m 0755 $f $(prefix)/bin && ):

clean:
	@rm -rf $(ROOT_BUILD_DIR)
