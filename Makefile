# The Makefile for libdistmesh
PROJECT := distmesh

# Directories
BUILD_DIR := build
PREFIX ?= /usr/local

# Target build architecture
ARM ?= 0
ifeq ($(ARM), 1)
	TARGET_ARCH := armv7-linux-gnueabihf
else
	TARGET_ARCH := x86_64-linux
endif
BUILD_DIR := $(BUILD_DIR)/$(TARGET_ARCH)

# The target shared library and static library names
LIB_BUILD_DIR := $(BUILD_DIR)/lib
NAME := $(LIB_BUILD_DIR)/lib$(PROJECT).so
STATIC_NAME := $(LIB_BUILD_DIR)/lib$(PROJECT)_static.a

# Compiler
AR := ar rcs
ifeq ($(ARM), 1)
	CXX := arm-linux-gnueabihf-g++
else
	CXX := clang++
endif

# Version Define
GIT_VERSION := $(shell git describe --tags --long)

# Includes and libraries
LIBRARIES := qhull
LIBRARY_DIRS :=
INCLUDE_DIRS := /usr/local/include ./include

# Compiler Flags
COMMON_FLAGS := $(addprefix -I, $(INCLUDE_DIRS)) -DGIT_VERSION=\"$(GIT_VERSION)\" -O3
CFLAGS := -std=c++11 -fPIC
LDFLAGS := $(addprefix -l, $(LIBRARIES)) $(addprefix -L, $(LIBRARY_DIRS))

# Source Files
CXX_SRCS := $(shell find src -name "*.cpp")
HXX_SRCS := $(shell find include -name "*.h")
EXAMPLE_SRCS := $(shell find examples -name "*.cpp")

# Object files
CXX_OBJS := $(addprefix $(BUILD_DIR)/, ${CXX_SRCS:.cpp=.o})
EXAMPLE_OBJS := $(addprefix $(BUILD_DIR)/, ${EXAMPLE_SRCS:.cpp=.o})
EXAMPLE_BINS := ${EXAMPLE_OBJS:.o=}

# Build targets
.PHONY: all install clean examples

all: $(NAME) $(STATIC_NAME)

examples: $(EXAMPLE_BINS)

$(EXAMPLE_BINS): % : %.o $(STATIC_NAME)
	@mkdir -p $(BUILD_DIR)/examples
	@cp examples/plot_mesh.py $(BUILD_DIR)/examples
	$(CXX) $< $(STATIC_NAME) -o $@ $(LDFLAGS)

$(NAME): $(CXX_OBJS) $(CU_OBJS)
	@mkdir -p $(LIB_BUILD_DIR)
	$(CXX) -shared -o $@ $(CXX_OBJS) $(CU_OBJS)

$(STATIC_NAME): $(CXX_OBJS) $(CU_OBJS)
	@mkdir -p $(LIB_BUILD_DIR)
	$(AR) $@ $(CXX_OBJS) $(CU_OBJS)

$(BUILD_DIR)/%.o: %.cpp $(HXX_SRCS)
	@$(foreach d, $(subst /, ,${@D}), mkdir -p $d && cd $d && ):
	$(CXX) $(CFLAGS) $(COMMON_FLAGS) -c -o $@ $<

install: $(NAME) $(STATIC_NAME) $(HXX_SRCS)
	install -m 0644 $(NAME) $(PREFIX)/lib
	install -m 0644 $(STATIC_NAME) $(PREFIX)/lib
	$(foreach f, $(HXX_SRCS), install -D -m 0644 $f $(PREFIX)/$f && ):

clean:
	@rm -rf $(BUILD_DIR)