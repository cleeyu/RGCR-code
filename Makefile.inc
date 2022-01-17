CXX = g++
# Need to change this for your location of SNAP
SNAP_DIR = /home/yinh/software/Snap-4.1

SNAP_INCLUDES = -I$(SNAP_DIR)/snap-core -I$(SNAP_DIR)/snap-adv -I$(SNAP_DIR)/snap-exp -I$(SNAP_DIR)/glib-core
SNAP_OBJ = $(SNAP_DIR)/snap-core/Snap.o
INCLUDES := -I. $(SNAP_INCLUDES)
DEFINES := -DNDEBUG

OPT := -O3
CFLAGS := $(OPT) $(INCLUDES) $(DEFINES) -std=c++11 $(DEFINES)

# LDFLAGS := -O3
LDFLAGS := -O3 -lgomp -lrt -fopenmp

vpath %.cpp src

