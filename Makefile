# NOTE change these things to match your setup!!!
OPENMC = /home/gavin/code/openmc
OPENMC_BUILD = $(OPENMC)/build_opt
CXX ='g++'

# These should be set automatically
INCLUDE = -I$(OPENMC)/include -I$(OPENMC)/vendor/pugixml -I$(OPENMC)/vendor/gsl/include
LIBS = -L$(OPENMC_BUILD)/lib
HDF_INCLUDES = -I/usr/include/hdf5/serial
HDF_LINKS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5
DEBUG = -O0 -g -fstack-protector -fopenmp
OPT = -O3 -ffast-math -march=native -fopenmp

# set whether its in debug or opt mode plus other options that regardless apply
FLAGS = $(DEBUG) $(INCLUDE) $(LIBS) $(HDF_LINKS) $(HDF_INCLUDES)

openmc-explorer: openmc-explorer.cc
	$(CXX) openmc-explorer.cc $(FLAGS) -lopenmc -lSDL2 -o openmc-explorer

clean:
	rm openmc-explorer
