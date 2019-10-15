#HDF5_DIR = $(HOME)/wsa/tools/lassen/hdf5-1.10.5
HDF5_DIR = /g/g20/emccarth/app/power9/spack/opt/spack/linux-rhel7-ppc64le/gcc-7.3.1/hdf5-1.10.5-w73fzgh7oj5tkzbbks44fj6zgm3odtzo/
FLAGS = -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5 -ldl -fopenmp

all: readFileHDF5_dat_plain transpose 
#partition_to_raw read_partitions

readFileHDF5_dat_plain: readFileHDF5.cpp common.hpp
	mpic++ -O2 -Wall $< -o $@ $(FLAGS)

transpose: transpose.cpp common.hpp
	mpic++ -O2 -Wall $< -o $@ $(FLAGS)

#partition_to_raw: partition_to_raw.cpp common.hpp
#	mpic++ -O2 -Wall $< -o $@ $(FLAGS)

#read_partitions: read_partitions.cpp common.hpp
#	mpic++ -O2 -Wall $< -o $@ $(FLAGS)
