#include <iostream>
#include "mpi.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <string>
#include <cassert>

#include "common.hpp"

void transpose_matrix(std::vector<short> &buf,
                      hsize_t *dims) {
  auto transposed = buf;
  auto get_offset = [&](hsize_t i, hsize_t j, hsize_t k, hsize_t l) {
                      return i * dims[1] * dims[2] * dims[3]
                          + j * dims[2] * dims[3]
                          + k * dims[3]
                          + l;};
  auto get_transposed_offset = [&](hsize_t i, hsize_t j, hsize_t k, hsize_t l) {

                                 return l * dims[0] * dims[1] * dims[2]
                                     + k * dims[0] * dims[1]
                                     + j * dims[0]
                                     + i;};
#pragma omp parallel for
  for (hsize_t i = 0; i < dims[0]; ++i) {
    for (hsize_t j = 0; j < dims[1]; ++j) {
      for (hsize_t k = 0; k < dims[2]; ++k) {
        for (hsize_t l = 0; l < dims[3]; ++l) {
          transposed.data()[get_transposed_offset(i, j, k, l)] =
              buf.data()[get_offset(i, j, k, l)];
        }
      }
    }
  }
  buf = transposed;
}

void transpose(const std::string &src_path,
               std::vector<short> &buf) {
  std::cout << "Reading " << src_path << std::endl;
  auto src_file = H5Fopen(src_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  auto dst_path = get_transposed_path(src_path);
  auto dst_file = CHECK_HDF5(H5Fcreate(dst_path.c_str(), H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT));

  std::vector<std::string> names = {"full", "physPar"};
  std::vector<hid_t> types = {H5T_STD_I16LE, H5T_IEEE_F32LE};
  hsize_t dims[RANK] = {spatial_dim, spatial_dim, spatial_dim,
                        channel_dim};
  hsize_t transposed_dims[RANK] = {channel_dim, spatial_dim,
                                   spatial_dim, spatial_dim};
  const hsize_t response_dim = response_nums;
  std::vector<hid_t> dataspaces = {
    CHECK_HDF5(H5Screate_simple(RANK, dims, NULL)),
    CHECK_HDF5(H5Screate_simple(1, &response_dim, NULL))};
  std::vector<hid_t> transposed_dataspaces = {
    CHECK_HDF5(H5Screate_simple(RANK, transposed_dims, NULL)),
    CHECK_HDF5(H5Screate_simple(1, &response_dim, NULL))};

  for (int i = 0; i < (int)names.size(); ++i) {
    auto name = names[i];
    auto type = types[i];
    auto dataset = CHECK_HDF5(H5Dopen(src_file, name.c_str(),
                                      H5P_DEFAULT));

    // Use buf even for redshifts as it's large enough
    void *read_buf = buf.data();
    CHECK_HDF5(H5Dread(dataset, type, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, read_buf));
    CHECK_HDF5(H5Dclose(dataset));
    auto transposed_dataset = CHECK_HDF5(
        H5Dcreate2(dst_file, name.c_str(), type,
                   transposed_dataspaces[i], H5P_DEFAULT,
                   H5P_DEFAULT, H5P_DEFAULT));
    if (name == "full") {
      transpose_matrix(buf, dims);
    }
    CHECK_HDF5(H5Dwrite(transposed_dataset, type, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, read_buf));
    CHECK_HDF5(H5Dclose(transposed_dataset));
  }
  CHECK_HDF5(H5Fclose(src_file));
  CHECK_HDF5(H5Fclose(dst_file));
}

int main(int argc, char *argv[]) {
  const size_t sample_size = spatial_dim * spatial_dim * spatial_dim *
      channel_dim;
  std::vector<std::string> dirs;
  for (int i = 1; i < argc; ++i) {
    dirs.push_back(argv[i]);
  }

  auto &&files = locate_sample_files(dirs, false, 1);

  const int numsamples = files.size();

  MPI_Init(&argc,&argv);
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Info mpi_info = MPI_INFO_NULL;

  if (rank == 0) {
    std::cout << "Sample dimension: "
              << spatial_dim << "x"
              << "x" << spatial_dim
              << "x" << spatial_dim
              << "x" << channel_dim
              << std::endl;
  }

  MPI_Comm m_comm;
  int color = rank;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &m_comm);

  // Setting up HDF5 data structures
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, m_comm, mpi_info);

  // read buffer
  std::vector<short> data_out(sample_size);

  for(int file_idx = rank; file_idx < numsamples; file_idx += nprocs) {
    if (rank == 0) {
      //std::cout << "file_idx: " << file_idx << std::endl;
    }
    auto filename = files[file_idx];
    transpose(filename, data_out);
  }

  MPI_Finalize();
  return 0;
}
