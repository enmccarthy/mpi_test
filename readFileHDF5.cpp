
#include <iostream>
#include "mpi.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <string>
#include <algorithm>
#include <cassert>

#include "common.hpp"

int main(int argc, char *argv[]) {
  const size_t sample_size = spatial_dim * spatial_dim * spatial_dim *
      channel_dim;
  int num_ranks_per_sample = 8; // default partitioning
  if (argc > 1) {
    num_ranks_per_sample = std::atoi(argv[1]);
  }
  const size_t zPerNode = spatial_dim / num_ranks_per_sample;
  const size_t local_sample_size = sample_size / num_ranks_per_sample;
  bool trans = false;
  if (argc > 2) {
    trans = std::atoi(argv[2]);
  }

  bool chunked = false;
  if (argc > 3) {
    chunked = std::atoi(argv[3]);
  }
  
  std::vector<std::string> dirs;
  for (int i = 4; i < argc; ++i) {
    dirs.push_back(argv[i]);
  }
  if (dirs.size() == 0) {
    dirs.push_back("21688988");
  } else {
    if (dirs[0] == "all") {
      dirs = get_all_dirs();
    }
  }

  auto &&files = locate_sample_files(dirs, trans,
                                     chunked ? num_ranks_per_sample : 1);
  const int numsamples = files.size();

  MPI_Init(&argc, &argv);
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Info mpi_info = MPI_INFO_NULL;
  int spatial_rank = rank % num_ranks_per_sample;

  if (nprocs < num_ranks_per_sample) {
    if (rank == 0) {
      std::cerr << "Error: num_ranks_per_sample is larger than the number of MPI ranks" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  if (rank == 0) {
    std::cout << "Local dimension: "
              << spatial_dim
              << "x" << spatial_dim
              << "x" << spatial_dim / num_ranks_per_sample
              << "x" << channel_dim
              << std::endl;
    if (trans) {
      std::cout << "Read transposed sample files\n";
    }
    if (chunked) {
      std::cout << "Read chunked sample files\n";
    }
  }

  MPI_Comm m_comm;
  int color = rank/ num_ranks_per_sample;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &m_comm);

  // Setting up HDF5 data structures
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  int ierr;
#if 0
  ierr = H5Pset_sieve_buf_size(fapl_id, 52262144);
  if(num_ranks_per_sample==8) {
   // ierr = H5Pset_sieve_buf_size(fapl_id, 52262144);
    ierr = H5Pset_sieve_buf_size(fapl_id, 1000000000);
  } else if (num_ranks_per_sample==16) {
    ierr = H5Pset_sieve_buf_size(fapl_id, 2000000000);
  } else if (num_ranks_per_sample==32) {
    ierr = H5Pset_sieve_buf_size(fapl_id, 4000000000);
  }
#endif
  H5Pset_fapl_mpio(fapl_id, m_comm, mpi_info);
  hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);


  const hsize_t count[RANK] = {1, 1, 1, 1};
  std::vector<hsize_t> dims_local = {spatial_dim, spatial_dim,
                                     zPerNode, channel_dim};
  std::vector<hsize_t> offset = {0, 0, zPerNode * spatial_rank, 0};
  if (trans) {
    std::reverse(std::begin(dims_local), std::end(dims_local));
    std::reverse(std::begin(offset), std::end(offset));
  }
  auto memspace = H5Screate_simple(RANK, dims_local.data(), NULL);

  if (rank == 0) {
    std::cout << dims_local[0] << ", " << dims_local[1]
              << ", " << dims_local[2] << ", " << dims_local[3]
              << std::endl;
  }

  // read buffer
  std::vector<short> data_out(local_sample_size);

  int trial_count = 1;
  double start = 0;

  for (int trial_idx = 0; trial_idx < trial_count; ++trial_idx) {
    if (trial_idx == trial_count - 1) {
      MPI_Barrier(MPI_COMM_WORLD);
      start = MPI_Wtime();
    }
    for(int file_idx = rank/num_ranks_per_sample; file_idx < numsamples;
        file_idx += nprocs/num_ranks_per_sample) {
      herr_t status;
      auto filename = files[file_idx];
#if 0
      if (rank == 0) {
        std::cout << "file_idx: " << file_idx
                  << ", file name: " << filename << std::endl;
      }
#endif
      //auto file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY,
      //H5P_DEFAULT);
      auto file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY,
                          fapl_id);
      char name[100];

      H5Gget_objname_by_idx(file, 0, name, 100);
      //auto dataset = H5Dopen(file, name, dxpl_id);
      auto dataset = H5Dopen(file, name, H5P_DEFAULT);
      auto filespace = H5Dget_space(dataset);

      CHECK_HDF5(H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                     offset.data(),
                                     NULL, count, dims_local.data()));

      // Use MPI-IO to read the data. This doesn't seem to have
      // measurable impact, though.
#if 1
      CHECK_HDF5(H5Dread(dataset, H5T_NATIVE_SHORT, memspace, filespace,
                         dxpl_id, data_out.data()));
#else
      CHECK_HDF5(H5Dread(dataset, H5T_NATIVE_SHORT, memspace, filespace,
                         H5P_DEFAULT, data_out.data()));
#endif
      assert(status >= 0);

      //std::cout<<data_out[10]<<"\n";
      //MPI_Comm_free(&file_com);

      H5Dclose(dataset);
      H5Fclose(file);
    
    if (false) {
      std::ofstream output;
      std::string outname = "output";
      std::string out;
      std::stringstream ss;
      ss << rank;
      out = ss.str();
      outname.append(out);
      outname.append(".txt");
      output.open(outname.c_str());
      //    std::cout<<(xPerNode*yPerNode*zPerNode*sPerNode) << "\n";
      std::cout<<"Rank: "<< rank<<"here\n";
      std::cout<<"file " << filename << "\n";
      for (uint iter = 0; iter < data_out.size(); iter++) {
        
       output << data_out.data()[iter] << " ";
      }
      output.close();
    }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double end = MPI_Wtime();
  //MPI_Group_free(&world_group);
  //  MPI_Group_free(&file_group);
  //  MPI_Comm_free(&file_comm);
  MPI_Finalize();
  std::cout<< "Rank " << rank << " \n";
  std::cout<< "The process took " << end - start << " seconds to run. \n";
  return 0;
}
