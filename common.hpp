#pragma once
//comment
#include <vector>
#include <string>
#include <dirent.h>
#include <iostream>

#include "hdf5.h"

#define BASE_PATH "/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/hdf5"
#define TRANSPOSED_PATH "/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/hdf5_transposed"
#define RAW_PATH "/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/hdf5_transposed_raw"
#define TRAIN_PATH "/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/hdf5_transposed/train"
#define VAL_PATH "/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/hdf5_transposed/val"
#define TEST_PATH "/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/hdf5_transposed/test"
#define ALL_PATH  "/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/hdf5_transposed/all"
constexpr int RANK = 4;
constexpr size_t spatial_dim = 512;
constexpr size_t channel_dim = 4;
constexpr size_t response_nums = 4;

inline hid_t check_hdf5(hid_t hid, const char *file, int line) {
  if (hid < 0) {
    std::cerr << "HDF5 error" << std::endl;
    std::cerr << "Error at " << file << ":" << line << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return hid;
}

#define CHECK_HDF5(call) check_hdf5(call, __FILE__, __LINE__)

inline std::vector<std::string> get_all_dirs() {
  std::vector<std::string> dirs = {
    "21688988", "21922619", "21997469", "22059249", "22098324",
    "22309462", "21812950", "21929749", "22021490", "22074825",
    "22118427"};
  return dirs;
}

inline std::vector<std::string> locate_sample_files(
  const std::vector<std::string> &dirs, bool trans, int num_chunks) {
  std::vector<std::string> files;
  struct dirent *entry;
  std::string dirPath = trans ? TRANSPOSED_PATH : BASE_PATH;
  
  if (num_chunks > 1) {
    dirPath += "_chunk" + std::to_string(num_chunks);
  }
  dirPath += "/";
  for (auto dir_name: dirs) {
    dir_name = dirPath + dir_name;
    if (dir_name.back() != '/') {
      dir_name += "/";
    }
    //std::cout << "Using samples under " << dir_name << std::endl;
    DIR *dir = opendir(dir_name.c_str());
    while ((entry = readdir(dir)) != NULL) {
      std::string file_name = entry->d_name;
      if (file_name.find(".hdf5") == std::string::npos) {
        // not an HDF file
        continue;
      }
      files.push_back(dir_name + file_name);
    }
    closedir(dir);
  }
  std::cout<<files.size()<< " files \n";
  return files;
}

inline std::string get_transposed_path(std::string path) {
  std::string key = "hdf5/";
  std::string replaced = "hdf5_transposed/";
  return path.replace(path.find(key), key.size(), replaced);
}

inline std::vector<std::string> locate_raw_files(
    const std::vector<std::string> &dirs,
    int num_ranks_per_sample, int zidx) {
  std::vector<std::string> files;
  struct dirent *entry;
  std::string dirPath = std::string(RAW_PATH) +
      std::to_string(num_ranks_per_sample) + "/z" +
      std::to_string(zidx) + "/";
  for (auto dir_name: dirs) {
    dir_name = dirPath + dir_name;
    if (dir_name.back() != '/') {
      dir_name += "/";
    }
    //std::cout << "Using samples under " << dir_name << std::endl;
    DIR *dir = opendir(dir_name.c_str());
    while ((entry = readdir(dir)) != NULL) {
      std::string file_name = entry->d_name;
      if (file_name.find(".raw") == std::string::npos) {
        continue;
      }
      files.push_back(dir_name + file_name);
    }
    closedir(dir);
  }
  std::cout<<files.size()<< " files \n";
  return files;
}
