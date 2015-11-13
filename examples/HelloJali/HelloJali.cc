#include <iostream>
#include <vector>
#include <array>

#include "mpi.h"

// Most basic application code compiled with Jali but not using it at all

int main(int argc, char *argv[]) {

  // Jali depends on MPI 

  MPI_Init(&argc,&argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int rank, nprocs;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&nprocs);

  if (rank == 0)
    std::cerr << "Running on " << nprocs << " ranks" << std::endl;

  std::cerr << "Hello, from Jali running on rank " << rank << std::endl;

  MPI_Finalize();
}

