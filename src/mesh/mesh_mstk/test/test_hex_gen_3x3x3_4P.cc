#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

// Test for generation of hex mesh distributed over 4 processors

TEST(MSTK_HEX_GEN_3x3x3_4P) {

  int i, j, k, err, nc, nf, nv;
  std::vector<Jali::Entity_ID> faces(6), nodes(8);
  std::vector<JaliGeometry::Point> ccoords(8), fcoords(4);

			
  int rank, size;

  int initialized;

  MPI_Initialized(&initialized);

  if (!initialized)
    MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  CHECK_EQUAL(4, size);

  if (size != 4) {
    std::cerr << "Test must be run with 4 processors" << std::endl;
  }

  // Create a 3x3x3 cell hex mesh

  Jali::Mesh *mesh(new Jali::Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3,
                                       MPI_COMM_WORLD));

}

