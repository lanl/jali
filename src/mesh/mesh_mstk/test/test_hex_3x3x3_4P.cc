#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

TEST(MSTK_HEX_3x3x3_4P)
{

  int i, j, k, err, nc, nf, nv;
  std::vector<Jali::Entity_ID> faces(6), nodes(8);
  std::vector<JaliGeometry::Point> ccoords(8), fcoords(4);

  int rank, size;

  int initialized;
  MPI_Initialized(&initialized);
  
  if (!initialized)
    MPI_Init(NULL,NULL);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  CHECK_EQUAL(4,size);

  if (size != 4) {
    std::cerr << "Test must be run with 4 processors" << std::endl;
    //    return;
  }

  //  if (rank == 0) {
  int DebugWait = 0;
  while (DebugWait);
  //  }

  // Load a single hex from the hex1.exo file

  Jali::Mesh *mesh(new Jali::Mesh_MSTK("test/hex_3x3x3_sets.exo",
                                                  MPI_COMM_WORLD,3));

}

