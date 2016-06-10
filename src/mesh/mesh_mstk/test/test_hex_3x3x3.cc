#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

TEST(MSTK_HEX_3x3x3) {

  int i, j, k, err, nc, nf, nv;
  Jali::Entity_ID faces[6], nodes[8];

  int NV = 64;
  int NF = 108;
  int NC = 27;

  // Load a mesh consisting of 3x3x3 elements

  Jali::Mesh *mesh(new Jali::Mesh_MSTK("test/hex_3x3x3_sets.exo",
                                       MPI_COMM_WORLD, 3));

  nf = mesh->num_entities(Jali::Entity_kind::FACE,
                          Jali::Entity_type::PARALLEL_OWNED);
  CHECK_EQUAL(NF, nf);

  nc = mesh->num_entities(Jali::Entity_kind::CELL,
                          Jali::Entity_type::PARALLEL_OWNED);
  CHECK_EQUAL(NC, nc);

}

