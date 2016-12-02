#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

// Test generation of quad mesh in serial

TEST(MSTK_QUAD_GEN_3x3) {

  int i, j, k, err, nc, nf, nv;
  Jali::Set_ID faces[6], nodes[8];

  int NV = 16;
  int NF = 24;
  int NC = 9;

  // Load a mesh consisting of 3x3 elements

  Jali::Mesh *mesh(new Jali::Mesh_MSTK(0.0, 0.0, 1.0, 1.0, 3, 3,
                                       MPI_COMM_WORLD));

  nv = mesh->num_entities(Jali::Entity_kind::NODE,
                          Jali::Entity_type::PARALLEL_OWNED);
  CHECK_EQUAL(NV, nv);

  nf = mesh->num_entities(Jali::Entity_kind::FACE,
                          Jali::Entity_type::PARALLEL_OWNED);
  CHECK_EQUAL(NF, nf);

  nc = mesh->num_entities(Jali::Entity_kind::CELL,
                          Jali::Entity_type::PARALLEL_OWNED);
  CHECK_EQUAL(NC, nc);

}

