#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

// Test generation of hex mesh in serial

TEST(MSTK_HEX_GEN_3x3x3)
{

  int i, j, k, err, nc, nf, nv;
  Jali::Set_ID faces[6], nodes[8];

  int NV = 64;
  int NF = 108;
  int NC = 27;

  // Generate a mesh consisting of 3x3x3 elements

  Jali::Mesh *mesh(new Jali::Mesh_MSTK(0.0,0.0,0.0,1.0,1.0,1.0,3,3,3,MPI_COMM_WORLD));

  nv = mesh->num_entities(Jali::NODE,Jali::OWNED);
  CHECK_EQUAL(NV,nv);

  nf = mesh->num_entities(Jali::FACE,Jali::OWNED);
  CHECK_EQUAL(NF,nf);

  nc = mesh->num_entities(Jali::CELL,Jali::OWNED);
  CHECK_EQUAL(NC,nc);

}

