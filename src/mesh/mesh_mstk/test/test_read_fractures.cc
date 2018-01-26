#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

TEST(MSTK_NON_MANIFOLD_SURFS) {

  int i, j, k, err, nc, nf, nv;
  Jali::Entity_ID faces[6], nodes[8];

  // Load a mesh that has two surfaces intersecting in a non-manifold fashion

  Jali::Mesh *mesh(new Jali::Mesh_MSTK("test/fractures.exo",
                                       MPI_COMM_WORLD));

  nf = mesh->num_entities(Jali::Entity_kind::FACE,
                          Jali::Entity_type::PARALLEL_OWNED);

  nc = mesh->num_entities(Jali::Entity_kind::CELL,
                          Jali::Entity_type::PARALLEL_OWNED);

}

