/*
 Copyright (c) 2019, Triad National Security, LLC
 All rights reserved.

 Copyright 2019. Triad National Security, LLC. This software was
 produced under U.S. Government contract 89233218CNA000001 for Los
 Alamos National Laboratory (LANL), which is operated by Triad
 National Security, LLC for the U.S. Department of Energy. 
 All rights in the program are reserved by Triad National Security,
 LLC, and the U.S. Department of Energy/National Nuclear Security
 Administration. The Government is granted for itself and others acting
 on its behalf a nonexclusive, paid-up, irrevocable worldwide license
 in this material to reproduce, prepare derivative works, distribute
 copies to the public, perform publicly and display publicly, and to
 permit others to do so

 
 This is open source software distributed under the 3-clause BSD license.
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:
 
 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
 3. Neither the name of Triad National Security, LLC, Los Alamos
    National Laboratory, LANL, the U.S. Government, nor the names of its
    contributors may be used to endorse or promote products derived from this
    software without specific prior written permission.

 
 THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND
 CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <iostream>
#include <cstdint>

#include "mpi.h"


#include <UnitTest++.h>

#include "../Mesh_MSTK.hh"

TEST(MSTK_HEX1)
{

  int i, j, k, err, nc, nv;
  std::vector<Jali::Entity_ID> faces(6), facenodes(4), cellnodes(8),
      expfacenodes(4);
  std::vector<Jali::dir_t> facedirs(6);
  std::vector<JaliGeometry::Point> ccoords(8), fcoords(4);

  int NV = 8;
  int NF = 6;
  int NC = 1;
  double xyz[12][3] = {{0, 0, 0},
                       {1, 0, 0},
                       {0, 1, 0},
                       {1, 1, 0},
                       {0, 0, 1},
                       {1, 0, 1},
                       {0, 1, 1},
                       {1, 1, 1}};
  Jali::Entity_ID local_cellnodes[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  Jali::Entity_ID local_facenodes[6][4] = {{0, 1, 5, 4},
                                           {1, 2, 6, 5},
                                           {2, 3, 7, 6},
                                           {3, 0, 4, 7},
                                           {0, 3, 2, 1},
                                           {4, 5, 6, 7}};


  // Load a single hex from the hex1.exo file

  Jali::Mesh *mesh(new Jali::Mesh_MSTK("test/hex_1x1x1_ss.exo",
                                       MPI_COMM_WORLD));


  // Check number of nodes and their coordinates

  nv = mesh->num_entities(Jali::Entity_kind::NODE,
                          Jali::Entity_type::PARALLEL_OWNED);
  CHECK_EQUAL(NV, nv);

  for (i = 0; i < nv; i++) {
    JaliGeometry::Point coords(mesh->space_dimension());

    //    coords.init(mesh->space_dimension());

    mesh->node_get_coordinates(i, &coords);
    CHECK_ARRAY_EQUAL(xyz[i], coords, 3);
  }


  // Check number of cells and their face nodes and their face coordinates

  nc = mesh->num_entities(Jali::Entity_kind::CELL,
                          Jali::Entity_type::PARALLEL_OWNED);
  CHECK_EQUAL(NC, nc);


  // Check cell coordinates directly

  mesh->cell_get_nodes(0, &cellnodes);
  mesh->cell_get_coordinates(0, &ccoords);

  for (j = 0; j < 8; j++) {
    CHECK_ARRAY_EQUAL(xyz[cellnodes[j]], ccoords[j], 3);
  }



  mesh->cell_get_faces_and_dirs(0, &faces, &facedirs, true);

  for (j = 0; j < 6; j++) {
    mesh->face_get_nodes(faces[j], &facenodes);
    mesh->face_get_coordinates(faces[j], &fcoords);


    for (k = 0; k < 4; k++)
      expfacenodes[k] = cellnodes[local_facenodes[j][k]];

    // The order of nodes returned may be different from what we expected
    // So make sure we have a matching node to start with

    int k0 = -1;
    int found = 0;
    for (k = 0; k < 4; k++) {
      if (expfacenodes[k] == facenodes[0]) {
        k0 = k;
        found = 1;
        break;
      }
    }

    CHECK_EQUAL(found, 1);

    if (facedirs[j] == 1) {
      for (k = 0; k < 4; k++) {
        CHECK_EQUAL(expfacenodes[(k0+k)%4], facenodes[k]);
        CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+k)%4]], fcoords[k], 3);
      }
    } else {
      for (k = 0; k < 4; k++) {
        CHECK_EQUAL(expfacenodes[(k0+4-k)%4], facenodes[k]);
        CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+4-k)%4]], fcoords[k], 3);
      }
    }

  }

}

