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
#include <vector>
#include "stdlib.h"
#include "math.h"

#include "UnitTest++.h"
#include "../Mesh_simple.hh"

SUITE(MeshSimple) {
  TEST(MAPS) {
    Jali::Mesh_simple Mm(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                         1, 1, 1, MPI_COMM_WORLD);

    double xc[] = { 2.0, 2.0, 2.0 };
    Mm.node_set_coordinates(7, xc);

    int expcellnodes[8] = {0, 1, 3, 2, 4, 5, 7, 6};
    double expnodecoords[8][3] = {{0.0, 0.0, 0.0},
                                  {1.0, 0.0, 0.0},
                                  {0.0, 1.0, 0.0},
                                  {1.0, 1.0, 0.0},
                                  {0.0, 0.0, 1.0},
                                  {1.0, 0.0, 1.0},
                                  {0.0, 1.0, 1.0},
                                  {2.0, 2.0, 2.0}};
    int expfacenodes[6][4] = {{0, 1, 3, 2},
                              {4, 5, 7, 6},
                              {0, 1, 5, 4},
                              {2, 3, 7, 6},
                              {0, 2, 6, 4},
                              {1, 3, 7, 5}};

    CHECK_EQUAL(1, Mm.num_entities(Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED));
    CHECK_EQUAL(6, Mm.num_entities(Jali::Entity_kind::FACE,
                                   Jali::Entity_type::PARALLEL_OWNED));
    CHECK_EQUAL(8, Mm.num_entities(Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED));

    std::vector<JaliGeometry::Point> x(8);
    std::vector<Jali::Entity_ID> nodes(8);
    std::vector<Jali::Entity_ID> faces(6);

    for (Jali::Entity_ID i = 0;
         i < Mm.num_entities(Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_OWNED);
         ++i) {
      Mm.cell_get_nodes(i, &nodes);

      CHECK_EQUAL(8, nodes.size());
      CHECK_ARRAY_EQUAL(expcellnodes, nodes, 8);

      for (int j = 0; j < 8; ++j) {
        Mm.node_get_coordinates(nodes[j], &(x[j]));
        CHECK_ARRAY_EQUAL(expnodecoords[expcellnodes[j]], x[j], 3);
      }

      Mm.cell_get_faces(i, &faces, true);
      for (int j = 0; j < 6; ++j) {
        Jali::Entity_ID_List fnodes;

        Mm.face_get_nodes(faces[j], &fnodes);
        CHECK_ARRAY_EQUAL(expfacenodes[faces[j]], fnodes, 4);

        Mm.face_get_coordinates(faces[j], &x);

        for (int k = 0; k < 4; ++k) {
          CHECK_ARRAY_EQUAL(expnodecoords[expfacenodes[faces[j]][k]], x[k], 3);
        }
      }


      Mm.cell_get_coordinates(i, &x);
      CHECK_EQUAL(8, x.size());
      for (int k = 0; k < 8; k++)
        CHECK_ARRAY_EQUAL(expnodecoords[expcellnodes[k]], x[k], 3);
    }
  }

  TEST(MAPS_1D) {
    // Make a 1 cell, 1d cartesian mesh
    std::vector<double> node_pts = {0.0, 1.0};
    Jali::Mesh_simple Mm(node_pts, MPI_COMM_WORLD, NULL,
                         true, true, true, true, true);

    // Move one of the nodes
    double xc[] = { 2.0 };
    Mm.node_set_coordinates(1, xc);

    // Indices of the nodes of this cell
    int expcellnodes[2] = {0, 1};
    // Coordinates of the nodes of this cell
    double expnodecoords[2][1] = {{0.0},
                                  {2.0}};
    // Indices for each face of this node
    int expfacenodes[2][1] = {{0},
                              {1}};

    // Make sure we have the right number of entities
    CHECK_EQUAL(1, Mm.num_entities(Jali::Entity_kind::CELL,
                                   Jali::Entity_type::PARALLEL_OWNED));
    CHECK_EQUAL(2, Mm.num_entities(Jali::Entity_kind::EDGE,
                                   Jali::Entity_type::PARALLEL_OWNED));
    CHECK_EQUAL(2, Mm.num_entities(Jali::Entity_kind::FACE,
                                   Jali::Entity_type::PARALLEL_OWNED));
    CHECK_EQUAL(2, Mm.num_entities(Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED));
    CHECK_EQUAL(4, Mm.num_entities(Jali::Entity_kind::WEDGE,
                                   Jali::Entity_type::PARALLEL_OWNED));
    CHECK_EQUAL(2, Mm.num_entities(Jali::Entity_kind::CORNER,
                                   Jali::Entity_type::PARALLEL_OWNED));

    std::vector<JaliGeometry::Point> x(2);
    std::vector<Jali::Entity_ID> nodes(2);
    std::vector<Jali::Entity_ID> faces(2);

    for (Jali::Entity_ID i = 0;
         i < Mm.num_entities(Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_OWNED);
         ++i) {
      Mm.cell_get_nodes(i, &nodes);

      CHECK_EQUAL(2, nodes.size());
      CHECK_ARRAY_EQUAL(expcellnodes, nodes, 2);

      // Check node coordinates
      for (int j = 0; j < 2; ++j) {
        Mm.node_get_coordinates(nodes[j], &(x[j]));
        CHECK_ARRAY_EQUAL(expnodecoords[expcellnodes[j]], x[j], 1);
      }

      // Check face coordinates
      Mm.cell_get_faces(i, &faces, true);
      for (int j = 0; j < 2; ++j) {
        Jali::Entity_ID_List fnodes;

        Mm.face_get_nodes(faces[j], &fnodes);
        CHECK_ARRAY_EQUAL(expfacenodes[faces[j]], fnodes, 1);

        Mm.face_get_coordinates(faces[j], &x);

        for (int k = 0; k < 1; ++k) {
          CHECK_ARRAY_EQUAL(expnodecoords[expfacenodes[faces[j]][k]], x[k], 1);
        }
      }


      Mm.cell_get_coordinates(i, &x);
      CHECK_EQUAL(2, x.size());
      for (int k = 0; k < 2; k++)
        CHECK_ARRAY_EQUAL(expnodecoords[expcellnodes[k]], x[k], 1);
    }
  }
}
