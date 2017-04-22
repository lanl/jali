/*
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.

Copyright 2017. Los Alamos National Security, LLC. This software was
produced under U.S. Government contract DE-AC52-06NA25396 for Los
Alamos National Laboratory (LANL), which is operated by Los Alamos
National Security, LLC for the U.S. Department of Energy. The
U.S. Government has rights to use, reproduce, and distribute this
software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so
as not to confuse it with the version available from LANL.
 
Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

1.  Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3.  Neither the name of Los Alamos National Security, LLC, Los Alamos
National Laboratory, LANL, the U.S. Government, nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
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

TEST(NODE_CELL_FACES) {
  const unsigned int exp_nnode = 27;

  Jali::Mesh_simple Mm(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2, MPI_COMM_WORLD);


  for (int i = 0; i < exp_nnode; ++i) {
      Jali::Entity_ID node = i;

      Jali::Entity_ID_List cells;

      Mm.node_get_cells(node, Jali::Entity_type::PARALLEL_OWNED, &cells);

      unsigned int ncells = cells.size();

      for (int j = 0; j < ncells; ++j) {
        Jali::Entity_ID cell = cells[j];

        Jali::Entity_ID_List faces;

        Mm.node_get_cell_faces(node, cell, Jali::Entity_type::PARALLEL_OWNED, &faces);

        // This is a hex mesh. In any given cell, number of faces
        // connected to a node should be 3

        CHECK_EQUAL(3, faces.size());

        for (int k = 0; k < 3; ++k) {
          Jali::Entity_ID face = faces[k];

          Jali::Entity_ID_List fnodes;

          Mm.face_get_nodes(face, &fnodes);

          unsigned int nfnodes = fnodes.size();

          unsigned int found = 0;

          for (int n = 0; n < nfnodes; ++n) {
            if (fnodes[n] == node) {
              found = 1;
              break;
            }
          }

          CHECK_EQUAL(1, found);
        }
      }
  }
}

TEST(NODE_CELL_FACES_1D) {
  // Create a 1d mesh with 2 cells, 3 nodes and 4 faces
  // Number of nodes
  const unsigned int exp_nnode = 2;

  std::vector<double> node_pts = {0.0, 0.5, 1.0};
  Jali::Mesh_simple Mm(node_pts, MPI_COMM_WORLD);


  for (int i = 0; i < exp_nnode; ++i) {
      Jali::Entity_ID node = i;

      Jali::Entity_ID_List cells;

      Mm.node_get_cells(node, Jali::Entity_type::PARALLEL_OWNED, &cells);

      unsigned int ncells = cells.size();

      for (int j = 0; j < ncells; ++j) {
        Jali::Entity_ID cell = cells[j];

        Jali::Entity_ID_List faces;

        Mm.node_get_cell_faces(node, cell, Jali::Entity_type::PARALLEL_OWNED, &faces);

        // In 1d, a node and a face are the same, number of faces connected to a
        // node should be 1

        CHECK_EQUAL(1, faces.size());

        for (int k = 0; k < 1; ++k) {
          Jali::Entity_ID face = faces[k];

          Jali::Entity_ID_List fnodes;

          Mm.face_get_nodes(face, &fnodes);

          unsigned int nfnodes = fnodes.size();

          unsigned int found = 0;

          for (int n = 0; n < nfnodes; ++n) {
            if (fnodes[n] == node) {
              found = 1;
              break;
            }
          }

          CHECK_EQUAL(1, found);
        }
      }
  }
}

