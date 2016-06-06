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

