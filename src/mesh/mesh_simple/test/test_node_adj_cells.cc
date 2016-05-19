#include <iostream>
#include <vector>
#include "stdlib.h"
#include "math.h"

#include "UnitTest++.h"
#include "../Mesh_simple.hh"

TEST(NODE_ADJ_CELLS) {
  const unsigned int exp_ncell = 27;

  // Test 4 representative cells - one at a boundary corner, one at a
  // boundary edge, one on a boundary face and one in the interior
  const unsigned int ntestcells = 4;
  const int test_cells[4] = {6, 3, 14, 13};
  const unsigned int exp_nadj[4] = {7, 11, 17, 26};
  const int exp_adjcells[4][26] = {
      {3, 12, 15, 16, 7, 13, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      { 0,  9, 12, 15,  6,  1, 10,  4, 13,  7, 16, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      { 2,  5,  8, 11, 17, 20, 23, 26,  1,  4,  7, 10, 13, 16, 19, 22, 25, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26}};

  Jali::Mesh_simple Mm(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                       3, 3, 3, MPI_COMM_WORLD);

  CHECK_EQUAL(exp_ncell, Mm.num_entities(Jali::Entity_kind::CELL,
                                         Jali::Parallel_type::OWNED));

  for (int i = 0; i < ntestcells; ++i) {
    Jali::Entity_ID_List adjcells;
    
    Mm.cell_get_node_adj_cells(test_cells[i], Jali::Parallel_type::OWNED,
                               &adjcells);
    
    unsigned int nadj = adjcells.size();
    CHECK_EQUAL(exp_nadj[i], nadj);
    
    for (int j = 0; j < nadj; ++j) {
      bool found = false;
      for (int k = 0; k < nadj; ++k)
        if (exp_adjcells[i][k] == adjcells[j]) {
          found = true;
          break;
        }
      CHECK(found);
    }
  }
}


