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
                                         Jali::Entity_type::PARALLEL_OWNED));

  for (int i = 0; i < ntestcells; ++i) {
    Jali::Entity_ID_List adjcells;
    
    Mm.cell_get_node_adj_cells(test_cells[i], Jali::Entity_type::PARALLEL_OWNED,
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


