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

TEST(FACE_ADJ_CELLS) {
  const unsigned int exp_ncell = 27, exp_nface = 108, exp_nnode = 64;
  const unsigned int exp_nadj[27] = {3, 4, 3,  4, 5, 4,  3, 4, 3,
                                     4, 5, 4,  5, 6, 5,  4, 5, 4,
                                     3, 4, 3,  4, 5, 4,  3, 4, 3};
  const int exp_adjcells[27][6] = {{  1,  3,  9, -1, -1, -1},
                                   {  2,  4,  0, 10, -1, -1},
                                   {  5,  1, 11, -1, -1, -1},
                                   {  0,  4,  6, 12, -1, -1},
                                   {  1,  5,  7,  3, 13, -1},
                                   {  2,  8,  4, 14, -1, -1},
                                   {  3,  7, 15, -1, -1, -1},
                                   {  4,  8,  6, 16, -1, -1},
                                   {  5,  7, 17, -1, -1, -1},

                                   { 10, 12,  0, 18, -1, -1},
                                   { 11, 13,  9,  1, 19, -1},
                                   { 14, 10,  2, 20, -1, -1},
                                   {  9, 13, 15,  3, 21, -1},
                                   { 10, 14, 16, 12,  4, 22},
                                   { 11, 17, 13,  5, 23, -1},
                                   { 12, 16,  6, 24, -1, -1},
                                   { 13, 17, 15,  7, 25, -1},
                                   { 14, 16,  8, 26, -1, -1},

                                   { 19, 21,  9, -1, -1, -1},
                                   { 20, 22, 18, 10, -1, -1},
                                   { 23, 19, 11, -1, -1, -1},
                                   { 18, 22, 24, 12, -1, -1},
                                   { 19, 23, 25, 21, 13, -1},
                                   { 20, 26, 22, 14, -1, -1},
                                   { 21, 25, 15, -1, -1, -1},
                                   { 22, 26, 24, 16, -1, -1},
                                   { 23, 25, 17, -1, -1, -1}};


  Jali::Mesh_simple Mm(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                       3, 3, 3, MPI_COMM_WORLD);

  CHECK_EQUAL(exp_ncell, Mm.num_entities(Jali::Entity_kind::CELL,
                                         Jali::Entity_type::PARALLEL_OWNED));
  CHECK_EQUAL(exp_nface, Mm.num_entities(Jali::Entity_kind::FACE,
                                         Jali::Entity_type::PARALLEL_OWNED));
  CHECK_EQUAL(exp_nnode, Mm.num_entities(Jali::Entity_kind::NODE,
                                         Jali::Entity_type::PARALLEL_OWNED));


  for (int i = 0; i < exp_ncell; ++i) {
      Jali::Entity_ID_List adjcells;

      Mm.cell_get_face_adj_cells(i, Jali::Entity_type::PARALLEL_OWNED, &adjcells);

      unsigned int nadj = adjcells.size();
      CHECK_EQUAL(exp_nadj[i], nadj);

      for (int j = 0; j < nadj; ++j)
        CHECK_EQUAL(exp_adjcells[i][j], adjcells[j]);
    }
}

TEST(FACE_ADJ_CELLS_1D) {
  // Make a mesh with 4 cells, which is 5 faces and 5 nodes
  const unsigned int exp_ncell = 4, exp_nface = 5, exp_nnode = 5;
  // Number of cells adjacent to each cell
  const unsigned int exp_nadj[exp_ncell] = {1, 2, 2, 1};
  // Indices of adjacent cells, for each cell
  const int exp_adjcells[exp_ncell][2] = {{  1,  -1},
                                          {  0,  2},
                                          {  1,  3},
                                          {  2, -1}};


  std::vector<double> node_pts = {0.0, 0.25, 0.5, 0.75, 1.0};
  Jali::Mesh_simple Mm(node_pts, MPI_COMM_WORLD);

  CHECK_EQUAL(exp_ncell, Mm.num_entities(Jali::Entity_kind::CELL,
                                         Jali::Entity_type::PARALLEL_OWNED));
  CHECK_EQUAL(exp_nface, Mm.num_entities(Jali::Entity_kind::FACE,
                                         Jali::Entity_type::PARALLEL_OWNED));
  CHECK_EQUAL(exp_nnode, Mm.num_entities(Jali::Entity_kind::NODE,
                                         Jali::Entity_type::PARALLEL_OWNED));


  for (int i = 0; i < exp_ncell; ++i) {
      Jali::Entity_ID_List adjcells;

      Mm.cell_get_face_adj_cells(i, Jali::Entity_type::PARALLEL_OWNED, &adjcells);

      unsigned int nadj = adjcells.size();
      CHECK_EQUAL(exp_nadj[i], nadj);

      for (int j = 0; j < nadj; ++j)
        CHECK_EQUAL(exp_adjcells[i][j], adjcells[j]);
    }
}

