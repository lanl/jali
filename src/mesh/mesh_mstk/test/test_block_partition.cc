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

#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

// Special class just for testing protected function of the base mesh
// class - block_partition_regular_mesh

namespace Jali {
class Mesh_MSTK_Test_Protected : public Mesh_MSTK {
 public:

  // Simplified constructor that generates a 3D mesh

  Mesh_MSTK_Test_Protected(double const x0, double const y0, double const z0,
                           double const x1, double const y1, double const z1,
                           int const nx, int const ny, int const nz,
                           MPI_Comm const& incomm) :
      Mesh_MSTK(x0, y0, z0, x1, y1, z1, nx, ny, nz, incomm) {}

  // Destructor

  ~Mesh_MSTK_Test_Protected() {}

  // function to test
  int block_partition_regular_mesh(int dim, double *domain,
                                   int *num_cells_in_dir,
                                   int num_blocks_requested,
                                   std::vector<std::array<double, 6>> *blocklimits,
                                   std::vector<std::array<int, 3>> *blocknumcells) {
    return Mesh_MSTK::block_partition_regular_mesh(dim, domain,
                                                   num_cells_in_dir,
                                                   num_blocks_requested,
                                                   blocklimits, blocknumcells);
  }

};  // end class Mesh_MSTK_Test_Protected
}  // end namespace Jali

int check_block_partitioning(int const dim, double const * const domain,
                             int const * const num_cells_in_dir,
                             int const num_blocks,
                             std::vector<std::array<double, 6>> const& blocklimits,
                             std::vector<std::array<int, 3>> const& blocknumcells) {

  if (blocklimits.size() != num_blocks) return 0;
  if (blocknumcells.size() != num_blocks) return 0;

  // Minimum sanity check - sum of volumes of blocks should be the
  // volume of the original domain - doesn't preclude blocks on top of
  // each other

  double total_volume = 1.0;
  for (int dir = 0; dir < dim; dir++)
    total_volume *= (domain[2*dir+1]-domain[2*dir]);
  double volume_sum = 0.0;
  for (int ib = 0; ib < num_blocks; ib++) {
    double block_volume = 1.0;
    for (int dir = 0; dir < dim; dir++)
      block_volume *= (blocklimits[ib][2*dir+1] - blocklimits[ib][2*dir]);
    volume_sum += block_volume;
  }
  double reldiff = fabs(volume_sum-total_volume)/total_volume;
  if (reldiff > 1.0e-09) {
    CHECK(reldiff > 1.0e-09); // so it prints a message
    return 0;
  }

  // Sum of the number of cells in blocks should be equal to each other
  int ncells = 1;
  for (int dir = 0; dir < dim; dir++)
    ncells *= num_cells_in_dir[dir];
  int ncells_sum = 0.0;
  for (int ib = 0; ib < num_blocks; ib++) {
    int ncells_block = 1;
    for (int dir = 0; dir < dim; dir++)
      ncells_block *= blocknumcells[ib][dir];
    ncells_sum += ncells_block;
  }
  if (ncells_sum != ncells) {
    CHECK_EQUAL(ncells_sum, ncells);  // so it prints out a message
    return 0;
  }



  // Extended checks

  // Make sure that no block overlaps another

  for (int ib = 0; ib < num_blocks; ib++) {
    bool found_overlap = false;
    for (int ib2 = 0; ib2 < num_blocks; ib2++) {
      if (ib2 == ib) continue;
      int noverlapdir = 0;
      for (int dir = 0; dir < dim; dir++)
        if ((blocklimits[ib][2*dir] < blocklimits[ib2][2*dir+1] &&
             blocklimits[ib][2*dir] > blocklimits[ib2][2*dir]) ||
            (blocklimits[ib][2*dir+1] > blocklimits[ib2][2*dir] &&
             blocklimits[ib][2*dir+1] < blocklimits[ib2][2*dir]))
          noverlapdir++;
      if (noverlapdir == dim) {
        found_overlap = true;
        std::cerr << "Blocks " << ib << " and " << ib2 << " overlap \n";
        std::cerr << " Block " << ib << " extents are (";
        for (int dir1 = 0; dir1 < dim-1; dir1++)
          std::cerr << blocklimits[ib][2*dir1] << ",";
        std::cerr << blocklimits[ib][2*(dim-1)] << ") ";
        for (int dir1 = 0; dir1 < dim-1; dir1++)
          std::cerr << blocklimits[ib][2*dir1+1] << ",";
        std::cerr << blocklimits[ib][2*dim-1] << ") ";
        break;
      }
      CHECK(!found_overlap);
    }
  }


  // Make sure that each block has neighbors on the interior side

  for (int ib = 0; ib < num_blocks; ib++) {
    for (int dir = 0; dir < dim; dir++) {
      for (int k = 0; k < 2; k++) {

        // If it is external boundary face, there will be no neighbors

        if (blocklimits[ib][2*dir+k] == domain[2*dir] ||
            blocklimits[ib][2*dir+k] == domain[2*dir+1]) continue;

        double block_face_coords[3][3];
        block_face_coords[0][0] = blocklimits[ib][2*dir+k];
        block_face_coords[0][1] = blocklimits[ib][2*((dir+1)%3)];
        block_face_coords[0][2] = blocklimits[ib][2*((dir+2)%3)];
        block_face_coords[1][0] = blocklimits[ib][2*dir+k];
        block_face_coords[1][1] = blocklimits[ib][2*((dir+1)%3)+1];
        block_face_coords[1][2] = blocklimits[ib][2*((dir+2)%3)+1];

        bool neighbor_found = false;
        for (int ib1 = 0; ib1 < num_blocks; ib1++) {
          if (ib == ib1) continue;

          double nbr_block_face_coords[3][3];
          nbr_block_face_coords[0][0] = blocklimits[ib][2*dir+k%2];
          nbr_block_face_coords[0][1] = blocklimits[ib][2*((dir+1)%3)];
          nbr_block_face_coords[0][2] = blocklimits[ib][2*((dir+2)%3)];
          nbr_block_face_coords[1][0] = blocklimits[ib][2*dir+k%2];
          nbr_block_face_coords[1][1] = blocklimits[ib][2*((dir+1)%3)+1];
          nbr_block_face_coords[1][2] = blocklimits[ib][2*((dir+2)%3)+1];

          if (block_face_coords[0][0] == nbr_block_face_coords[0][0] &&
              block_face_coords[0][1] == nbr_block_face_coords[0][1] &&
              block_face_coords[0][2] == nbr_block_face_coords[0][2] &&
              block_face_coords[1][0] == nbr_block_face_coords[1][0] &&
              block_face_coords[1][1] == nbr_block_face_coords[1][1] &&
              block_face_coords[1][2] == nbr_block_face_coords[1][2]) {
            neighbor_found = true;
            break;
          }
        }

        CHECK(neighbor_found);
        if (!neighbor_found) {
          std::cerr << "Could not find neighbor for block " << "[(" <<
              blocklimits[ib][0] << "," << blocklimits[ib][2] << "," <<
              blocklimits[ib][4] << "), (" <<
              blocklimits[ib][1] << "," << blocklimits[ib][3] << "," <<
              blocklimits[ib][5] << ")] across face [(" <<
              block_face_coords[0][0] << "," << block_face_coords[0][1] <<
              "," << block_face_coords[0][2] << "), (" <<
              block_face_coords[1][0] << "," << block_face_coords[1][1] <<
              "," << block_face_coords[1][2] << ")]\n";
          return 0;
        }
      }
    }
  }

  return 1;
}

TEST(BLOCK_PARTITION) {

  // Generate a mesh consisting of 3x3x3 elements or whatever -
  // doesn't matter since we won't get a partitioning based on the
  // actual mesh - its just a way to access the method

  Jali::Mesh_MSTK_Test_Protected *mesh =
      new Jali::Mesh_MSTK_Test_Protected(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                         3, 3, 3, MPI_COMM_WORLD);

  {
    // Partition a 1D domain
    int dim = 1;

    std::vector<std::array<double, 6>> blocklimits;
    std::vector<std::array<int, 3>> blocknumcells;

    // Domain is [0.0,1.0]
    double domain[6] = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

    // 16 cells to force regular partition
    int num_cells_in_dir[3] = {16, 0, 0};

    int nparts = 4;

    int ok = mesh->block_partition_regular_mesh(dim, domain, num_cells_in_dir,
                                               nparts,
                                               &blocklimits, &blocknumcells);
    CHECK(ok);

    ok = check_block_partitioning(dim, domain, num_cells_in_dir, nparts,
                                  blocklimits, blocknumcells);
    CHECK(ok);

    // 15 cells to force irregular partition
    num_cells_in_dir[0] = 15;

    ok = mesh->block_partition_regular_mesh(dim, domain, num_cells_in_dir,
                                           nparts,
                                           &blocklimits, &blocknumcells);
    CHECK(ok);

    ok = check_block_partitioning(dim, domain, num_cells_in_dir, nparts,
                                  blocklimits, blocknumcells);
    CHECK(ok);
  }

  {
    // Partition a 2D domain
    int dim = 2;

    std::vector<std::array<double, 6>> blocklimits;
    std::vector<std::array<int, 3>> blocknumcells;

    // Domain is [0.0,1.0]x[0.0,1.0]
    double domain[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 0.0};

    // 10 cells to force regular partition
    int num_cells_in_dir[3] = {10, 10, 0};

    int nparts = 4;

    int ok = mesh->block_partition_regular_mesh(dim, domain, num_cells_in_dir,
                                                nparts,
                                                &blocklimits, &blocknumcells);
    CHECK(ok);

    ok = check_block_partitioning(dim, domain, num_cells_in_dir, nparts,
                                  blocklimits, blocknumcells);
    CHECK(ok);

    // larger number of partitions to force irregular partitioning
    nparts = 16;

    ok = mesh->block_partition_regular_mesh(dim, domain, num_cells_in_dir,
                                            nparts,
                                            &blocklimits, &blocknumcells);
    CHECK(ok);

    ok = check_block_partitioning(dim, domain, num_cells_in_dir, nparts,
                                  blocklimits, blocknumcells);
    CHECK(ok);
  }

  {
    // Partition a 3D domain
    int dim = 3;

    std::vector<std::array<double, 6>> blocklimits;
    std::vector<std::array<int, 3>> blocknumcells;

    // Domain is [0.0,1.0]x[0.0,1.0]x[0.0,2.0]
    double domain[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 2.0};

    // 10 cells to force regular partition
    int num_cells_in_dir[3] = {10, 10, 10};

    int nparts = 8;

    int ok = mesh->block_partition_regular_mesh(dim, domain, num_cells_in_dir,
                                                nparts,
                                                &blocklimits, &blocknumcells);
    CHECK(ok);

    ok = check_block_partitioning(dim, domain, num_cells_in_dir, nparts,
                                  blocklimits, blocknumcells);
    CHECK(ok);

    // More partitions and irregular cells to force irregular partition
    num_cells_in_dir[0] = 15;
    nparts = 32;

    ok = mesh->block_partition_regular_mesh(dim, domain, num_cells_in_dir,
                                            nparts,
                                            &blocklimits, &blocknumcells);
    CHECK(ok);

    ok = check_block_partitioning(dim, domain, num_cells_in_dir, nparts,
                                  blocklimits, blocknumcells);
    CHECK(ok);
  }

  delete(mesh);
}

