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

#include <vector>
#include <array>
#include <algorithm>
#include <cassert>


namespace Jali {


/*
  Get the partitioning of a regular mesh such that each
  partition is a rectangular block

  dim                   Dimension of problem - 1, 2 or 3
  domain                2*dim values for min/max of domain
                        (xmin, xmax, ymin, ymax, zmin, zmax)
  num_cells_in_dir      number of cells in each direction
  num_blocks_requested  number of blocks requested
  block_start_indices   start node/cell indices for each block
  block_num_cells num   cells in each direction for blocks

  Returns 1 if successful, 0 otherwise
*/

int block_partition_regular_mesh(int const dim,
                                 double const * const domain,
                                 int const * const num_cells_in_dir,
                                 int const num_blocks_requested,
                                 std::vector<std::array<int, 3>> *block_start_indices,
                                 std::vector<std::array<int, 3>> *block_num_cells) {

  // Create local block arrays that are larger than the requested number of
  // blocks. The extra storage is for temporary blocks while we are subdividing

  int nblocks = 1;
  int nblocks_in_dir[3] = {1, 1, 1};
  int nblockcells_in_dir[3] = {0, 0, 0};
  for (int i = 0; i < dim; i++)
    nblockcells_in_dir[i] = num_cells_in_dir[i];

  int ncells_total = 1;
  for (int i = 0; i < dim; i++)
    ncells_total *= num_cells_in_dir[i];
  assert (num_blocks_requested <= ncells_total);
    
  
  if (num_blocks_requested > 1) {
    // First try bisection
    
    bool done = false;
    while (!done) {
      bool bisected = false;
      if (nblockcells_in_dir[0]%2 == 0) {  // even number of cells in x
        if (2*nblocks <= num_blocks_requested) {
          nblockcells_in_dir[0] /= 2;
          nblocks *= 2;
          nblocks_in_dir[0] *= 2;
          bisected = true;
          
          if (nblocks == num_blocks_requested) {
            done = 1;
            continue;
          }
        }
      }
      if (dim > 1 && nblockcells_in_dir[1]%2 == 0) {  // even num of cells in y
        if (2*nblocks <= num_blocks_requested) { 
          nblockcells_in_dir[1] /= 2;
          nblocks *= 2;
          nblocks_in_dir[1] *= 2;
          
          bisected = true;
          if (nblocks == num_blocks_requested) {
            done = 1;
            continue;
          }
        }
      }
      if (dim > 2 && nblockcells_in_dir[2]%2 == 0) {  // even num of cells in z
        if (2*nblocks <= num_blocks_requested) {
          nblockcells_in_dir[2] /= 2;
          nblocks *= 2;
          nblocks_in_dir[2] *= 2;
          bisected = true;
          
          if (nblocks == num_blocks_requested) {
            done = 1;
            continue;
          }
        }
      }
      if (!bisected) {
        // Reached a state where we cannot evenly subdivide the number
        // of elements in any one direction
        done = 1;
      }
    }
  }

  std::array<int, 3> iarray0 = {0, 0, 0};
  std::vector<std::array<int, 3>>
      bnumcells(8*nblocks_in_dir[0]*nblocks_in_dir[1]*nblocks_in_dir[2],
                iarray0);
  std::vector<std::array<int, 3>>
      blimits(8*nblocks_in_dir[0]*nblocks_in_dir[1]*nblocks_in_dir[2],
              iarray0);

  // Populate the block details
  int partnum = 0;
  for (int i = 0; i < nblocks_in_dir[0]; i++) {
    for (int j = 0; j < nblocks_in_dir[1]; j++) {
      for (int k = 0; k < nblocks_in_dir[2]; k++) {
        blimits[partnum][0] = i*nblockcells_in_dir[0];
        bnumcells[partnum][0] = nblockcells_in_dir[0];

        if (dim > 1) {
          blimits[partnum][1] = j*nblockcells_in_dir[1];
          bnumcells[partnum][1] = nblockcells_in_dir[1];

          if (dim > 2) {
            blimits[partnum][2] = k*nblockcells_in_dir[2];
            bnumcells[partnum][2] = nblockcells_in_dir[2];
          }
        }
        partnum++;
      }
    }
  }

  if (nblocks != num_blocks_requested) {  // assuming that nblocks <= num_blocks_requested

    // Start dividing the blocks as unevenly to get the number of
    // partitions we need
    
    bool done = false;
    while (!done) {
      for (int dir = 0; dir < dim; dir++) {  // split blocks in x, then y etc
        int nnewblocks = 0;
        for (int ib = 0; ib < nblocks; ib++) {
          if (bnumcells[ib][0] == 0 && bnumcells[ib][1] == 0 &&
              bnumcells[ib][2] == 0) continue;  // Blanked out block
          
          int ncells1 = bnumcells[ib][dir]/2;
          int ncells2 = bnumcells[ib][dir] - ncells1;
          if (ncells1 == 0 || ncells2 == 0) continue;
          
          // Make two new partitions at the end of the list
          // But first make sure there is enough room to hold them

          int nalloc = bnumcells.size();
          if (nalloc < nblocks + nnewblocks + 2) {
            nalloc *= 2;
            bnumcells.resize(nalloc, iarray0);
            blimits.resize(nalloc, iarray0);
          }
          
          // Copy the original block details over to the first child block
          int ib1 = nblocks + nnewblocks;
          for (int dir1 = 0; dir1 < dim; dir1++) {
            bnumcells[ib1][dir1] = bnumcells[ib][dir1];
            blimits[ib1][dir1] = blimits[ib][dir1];
          }
          /* overwrite the data in the direction of the refinement */
          bnumcells[ib1][dir] = ncells1;
          blimits[ib1][dir] = blimits[ib][dir];
          
          // copy the initial block details over to the second child block
          int ib2 = nblocks + nnewblocks + 1;
          for (int dir1 = 0; dir1 < dim; dir1++) {
            bnumcells[ib2][dir1] = bnumcells[ib][dir1];
            blimits[ib2][dir1] = blimits[ib][dir1];
          }
          /* overwrite the data in the direction of the refinement */
          bnumcells[ib2][dir] = ncells2;
          blimits[ib2][dir] = blimits[ib][dir] + ncells1;
          
          // Blank out the parent block
          bnumcells[ib][0] = bnumcells[ib][1] = bnumcells[ib][2] = 0;
          nnewblocks += 2;
          
          // Check if we reached the requested number of blocks. Each
          // block that was split into two will cause the loss of one
          // block and gain of two new blocks, so the net gain is just
          // nnewblocks/2
          
          if (nblocks + nnewblocks/2 >= num_blocks_requested) {
            done = 1;
            break;
          }
        }
        
        // Squeeze out the blocks that were split and dummied out
        for (int ib = nblocks-1; ib >= 0; ib--) {
          if (bnumcells[ib][0] == 0 && bnumcells[ib][1] == 0 &&
              bnumcells[ib][2] == 0) {  // dummy block
            for (int ib1 = ib; ib1 < nblocks+nnewblocks-1; ib1++) {
              for (int dir1 = 0; dir1 < 3; dir1++) {
                bnumcells[ib1][dir1] = bnumcells[ib1+1][dir1];
                blimits[ib1][dir1] = blimits[ib1+1][dir1];
              }
            }
            // Blank out the last block
            bnumcells[nblocks+nnewblocks-1][0] =
                bnumcells[nblocks+nnewblocks-1][1] =
                bnumcells[nblocks+nnewblocks-1][2] = 0;
            nblocks--;
          }
        }
        nblocks += nnewblocks;
        if (done)
          break;
      }
    }
  }  // if (nblocks != num_blocks_requested)

  blimits.resize(nblocks);
  bnumcells.resize(nblocks);

  block_start_indices->resize(nblocks);
  block_num_cells->resize(nblocks);
  std::copy(blimits.begin(), blimits.end(), block_start_indices->begin());
  std::copy(bnumcells.begin(), bnumcells.end(), block_num_cells->begin());

  return 1;
}

}  // close namespace Jali
