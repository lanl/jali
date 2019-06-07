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

// Test for generation of hex mesh distributed over 4 processors

TEST(MSTK_HEX_GEN_3x3x3_4P) {
  int rank, size;

  int initialized;
  MPI_Initialized(&initialized);

  if (!initialized)
    MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  CHECK_EQUAL(4, size);

  if (size != 4) {
    std::cerr << "Test must be run with 4 processors" << std::endl;
  }

  // Create a 3x3x3 cell hex mesh distributed over 4 processors

  Jali::Mesh_MSTK parmesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, MPI_COMM_WORLD);

  // Create a serial version of the same mesh
  MPI_Comm serialcomm;
  MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &serialcomm);
  Jali::Mesh_MSTK serialmesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, serialcomm);


  // For each mesh vertex, find a vertex with the same GID in the
  // serial mesh and compare its coordinates.

  int NN = serialmesh.num_nodes();
  std::vector<JaliGeometry::Point> serialpnts(NN);
  for (int i = 0; i < NN; i++) {
    int gid = serialmesh.GID(i, Jali::Entity_kind::NODE);
    JaliGeometry::Point pnt;
    serialmesh.node_get_coordinates(i, &pnt);
    serialpnts[gid] = pnt;
  }

  int nn = parmesh.num_nodes();
  for (int i = 0; i < nn; i++) {
    int gid = parmesh.GID(i, Jali::Entity_kind::NODE);
    JaliGeometry::Point pnt;
    parmesh.node_get_coordinates(i, &pnt);
    JaliGeometry::Point diff = pnt - serialpnts[gid];
    if (JaliGeometry::L22(diff) > 1.0e-16) {
      std::cout << "Point " << gid << " in local mesh has coordinates " << pnt << "\n";
      std::cout << "Point " << gid << " in GLOBAL mesh has coordinates " << serialpnts[gid] << "\n";
    }
    //    CHECK(JaliGeometry::L22(diff) < 1.0e-16);
  }
  
  // For each mesh cell, find a cell with the same GID in the
  // serial mesh and compare its centroid.

  int NC = serialmesh.num_cells();
  serialpnts.resize(NC);
  for (int i = 0; i < NC; i++) {
    int gid = serialmesh.GID(i, Jali::Entity_kind::CELL);
    serialpnts[gid] = serialmesh.cell_centroid(i);
  }

  int nc = parmesh.num_cells();
  for (int i = 0; i < nc; i++) {
    int gid = parmesh.GID(i, Jali::Entity_kind::CELL);
    JaliGeometry::Point pnt = parmesh.cell_centroid(i);
    JaliGeometry::Point diff = pnt - serialpnts[gid];
    CHECK(JaliGeometry::L22(diff) < 1.0e-16);
  }
  
}

