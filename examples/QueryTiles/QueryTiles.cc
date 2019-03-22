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

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

using namespace Jali;

// Fire up Jali, create a mesh, iterate through tiles (groups of
// elements) of the mesh and ask each tile about the local mesh
// topology.

// This program can be run in serial or parallel


int main(int argc, char *argv[]) {

  // Jali depends on MPI

  MPI_Init(&argc, &argv);

  // Create a mesh factory object - this object has methods for
  // specifying the preference of mesh frameworks and unified
  // interfaces for instantiating a mesh object of a particular
  // framework type

  MPI_Comm comm = MPI_COMM_WORLD;
  MeshFactory mesh_factory(comm);

  // Specify that MSTK is the preferred mesh framework. Currently Jali is
  // compiled only with MSTK and SimpleMesh support

  std::shared_ptr<Mesh> mymesh;  // Pointer to a mesh object
  bool parallel_mesh = false;
  int mesh_dimension = 3;
  if (framework_available(MSTK) &&
      framework_generates(MSTK, parallel_mesh, mesh_dimension)) {

    mesh_factory.framework(MSTK);
    
    // Create a 3D mesh from (0.0, 0.0, 0.0) to (1.0, 1.0, 1.0)
    // with 3, 3 and 3 elements in the X, Y and Z directions.
    // request faces, edges, wedges and corners and tiles

    mesh_factory.included_entities(Entity_kind::ALL_KIND);

    mesh_factory.num_tiles(4);
    mesh_factory.num_ghost_layers_tile(1);

    mymesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
  }

  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  // Iterate through tiles of the mesh and query the mesh topology

  int ntiles = 0;
  for (auto const& t : mymesh->tiles()) {

    std::cerr << "Processor " << rank << "    Tile " << ntiles << std::endl;

    // Cells of the mesh tile

    int ncells_owned = t->num_cells<Entity_type::PARALLEL_OWNED>();
    int ncells_ghost = t->num_cells<Entity_type::PARALLEL_GHOST>();

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Owned Cells " << ncells_owned << std::endl;
    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Ghost Cells " << ncells_ghost << std::endl;

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Cells:  ";
    for (auto const& c : t->cells())
      std::cerr << " " << c << " ";
    std::cerr << std::endl << std::endl;

    // Nodes of the mesh tile

    int nnodes_owned = t->num_nodes<Entity_type::PARALLEL_OWNED>();
    int nnodes_ghost = t->num_nodes<Entity_type::PARALLEL_GHOST>();

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Owned Nodes " << nnodes_owned << std::endl;
    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Ghost Nodes " << nnodes_ghost << std::endl;

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Nodes:  ";
    for (auto const& n : t->nodes())
      std::cerr << " " << n << " ";
    std::cerr << std::endl << std::endl;

    // Edges of the mesh tile

    int nedges_owned = t->num_edges<Entity_type::PARALLEL_OWNED>();
    int nedges_ghost = t->num_edges<Entity_type::PARALLEL_GHOST>();

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Owned Edges " << nedges_owned << std::endl;
    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Ghost Edges " << nedges_ghost << std::endl;
    std::cerr << std::endl;

    // Faces of the mesh tile

    int nfaces_owned = t->num_faces<Entity_type::PARALLEL_OWNED>();
    int nfaces_ghost = t->num_faces<Entity_type::PARALLEL_GHOST>();

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Owned Faces " << nfaces_owned << std::endl;
    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Ghost Faces " << nfaces_ghost << std::endl;
    std::cerr << std::endl;

    // Wedges of the mesh tile

    int nwedges_owned = t->num_wedges<Entity_type::PARALLEL_OWNED>();
    int nwedges_ghost = t->num_wedges<Entity_type::PARALLEL_GHOST>();

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Owned Wedges " << nwedges_owned << std::endl;
    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Ghost Wedges " << nwedges_ghost << std::endl;
    std::cerr << std::endl;

    // Corners of the mesh tile

    int ncorners_owned = t->num_corners<Entity_type::PARALLEL_OWNED>();
    int ncorners_ghost = t->num_corners<Entity_type::PARALLEL_GHOST>();

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Owned Corners " << ncorners_owned << std::endl;
    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Ghost Corners " << ncorners_ghost << std::endl;
    std::cerr << std::endl;

    std::cerr << std::endl;

    ntiles++;
  }

  // Clean up and exit

  MPI_Finalize();

}

