//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#include <iostream>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"

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

  FrameworkPreference pref;
  pref.push_back(MSTK);

  std::shared_ptr<Mesh> mymesh;  // Pointer to a mesh object
  if (framework_available(MSTK)) {  // check if framework is available
    mesh_factory.preference(pref);
  
    // Create a 3D mesh from (0.0, 0.0, 0.0) to (1.0, 1.0, 1.0)
    // with 3, 3 and 3 elements in the X, Y and Z directions. Specify
    // that we did not instantiate a geometric model (NULL). Also,
    // request faces, edges, wedges and corners (true, true, true,
    // true)

    int num_tiles_requested = 4;
    mymesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, NULL,
                          true, true, true, true, num_tiles_requested);
  }

  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  // Iterate through tiles of the mesh and query the mesh topology

  int ntiles = 0;
  for (auto const& t : mymesh->tiles()) {

    std::cerr << "Processor " << rank << "    Tile " << ntiles << std::endl;

    // Cells of the mesh tile

    int ncells_owned = t->num_cells<Parallel_type::OWNED>();
    int ncells_ghost = t->num_cells<Parallel_type::GHOST>();

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

    int nnodes_owned = t->num_nodes<Parallel_type::OWNED>();
    int nnodes_ghost = t->num_nodes<Parallel_type::GHOST>();

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

    int nedges_owned = t->num_edges<Parallel_type::OWNED>();
    int nedges_ghost = t->num_edges<Parallel_type::GHOST>();

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Owned Edges " << nedges_owned << std::endl;
    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Ghost Edges " << nedges_ghost << std::endl;
    std::cerr << std::endl;

    // Faces of the mesh tile

    int nfaces_owned = t->num_faces<Parallel_type::OWNED>();
    int nfaces_ghost = t->num_faces<Parallel_type::GHOST>();

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Owned Faces " << nfaces_owned << std::endl;
    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Ghost Faces " << nfaces_ghost << std::endl;
    std::cerr << std::endl;

    // Wedges of the mesh tile

    int nwedges_owned = t->num_wedges<Parallel_type::OWNED>();
    int nwedges_ghost = t->num_wedges<Parallel_type::GHOST>();

    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Owned Wedges " << nwedges_owned << std::endl;
    std::cerr << "Processor " << rank << "    Tile " << ntiles <<
        "      Num Ghost Wedges " << nwedges_ghost << std::endl;
    std::cerr << std::endl;

    // Corners of the mesh tile

    int ncorners_owned = t->num_corners<Parallel_type::OWNED>();
    int ncorners_ghost = t->num_corners<Parallel_type::GHOST>();

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

