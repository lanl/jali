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

// Fire up Jali, create a mesh and ask the mesh a very basic question

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
  // compiled only with MSTK support

  FrameworkPreference pref;
  pref.push_back(MSTK);

  std::shared_ptr<Mesh> mymesh;  // Pointer to a mesh object
  if (framework_available(MSTK)) {  // check if framework is available
    mesh_factory.preference(pref);
  
    // Read in an exodus file.
    // request faces, edges, wedges and corners

    mesh_factory.included_entities(Entity_kind::ALL_KIND);

    mymesh = mesh_factory("quadtri.exo");
  }


  // Print out the spatial dimension of the mesh

  std::cerr << "Spatial dimension of mesh: " << mymesh->space_dimension() <<
      std::endl;
  
  // Print out the number of cells in the mesh

  std::cerr << "Number of mesh cells: " <<
    mymesh->num_entities(Entity_kind::CELL, Entity_type::PARALLEL_ALL) << std::endl;

  // Print out the number of nodes in the mesh

  std::cerr << "Number of mesh nodes: " <<
    mymesh->num_entities(Entity_kind::NODE, Entity_type::PARALLEL_ALL) << std::endl;


  // Clean up and exit

  MPI_Finalize();

}

