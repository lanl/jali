//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#include <iostream>
#include <cstdint>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

using namespace Jali;

// Fire up Jali, create a mesh and ask the mesh about cell topology

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

  std::shared_ptr<Mesh> mymesh;  // Pointer to a mesh object
  bool parallel_mesh = false;
  int mesh_dimension = 3;
  if (framework_available(MSTK) &&
      framework_generates(MSTK, parallel_mesh, mesh_dimension)) {

    mesh_factory.framework(MSTK);
  
    // Create a 3D mesh from (0.0, 0.0, 0.0) to (1.0, 1.0, 1.0)
    // with 3, 3 and 3 elements in the X, Y and Z directions.
    // Request faces, edges, wedges and corners

    mesh_factory.included_entities(Entity_kind::ALL_KIND);

    mymesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
  }


  // Iterate through the cells of the mesh and get the faces and of the cell

  for (auto c : mymesh->cells<Entity_type::PARALLEL_OWNED>()) {

    Entity_ID_List cfaces;   // Entity_ID_List is just a std::vector<Entity_ID>
    std::vector<dir_t> cfdirs;

    mymesh->cell_get_faces_and_dirs(c, &cfaces, &cfdirs);

    std::cerr << "Cell " << c << ":" << std::endl;
    std::cerr << "  Faces:";
    for (auto const & f : cfaces)
      std::cerr << " " << c;
    std::cerr << std::endl;

    std::cerr << "  Dirs:";
    for (auto const & fdir : cfdirs)
      std::cerr << " " << fdir;
    std::cerr << std::endl;

    std::cerr << std::endl;
  }

  // Clean up and exit

  MPI_Finalize();

}

