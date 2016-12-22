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

using namespace Jali;
using namespace JaliGeometry;

// Fire up Jali, create a mesh and ask the mesh about cell geometry

int main(int argc, char *argv[]) {

  // Jali depends on MPI

  MPI_Init(&argc, &argv);

  // Create a mesh factory object - this object has methods for
  // specifying the preference of mesh frameworks and unified
  // interfaces for instantiating a mesh object of a particular
  // framework type

  MPI_Comm comm = MPI_COMM_WORLD;
  MeshFactory mesh_factory(comm);

  // For this example we will use the in-built Simple mesh framework


  std::shared_ptr<Mesh> mymesh;  // Pointer to a mesh object
  bool parallel_mesh = false;
  int mesh_dimension = 3;
  if (framework_available(Simple) &&
      framework_generates(Simple, parallel_mesh, mesh_dimension)) {

    mesh_factory.framework(Simple);
  
    // Create a 3D mesh from (0.0,0.0,0.0) to (1.0,1.0,1.0)
    // with 3, 3 and 3 elements in the X, Y and Z directions.
    // Request faces, edges, wedges and corners

    mesh_factory.included_entities(Entity_kind::ALL_KIND);

    mymesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
  }


  // Iterate through the cells of the mesh and get their node vertices and
  // their volumes+centroids
  
  for (auto c : mymesh->cells<Entity_type::PARALLEL_OWNED>()) {
    std::cerr << "Cell " << c << ":" << std::endl;
    
    // Get coordinates of nodes of cell

    std::vector<Point> ccoords;
    mymesh->cell_get_coordinates(c, &ccoords);

    std::cerr << "  Node Coordinates:" << std::endl;
    for (auto const& xyz : ccoords)
      std::cerr << "     " << xyz << std::endl;

    // Get volume of cell
    
    double cellvol = mymesh->cell_volume(c);

    std::cerr << "  Cell Volume/Area: " << cellvol << std::endl;

    // Get centroid of cell

    JaliGeometry::Point ccen;
    ccen = mymesh->cell_centroid(c);

    std::cerr << "  Cell Centroid: " << ccen << std::endl;

    std::cerr << std::endl;
  }

  // Clean up and exit
  
  MPI_Finalize();

}

