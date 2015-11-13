#include <iostream>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"

using namespace Jali;

// Fire up Jali, create a mesh and print out number of CELLS and NODES

int main(int argc, char *argv[]) {

  // Jali depends on MPI 

  MPI_Init(&argc,&argv);

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

  Mesh *mymesh;  // Pointer to a mesh object
  if (framework_available(MSTK)) {  // check if framework is available
    mesh_factory.preference(pref);  
  
    // Create a 3D mesh from (0.0,0.0,0.0) to (1.0,1.0,1.0)
    // with 3, 3 and 3 elements in the X, Y and Z directions. Specify
    // that we did not instantiate a geometric model (NULL). Also,
    // request faces, edges, wedges and corners (true, true, true,
    // true)

    mymesh = mesh_factory(0.0,0.0,0.0,1.0,1.0,1.0,3,3,3,NULL,
                          true,true,true,true);
  }


  // Print out the topological dimension of cells in the mesh

  std::cerr << "Cells are of dimension: " << mymesh->cell_dimension() << std::endl;
  
  // Print out the number of cells in the mesh

  std::cerr << "Number of mesh cells: " << mymesh->num_entities(CELL,ALL) << std::endl;

  // Print out the number of nodes in the mesh

  std::cerr << "Number of mesh nodes: " << mymesh->num_entities(NODE,ALL) << std::endl;



  // Clean up and exit

  delete mymesh;
  
  MPI_Finalize();

}

