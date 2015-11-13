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


  // Iterate through the cells of the mesh and get their node vertices and
  // their volumes+centroids
  
  Mesh::cell_iterator itc = mymesh->begin_cells();
  while (itc != mymesh->end_cells()) {
    Entity_ID c = *itc;

    std::cerr << "Cell " << c << ":" << std::endl;
    
    // Get coordinates of nodes of cell

    std::vector<JaliGeometry::Point> ccoords;
    mymesh->cell_get_coordinates(c,&ccoords);

    std::cerr << "  Node Coordinates:" << std::endl;
    int nnodes = ccoords.size();
    for (int i = 0; i < nnodes; i++)
      std::cerr << "     " << ccoords[i] << std::endl;

    // Get volume of cell
    
    double cellvol = mymesh->cell_volume(c);

    std::cerr << "  Cell Volume/Area: " << cellvol << std::endl;

    // Get centroid of cell

    JaliGeometry::Point ccen;
    ccen = mymesh->cell_centroid(c);

    std::cerr << "  Cell Centroid: " << ccen << std::endl; 

    std::cerr << std::endl;

    ++itc;
  }

  // Clean up and exit

  delete mymesh;
  
  MPI_Finalize();

}

