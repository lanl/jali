#include <iostream>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"

using namespace Jali;

// Fire up Jali, create a mesh and ask the mesh about cell topology

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


  
  // Iterate through the cells of the mesh and get the faces and of the cell
  
  Mesh::cell_iterator itc = mymesh->begin_cells();
  while (itc != mymesh->end_cells()) {
    Entity_ID c = *itc;
    
    Entity_ID_List cfaces;   // Entity_ID_List is just a std::vector<Entity_ID>
    std::vector<int> cfdirs;

    mymesh->cell_get_faces_and_dirs(c,&cfaces,&cfdirs);

    std::cerr << "Cell " << c << ":" << std::endl;
    std::cerr << "  Faces:";
    int nfaces = cfaces.size();
    for (int i = 0; i < nfaces; i++) 
      std::cerr << " " << cfaces[i];
    std::cerr << std::endl;

    std::cerr << "  Dirs:";
    for (int i = 0; i < nfaces; i++)
      std::cerr << " " << cfdirs[i];
    std::cerr << std::endl;

    std::cerr << std::endl;


    ++itc;
  }

  // Clean up and exit

  delete mymesh;
  
  MPI_Finalize();

}

