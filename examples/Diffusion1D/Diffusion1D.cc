#include <iostream>
#include <memory>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "Geometry.hh"

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
  pref.push_back(Simple);

  std::shared_ptr<Mesh> mymesh;  // Pointer to a mesh object
  if (framework_available(Simple)) {  // check if framework is available
    mesh_factory.preference(pref);

    // Create a 1D mesh from 0.0 to 1.0
    // with 10 elements in the X direction. Specify
    // that we did not instantiate a geometric model (NULL). Also,
    // request faces, edges, wedges and corners (true, true, true,
    // true)
    // mymesh = mesh_factory(0.0, 1.0, 11, NULL, true, true, true, true);
    std::vector<double> node_coords = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0};
    mymesh = mesh_factory(node_coords, NULL, true, true, true, true,
                          JaliGeometry::SPHERICAL);
    //    mymesh->set_geom_type(SPHERICAL);
  }


  // Print out the topological dimension of cells in the mesh

  std::cerr << "Cells are of dimension: " << mymesh->cell_dimension()
            << std::endl;

  // Print out the number of cells in the mesh

  Entity_ID numcells = mymesh->num_entities(CELL,ALL);
  std::cerr << "Number of mesh cells: " << numcells << std::endl;

  // Print out the number of nodes in the mesh

  std::cerr << "Number of mesh nodes: " << mymesh->num_entities(NODE,ALL)
            << std::endl;
  std::cerr << "Number of mesh edges: " << mymesh->num_entities(EDGE,ALL)
            << std::endl;
  Entity_ID numwedges = mymesh->num_entities(WEDGE,ALL);
  std::cerr << "Number of mesh wedges: " << numwedges << std::endl;
  std::cerr << "Number of mesh corners: " << mymesh->num_entities(CORNER,ALL)
            << std::endl;
  std::cerr << "Number of mesh faces: " << mymesh->num_entities(FACE,ALL)
            << std::endl;

  for (Entity_ID i = 0; i < numcells; i++) {
    std::cout << "Cell " << i << std::endl;
    std::cout << " volume: " << mymesh->cell_volume(i)
              << std::endl;
    Entity_ID_List faceids;
    std::vector<int> facedirs;
    mymesh->cell_get_faces_and_dirs(i, &faceids, &facedirs);
    for (Entity_ID f : faceids) {
      std::cout << " face " << f << " area: "
                << mymesh->face_normal(f, false, i)
                << std::endl;
    }
  }

  for (Entity_ID i = 0; i < numwedges; i++) {
    std::cout << "Wedge " << i << std::endl;
    std::cout << " volume: " << mymesh->wedge_volume(i)
              << std::endl;
  }

  // Clean up and exit

  MPI_Finalize();

}

