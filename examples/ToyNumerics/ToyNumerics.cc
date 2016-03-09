#include <iostream>
#include <vector>
#include <array>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

// Code that attempts to _mock_ what may be done in a staggered grid 
// numerical scheme with some fields on cells and others on nodes
// The "numerics" have no relation to any real numerical algorithm but
// but illustrate how a real algorithm may be written

using namespace Jali;

// Define two variable names that will be used to define state data in
// the "initialize_data" routine and retrieve it in "main"

std::string density_name("rhoMetal");
std::string velocity_name("nodevel");

// Forward declaration of routine to initialize data - meant to
// illustrate definition/initialization of state data in one place and
// retrieval of the data by name in another

void initialize_data(Mesh & mesh, State & state);



// Start main routine


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

  std::shared_ptr<Mesh> mymesh;  // Pointer to a mesh object
  if (framework_available(MSTK)) {  // check if framework is available
    mesh_factory.preference(pref);  
  
    // Create a 3D mesh from (0.0,0.0,0.0) to (1.0,1.0,1.0) with 10, 5
    // and 5 elements in the X, Y and Z directions. Specify that we
    // did not instantiate a geometric model (NULL). Also, request
    // faces, edges, wedges and corners (true, true, true, true)

    mymesh = mesh_factory(0.0,0.0,0.0,1.0,1.0,1.0,10,5,5,NULL,
                          true,true,true,true);
  }


  // Spatial dimension of points in the mesh - since this is a 3D mesh, this
  // should be 3

  int spdim = mymesh->space_dimension();


  // Topological dimension of cells in the mesh - since this is a 3D
  // mesh, this should be 3

  int celldim = mymesh->cell_dimension();


  // Create a state manager for handling data associated with mesh entities

  State mystate(mymesh.get());


  // Initialize the state manager with some data - see routine at end
  // of file.  It sets up a scalar array on cells using the string in
  // density_name and a vector array on nodes using the string in
  // velocity_name

  initialize_data(*mymesh,mystate);


  // Retrieve the density data on cells. This is indirect. The find
  // routine returns an iterator value and one has to dereference it
  // to get at the StateVector. This means one needs to know the template
  // argument for StateVector (in this case it is "double")

  StateVector<double> rhovec;
  bool found = mystate.get(density_name,CELL,&rhovec);
  if (!found) {
    std::cerr << "Could not find state vector on cells with name " << density_name << std::endl;
    exit(-1);
  }

  // Compute the average density on cells using all node connected neighbors

  int nc = mymesh->num_entities(CELL,ALL);
  double *ave_density = new double[nc];
  for (int i = 0; i < nc; ++i) ave_density[i] = 0.0;

  auto itc = mymesh->begin_cells();
  while (itc != mymesh->end_cells()) {
    auto c = *itc;

    // Get all (owned or ghost) node connected/adjacent neighbors of a cell 

    Entity_ID_List nbrs;
    mymesh->cell_get_node_adj_cells(c,ALL,&nbrs);

    auto itc2 = nbrs.begin();
    while (itc2 != nbrs.end()) {
      auto c = *itc2;
      ave_density[c] += rhovec[c];
      ++itc2;
    }
    ave_density[c] /= nbrs.size();

    ++itc;
  }

  // Add the average density data as a new state vector

  StateVector<double> &rhobarvec = mystate.add("rhobar",CELL,ave_density);
  
  delete [] ave_density;



  // Retrieve the vector of velocities

  StateVector<std::array<double,3>> vels;
  found = mystate.get(velocity_name,NODE,&vels);
  if (!found) {
    std::cerr << "Could not find state vector on nodes with name " << 
        velocity_name << std::endl;
    exit(-1);
  }

  // Update them to be the average of the centroids of the connected
  // cells weighted by the average cell density

  auto itn = mymesh->begin_nodes();
  while (itn != mymesh->end_nodes()) {
    auto n = *itn;

    // Get the cells using (connected to) this node

    Entity_ID_List nodecells;
    mymesh->node_get_cells(n,ALL,&nodecells);
  
    std::array<double,3> tmpvels;
    for (int i = 0; i < 3; ++i) tmpvels[i] = 0.0;

    auto itc2 = nodecells.begin();
    while (itc2 != nodecells.end()) {
      auto c = *itc2;

      // Get cell centroid - this is computed once and cached unless the 
      // mesh changes or the routine is explicitly asked to recompute it
      // to help understand the effect of a temporary change

      JaliGeometry::Point ccen = mymesh->cell_centroid(c);
      for (int i = 0; i < 3; ++i) tmpvels[i] += rhobarvec[c]*ccen[i];
      ++itc2;
    }

    vels[n] = tmpvels;
    ++itn;
  }

  

  // Print out the average densities at cell centers
  
  std::cerr << "Average densities at cell centers:" << std::endl;

  itc = mymesh->begin_cells();
  while (itc != mymesh->end_cells()) {
    auto c = *itc;

    JaliGeometry::Point ccen = mymesh->cell_centroid(c);
    
    std::cerr << "Cell " << c << "    Centroid (" <<
        ccen[0] << "," << ccen[1] << "," << ccen[2] <<
        ")    Ave density " << rhobarvec[c] << std::endl;
    ++itc;
  }
  std::cerr << std::endl << std::endl;

  
  // Print out computed velocities at nodes

  std::cerr << "Computed velocities at nodes:" << std::endl;

  itn = mymesh->begin_nodes();
  while (itn != mymesh->end_nodes()) {
    auto n = *itn;

    JaliGeometry::Point npnt;
    mymesh->node_get_coordinates(n,&npnt);

    std::cerr << "Node " << n << "    Coord (" <<
        npnt[0] << "," << npnt[1] << "," << npnt[2] <<
        ")    Velocity (" << vels[n][0] << "," << vels[n][1] <<
        "," << vels[n][2] << ")" << std::endl << std::endl;
    ++itn;
  }


  // Wrap up

  MPI_Finalize();
}




// Routine for initialization of state data


void initialize_data(Mesh & mesh, State & state) {

  // number of cells in the mesh - ALL means OWNED+GHOST
  int nc = mesh.num_entities(CELL,ALL);

  // Create a density vector that will be used to initialize a state
  // variable called 'rho99' on cells

  std::vector<double> density(nc);
  auto itc = mesh.begin_cells();
  while (itc != mesh.end_cells()) {
    auto c = *itc;
    JaliGeometry::Point ccen = mesh.cell_centroid(c);  
    density[c] = ccen[0]+ccen[1]+ccen[2];
    ++itc;
  }

  // Create a state vector of densities on cells using the
  // initialization data.  This says to create a vector named
  // "rhoMetal" on each cell and populate it with the given
  // data. Since density is a std::vector<double> we have to send in
  // the address of the first element.  as &(density[0]).

  state.add(density_name,CELL,&(density[0]));


  // Create a velocity vector

  int dim = mesh.space_dimension();
  int nn = mesh.num_entities(NODE,ALL);

  // Initialize to zero

  std::array<double,3> initarray;
  for (int i = 0; i < dim; ++i) initarray[i] = 0.0;

  std::vector<std::array<double,3>> vels(nn,initarray);

  // Add it to the state manager

  state.add(velocity_name,NODE,&(vels[0]));
}

