#include <iostream>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

using namespace Jali;

// Fire up Jali, create a mesh and ask the mesh a very basic question

//**************************************************************************
//
// NOTE: One must be very careful with assigning and copying state
// vectors Some operations do a shallow copy - they copy the name,
// entity_kind and the pointer to the data which means both vectors
// point to the same data.
//
// Other operations do a copy construction which is a deep copy - in
// such an operation, not only is the meta data copied but new memory
// is allocated in the new vector and the data from the old vector
// copied to the new vector
//
// So, as an example:
//
// Jali::StateVector<int> v1("vec1",Jali::CELL,data);
//
// Jali::StateVector<int> & v1_ref = v1; - both point to SAME data
//
// Jali::StateVector<int> v2;            - default construction
// v2 = v1;                              - v1 and v2 point to SAME data
//
// Jali::StateVector<int> v3 = v1;       - copy construction
//                                       - v1 and v3 point to DIFFERENT data
//
// Jali::StateVector<int> v4(v1);        - copy construction
//                                       - v1 and v4 point to DIFFERENT data
//
// IN PARTICULAR NOTE THE DIFFERENCE BETWEEN v2 AND v3
//




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
  
    // Create a 2D mesh from (0.0,0.0) to (1.0,1.0)
    // with 3, 3 and 3 elements in the X, Y directions. Specify
    // that we did not instantiate a geometric model (NULL). Also,
    // request faces, edges, wedges and corners (true, true, true,
    // true)

    mymesh = mesh_factory(0.0,0.0,1.0,1.0,2,2,NULL,true,true,true,true);
  }


  
  // Print out the number of cells in the mesh

  std::cerr << "Number of mesh cells: " << mymesh->num_entities(CELL,ALL) 
            << std::endl;

  // Print out the number of nodes in the mesh

  std::cerr << "Number of mesh nodes: " << mymesh->num_entities(NODE,ALL)
            << std::endl;


  // Create a Jali State Manager

  Jali::State mystate(mymesh);


  // Create some CELL-based data

  double data[4] = {0.0,1.0,2.0,3.0};

  // Add it to state and get back a reference
  //
  // DON'T DO
  //   Jali::StateVector<double> myvec = mystate.add()
  //
  // This will do a copy construction and myvec data space will be different
  // from the state vector data space!!

  Jali::StateVector<double> & myvec = mystate.add("myzonevar",Jali::CELL,data);


  // Try to retrieve it through a get function

  Jali::StateVector<double> myvec_copy;
  bool found = mystate.get("myzonevar",Jali::CELL,&myvec_copy);

  int ndata = myvec_copy.size();
  if (myvec.size() != myvec_copy.size()) {
    std::cerr << "Stored and retrieved vectors have different sizes?" << std::endl;
    exit(-1);
  }
   
  for (int i = 0; i < ndata; i++) {
    if (myvec[i] != myvec_copy[i]) {
      std::cerr << "Stored and retrieved vectors differ at element " << i << std::endl;
    }
  }


  // Assign to another state vector - should be a shallow copy

  Jali::StateVector<double> myvec_copy2;
  myvec_copy2 = myvec;



  // Modify a value in the state vector

  myvec_copy[3] = 25;


  // It should get reflected in myvec since it points to the same data

  if (myvec[3] != myvec_copy[3]) {
    std::cerr << "Stored and retrieved vectors don't point to the same data?" <<
        std::endl;
    std::cerr << "A change in one was not reflected in the other" << std::endl;
    exit(-1);
  }


  if (myvec[3] != myvec_copy2[3]) {
    std::cerr << "Original and assigned vectors don't point to the same data?" <<
        std::endl;
    std::cerr << "A change in one was not reflected in the other" << std::endl;
    exit(-1);
  }


  // Make a new state vector using a copy constructor - DEEP COPY OF DATA

  Jali::StateVector<double> newvec = myvec; // Also note: newvec is not in mystate!



  // changing an element in myvec will not change newvec since new vec
  // has its own data space

  double old_val = myvec[2];
  myvec[2] = 33;

  if (newvec[2] == myvec[2] || newvec[2] != old_val) {
    std::cerr << "Original and copy constructed vectors share data space?" <<
        std::endl;
    exit(-1);
  }


  // Print out myvec
  
  std::cerr << "StateVector " << myvec.name() << ":" << std::endl;
  std::cerr << myvec << std::endl;


  // Define a more complicated StateVector and print it
  // Note: This vector is not added to the state

  // 

  std::array<double,3> arrdata[9] = {{0.0,1.0,2.0},
                                     {3.0,4.0,5.0},
                                     {6.0,7.0,8.0},
                                     {9.0,10.0,11.0},
                                     {12.0,13.0,14.0},
                                     {15.0,16.0,17.0},
                                     {18.0,19.0,20.0},
                                     {21.0,22.0,23.0},
                                     {24.0,25.0,26.0}};


  Jali::StateVector<std::array<double,3>> pntvec("vec3",Jali::NODE,mymesh,&(arrdata[0]));

  std::cerr << pntvec << std::endl;


  // Clean up and exit

  delete mymesh;
  
  MPI_Finalize();

}

