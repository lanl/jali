/*
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.

Copyright 2017. Los Alamos National Security, LLC. This software was
produced under U.S. Government contract DE-AC52-06NA25396 for Los
Alamos National Laboratory (LANL), which is operated by Los Alamos
National Security, LLC for the U.S. Department of Energy. The
U.S. Government has rights to use, reproduce, and distribute this
software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so
as not to confuse it with the version available from LANL.
 
Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

1.  Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3.  Neither the name of Los Alamos National Security, LLC, Los Alamos
National Laboratory, LANL, the U.S. Government, nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <iostream>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

using namespace Jali;

// Fire up Jali, create a mesh and a state manager
// Add and retrieve state data on meshes and meshtiles

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
// Jali::StateVector<int> v1("vec1", Jali::CELL, data);
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
  int mesh_dimension = 2;
  if (framework_available(MSTK) &&
      framework_generates(MSTK, parallel_mesh, mesh_dimension)) {

    mesh_factory.framework(MSTK);

    // Create a 2D mesh from (0.0, 0.0) to (1.0, 1.0)
    // with 4 and 4 elements in the X, Y directions. Also,
    // request faces, edges, wedges and corners. Finally, request
    // four tiles in the mesh (with 4 elements each)

    int num_tiles = 4;

    mesh_factory.included_entities(Entity_kind::ALL_KIND);

    mesh_factory.num_tiles(num_tiles);

    mymesh = mesh_factory(0.0, 0.0, 1.0, 1.0, 4, 4);
  }


  // Print out the number of cells in the mesh

  std::cerr << "Number of mesh cells: " <<
    mymesh->num_cells<Entity_type::ALL>() << std::endl;

  // Print out the number of nodes in the mesh

  std::cerr << "Number of mesh nodes: " <<
    mymesh->num_nodes<Entity_type::ALL>() << std::endl;


  // Create a Jali State Manager

  Jali::State mystate(mymesh);


  // Create some CELL-based data

  double data[16] = {0.0, 1.0, 2.0, 3.0, 5.5, -1.0, 2.0, 3.0,
                     2.2, 3.0, 9.0, 0.0, -7, -10, 4.2, 1.1};

  // Add it to state and get back a reference
  //
  // DO
  //   StateVector<double>& myvec = mystate.add(...)
  //
  // DON'T DO
  //   StateVector<double> myvec = mystate.add(...)
  // This will do a copy construction and myvec data space will be different
  // from the state vector data space!!

  // NOTE that there is an implicit second template parameter "Mesh" here
  // that is figured out from the arguments to the add function. 

  StateVector<double> & myvec = mystate.add("myzonevar", mymesh,
                                             Entity_kind::CELL,
                                             Entity_type::ALL, data);


  // Try to retrieve it through a get function

  StateVector<double, Jali::Mesh> myvec_copy1;
  bool found = mystate.get("myzonevar", mymesh, Entity_kind::CELL,
                           Entity_type::ALL, &myvec_copy1);

  int ndata = myvec_copy1.size();
  if (myvec.size() != myvec_copy1.size()) {
    std::cerr << "Stored and retrieved vectors have different sizes?" <<
      std::endl;
    exit(-1);
  }
   
  for (int i = 0; i < ndata; i++) {
    if (myvec[i] != myvec_copy1[i]) {
      std::cerr << "Stored and retrieved vectors differ at element " << i <<
        std::endl;
    }
  }


  // Assign to another state vector AFTER creating the vector as a
  // default vector. Should be a shallow copy

  StateVector<double> myvec_copy2;
  myvec_copy2 = myvec;



  // Modify a value in the state vector

  myvec_copy1[3] = 25;


  // It should get reflected in myvec since it points to the same data

  if (myvec[3] != myvec_copy1[3]) {
    std::cerr << "Stored and retrieved vectors don't point to the same data?" <<
        std::endl;
    std::cerr << "A change in one was not reflected in the other" << std::endl;
    exit(-1);
  }


  if (myvec[3] != myvec_copy2[3]) {
    std::cerr <<
        "Original and assigned vectors don't point to the same data?" <<
        std::endl;
    std::cerr << "A change in one was not reflected in the other" << std::endl;
    exit(-1);
  }


  // Add a new uninitialized state vector - Note that we have to
  // explicitly tell the state manager the data type (double) since
  // there is no input data for it to infer it from

  StateVector<double>& newvec = mystate.add<double>("vector2", mymesh, 
                                                   Entity_kind::CELL,
                                                   Entity_type::ALL);

  // Modify newvec 

  for (auto const& c : mymesh->cells())
    newvec[c] = 2.0*c;
    

  // Retrieve the vector from the state manager separately and make sure
  // that the changes in newvec were reflected in the state manager

  StateVector<double, Mesh> newvec_copy;
  found = mystate.get("vector2", mymesh, Entity_kind::CELL,
                      Entity_type::ALL, &newvec_copy);


  ndata = newvec.size();
  for (int i = 0; i < ndata; ++i) {
    if (newvec[i] != newvec_copy[i]) {
      std::cerr <<
          "Changes to state vector reference not reflected in State Manager?" <<
          std::endl;
      break;
    }
  }

  // Print out myvec

  std::cerr << "StateVector " << myvec.name() << ":" << std::endl;
  std::cerr << myvec << std::endl;


  // Define a more complicated StateVector and print it.

  // Note: This is a standalone vector that is not added to the
  // state. Such vectors could be used as temporaries in a calculation

  std::array<double, 3> arrdata[25] = {{0.0, 1.0, 2.0},    {3.0, 4.0, 5.0},
                                       {6.0, 7.0, 8.0},    {9.0, 10.0, 11.0},
                                       {12.0, 13.0, 14.0}, {15.0, 16.0, 17.0},
                                       {18.0, 19.0, 20.0}, {21.0, 22.0, 23.0},
                                       {24.0, 25.0, 26.0}, {3.0, 4.0, 2.2},
                                       {18.0, 19.0, 20.0}, {21.0, 22.0, 23.0},
                                       {-1.0, 2.2, 9.8},   {99, 3.7, -1.0},
                                       {0.0, 0.0, 8.9},    {-9.2, -1.0, 4.2},
                                       {12.0, 13.0, 14.0}, {15.0, 16.0, 17.0},
                                       {21.5, -2.4, -1},   {30.0, 21.1, 2.2},
                                       {24.0, 25.0, 26.0}, {3.0, 4.0, 2.2},
                                       {-1.0, 2.2, 9.8},   {99, 3.7, -1.0},
                                       {21.5, -2.4, -1}};


  StateVector<std::array<double, 3>> vec2d("vec3", mymesh, Entity_kind::NODE,
                                           Entity_type::PARALLEL_OWNED,
                                           &(arrdata[0]));

  std::cerr << vec2d << std::endl;



  // Define data on meshtiles and later print it

  int i = 0;
  for (auto const& meshtile : mymesh->tiles()) {
    auto tilevec_put = mystate.add("vector_on_tile", meshtile,
                                   Entity_kind::CELL, Entity_type::PARALLEL_OWNED,
                                   &(data[4*i]));
    
    auto tilevec2d_put = mystate.add("2Darray_on_tile", meshtile,
                                     Entity_kind::CORNER, Entity_type::PARALLEL_OWNED,
                                     &(arrdata[0]));
    i++;
  }


  // Retrieve it and print it

  i = 0;
  for (auto const& meshtile : mymesh->tiles()) {
    StateVector<double, MeshTile> tilevec_get;
    found = mystate.get("vector_on_tile", meshtile, Entity_kind::CELL,
                        Entity_type::PARALLEL_OWNED, &tilevec_get);

    std::cerr << "1D vector on tile cells: " << std::endl;
    std::cerr << tilevec_get << std::endl;
    std::cerr << std::endl;

    StateVector<std::array<double, 3>, MeshTile> tilevec2d_get;
    found = mystate.get("2Darray_on_tile", meshtile, Entity_kind::CORNER,
                        Entity_type::PARALLEL_OWNED, &tilevec2d_get);


    std::cerr << "2D vector on tile cells: " << std::endl;
    std::cerr << tilevec2d_get << std::endl;
    std::cerr << std::endl;
  }

  // Clean up and exit

  MPI_Finalize();

}

