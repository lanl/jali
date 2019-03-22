/*
 Copyright (c) 2019, Triad National Security, LLC
 All rights reserved.

 Copyright 2019. Triad National Security, LLC. This software was
 produced under U.S. Government contract 89233218CNA000001 for Los
 Alamos National Laboratory (LANL), which is operated by Triad
 National Security, LLC for the U.S. Department of Energy. 
 All rights in the program are reserved by Triad National Security,
 LLC, and the U.S. Department of Energy/National Nuclear Security
 Administration. The Government is granted for itself and others acting
 on its behalf a nonexclusive, paid-up, irrevocable worldwide license
 in this material to reproduce, prepare derivative works, distribute
 copies to the public, perform publicly and display publicly, and to
 permit others to do so

 
 This is open source software distributed under the 3-clause BSD license.
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:
 
 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
 3. Neither the name of Triad National Security, LLC, Los Alamos
    National Laboratory, LANL, the U.S. Government, nor the names of its
    contributors may be used to endorse or promote products derived from this
    software without specific prior written permission.

 
 THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND
 CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <iostream>
#include <vector>
#include <array>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

// Code that attempts to _mock_ what may be done in a staggered grid
// numerical scheme with some fields on cells and others on nodes
// The "numerics" have no relation to any real numerical algorithm but
// but illustrate how a real algorithm may be written

using namespace Jali;
using namespace JaliGeometry;

// Define two variable names that will be used to define state data in
// the "initialize_data" routine and retrieve it in "main"

std::string density_name("rhoMetal");
std::string velocity_name("nodevel");

// Forward declaration of routine to initialize data - meant to
// illustrate definition/initialization of state data in one place and
// retrieval of the data by name in another

void initialize_data(const std::shared_ptr<Mesh> mesh,
                     std::shared_ptr<State> state);


// Start main routine


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
  
    // Create a 3D mesh from (0.0,0.0,0.0) to (1.0,1.0,1.0) with 10, 5
    // and 5 elements in the X, Y and Z directions. Specify that we
    // want all kinds of entities (faces, edges, wedges and corners)
    // to be present. Finally, request that the mesh be divided into 10
    // tiles.

    int num_tiles_requested = 10;
    mesh_factory.included_entities(Entity_kind::ALL_KIND);
    mesh_factory.num_tiles(num_tiles_requested);
    mesh_factory.num_ghost_layers_tile(1);
    mymesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 5, 5);
  }


  // Spatial dimension of points in the mesh - since this is a 3D mesh, this
  // should be 3

  int spdim = mymesh->space_dimension();


  // Topological dimension of cells in the mesh - since this is a 3D
  // mesh, this should be 3

  int celldim = mymesh->manifold_dimension();


  // Create a state manager for handling data associated with mesh entities

  std::shared_ptr<State> mystate = Jali::State::create(mymesh);


  // Initialize the state manager with some data - see routine at end
  // of file.  It sets up a scalar array on cells using the string in
  // density_name and a vector array on nodes using the string in
  // velocity_name

  initialize_data(mymesh, mystate);


  // Retrieve the density data on cells. This is indirect. The find
  // routine returns an iterator value and one has to dereference it
  // to get at the StateVector. This means one needs to know the template
  // argument for StateVector (in this case it is "double")

  UniStateVector<double> rhovec;
  bool found = mystate->get(density_name, mymesh, Entity_kind::CELL,
                            Entity_type::ALL, &rhovec);
  if (!found) {
    std::cerr << "Could not find state vector on cells with name " <<
        density_name << std::endl;
    exit(-1);
  }

  // Compute the average density on cells using all node connected neighbors

  int nc = mymesh->num_entities(Entity_kind::CELL, Entity_type::ALL);
  double *ave_density = new double[nc];
  for (int i = 0; i < nc; ++i) ave_density[i] = 0.0;

  for (auto const& t : mymesh->tiles()) {
    for (auto const& c : t->cells()) {

      // Get all (owned or ghost) node connected/adjacent neighbors of a cell

      Entity_ID_List nbrs;
      mymesh->cell_get_node_adj_cells(c, Entity_type::ALL, &nbrs);

      for (auto const & nc : nbrs)
        ave_density[c] += rhovec[nc];
      ave_density[c] /= nbrs.size();
    }
  }

  // Add the average density data as a new state vector

  UniStateVector<double>& rhobarvec = mystate->add("rhobar", mymesh,
                                                Entity_kind::CELL,
                                                Entity_type::ALL,
                                                ave_density);

  delete [] ave_density;



  // Retrieve the vector of velocities

  UniStateVector<std::array<double, 3>, Mesh> vels;
  found = mystate->get(velocity_name, mymesh, Entity_kind::NODE,
                       Entity_type::ALL, &vels);
  if (!found) {
    std::cerr << "Could not find state vector on nodes with name " <<
        velocity_name << std::endl;
    exit(-1);
  }

  // Update them to be the average of the centroids of the connected
  // cells weighted by the average cell density

  for (auto const& t : mymesh->tiles()) {
    for (auto const n : t->nodes()) {

      // Get the cells using (connected to) this node

      Entity_ID_List nodecells;
      mymesh->node_get_cells(n, Entity_type::ALL, &nodecells);

      std::array<double, 3> tmpvels;
      for (int i = 0; i < 3; ++i) tmpvels[i] = 0.0;

      for (auto c : nodecells) {

        // Get cell centroid - this is computed once and cached unless the
        // mesh changes or the routine is explicitly asked to recompute it
        // to help understand the effect of a temporary change

        JaliGeometry::Point ccen = mymesh->cell_centroid(c);
        for (int i = 0; i < 3; ++i) tmpvels[i] += rhobarvec[c]*ccen[i];

      }

      vels[n] = tmpvels;
    }
  }



  // Print out the average densities at cell centers

  std::cerr << "Average densities at cell centers:" << std::endl;

  for (auto c : mymesh->cells<Entity_type::PARALLEL_OWNED>()) {
    Point ccen = mymesh->cell_centroid(c);
    
    std::cerr << "Cell " << c << "    Centroid (" <<
        ccen[0] << ", " << ccen[1] << ", " << ccen[2] <<
        ")    Ave density " << rhobarvec[c] << std::endl;
  }
  std::cerr << std::endl << std::endl;


  // Print out computed velocities at nodes

  std::cerr << "Computed velocities at nodes:" << std::endl;

  for (auto n : mymesh->nodes<Entity_type::PARALLEL_OWNED>()) {
    Point npnt;
    mymesh->node_get_coordinates(n, &npnt);

    std::cerr << "Node " << n << "    Coord (" <<
        npnt[0] << ", " << npnt[1] << ", " << npnt[2] <<
        ")    Velocity (" << vels[n][0] << ", " << vels[n][1] <<
        ", " << vels[n][2] << ")" << std::endl << std::endl;
  }


  // Wrap up

  MPI_Finalize();
}




// Routine for initialization of state data


void initialize_data(const std::shared_ptr<Mesh> mesh,
                     std::shared_ptr<State> state) {

  // number of cells in the mesh - ALL means OWNED+GHOST
  int nc = mesh->num_cells<Entity_type::ALL>();

  // Add a state vector, called "rho99", of densities on cells.
  // This says to create a vector named "rhoMetal" on each cell
  // initialized to 0.0

  UniStateVector<double, Mesh>& density =
      state->add<double, Mesh, UniStateVector>(density_name, mesh,
                                               Entity_kind::CELL,
                                               Entity_type::ALL, 0.0);


  // Create a density vector that will be used to initialize a state
  // variable called 'rho99' on cells

  for (auto c : mesh->cells<Entity_type::PARALLEL_OWNED>()) {
    Point ccen = mesh->cell_centroid(c);
    density[c] = ccen[0]+ccen[1]+ccen[2];
  }


  // Add a state vector for velocities

  int dim = mesh->space_dimension();
  int nn = mesh->num_nodes<Entity_type::ALL>();

  std::array<double, 3> initarray;
  for (int i = 0; i < dim; ++i) initarray[i] = 0.0;

  // Note that we did not send in the second template parameter Mesh -
  // it is the default

  UniStateVector<std::array<double, 3>>& velocity =
      state->add<std::array<double, 3>, Mesh, UniStateVector>(velocity_name,
                                                              mesh,
                                                              Entity_kind::NODE,
                                                              Entity_type::ALL,
                                                              initarray);
}

