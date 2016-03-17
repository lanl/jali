/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include <mpi.h>
#include <stdlib.h>

#include <iostream>

#include "JaliState.h"
#include "JaliStateVector.h"
#include "Mesh.hh"
#include "MeshFactory.hh"

#include "UnitTest++.h"

// Vector type for 2d doubles
struct Vec2d {
  double x;
  double y;

  void set(double xvalue, double yvalue) {
    x = xvalue;  y = yvalue;
  }

  friend std::ostream &operator<<(std::ostream &output, const Vec2d &v) {
    output << "(" << v.x << ", " << v.y << ")";
    return output;
  }
};


TEST(Jali_State_Define_Mesh) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh1 = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh1);

  // Define two state vectors

  std::vector<double> data1 = {1.0, 3.0, 2.5, 4.5};
  Jali::StateVector<double, Jali::Mesh> myvec1("cellvars", mesh1,
                                               Jali::Entity_kind::CELL,
                                               Jali::Parallel_type::ALL,
                                               &(data1[0]));

  std::vector<double> data2 = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  Jali::StateVector<double, Jali::Mesh> myvec2("nodevars", mesh1,
                                               Jali::Entity_kind::NODE,
                                               Jali::Parallel_type::ALL,
                                               &(data2[0]));

  // Define another mesh and another statevector on that mesh

  std::shared_ptr<Jali::Mesh> mesh2 = mf(0.0, 0.0, 1.0, 1.0, 3, 3);

  std::vector<double> data3 = {1.0, 3.0, 2.5, 4.5, 1.0, 2.0};
  Jali::StateVector<double, Jali::Mesh> myvec3("cellvars2", mesh2,
                                               Jali::Entity_kind::CELL,
                                               Jali::Parallel_type::ALL,
                                               &(data3[0]));

  // Create a state object and add the first two vectors to it

  Jali::State mystate(mesh1);

  int add_status;
  Jali::StateVector<double, Jali::Mesh> &addvec1 = mystate.add(myvec1);
  CHECK_EQUAL(addvec1.size(), myvec1.size());
  for (int i = 0; i < addvec1.size(); ++i)
    CHECK_EQUAL(addvec1[i], myvec1[i]);

  Jali::StateVector<double, Jali::Mesh> &addvec2 =
      mystate.add("nodevars", mesh1, Jali::Entity_kind::NODE,
                  Jali::Parallel_type::ALL, &(data2[0]));
  CHECK_EQUAL(addvec2.size(), myvec2.size());
  for (int i = 0; i < addvec2.size(); ++i)
    CHECK_EQUAL(addvec2[i], myvec2[i]);
  

  // Try to add the third vector (defined on a different mesh) to it - it
  // should copy the data but be assigned to mesh1 instead of mesh2

  Jali::StateVector<double, Jali::Mesh> &addvec3 = mystate.add(myvec3);

  // The mesh() functions gives references to the Mesh object and the
  // Mesh object has no == or != operator (too expensive), so make
  // sure their addresses are the same

  CHECK(&(addvec3.mesh()) != &(myvec3.mesh()));


  // Now retrieve the state vectors from the state object in different ways

  Jali::State::const_iterator itc;

  // Make sure we can retrieve the object by name

  itc = mystate.find<double>("cellvars", mesh1, Jali::Entity_kind::CELL,
                             Jali::Parallel_type::ALL);
  CHECK(mystate.end() != itc);

  // Make sure the object we retrieved is identical to the one we put in

  Jali::StateVector<double> myvec1_copy =
      *(std::static_pointer_cast<Jali::StateVector<double, Jali::Mesh>>(*itc));

  CHECK_EQUAL(myvec1.size(), myvec1_copy.size());
  for (int i = 0; i < myvec1.size(); ++i)
    CHECK_EQUAL(myvec1[i], myvec1_copy[i]);

  // Retrieve the state vector more easily as a shared_ptr

  std::shared_ptr<Jali::StateVector<double, Jali::Mesh>> myvec1_ptr;
  bool found;
  found = mystate.get<double, Jali::Mesh>("cellvars", mesh1,
                                          Jali::Entity_kind::CELL,
                                          Jali::Parallel_type::ALL,
                                          &myvec1_ptr);

  CHECK(found);
  CHECK_EQUAL(myvec1.size(), myvec1_ptr->size());
  for (int i = 0; i < myvec1.size(); ++i)
    CHECK_EQUAL(myvec1[i], (*myvec1_ptr)[i]);

  // Retrieve the state vector even more easily

  found = mystate.get<double, Jali::Mesh>("cellvars", mesh1,
                                          Jali::Entity_kind::CELL,
                                          Jali::Parallel_type::ALL,
                                          &myvec1_copy);

  CHECK(found);
  CHECK_EQUAL(myvec1.size(), myvec1_copy.size());
  for (int i = 0; i < myvec1.size(); ++i)
    CHECK_EQUAL(myvec1[i], myvec1_copy[i]);



  // Make sure the code fails if we ask for the right name but wrong entity type

  itc = mystate.find<double>("cellvars", mesh1, Jali::Entity_kind::FACE,
                             Jali::Parallel_type::ALL);
  CHECK(mystate.end() == itc);


  // Try to retrieve a different vector by name

  itc = mystate.find<double>("nodevars", mesh1, Jali::Entity_kind::NODE,
                             Jali::Parallel_type::ALL);
  CHECK(mystate.end() != itc);

  // Make sure the object we retrieved is identical to the one we put in

  Jali::StateVector<double, Jali::Mesh> myvec2_copy =
      *(std::static_pointer_cast<Jali::StateVector<double, Jali::Mesh>>(*itc));

  CHECK_EQUAL(myvec2.size(), myvec2_copy.size());
  for (int i = 0; i < myvec2.size(); ++i)
    CHECK_EQUAL(myvec2[i], myvec2_copy[i]);


  // Try to retrieve the vector by name but without giving a specific type

  itc = mystate.find<double>("nodevars", mesh1);
  CHECK(mystate.end() != itc);

  // Make sure the object we retrieved is identical to the one we put in

  myvec2_copy =
      *(std::static_pointer_cast<Jali::StateVector<double, Jali::Mesh>>(*itc));

  CHECK_EQUAL(myvec2.size(), myvec2_copy.size());
  for (int i = 0; i < myvec2.size(); ++i)
    CHECK_EQUAL(myvec2[i], myvec2_copy[i]);


  // Retrieve state data through iterators and [] operators

  Jali::State::iterator it = mystate.begin();
  while (it != mystate.end()) {
    Jali::StateVector<double, Jali::Mesh> myvec4 =
        *(std::static_pointer_cast<Jali::StateVector<double, Jali::Mesh>>(*it));

    CHECK((myvec4.name() == "cellvars" &&
           myvec4.on_what() == Jali::Entity_kind::CELL &&
           myvec4.parallel_type() == Jali::Parallel_type::ALL)
          ||
          (myvec4.name() == "cellvars2" &&
           myvec4.on_what() == Jali::Entity_kind::CELL &&
           myvec4.parallel_type() == Jali::Parallel_type::ALL)
          ||
          (myvec4.name() == "nodevars" &&
           myvec4.on_what() == Jali::Entity_kind::NODE &&
           myvec4.parallel_type() == Jali::Parallel_type::ALL));

    ++it;
  }


  myvec1_copy =
      *(std::static_pointer_cast<Jali::StateVector<double, Jali::Mesh>>(mystate[0]));
  CHECK(myvec1_copy.name() == "cellvars" &&
        myvec1_copy.on_what() == Jali::Entity_kind::CELL &&
        myvec1_copy.parallel_type() == Jali::Parallel_type::ALL);

  myvec2_copy =
      *(std::static_pointer_cast<Jali::StateVector<double, Jali::Mesh>>(mystate[1]));
  CHECK(myvec2_copy.name() == "nodevars" &&
        myvec2_copy.on_what() == Jali::Entity_kind::NODE &&
        myvec2_copy.parallel_type() == Jali::Parallel_type::ALL);

  // Print out state

  std::cout << mystate;


  // Add state vectors of different data types

  const int n_cells = 4;
  int n_nodes = 9;
  float ftest[] = {1.1, 2.2, 3.3, 4.4};
  int itest[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Vec2d vtest[n_cells];
  for (unsigned int i = 0; i < n_cells; i++) vtest[i].set(1.0*i, 2.0*i);

  Jali::MeshFactory factory(MPI_COMM_WORLD);

  std::shared_ptr<Jali::Mesh> dataMesh = factory(0.0, 0.0, 1.0, 1.0, 2, 2);
  Jali::State dstate(dataMesh);

  dstate.add("f1", dataMesh, Jali::Entity_kind::CELL, Jali::Parallel_type::ALL,
             ftest);
  dstate.add("i1", dataMesh, Jali::Entity_kind::NODE, Jali::Parallel_type::ALL,
             itest);
  dstate.add("v1", dataMesh, Jali::Entity_kind::CELL, Jali::Parallel_type::ALL,
             vtest);

  // Iterate through all state vectors and count them

  int cnt = 0;
  for (Jali::State::iterator it = dstate.begin(); it != dstate.end(); it++)
    cnt++;
  CHECK_EQUAL(cnt, 3);

  // Iterate through all cell state vectors and count them

  cnt = 0;
  for (Jali::State::permutation_type it =
           dstate.entity_begin(Jali::Entity_kind::CELL);
       it != dstate.entity_end(Jali::Entity_kind::CELL);
       it++)
    cnt++;
  CHECK_EQUAL(cnt, 2);

  // Iterate through all node state vectors and count them

  cnt = 0;
  for (Jali::State::permutation_type it =
           dstate.entity_begin(Jali::Entity_kind::NODE);
       it != dstate.entity_end(Jali::Entity_kind::NODE);
       it++)
    cnt++;
  CHECK_EQUAL(cnt, 1);

  // Iterate through all state vectors and get their type

  int testCnt = 0;
  for (Jali::State::iterator it = dstate.begin(); it != dstate.end(); it++) {
    if (typeid(float) == (*it)->get_type())
      CHECK_EQUAL(testCnt, 0);
    else if (typeid(int) == (*it)->get_type())
      CHECK_EQUAL(testCnt, 1);
    else if (typeid(Vec2d) == (*it)->get_type())
      CHECK_EQUAL(testCnt, 2);
    else
      CHECK_EQUAL(0, 1);  // This else should never be reached in this test
    testCnt++;
  }
}


TEST(Jali_State_Define_MeshTiles) {

  // Create a 6x6 mesh and ask for 4 tiles on it so that each tile has 4 cells
  constexpr int NXY = 6;  // cells in any direction
  constexpr int NTILES = 4;
  constexpr int NCELLS_PER_TILE = (NXY*NXY)/NTILES;
  constexpr int NCORNERS_PER_TILE = NCELLS_PER_TILE*4;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mymesh = mf(0.0, 0.0, 1.0, 1.0, NXY, NXY, nullptr,
                                        true, true, true, true, NTILES);

  CHECK(mymesh);

  unsigned int seed = 27;

  // Create data for the  CELLS on each tile
  double data1[NTILES][NCELLS_PER_TILE];
  for (int i = 0; i < NTILES; i++)
    for (int j = 0; j < NCELLS_PER_TILE; j++)
      data1[i][j] = (static_cast<double>(rand_r(&seed)))/RAND_MAX;

  // Create data for the CORNERS on each tile
  double data2[NTILES][NCORNERS_PER_TILE];
  for (int i = 0; i < NTILES; i++)
    for (int j = 0; j < NCORNERS_PER_TILE; j++)
      data2[i][j] = (static_cast<double>(rand_r(&seed)))/RAND_MAX;

  std::array<double, 3> data3[4][4];
  for (int i = 0; i < NTILES; i++)
    for (int j = 0; j < NCELLS_PER_TILE; j++)
      for (int k = 0; k < 3; k++)
        data3[i][j][k] = (static_cast<double>(rand_r(&seed)))/RAND_MAX;
  

  // Create a state object
  Jali::State mystate(mymesh);

  // Iterate through tiles and add state vectors to it

  int i = 0;
  for (auto const& meshtile : mymesh->tiles()) {

    auto myvec1 = mystate.add("cellvars", meshtile,
                            Jali::Entity_kind::CELL, Jali::Parallel_type::OWNED,
                            &(data1[i][0]));

    Jali::StateVector<double, Jali::MeshTile> myvec2 =
        mystate.add("cornervars",
                    meshtile,
                    Jali::Entity_kind::CORNER, Jali::Parallel_type::OWNED,
                    &(data2[i][0]));

    auto myvec3 = mystate.add("cellarrays", meshtile,
                              Jali::Entity_kind::CELL,
                              Jali::Parallel_type::OWNED,
                              &(data3[i][0]));

    ++i;
  }




  // Now iterate through the tiles, retrieve the state vectors and make
  // sure the results are what we expected

  i = 0;
  for (auto const& meshtile : mymesh->tiles()) {

    Jali::StateVector<double, Jali::MeshTile> svec1;
    bool found = mystate.get("cellvars", meshtile,
                        Jali::Entity_kind::CELL, Jali::Parallel_type::OWNED,
                        &svec1);
    CHECK(found);

    if (found) {
      CHECK_EQUAL(NCELLS_PER_TILE, svec1.size());
      for (int j = 0; j < NCELLS_PER_TILE; ++j)
        CHECK_EQUAL(data1[i][j], svec1[j]);
    }
      
    found = mystate.get("cornervars", meshtile,
                        Jali::Entity_kind::CORNER, Jali::Parallel_type::OWNED,
                        &svec1);
    CHECK(found);

    if (found) {
      CHECK_EQUAL(NCORNERS_PER_TILE, svec1.size());
      for (int j = 0; j < NCORNERS_PER_TILE; ++j)
        CHECK_EQUAL(data2[i][j], svec1[j]);
    }


    Jali::StateVector<std::array<double, 3>, Jali::MeshTile> svec2;
    found = mystate.get("cellarrays", meshtile,
                        Jali::Entity_kind::CELL, Jali::Parallel_type::OWNED,
                        &svec2);
    CHECK(found);

    if (found) {
      CHECK_EQUAL(NCELLS_PER_TILE, svec2.size());
      for (int j = 0; j < NCELLS_PER_TILE; ++j)
        for (int k = 0; k < 3; ++k)
          CHECK_EQUAL(data3[i][j][k], svec2[j][k]);
    }

    ++i;
  }

}


TEST(State_Write_Read_With_Mesh) {

  // Define mesh with 4 cells and 9 nodes

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh1 = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh1);

  // Create a state object associated with this mesh

  Jali::State mystate1(mesh1);

  // Add a state vector of scalars on cells

  std::vector<double> data1 = {1.0, 3.0, 2.5, 4.5};
  Jali::StateVector<double> & outvec1 =
      mystate1.add("cellvars", mesh1, Jali::Entity_kind::CELL,
                   Jali::Parallel_type::ALL, &(data1[0]));

  std::vector<std::array<double, 2>> data2(9);
  // Intel compiler 15.0.3 is not allowing me to initialize with curly brace list
  for (int i = 0; i < 9; i++)
    for (int j = 0; j < 3; j++)
      data2[i][j] = 0.4*i+0.1*j;

  Jali::StateVector<std::array<double, 2>> & outvec2 =
      mystate1.add("nodevars", mesh1, Jali::Entity_kind::NODE,
                   Jali::Parallel_type::ALL, &(data2[0]));

  // Export the fields to the mesh

  mystate1.export_to_mesh();

  // Export the mesh (along with the state fields) to a temporary exodus file

  bool with_fields = true;
  mesh1->write_to_exodus_file("temp.exo", with_fields);




  // Now read the mesh back in - No need of wedges, corners, faces, etc
  // so just the filename is needed

  std::shared_ptr<Jali::Mesh> mesh2 = mf("temp.exo");

  // Create a state object associated with this mesh

  Jali::State mystate2(mesh2);

  // Initialize the state object from the mesh

  mystate2.init_from_mesh();

  // Retrieve the cell field and make sure we got back what we put in

  Jali::StateVector<double, Jali::Mesh> invec1;
  bool status = mystate2.get("cellvars", mesh2, Jali::Entity_kind::CELL,
                             Jali::Parallel_type::ALL, &invec1);
  CHECK(status);

  CHECK_EQUAL(outvec1.size(), invec1.size());
  for (int i = 0; i < outvec1.size(); i++)
    CHECK_EQUAL(outvec1[i], invec1[i]);


  // Retrieve the node field and make sure we got back what we put in

  Jali::StateVector<std::array<double, 2>> invec2;
  status = mystate2.get("nodevars", mesh2, Jali::Entity_kind::NODE,
                        Jali::Parallel_type::ALL, &invec2);
  CHECK(status);

  CHECK_EQUAL(outvec2.size(), invec2.size());
  for (int i = 0; i < outvec2.size(); i++)
    for (int j = 0; j < 2; j++)
      CHECK_EQUAL(outvec2[i][j], invec2[i][j]);
}
