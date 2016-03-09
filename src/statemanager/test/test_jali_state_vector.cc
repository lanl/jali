/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "JaliStateVector.h"
#include "Mesh.hh"
#include "MeshFactory.hh"

#include <iostream>

#include "UnitTest++.h"
#include "mpi.h"

TEST(JaliStateVectorCells) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0,0.0,1.0,1.0,2,2);

  CHECK(mesh != NULL);

  std::vector<double> data1 = {1.0,3.0,2.5,4.5}; 
  Jali::StateVector<double> myvec1("var1",Jali::CELL,mesh.get(),&(data1[0]));

  int ncells = mesh->num_entities(Jali::CELL,Jali::ALL);
  CHECK_EQUAL(ncells,myvec1.size());
  CHECK_EQUAL(data1[0],myvec1[0]);
  CHECK_EQUAL(data1[1],myvec1[1]);
  CHECK_EQUAL(data1[2],myvec1[2]);
  CHECK_EQUAL(data1[3],myvec1[3]);

  std::cout << myvec1 << std::endl;
}
 

TEST(JaliStateVectorNodes) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0,0.0,1.0,1.0,2,2);

  CHECK(mesh != NULL);

  std::vector<double> data1 = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0}; 
  Jali::StateVector<double> myvec1("var1",Jali::NODE,mesh.get(),&(data1[0]));

  int nnodes = mesh->num_entities(Jali::NODE,Jali::ALL);
  CHECK_EQUAL(nnodes,myvec1.size());
  for (int i = 0; i < nnodes; ++i)
    CHECK_EQUAL(data1[i],myvec1[i]);
}


TEST(JaliStateVectorAssignCopy) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0,0.0,1.0,1.0,2,2);

  CHECK(mesh != NULL);

  std::vector<double> data1 = {1.0,3.0,2.5,4.5}; 
  Jali::StateVector<double> myvec1("var1",Jali::CELL,mesh.get(),&(data1[0]));
  Jali::StateVector<double> myvec2;

  myvec2 = myvec1;

  // Pick a particular element from both vectors to compare

  CHECK_EQUAL(myvec1[2],myvec2[2]);

  // Modify that element in myvec2 and verify that the element in
  // myvec1 changed too

  myvec2[2] = 99;
  CHECK_EQUAL(myvec2[2],myvec1[2]);

  // Copy construct two ways

  Jali::StateVector<double> myvec3 = myvec1;
  Jali::StateVector<double> myvec4(myvec1);

  // Verify that one of the elements is equal

  CHECK_EQUAL(myvec1[1],myvec3[1]);
  CHECK_EQUAL(myvec1[1],myvec4[1]);

  // Modify the element in myvec1 and verify that the element did not
  // change in myvec3 and myvec4

  double old_myvec3_1 = myvec3[1];
  double old_myvec4_1 = myvec4[1];

  myvec1[1] = -99;

  CHECK_EQUAL(old_myvec3_1,myvec3[1]);
  CHECK_EQUAL(old_myvec4_1,myvec4[1]);

  // Verify output

  std::cout << myvec1 << std::endl;
}

TEST(JaliStateVectorArray) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0,0.0,1.0,1.0,2,2);

  CHECK(mesh != NULL);

  // Intel compiler 15.0.3 on moonlight is not letting me code this as
  // std::vector<std::array<double,2>> data1 = {{-1.0,1.0},
  //                                            {-2.0,2.0},
  //                                            {3.0,-3.0},
  //                                            {4.0,-4.0}};


  int ncells = mesh->num_entities(Jali::CELL,Jali::ALL);

  std::vector<std::array<double,2>> data1(ncells);
  data1[0][0] = -1.0; data1[0][1] = 1.0;
  data1[1][0] = -2.0; data1[1][1] = 2.0;
  data1[2][0] = 3.0; data1[2][1] = -3.0;
  data1[3][0] = 4.0; data1[3][1] = -4.0;

  Jali::StateVector<std::array<double,2>> myvec1("var1",Jali::CELL,mesh.get(),&(data1[0]));

  // Verify we can retrieve the data as expected

  CHECK_EQUAL(data1[3][0],myvec1[3][0]);
  CHECK_EQUAL(data1[2][1],myvec1[2][1]);

  std::array<double,2> myarray = myvec1[2];
  CHECK_EQUAL(myarray[0],data1[2][0]);
  CHECK_EQUAL(myarray[1],data1[2][1]);

  std::cout << myvec1 << std::endl;
}
