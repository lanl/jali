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
  Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2);

  CHECK(mesh != NULL);

  std::vector<double> data1 = {1.0,3.0,2.5,4.5}; 
  Jali::StateVector<double> myvec1("var1",Jali::CELL,mesh,&(data1[0]));

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
  Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2);

  CHECK(mesh != NULL);

  std::vector<double> data1 = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0}; 
  Jali::StateVector<double> myvec1("var1",Jali::NODE,mesh,&(data1[0]));

  int nnodes = mesh->num_entities(Jali::NODE,Jali::ALL);
  CHECK_EQUAL(nnodes,myvec1.size());
  for (int i = 0; i < nnodes; ++i)
    CHECK_EQUAL(data1[i],myvec1[i]);
  
}


