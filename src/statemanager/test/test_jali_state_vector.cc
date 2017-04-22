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


#include "mpi.h"

#include <iostream>

#include "JaliStateVector.h"
#include "Mesh.hh"
#include "MeshFactory.hh"

#include "UnitTest++.h"

TEST(JaliStateVector_Cells_Mesh) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh);

  // Check array initialization

  std::vector<double> data1 = {1.0, 3.0, 2.5, 4.5};
  Jali::StateVector<double> myvec1("var1", mesh, Jali::Entity_kind::CELL,
                                   Jali::Entity_type::ALL, &(data1[0]));

  int ncells = mesh->num_entities(Jali::Entity_kind::CELL,
                                  Jali::Entity_type::ALL);
  CHECK_EQUAL(ncells, myvec1.size());
  CHECK_EQUAL(data1[0], myvec1[0]);
  CHECK_EQUAL(data1[1], myvec1[1]);
  CHECK_EQUAL(data1[2], myvec1[2]);
  CHECK_EQUAL(data1[3], myvec1[3]);

  std::cout << myvec1 << std::endl;  // Just to check ability to print

  // Check constant initialization

  double data2 = -1.33;
  Jali::StateVector<double> myvec2("var2", mesh, Jali::Entity_kind::CELL,
                                   Jali::Entity_type::ALL, data2);

  CHECK_EQUAL(ncells, myvec2.size());
  CHECK_EQUAL(data2, myvec2[0]);
  CHECK_EQUAL(data2, myvec2[1]);
  CHECK_EQUAL(data2, myvec2[2]);
  CHECK_EQUAL(data2, myvec2[3]);

  // Check no initialization

  Jali::StateVector<double> myvec3("var3", mesh,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::ALL);                                               
  CHECK_EQUAL(ncells, myvec3.size());
}
 

TEST(JaliStateVectorNodes) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh);

  std::vector<double> data1 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  Jali::StateVector<double> myvec1("var1", mesh, Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data1[0]));

  int nnodes = mesh->num_entities(Jali::Entity_kind::NODE,
                                  Jali::Entity_type::PARALLEL_OWNED);
  CHECK_EQUAL(nnodes, myvec1.size());
  for (int i = 0; i < nnodes; ++i)
    CHECK_EQUAL(data1[i], myvec1[i]);
}


TEST(JaliStateVectorAssignCopy) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh);

  std::vector<double> data1 = {1.0, 3.0, 2.5, 4.5};
  Jali::StateVector<double> myvec1("var1", mesh, Jali::Entity_kind::CELL,
                                   Jali::Entity_type::ALL,
                                   &(data1[0]));

  // Assignment (NOTE that doing myvec2 = myvec1 is a copy not an assignment)

  Jali::StateVector<double>& myvec2 = myvec1;

  // Pick a particular element from both vectors to compare

  CHECK_EQUAL(myvec1[2], myvec2[2]);

  // Modify that element in myvec2 and verify that the element in
  // myvec1 changed too

  myvec2[2] = 99;
  CHECK_EQUAL(myvec2[2], myvec1[2]);

  // Copy construct two ways

  Jali::StateVector<double> myvec3 = myvec1;
  Jali::StateVector<double> myvec4(myvec1);

  // Verify that one of the elements is equal

  CHECK_EQUAL(myvec1[1], myvec3[1]);
  CHECK_EQUAL(myvec1[1], myvec4[1]);

  // Modify the element in myvec1 and verify that the element did not
  // change in myvec3 and myvec4

  double old_myvec3_1 = myvec3[1];
  double old_myvec4_1 = myvec4[1];

  myvec1[1] = -99;

  CHECK_EQUAL(old_myvec3_1, myvec3[1]);
  CHECK_EQUAL(old_myvec4_1, myvec4[1]);

  // Verify output

  std::cout << myvec1 << std::endl;
}

TEST(JaliStateVectorArray) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh != NULL);

  // Intel compiler 15.0.3 on moonlight is not letting me code this as
  // std::vector<std::array<double, 2>> data1 = {{-1.0, 1.0},
  //                                            {-2.0, 2.0},
  //                                            {3.0, -3.0},
  //                                            {4.0, -4.0}};


  int ncells = mesh->num_entities(Jali::Entity_kind::CELL,
                                  Jali::Entity_type::ALL);

  std::vector<std::array<double, 2>> data1(ncells);
  data1[0][0] = -1.0; data1[0][1] = 1.0;
  data1[1][0] = -2.0; data1[1][1] = 2.0;
  data1[2][0] = 3.0; data1[2][1] = -3.0;
  data1[3][0] = 4.0; data1[3][1] = -4.0;

  Jali::StateVector<std::array<double, 2>> myvec1("var1", mesh,
                                                  Jali::Entity_kind::CELL,
                                                  Jali::Entity_type::ALL,
                                                  &(data1[0]));

  // Verify we can retrieve the data as expected

  CHECK_EQUAL(data1[3][0], myvec1[3][0]);
  CHECK_EQUAL(data1[2][1], myvec1[2][1]);

  std::array<double, 2> myarray = myvec1[2];
  CHECK_EQUAL(myarray[0], data1[2][0]);
  CHECK_EQUAL(myarray[1], data1[2][1]);

  std::cout << myvec1 << std::endl;
}
