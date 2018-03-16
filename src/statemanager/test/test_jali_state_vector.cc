/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "mpi.h"

#include <iostream>

#include "JaliStateVector.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"

#include "UnitTest++.h"

TEST(JaliUniStateVector_Cells_Mesh) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh);

  // Check array initialization

  std::vector<double> data1 = {1.0, 3.0, 2.5, 4.5};
  Jali::UniStateVector<double> myvec1("var1", mesh, nullptr,
                                   Jali::Entity_kind::CELL,
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
  Jali::UniStateVector<double> myvec2("var2", mesh, nullptr,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::ALL, data2);

  CHECK_EQUAL(ncells, myvec2.size());
  CHECK_EQUAL(data2, myvec2[0]);
  CHECK_EQUAL(data2, myvec2[1]);
  CHECK_EQUAL(data2, myvec2[2]);
  CHECK_EQUAL(data2, myvec2[3]);

  // Check no initialization

  Jali::UniStateVector<double> myvec3("var3", mesh, nullptr,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::ALL);

  CHECK_EQUAL(ncells, myvec3.size());
}
 

TEST(JaliUniStateVectorNodes) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh);

  std::vector<double> data1 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  Jali::UniStateVector<double> myvec1("var1", mesh, nullptr,
                                   Jali::Entity_kind::NODE,
                                   Jali::Entity_type::PARALLEL_OWNED,
                                   &(data1[0]));

  int nnodes = mesh->num_entities(Jali::Entity_kind::NODE,
                                  Jali::Entity_type::PARALLEL_OWNED);
  CHECK_EQUAL(nnodes, myvec1.size());
  for (int i = 0; i < nnodes; ++i)
    CHECK_EQUAL(data1[i], myvec1[i]);
}


TEST(JaliUniStateVectorAssignCopy) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh);

  std::vector<double> data1 = {1.0, 3.0, 2.5, 4.5};
  Jali::UniStateVector<double> myvec1("var1", mesh, nullptr,
                                   Jali::Entity_kind::CELL,
                                   Jali::Entity_type::ALL,
                                   &(data1[0]));

  // Assignment (NOTE that doing myvec2 = myvec1 is a copy not an assignment)

  Jali::UniStateVector<double>& myvec2 = myvec1;

  // Pick a particular element from both vectors to compare

  CHECK_EQUAL(myvec1[2], myvec2[2]);

  // Modify that element in myvec2 and verify that the element in
  // myvec1 changed too

  myvec2[2] = 99;
  CHECK_EQUAL(myvec2[2], myvec1[2]);

  // Copy construct two ways

  Jali::UniStateVector<double> myvec3 = myvec1;
  Jali::UniStateVector<double> myvec4(myvec1);

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

TEST(JaliUniStateVectorArray) {

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

  Jali::UniStateVector<std::array<double, 2>> myvec1("var1", mesh, nullptr,
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


TEST(Jali_MultiStateVector_Cells_Mesh) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 3.0, 3.0, 3, 3);

  CHECK(mesh);

  // MMState vectors need Jali State

  std::shared_ptr<Jali::State> state = Jali::State::create(mesh);

  // Create 3 material sets in the state corresponding to a T-junction
  // configuration
  //
  //     *--------*----:---*--------*
  //     |        |    :   |        |
  //     |    0   | 0  : 2 |    2   |
  //     |        |    :   |        |
  //     *--------*----:---*--------*
  //     |        |    : 2 |    2   |
  //     |        |    +............|
  //     |    0   |  0 : 1 |    1   |
  //     *--------*----:---*--------*
  //     |        |    :   |        |
  //     |    0   |  0 : 1 |    1   |
  //     |        |    :   |        |
  //     *--------*----:---*--------*

  std::vector<std::vector<int>> matcells = {{0, 1, 3, 4, 6, 7},
                                            {1, 2, 4, 5},
                                            {4, 5, 7, 8}};
  std::vector<std::vector<int>> cell_matindex = {{0, 1, -1, 2, 3, -1, 4, 5, -1},
                                                 {-1, 0, 1, -1, 2, 3, -1, -1, -1},
                                                 {-1, -1, -1, -1, 0, 1, -1, 2, 3}};

  state->add_material("steel", matcells[0]);
  state->add_material("aluminum", matcells[1]);
  state->add_material("air", matcells[2]);

  // Create a multi-material state vector corresponding to volume
  // fractions of materials as shown in fig above. 'vfm' is laid out in
  // a material centric manner

  double **vfm = new double*[3];
  for (int i = 0; i < 3; i++)
    vfm[i] = new double[9];

  double vfarr[3][9] = {{1.0, 0.5, 0.0, 1.0, 0.5, 0.0, 1.0, 0.5, 0.0},
                        {0.0, 0.5, 1.0, 0.0, 0.25, 0.5, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.0, 0.5, 1.0}};
  for (int m = 0; m < 3; m++)
    for (int c = 0; c < 9; c++)
      vfm[m][c] = vfarr[m][c];

  Jali::MultiStateVector<double> myvec1("volfrac1", mesh, state,
                                     Jali::Entity_kind::CELL,
                                     Jali::Entity_type::ALL,
                                     Jali::Data_layout::MATERIAL_CENTRIC,
                                     (double const **) vfm);
  
  CHECK(Jali::StateVector_type::MULTIVAL == myvec1.type());
  
  // Check the operator() for MultiStateVector
  
  for (int m = 0; m < 3; m++)
    for (int c = 0; c < 9; c++)
      if (vfarr[m][c] > 0.0)
        CHECK_EQUAL(vfm[m][c], myvec1(m, c));

  for (int i = 0; i < 3; i++)
    delete [] vfm[i];
  delete [] vfm;

  // Create a multi-material state vector corresponding to volume
  // fractions of materials as shown in fig above. 'vfc' is laid out in
  // a material centric manner

  double **vfc = new double*[9];
  for (int i = 0; i < 9; i++)
    vfc[i] = new double[3];

  for (int m = 0; m < 3; m++)
    for (int c = 0; c < 9; c++)
      vfc[c][m] = vfarr[m][c];

  Jali::MultiStateVector<double> myvec2("volfrac2", mesh, state,
                                     Jali::Entity_kind::CELL,
                                     Jali::Entity_type::ALL,
                                     Jali::Data_layout::CELL_CENTRIC,
                                     (double const **) vfc);
  
  CHECK(Jali::StateVector_type::MULTIVAL == myvec2.type());
  
  // Check the operator() for MultiStateVector - Note that the indices used
  // are switched for the call to operator() for myvec2 and so we have to
  // tell the code we are requesting the data in a non-default manner
  
  for (int m = 0; m < 3; m++)
    for (int c = 0; c < 9; c++)
      if (vfarr[m][c] > 0.0)
        CHECK_EQUAL(vfc[c][m],
                    myvec2(c, m, Jali::Data_layout::CELL_CENTRIC));

  for (int i = 0; i < 9; i++)
    delete [] vfc[i];
  delete [] vfc;


  // Also the two state vectors should be equivalent

  for (int m = 0; m < 3; m++)
    for (int c = 0; c < 9; c++) 
      if (vfarr[m][c] > 0.0)
        CHECK_EQUAL(myvec1(m, c), myvec2(c, m,
                                         Jali::Data_layout::CELL_CENTRIC));


  // Get data for materials one at a time and verify it
  std::vector<std::vector<double>> matvf = {{1.0, 0.5, 1.0, 0.5, 1.0, 0.5},
                                            {0.5, 1.0, 0.25, 0.5},
                                            {0.25, 0.5, 0.5, 1.0}};

  double vol = 0.0;  // sum of material volumes over the mesh
  for (int m = 0; m < 3; m++) {
    double matvol = 0.0;  // material volume
    std::vector<double> & matvec = myvec1.get_matdata(m);
    for (int c = 0; c < matvec.size(); c++) {
      CHECK_EQUAL(matvf[m][c], matvec[c]);
      double cellvol = mesh->cell_volume(matcells[m][c]);
      matvol += matvec[c]*cellvol;
    }
    if (m == 0)
      CHECK_CLOSE(4.5, matvol, 1.0e-8);
    else
      CHECK_CLOSE(2.25, matvol, 1.0e-8);
    vol += matvol;
  }
  CHECK_CLOSE(vol, 9.0, 1.0e-8);

  // Shift the vertical interface left so that material 0 is occupying a 1/3rd
  // of cells 1, 4 and 7.

  // Modify the local matvf vectors to be what we know they should be
  // (remember the second index is a local cell index in the material)
  matvf[0][cell_matindex[0][1]] = 1.0/3.0;
  matvf[0][cell_matindex[0][4]] = 1.0/3.0;
  matvf[0][cell_matindex[0][7]] = 1.0/3.0;
  matvf[1][cell_matindex[1][1]] = 2.0/3.0;
  matvf[1][cell_matindex[1][4]] = 1.0/3.0;
  matvf[2][cell_matindex[2][4]] = 1.0/3.0;
  matvf[2][cell_matindex[2][7]] = 2.0/3.0;

  // Then modify the state vector 'myvec1' in two different ways
  // If we get the cells of the material we have to use a local
  // indexing scheme to address the entries

  std::vector<double> & matvec = myvec1.get_matdata(0);
  int cloc = cell_matindex[0][4];
  matvec[cell_matindex[0][1]] = 1.0/3.0;
  matvec[cell_matindex[0][4]] = 1.0/3.0;
  matvec[cell_matindex[0][7]] = 1.0/3.0;

  // Alternate way of modifying values in myvec1
  // If we use the operator() then we can use mesh cell indices directly
  myvec1(1, 1) = 2.0/3.0;
  myvec1(1, 4) = 1.0/3.0;

  myvec1(2, 4) = 1.0/3.0;
  myvec1(2, 7) = 2.0/3.0;
   
  vol = 0.0;  // sum of material volumes over the mesh
  for (int m = 0; m < 3; m++) {
    double matvol = 0.0;  // material volume
    std::vector<double> & matvec = myvec1.get_matdata(m);
    for (int c = 0; c < matvec.size(); c++) {
      CHECK_EQUAL(matvf[m][c], matvec[c]);
      double cellvol = mesh->cell_volume(matcells[m][c]);
      matvol += matvec[c]*cellvol;
    }
    if (m == 0)
      CHECK_CLOSE(4.0, matvol, 1.0e-8);
    else
      CHECK_CLOSE(2.5, matvol, 1.0e-8);
    vol += matvol;
  }
  CHECK_CLOSE(vol, 9.0, 1.0e-8);


  // Add a complex data type (type representing material centroids)
  // and retrieve it
  class Point2 {
   public:
    Point2(double xin, double yin) : x(xin), y(yin) {};
    Point2() : x(0.0), y(0.0) {};
    double x;
    double y;
  };


  std::vector<std::vector<Point2>> matcenarr =
      {{Point2(0.5, 0.5), Point2(1.25, 0.5), Point2(0.0, 0.0),
        Point2(0.5, 1.5), Point2(1.25, 1.5), Point2(0.0, 0.0),
        Point2(0.5, 2.5), Point2(1.25, 1.5), Point2(0.0, 0.0)},
       {Point2(0.0, 0.0), Point2(1.75, 0.5), Point2(2.5, 0.5),
        Point2(0.0, 0.0), Point2(1.75, 1.25), Point2(2.5, 1.25),
        Point2(0.0, 0.0), Point2(0.0, 0.0), Point2(0.0, 0.0)},
       {Point2(0.0, 0.0), Point2(0.0, 0.0), Point2(0.0, 0.0),
        Point2(0.0, 0.0), Point2(1.75, 1.75), Point2(2.5, 1.75),
        Point2(0.0, 0.0), Point2(1.75, 2.5), Point2(2.5, 2.5)}};

  Point2 **matcen = new Point2 *[3];
  for (int m = 0; m < 3; m++)
    matcen[m] = new Point2[9];
  for (int m = 0; m < 3; m++)
    for (int c = 0; c < 9; c++)
      matcen[m][c] = matcenarr[m][c];
  
  Jali::MultiStateVector<Point2>& cenvec =
      state->add("mat_centroids", mesh, Jali::Entity_kind::CELL,
                 Jali::Entity_type::ALL, Jali::Data_layout::MATERIAL_CENTRIC,
                 (Point2 const * const *) matcen);

  for (int m = 0; m < 3; m++) {
    std::vector<int> const& matcells = state->material_cells(m);

    std::vector<Point2> matcen_out = cenvec.get_matdata(m);

    int nmatcells = matcells.size();
    for (int ic = 0; ic < nmatcells; ic++) {
      int c = matcells[ic];
      CHECK_CLOSE(matcen[m][c].x, matcen_out[ic].x, 1.0e-12);
      CHECK_CLOSE(matcen[m][c].y, matcen_out[ic].y, 1.0e-12);
    }
  }
  
  for (int m = 0; m < 3; m++)
    delete [] matcen[m];
  delete [] matcen;
}
 
