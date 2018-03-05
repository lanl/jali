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



// Test to add state vectors of different data types

TEST(Jali_State_Var_Types) {

  const int n_cells = 4;
  int n_nodes = 9;
  float ftest[] = {1.1, 2.2, 3.3, 4.4};
  int itest[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Vec2d vtest[n_cells];
  for (unsigned int i = 0; i < n_cells; i++) vtest[i].set(1.0*i, 2.0*i);

  Jali::MeshFactory factory(MPI_COMM_WORLD);

  std::shared_ptr<Jali::Mesh> dataMesh = factory(0.0, 0.0, 1.0, 1.0, 2, 2);
  std::shared_ptr<Jali::State> dstate = Jali::State::create(dataMesh);

  dstate->add("fc1", dataMesh, Jali::Entity_kind::CELL, Jali::Entity_type::ALL,
              ftest);
  dstate->add("in1", dataMesh, Jali::Entity_kind::NODE,
              Jali::Entity_type::ALL, itest);
  dstate->add("vc1", dataMesh, Jali::Entity_kind::CELL, Jali::Entity_type::ALL,
              vtest);

  // Iterate through all state vectors and count them

  int cnt = 0;
  for (Jali::State::iterator it = dstate->begin(); it != dstate->end(); it++)
    cnt++;
  CHECK_EQUAL(cnt, 3);

  // Iterate through all cell state vectors and count them

  cnt = 0;
  for (Jali::State::permutation_type it =
           dstate->entity_begin(Jali::Entity_kind::CELL);
       it != dstate->entity_end(Jali::Entity_kind::CELL);
       it++)
    cnt++;
  CHECK_EQUAL(cnt, 2);

  // Iterate through all node state vectors and count them

  cnt = 0;
  for (Jali::State::permutation_type it =
           dstate->entity_begin(Jali::Entity_kind::NODE);
       it != dstate->entity_end(Jali::Entity_kind::NODE);
       it++)
    cnt++;
  CHECK_EQUAL(cnt, 1);

  // Iterate through all state vectors and get their type

  int testCnt = 0;
  for (Jali::State::iterator it = dstate->begin(); it != dstate->end(); it++) {
    std::shared_ptr<Jali::BaseStateVector> sv = *it;
    if (typeid(float) == sv->get_data_type())
      CHECK_EQUAL(testCnt, 0);
    else if (typeid(int) == sv->get_data_type())
      CHECK_EQUAL(testCnt, 1);
    else if (typeid(Vec2d) == sv->get_data_type())
      CHECK_EQUAL(testCnt, 2);
    else
      CHECK_EQUAL(0, 1);  // This else should never be reached in this test
    testCnt++;
  }
}


TEST(Jali_State_On_Mesh) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.boundary_ghosts_requested(true);
  std::shared_ptr<Jali::Mesh> mesh1 = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh1 != nullptr);

  // Create a state object

  std::shared_ptr<Jali::State> mystate = Jali::State::create(mesh1);

  // Define state vector on cells with string ids and initialized from
  // array. We requested boundary ghosts (have 8 of them in this mesh)
  // and they have to be initialized too

  std::vector<double> data1 = {1.0, 3.0, 2.5, 4.5, -1.0, -1.0, -1.0, -1.0,
                               -1.0, -1.0, -1.0, -1.0};
  Jali::StateVector<double, Jali::Mesh> myvec1("cellvars", mesh1, nullptr,
                                               Jali::Entity_kind::CELL,
                                               Jali::Entity_type::ALL,
                                               &(data1[0]));

  // Add the first vector to state using the state vector object (deep copy)

  int add_status;
  Jali::StateVector<double, Jali::Mesh>& addvec1 = mystate->add(myvec1);
  CHECK_EQUAL(addvec1.size(), myvec1.size());
  for (int i = 0; i < addvec1.size(); ++i)
    CHECK_EQUAL(addvec1[i], myvec1[i]);




  // Define a state vector on nodes with string ids and initialized from
  // array. Even if we ask for boundary ghost cells, there will be no
  // boundary ghost nodes

  double data2[9] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  Jali::StateVector<double, Jali::Mesh> myvec2("nodevars", mesh1, nullptr,
                                               Jali::Entity_kind::NODE,
                                               Jali::Entity_type::ALL,
                                               &(data2[0]));

  // Add a second vector to state using the pointer to the array
  // data

  Jali::StateVector<double, Jali::Mesh>& addvec2 =
      mystate->add("nodevars", mesh1, Jali::Entity_kind::NODE,
                  Jali::Entity_type::ALL, data2);
  CHECK_EQUAL(9, addvec2.size());
  for (int i = 0; i < addvec2.size(); ++i)
    CHECK_EQUAL(data2[i], addvec2[i]);




  // Define a state vector on cells initialized from constant value

  double constval = 5.5;
  Jali::StateVector<double, Jali::Mesh> myvec3("cellvars_const", mesh1,
                                               nullptr,
                                               Jali::Entity_kind::CELL,
                                               Jali::Entity_type::ALL,
                                               constval);


  // Add a third vector to state using the same constant value

  Jali::StateVector<double, Jali::Mesh>& addvec3 =
      mystate->add<double, Jali::Mesh, Jali::StateVector>("cellvars_const",
                   mesh1, Jali::Entity_kind::CELL, Jali::Entity_type::ALL,
                   constval);
  CHECK_EQUAL(mesh1->num_cells(), addvec3.size());
  for (int i = 0; i < addvec3.size(); ++i)
    CHECK_EQUAL(constval, addvec3[i]);



  // Add a state vector on nodes not initialized
  // to anything - NOTE THAT WE HAVE TO TELL IT THAT IT IS TYPE 'int'
  // SINCE THERE IS NO INPUT DATA TO INFER THIS FROM

  Jali::StateVector<int, Jali::Mesh>& addvec4 =
      mystate->add<int, Jali::Mesh, Jali::StateVector>("nodevars_postinit",
                                                       mesh1,
                                                       Jali::Entity_kind::NODE,
                                                       Jali::Entity_type::ALL);
  
  // Then set the entries (this should set it in the memory in the
  // state manager)

  for (int i = 0; i < mesh1->num_nodes(); ++i)
    addvec4[i] = (i+2);

  

  // Define another mesh and another statevector on that mesh

  std::shared_ptr<Jali::Mesh> mesh2 = mf(0.0, 0.0, 1.0, 1.0, 3, 3);

  std::vector<double> data5 = {1.0, 3.0, 2.5, 4.5, 1.0, 2.0, 7.0, 2.0, 9.0};
  Jali::StateVector<double, Jali::Mesh> myvec5("cellvars2", mesh2, nullptr,
                    Jali::Entity_kind::CELL, Jali::Entity_type::PARALLEL_OWNED,
                    &(data5[0]));

  // Try to add the fifth vector (defined on a different mesh) to it - it
  // should copy the data but be assigned to mesh1 instead of mesh2 (deep copy)

  Jali::StateVector<double, Jali::Mesh> &addvec5 = mystate->add(myvec5);

  // The mesh() functions gives references to the Mesh object and the
  // Mesh object has no == or != operator (too expensive), so make
  // sure their addresses are the same

  CHECK(&(addvec5.mesh()) != &(myvec5.mesh()));






  // Now retrieve the state vectors from the state object in different ways

  Jali::State::const_iterator itc;


  // Make sure retrieve first state vector by name

  itc = mystate->find<double, Jali::Mesh, Jali::StateVector>("cellvars", mesh1,
                      Jali::Entity_kind::CELL, Jali::Entity_type::ALL);
  CHECK(mystate->end() != itc);

  // Make sure the object we retrieved is identical to the one we put in

  Jali::StateVector<double> myvec1_copy =
      *(std::static_pointer_cast<Jali::StateVector<double, Jali::Mesh>>(*itc));

  CHECK_EQUAL(myvec1.size(), myvec1_copy.size());
  for (int i = 0; i < myvec1.size(); ++i)
    CHECK_EQUAL(myvec1[i], myvec1_copy[i]);


  // Retrieve the second state vector more directly

  Jali::StateVector<double, Jali::Mesh> myvec2_copy;
  bool found;
  found = mystate->get<double, Jali::Mesh>("nodevars", mesh1,
                                           Jali::Entity_kind::NODE,
                                           Jali::Entity_type::ALL,
                                           &myvec2_copy);

  CHECK(found);
  CHECK_EQUAL(mesh1->num_nodes(), myvec2_copy.size());
  for (int i = 0; i < myvec2_copy.size(); ++i)
    CHECK_EQUAL(data2[i], myvec2_copy[i]);



  // Retrieve the third state vector and check its contents

  found = mystate->get<double, Jali::Mesh>("cellvars_const", mesh1,
                                           Jali::Entity_kind::CELL,
                                           Jali::Entity_type::ALL,
                                           &myvec2_copy);
  
  CHECK(found);
  for (int i = 0; i < myvec2_copy.size(); ++i)
    CHECK_EQUAL(constval, myvec2_copy[i]);


  

  // Retrieve the fourth state vector and check its contents

  Jali::StateVector<int, Jali::Mesh> myintvec_copy;
  found = mystate->get<int, Jali::Mesh>("nodevars_postinit", mesh1,
                                        Jali::Entity_kind::NODE,
                                        Jali::Entity_type::ALL,
                                        &myintvec_copy);

  CHECK(found);
  CHECK_EQUAL(addvec4.size(), myintvec_copy.size());
  for (int i = 0; i < myintvec_copy.size(); ++i)
    CHECK_EQUAL(addvec4[i], myintvec_copy[i]);





  // Try to retrieve a vector without explicitly specifying domain type

  found = mystate->get<double>("cellvars", mesh1, Jali::Entity_kind::CELL,
                               Jali::Entity_type::ALL, &myvec2_copy);
  CHECK(found);
  CHECK_EQUAL(myvec1.size(), myvec2_copy.size());
  for (int i = 0; i < myvec2_copy.size(); ++i)
    CHECK_EQUAL(myvec1[i], myvec2_copy[i]);




  // Try to retrieve the vector by name and mesh but without the
  // kind/type of entity it lives on

  itc = mystate->find<int, Jali::Mesh, Jali::StateVector>("nodevars_postinit",
                                                          mesh1);
  CHECK(mystate->end() != itc);

  // Make sure the object we retrieved is identical to the one we put in

  myintvec_copy =
      *(std::static_pointer_cast<Jali::StateVector<int, Jali::Mesh>>(*itc));

  CHECK_EQUAL(addvec4.size(), myintvec_copy.size());
  for (int i = 0; i < addvec4.size(); ++i)
    CHECK_EQUAL(addvec4[i], myintvec_copy[i]);




  // Retrieve a shared pointer to a state vector and non-const data that
  // it is pointing to

  std::shared_ptr<Jali::StateVector<double, Jali::Mesh>> myvec2_sptr;
  found = mystate->get("cellvars", mesh1, Jali::Entity_kind::CELL,
                       Jali::Entity_type::ALL, &myvec2_sptr);
  CHECK(found);
  CHECK_EQUAL(myvec1.size(), myvec2_sptr->size());

  double *myarray = myvec2_sptr->get_raw_data();

  for (int i = 0; i < myvec1.size(); ++i)
    CHECK_EQUAL(myvec1[i], myarray[i]);


  // Retrieve a shared pointer to a state vector and const data that
  // it is pointing to

  std::shared_ptr<Jali::StateVector<double, Jali::Mesh> const> myvec2_sptr2;
  found = mystate->get("cellvars", mesh1, Jali::Entity_kind::CELL,
                       Jali::Entity_type::ALL, &myvec2_sptr2);
  CHECK(found);
  CHECK_EQUAL(myvec1.size(), myvec2_sptr->size());

  double const *myarray2 = myvec2_sptr2->get_raw_data();

  for (int i = 0; i < myvec1.size(); ++i)
    CHECK_EQUAL(myvec1[i], myarray2[i]);



  // Retrieve state vectors (not contents of a state vector) through iterators

  Jali::State::iterator it = mystate->begin();
  while (it != mystate->end()) {
    std::shared_ptr<Jali::BaseStateVector> sv = *it;
    if (sv->name() == "nodevars_postinit") {
      Jali::StateVector<int, Jali::Mesh> myvec6 =
          *(std::dynamic_pointer_cast<Jali::StateVector<int, Jali::Mesh>>(sv));

      CHECK(myvec6.entity_kind() == Jali::Entity_kind::NODE &&
            myvec6.entity_type() == Jali::Entity_type::ALL);
    } 
    else if (sv->name() == "cellvars_const") {
      Jali::StateVector<double, Jali::Mesh> myvec6 =
          *(std::dynamic_pointer_cast<Jali::StateVector<double, Jali::Mesh>>(sv));

      CHECK(myvec6.entity_kind() == Jali::Entity_kind::CELL &&
            myvec6.entity_type() == Jali::Entity_type::ALL);
    } 
    else {
      Jali::StateVector<double, Jali::Mesh> myvec6 =
          *(std::dynamic_pointer_cast<Jali::StateVector<double, Jali::Mesh>>(sv));
      
      CHECK((myvec6.name() == "cellvars" &&
             myvec6.entity_kind() == Jali::Entity_kind::CELL &&
             myvec6.entity_type() == Jali::Entity_type::ALL)
            ||
            (myvec6.name() == "cellvars2" &&
             myvec6.entity_kind() == Jali::Entity_kind::CELL &&
             myvec6.entity_type() == Jali::Entity_type::PARALLEL_OWNED)
            ||
            (myvec6.name() == "nodevars" &&
             myvec6.entity_kind() == Jali::Entity_kind::NODE &&
             myvec6.entity_type() == Jali::Entity_type::ALL));

    }
    ++it;
  }


  // Retrieve state vectors (not contents of a state vector) through [] operator

  myvec1_copy =
      *(std::dynamic_pointer_cast<Jali::StateVector<double, Jali::Mesh>>((*mystate)[0]));
  CHECK(myvec1_copy.name() == "cellvars" &&
        myvec1_copy.entity_kind() == Jali::Entity_kind::CELL &&
        myvec1_copy.entity_type() == Jali::Entity_type::ALL);


  myvec1_copy =
      *(std::dynamic_pointer_cast<Jali::StateVector<double, Jali::Mesh>>((*mystate)[1]));
  CHECK(myvec1_copy.name() == "nodevars" &&
        myvec1_copy.entity_kind() == Jali::Entity_kind::NODE &&
        myvec1_copy.entity_type() == Jali::Entity_type::ALL);





  // Make sure the code fails if we ask for the right name but wrong entity type

  itc = mystate->find<double, Jali::Mesh, Jali::StateVector>("cellvars",
                      mesh1, Jali::Entity_kind::FACE, Jali::Entity_type::ALL);
  CHECK(mystate->end() == itc);



  // Print out state

  std::cout << mystate;
}


TEST(Jali_MMState_On_Mesh) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 3.0, 3.0, 3, 3);

  int ncells = mesh->num_entities(Jali::Entity_kind::CELL,
                                  Jali::Entity_type::ALL);

  CHECK(mesh != nullptr);

  // Create a state object

  std::shared_ptr<Jali::State> mystate = Jali::State::create(mesh);


  // Define a single valued state vector on cells with string ids and
  // initialized to 0

  Jali::StateVector<double, Jali::Mesh>& rhocell =
      mystate->add<double, Jali::Mesh, Jali::StateVector>("cell_density",
                   mesh, Jali::Entity_kind::CELL, Jali::Entity_type::ALL, 0.0);


  // Define a multi-material state vector on cells to store volume fractions
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
  std::vector<std::vector<int>> cellmats = {{0}, {0, 1}, {1},
                                            {0}, {0, 1, 2}, {1, 2},
                                            {0}, {0, 2}, {2}};

  mystate->add_material("steel1", matcells[0]);
  mystate->add_material("aluminum1", matcells[1]);

  int nmats = mystate->num_materials();
  CHECK_EQUAL(2, nmats);

  mystate->add_material("aluminum1", matcells[2]);  // Name reused
  CHECK_EQUAL(nmats, mystate->num_materials());  // no change expected

  mystate->add_material("aluminum2", matcells[2]);  // Unique ID and name
  nmats++;
  CHECK_EQUAL(nmats, mystate->num_materials());  // should have incremented

  

  CHECK_EQUAL(mystate->material_name(0), "steel1");
  CHECK_EQUAL(mystate->material_name(1), "aluminum1");
  CHECK_EQUAL(mystate->material_name(2), "aluminum2");

  for (int i = 0; i < 9; i++) {
    CHECK_EQUAL(cellmats[i].size(), mystate->num_cell_materials(i));
    std::vector<int> const& cellmats2 = mystate->cell_materials(i);
    for (int j = 0; j < cellmats[i].size(); j++)
      CHECK_EQUAL(cellmats[i][j], cellmats2[j]);
  }

  // Create a multi-material state vector corresponding to volume fractions
  // of materials as shown in fig above. Because we will add another material
  // later in the test, we will make space for 4 materials

  double **vf_in = new double*[4];
  for (int i = 0; i < 4; i++)
    vf_in[i] = new double[9];

  double vfarr[4][9] = {{1.0, 0.5, 0.0, 1.0, 0.5,  0.0, 1.0, 0.5, 0.0},
                        {0.0, 0.5, 1.0, 0.0, 0.25, 0.5, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.0, 0.5, 1.0}};
  for (int m = 0; m < 3; m++)
    for (int c = 0; c < 9; c++)
      vf_in[m][c] = vfarr[m][c];
  for (int c = 0; c < 9; c++)
    vf_in[3][c] = 0.0;

  Jali::MMStateVector<double>& vf =
      mystate->add<double, Jali::Mesh>("volfrac", mesh,
                                       Jali::Entity_kind::CELL,
                                       Jali::Entity_type::ALL,
                                       Jali::Data_layout::MATERIAL_CENTRIC,
                                       (double const **) vf_in);
  
  // Create another UNINITIALIZED multimaterial vector and then assign it data
  Jali::MMStateVector<double>& vf_alt =
      mystate->add<double, Jali::Mesh, Jali::MMStateVector>("volfrac", mesh,
                   Jali::Entity_kind::CELL, Jali::Entity_type::ALL);
  vf_alt.assign(Jali::Data_layout::MATERIAL_CENTRIC, (double const **)vf_in);

  // Check that the multimaterial vectors created differently are equivalent

  for (int m = 0; m < nmats; m++) {
    std::vector<double>& vfmat = vf.get_matdata(m);
    std::vector<double>& vfmat_alt = vf_alt.get_matdata(m);
    CHECK_EQUAL(vfmat.size(), vfmat_alt.size());
    for (int c = 0; c < vfmat.size(); c++)
      CHECK_EQUAL(vfmat[c], vfmat_alt[c]);
  }
      


  // Create a multi-material state vector corresponding to density of the
  // materials (rho_0 = 10.0; rho_1 = 2.0; rho_2 = 2.0). Similar to the 
  // volume fraction array, we will keep a spot for the fourth material
  double rho_in[4] = {10.0, 2.0, 2.0, 0.0};

  Jali::MMStateVector<double, Jali::Mesh> & rhomat =
      mystate->add<double, Jali::Mesh, Jali::MMStateVector>("mat_density",
                   mesh, Jali::Entity_kind::CELL, Jali::Entity_type::ALL);

  for (int m = 0; m < nmats; m++) {
    std::vector<double>& matvec = rhomat.get_matdata(m);
    for (int c = 0; c < matvec.size(); c++) matvec[c] = rho_in[m];
  }

  
  // Expected cell density calculated as the volume-weighted average
  // of material densities in the cell

  std::vector<double> rhocell_exp(ncells, 0.0);
  for (int c = 0; c < ncells; c++)
    for (int m = 0; m < nmats; m++)
      rhocell_exp[c] += rho_in[m]*vf_in[m][c];

  // Verify against values from state vectors

  // Method 1: Easy to code but maybe poor cache performance

  for (int m = 0; m < nmats; m++)
    for (int c = 0; c < ncells; c++)
      if (vf_in[m][c] != 0.0)
        rhocell[c] += rhomat(m, c)*vf(m, c);

  for (int c = 0; c < ncells; c++)
    CHECK_CLOSE(rhocell_exp[c], rhocell[c], 1.0e-12);

  // Method 2: Less intuitive but better access patterns

  for (int c = 0; c < ncells; c++) rhocell[c] = 0.0;
  for (int m = 0; m < nmats; m++) {
    std::vector<double>& rhomatvec = rhomat.get_matdata(m);
    std::vector<double>& vfmatvec = vf.get_matdata(m);
    std::vector<int> const& matcells = mystate->material_cells(m);

    int nmatcells = matcells.size();
    for (int i = 0; i < nmatcells; i++) {
      int c = matcells[i];
      rhocell[c] += rhomatvec[i]*vfmatvec[i];
    }
  }
      
  for (int c = 0; c < ncells; c++)
    CHECK_CLOSE(rhocell_exp[c], rhocell[c], 1.0e-12);

  // Now lets introduce a 4th material into the state and have it
  // occupy the top left corner of the mesh, pushing material 1 down
  // to the bottom left corner
  //
  //     *--------*----:---*--------*
  //     |        |    :   |        |
  //     |    3   | 3  : 2 |    2   |
  //     |        |    :   |        |
  //     *--------*----:---*--------*
  //     |    3   | 3  : 2 |    2   |
  //     |.............+............|
  //     |    0   |  0 : 1 |    1   |
  //     *--------*----:---*--------*
  //     |        |    :   |        |
  //     |    0   |  0 : 1 |    1   |
  //     |        |    :   |        |
  //     *--------*----:---*--------*

  std::vector<int> matcells3 = {3, 4, 6, 7};
  mystate->add_material("steel2", matcells3);
  nmats = mystate->num_materials();

  // Material 3 is added to cells 3, 4, 6, 7 but material 0 is not removed
  // (only its volume fraction is zeroed out)
  cellmats[3].push_back(3);
  cellmats[4].push_back(3);
  cellmats[6].push_back(3);
  cellmats[7].push_back(3);
  
  for (int i = 0; i < 9; i++) {
    CHECK_EQUAL(cellmats[i].size(), mystate->num_cell_materials(i));
    std::vector<int> const& cellmats2 = mystate->cell_materials(i);
    for (int j = 0; j < cellmats[i].size(); j++)
      CHECK_EQUAL(cellmats[i][j], cellmats2[j]);
  }

  
  // Set the volume fractions material 3 in cells in the input vectors
  vf_in[3][3] = 0.5;
  vf_in[3][4] = 0.25;
  vf_in[3][6] = 1.0;
  vf_in[3][7] = 0.5;

  // Set the volume fractions material 3 in cells in the state vectors
  vf(3, 3) = 0.5;
  vf(3, 4) = 0.25;
  vf(3, 6) = 1.0;
  vf(3, 7) = 0.5;

  // Set the density of material 3 to be the same as material 0
  rho_in[3] = rho_in[0];

  std::vector<double>& rhomatvec = rhomat.get_matdata(3);
  for (auto & rho : rhomatvec)
    rho = rho_in[3];

  // Adjust the volume fractions of material 0 in the input vector
  vf_in[0][3] = 0.5;
  vf_in[0][4] = 0.25;
  vf_in[0][6] = vf_in[0][7] = 0.0;

  // Adjust the volume fractions of material 0 in the state vector (we
  // don't have rem_cells_from_material implemented yet)
  vf(0, 3) = 0.5;
  vf(0, 4) = 0.25;
  vf(0, 6) = vf(0, 7) = 0.0;


  
  // Since the new material that displaced the old material has the
  // same density, the expected cell densities should be the same as
  // before. Verify against values from state vectors

  // Method 1: Easy to code but maybe poor cache performance

  for (int c = 0; c < ncells; c++) rhocell[c] = 0.0;
  for (int m = 0; m < nmats; m++)
    for (int c = 0; c < ncells; c++)
      if (vf_in[m][c] != 0.0)
        rhocell[c] += rhomat(m, c)*vf(m, c);

  for (int c = 0; c < ncells; c++)
    CHECK_CLOSE(rhocell_exp[c], rhocell[c], 1.0e-12);

  // Method 2: Less intuitive but better access patterns

  for (int c = 0; c < ncells; c++) rhocell[c] = 0.0;
  for (int m = 0; m < nmats; m++) {
    std::vector<double>& rhomatvec = rhomat.get_matdata(m);
    std::vector<double>& vfmatvec = vf.get_matdata(m);
    std::vector<int> const& matcells = mystate->material_cells(m);

    int nmatcells = matcells.size();
    for (int i = 0; i < nmatcells; i++) {
      int c = matcells[i];
      rhocell[c] += rhomatvec[i]*vfmatvec[i];
    }
  }
      
  for (int c = 0; c < ncells; c++)
    CHECK_CLOSE(rhocell_exp[c], rhocell[c], 1.0e-12);
 
}




TEST(Jali_State_Define_MeshTiles) {

  // Create a 6x6 mesh and ask for 4 tiles on it so that each tile has
  // 9 owned cells. Right now since we are asking for 1 halo layer,
  // each tile will also have 16 ghost cells

  constexpr int NXY = 6;  // cells in any direction
  constexpr int NTILES = 4;
  constexpr int NCELLS_PER_TILE_OWNED = (NXY*NXY)/NTILES;
  constexpr int MAX_CELLS_PER_TILE_ALL = (NXY+2)*(NXY+2)/NTILES;  // upper bound
  constexpr int NCORNERS_PER_TILE_OWNED = NCELLS_PER_TILE_OWNED*4;
  constexpr int MAX_CORNERS_PER_TILE_ALL = MAX_CELLS_PER_TILE_ALL*4;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::vector<Jali::Entity_kind> entitylist = {Jali::Entity_kind::EDGE,
                                               Jali::Entity_kind::FACE,
                                               Jali::Entity_kind::WEDGE,
                                               Jali::Entity_kind::CORNER};
  mf.included_entities(entitylist);

  mf.num_tiles(NTILES);
  mf.num_ghost_layers_tile(1);
  mf.partitioner(Jali::Partitioner_type::BLOCK);

  std::shared_ptr<Jali::Mesh> mymesh = mf(0.0, 0.0, 1.0, 1.0, NXY, NXY);

  CHECK(mymesh);

  unsigned int seed = 27;

  // Create data for the  CELLS on each tile
  double data1[NTILES][NCELLS_PER_TILE_OWNED];
  for (int i = 0; i < NTILES; i++)
    for (int j = 0; j < NCELLS_PER_TILE_OWNED; j++)
      data1[i][j] = (static_cast<double>(rand_r(&seed)))/RAND_MAX;

  // Create data for the CORNERS on each tile
  double data2[NTILES][NCORNERS_PER_TILE_OWNED];
  for (int i = 0; i < NTILES; i++)
    for (int j = 0; j < NCORNERS_PER_TILE_OWNED; j++)
      data2[i][j] = (static_cast<double>(rand_r(&seed)))/RAND_MAX;

  std::array<double, 3> data3[NTILES][NCELLS_PER_TILE_OWNED];
  for (int i = 0; i < NTILES; i++)
    for (int j = 0; j < NCELLS_PER_TILE_OWNED; j++)
      for (int k = 0; k < 3; k++)
        data3[i][j][k] = (static_cast<double>(rand_r(&seed)))/RAND_MAX;
  

  std::vector<double> zeroscalar1(MAX_CELLS_PER_TILE_ALL, 0.0);
  std::vector<double> zeroscalar2(MAX_CORNERS_PER_TILE_ALL, 0.0);
  std::array<double, 3> temparray = {0.0, 0.0, 0.0};
  std::vector<std::array<double, 3>> zeroarray3(MAX_CELLS_PER_TILE_ALL,
                                                temparray);


  // Create a state object
  std::shared_ptr<Jali::State> mystate = Jali::State::create(mymesh);

  // Retrieve vector of tiles

  auto const& meshtiles = mymesh->tiles();

  // Iterate through tiles and add state vectors to it. Then
  // initialize the owned data from the 'data' vectors

  int tileID = 0;
  for (auto const& meshtile : meshtiles) {

    // Don't forget the & or else myvec1 will be a copy of the data and
    // changes to myvec1 will not be reflected in the state vector

    auto& myvec1 = mystate->add("cellvars", meshtile,
                                Jali::Entity_kind::CELL,
                                Jali::Entity_type::ALL,
                                &(zeroscalar1[0]));
    for (int j = 0; j < NCELLS_PER_TILE_OWNED; j++)
      myvec1[j] = data1[tileID][j];

    // Declare myvec2 more traditionally

    Jali::StateVector<double, Jali::MeshTile>& myvec2 =
        mystate->add("cornervars",
                     meshtile,
                     Jali::Entity_kind::CORNER,
                     Jali::Entity_type::ALL,
                     &(zeroscalar2[0]));
    for (int j = 0; j < NCORNERS_PER_TILE_OWNED; j++)
      myvec2[j] = data2[tileID][j];

    auto& myvec3 = mystate->add("cellarrays", meshtile,
                                Jali::Entity_kind::CELL,
                                Jali::Entity_type::ALL,
                                &(zeroarray3[0]));
    for (int j = 0; j < NCELLS_PER_TILE_OWNED; j++)
      myvec3[j] = data3[tileID][j];

    ++tileID;
  }


  // Now retrieve the data we stored to make sure we can find data on tiles

  Jali::StateVector<double, Jali::MeshTile> tilevec1[NTILES];
  Jali::StateVector<double, Jali::MeshTile> tilevec2[NTILES];
  Jali::StateVector<std::array<double, 3>, Jali::MeshTile> tilearray3[NTILES];

  tileID = 0;
  for (auto const& meshtile : meshtiles) {
    bool found = mystate->get("cellvars", meshtile,
                              Jali::Entity_kind::CELL, Jali::Entity_type::ALL,
                              &(tilevec1[tileID]));
    CHECK(found);
    
    found = mystate->get("cornervars", meshtile,
                         Jali::Entity_kind::CORNER, Jali::Entity_type::ALL,
                         &(tilevec2[tileID]));
    CHECK(found);
    
    found = mystate->get("cellarrays", meshtile,
                         Jali::Entity_kind::CELL, Jali::Entity_type::ALL,
                         &(tilearray3[tileID]));
    CHECK(found);
    ++tileID;
  }



  // Now update the ghost values on each tile using the state vector
  // values from owned cells on other tiles

  tileID = 0;
  for (auto const& meshtile : meshtiles) {
    
    int ncells_owned = meshtile->num_cells<Jali::Entity_type::PARALLEL_OWNED>();
    auto const& tilecells_ghost = meshtile->cells<Jali::Entity_type::PARALLEL_GHOST>();

    // Like the entity IDs in the list of owned+ghost or all entities
    // of the tile, the data for the ghost entities starts after the
    // data for the owned entities. So, start indexing into the data
    // from ncells_owned

    int j = ncells_owned;
    for (auto const& c : tilecells_ghost) {
      int tileID2 = mymesh->master_tile_ID_of_cell(c);
      auto const& meshtile2 = meshtiles[tileID2];
      auto const& tilecells2_owned = meshtile2->cells<Jali::Entity_type::PARALLEL_OWNED>();
      bool found2 = false;
      int k = 0;
      for (auto const& c2 : tilecells2_owned) {
        if (c == c2) {
          found2 = true;
          tilevec1[tileID][j] = tilevec1[tileID2][k];
          break;
        }
        else
          ++k;
      }
      ++j;
    }


    
    int ncorners_owned = meshtile->num_corners<Jali::Entity_type::PARALLEL_OWNED>();
    auto const& tilecorners_ghost = meshtile->corners<Jali::Entity_type::PARALLEL_GHOST>();

    // Like the entity IDs in the list of owned+ghost or all entities
    // of the tile, the data for the ghost entities starts after the
    // data for the owned entities. So, start indexing into the data
    // from ncorners_owned

    j = ncorners_owned;
    for (auto const& cn : tilecorners_ghost) {
      int tileID2 = mymesh->master_tile_ID_of_corner(cn);
      auto const& meshtile2 = meshtiles[tileID2];
      auto const& tilecorners2_owned = meshtile2->corners<Jali::Entity_type::PARALLEL_OWNED>();
      bool found2 = false;
      int k = 0;
      for (auto const& cn2 : tilecorners2_owned) {
        if (cn == cn2) {
          found2 = true;
          tilevec2[tileID][j] = tilevec2[tileID2][k];
          break;
        }
        else
          ++k;
      }
      ++j;
    }


    // Like the entity IDs in the list of owned+ghost or all entities
    // of the tile, the data for the ghost entities starts after the
    // data for the owned entities. So, start indexing into the data
    // from ncells_owned

    j = ncells_owned;
    for (auto const& c : tilecells_ghost) {
      int tileID2 = mymesh->master_tile_ID_of_cell(c);
      auto const& meshtile2 = meshtiles[tileID2];
      auto const& tilecells2_owned = meshtile2->cells<Jali::Entity_type::PARALLEL_OWNED>();
      bool found2 = false;
      int k = 0;
      for (auto const& c2 : tilecells2_owned) {
        if (c == c2) {
          found2 = true;
          tilearray3[tileID][j] = tilearray3[tileID2][k];
          break;
        }
        else
          ++k;
      }
      ++j;
    }

    ++tileID;
  }



  // Now retrieve the state vectors and make sure the values of the
  // owned entities match what we put in directly and the values of
  // the ghost entities are what they are supposed to be according to
  // the corresponding owned entities in adjacent tiles

  tileID = 0;
  for (auto const& meshtile : meshtiles) {

    Jali::StateVector<double, Jali::MeshTile> svec1;
    bool found = mystate->get("cellvars", meshtile,
                              Jali::Entity_kind::CELL, Jali::Entity_type::ALL,
                              &svec1);
    CHECK(found);

    if (found) {
      auto const& tilecells = meshtile->cells();
      for (int j = 0; j < tilecells.size(); ++j) {
        int c = tilecells[j];
        int tileID2 = mymesh->master_tile_ID_of_cell(c);
        if (tileID == tileID2)
          CHECK_EQUAL(data1[tileID][j], svec1[j]);
        else {
          auto const& meshtile2 = meshtiles[tileID2];
          auto const& tilecells2 = meshtile2->cells<Jali::Entity_type::PARALLEL_OWNED>();
          int k = 0;
          bool found2 = false;
          for (auto const& c2 : tilecells2) {
            if (c == c2) {
              found2 = true;
              break;
            }
            else
              ++k;
          }
          CHECK(found2);
          CHECK_EQUAL(data1[tileID2][k], svec1[j]);
        }
      }
    }
      
    found = mystate->get("cornervars", meshtile,
                         Jali::Entity_kind::CORNER, Jali::Entity_type::ALL,
                         &svec1);
    CHECK(found);

    if (found) {
      auto const& tilecorners = meshtile->corners();
      for (int j = 0; j < tilecorners.size(); ++j) {
        int cn = tilecorners[j];
        int tileID2 = mymesh->master_tile_ID_of_corner(cn);
        if (tileID == tileID2)
          CHECK_EQUAL(data2[tileID][j], svec1[j]);
        else {
          auto const& meshtile2 = meshtiles[tileID2];
          auto const& tilecorners2 = meshtile2->corners<Jali::Entity_type::PARALLEL_OWNED>();
          int k = 0;
          bool found2 = false;
          for (auto const& cn2 : tilecorners2) {
            if (cn == cn2) {
              found2 = true;
              break;
            }
            else
              ++k;
          }
          CHECK(found2);
          CHECK_EQUAL(data2[tileID2][k], svec1[j]);
        }
      }
    }


    Jali::StateVector<std::array<double, 3>, Jali::MeshTile> svec2;
    found = mystate->get("cellarrays", meshtile,
                         Jali::Entity_kind::CELL, Jali::Entity_type::ALL,
                         &svec2);
    CHECK(found);

    if (found) {
      auto const& tilecells = meshtile->cells();
      for (int j = 0; j < tilecells.size(); ++j) {
        int c = tilecells[j];
        int tileID2 = mymesh->master_tile_ID_of_cell(c);
        if (tileID == tileID2)
          for (int d = 0; d < 3; ++d)
            CHECK_EQUAL(data3[tileID][j][d], svec2[j][d]);
        else {
          auto const& meshtile2 = meshtiles[tileID2];
          auto const& tilecells2 = meshtile2->cells<Jali::Entity_type::PARALLEL_OWNED>();
          int k = 0;
          bool found2 = false;
          for (auto const& c2 : tilecells2) {
            if (c == c2) {
              found2 = true;
              break;
            }
            else
                ++k;
          }
          CHECK(found2);
          for (int d = 0; d < 3; ++d)
            CHECK_EQUAL(data3[tileID2][k][d], svec2[j][d]);
        }
      }
    }

    ++tileID;
  }
}


TEST(State_Write_Read_With_Mesh) {

  // Define mesh with 4 cells and 9 nodes

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.boundary_ghosts_requested(true);
  std::shared_ptr<Jali::Mesh> mesh1 = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  CHECK(mesh1);

  // Create a state object associated with this mesh

  std::shared_ptr<Jali::State> mystate1 = Jali::State::create(mesh1);

  // Add a state vector of scalars on cells - remember that we
  // requested boundary ghost cells - so there are 8 boundary ghosts.
  // The value assigned to the boundary ghosts has to be 0.0 because
  // when we read back the mesh (without boundary ghosts) and recreate
  // them, their values will be 0.0

  std::vector<double> data1 = {1.0, 3.0, 2.5, 4.5, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0};
  Jali::StateVector<double> & outvec1 =
      mystate1->add("cellvars", mesh1, Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, &(data1[0]));

  // Even though we requested boundary ghost cells, there are no
  // boundary ghost nodes - so ALL nodes implies OWNED nodes only
  std::vector<std::array<double, 2>> data2(9);
  // Intel compiler 15.0.3 is not allowing me to initialize with curly brace list
  for (int i = 0; i < 9; i++)
    for (int j = 0; j < 2; j++)
      data2[i][j] = 0.4*i+0.1*j;

  Jali::StateVector<std::array<double, 2>> & outvec2 =
      mystate1->add("nodevars", mesh1, Jali::Entity_kind::NODE,
                    Jali::Entity_type::ALL, &(data2[0]));

  // Export the fields to the mesh

  mystate1->export_to_mesh();

  // Export the mesh (along with the state fields) to a temporary exodus file

  bool with_fields = true;
  mesh1->write_to_exodus_file("temp.exo", with_fields);




  // Now read the mesh back in - No need of wedges, corners, faces, etc
  // so just the filename is needed

  std::shared_ptr<Jali::Mesh> mesh2 = mf("temp.exo");

  // Create a state object associated with this mesh

  std::shared_ptr<Jali::State> mystate2 = Jali::State::create(mesh2);

  // Initialize the state object from the mesh

  mystate2->init_from_mesh();

  // Retrieve the cell field and make sure we got back what we put in

  Jali::StateVector<double, Jali::Mesh> invec1;
  bool status = mystate2->get("cellvars", mesh2, Jali::Entity_kind::CELL,
                              Jali::Entity_type::ALL, &invec1);
  CHECK(status);

  CHECK_EQUAL(outvec1.size(), invec1.size());

  int nc = mesh1->num_entities(Jali::Entity_kind::CELL, Jali::Entity_type::ALL);
  for (int i = 0; i < nc; i++) {
    JaliGeometry::Point incen = mesh2->cell_centroid(i);

    // Search for another cell that has the same centroid and compare
    // the values that were output and the values that were read
    // in. We cannot rely on the cell numbering being the same
    bool found = false;
    for (int j = 0; j < nc; j++) {
      JaliGeometry::Point outcen = mesh1->cell_centroid(j);
      JaliGeometry::Point vec = incen - outcen;
      if (JaliGeometry::norm(vec) < 1.0e-12) {
        found = true;
        CHECK_EQUAL(outvec1[j], invec1[i]);
        break;
      }
    }
    CHECK(found);
  }

  // Retrieve the node field and make sure we got back what we put in

  Jali::StateVector<std::array<double, 2>> invec2;
  status = mystate2->get("nodevars", mesh2, Jali::Entity_kind::NODE,
                         Jali::Entity_type::ALL, &invec2);
  CHECK(status);

  CHECK_EQUAL(outvec2.size(), invec2.size());

  int nn = mesh1->num_entities(Jali::Entity_kind::NODE, Jali::Entity_type::ALL);
  for (int i = 0; i < nn; i++) {
    JaliGeometry::Point xyz_in;
    mesh2->node_get_coordinates(i, &xyz_in);

    // Search for another node that has the same coordinates and
    // compare the values that were output and the values that were
    // read in. We cannot rely on the node numbering being the same
    bool found = false;
    for (int j = 0; j < nn; j++) {
      JaliGeometry::Point xyz_out;
      mesh1->node_get_coordinates(j, &xyz_out);
      JaliGeometry::Point vec = xyz_in - xyz_out;
      if (JaliGeometry::norm(vec) < 1.0e-12) {
        found = true;
        CHECK_ARRAY_EQUAL(outvec2[j], invec2[i], 2);
        break;
      }
    }
    CHECK(found);
  }
}
