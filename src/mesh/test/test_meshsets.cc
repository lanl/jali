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

// -------------------------------------------------------------
/**
 * @file   test_meshsets.cc
 * @author Rao V. Garimella
 * @date   Nov 2, 2016
 *
 * @brief  Unit tests for retrieving sets in a mesh
 *
 * We are testing retrieval of cells and nodes in box region, retrieval of 
 * faces in a plane region, and cells in the four types of logical regions
 * We are not yet testing the labeled set region - we would need to read an
 * Exodus II file that supports such sets and exclude the Simple mesh in
 * that test
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Point.hh"
#include "BoxRegion.hh"
#include "PlaneRegion.hh"
#include "LogicalRegion.hh"
#include "LabeledSetRegion.hh"
#include "GeometricModel.hh"

TEST(MESH_SETS_3D) {
  return;

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  int dim = 3;

  // Create a box region which does not quite cover the whole
  // domain - only a central square that is a bit over half the size
  // of the domain in each direction
    
  std::vector<JaliGeometry::RegionPtr> gregions;
  JaliGeometry::Point box1lo(-0.51, -0.51, -0.51), box1hi(0.51, 0.51, 0.51);
  JaliGeometry::BoxRegion box1("box1", 1, box1lo, box1hi);
  gregions.push_back(&box1);
    
  // Create another box region which occupies the "hi" octant of the
  // domain in all three directions

  JaliGeometry::Point box2lo(-0.01, -0.01, -0.01), box2hi(1.01, 1.01, 1.01);
  JaliGeometry::BoxRegion box2("box2", 2, box2lo, box2hi);
  gregions.push_back(&box2);
    
  // Create a third box region which occupies the "lo" octant of the
  // domain in all three directions

  JaliGeometry::Point box3lo(-1.01, -1.01, -1.01), box3hi(0.01, 0.01, 0.01);
  JaliGeometry::BoxRegion box3("box3", 3, box3lo, box3hi);
  gregions.push_back(&box3);
    
  std::vector<std::string> regnames1 = {"box1", "box2", "box3"};
  std::vector<std::string> regnames2 = {"box1", "box2"};
  std::vector<std::string> regnames3 = {"box1", "box3"};

  // Create four logical regions for each of the logical operations possible
    
  // union of 3 regions
  JaliGeometry::LogicalRegion ureg("ureg", 4, JaliGeometry::Bool_type::UNION,
                                   regnames1);
  gregions.push_back(&ureg);

  // subtract box2 from box1
  JaliGeometry::LogicalRegion sreg("sreg", 5,
                                   JaliGeometry::Bool_type::SUBTRACT,
                                   regnames2);
  gregions.push_back(&sreg);

  // intersect the box1, box3 regions
  JaliGeometry::LogicalRegion ireg("ireg", 6,
                                   JaliGeometry::Bool_type::INTERSECT,
                                   regnames3);
  gregions.push_back(&ireg);

  // Complement of the 3 regions
  JaliGeometry::LogicalRegion creg("creg", 7,
                                   JaliGeometry::Bool_type::COMPLEMENT,
                                   regnames1);
  gregions.push_back(&creg);
    

  // Create a "plane" region consisting of the left boundary
    
  JaliGeometry::Point planepnt1(-1.0, 0.0, 0.0);
  JaliGeometry::Point planenormal1(-1.0, 0.0, 0.0);
  JaliGeometry::PlaneRegion plane1("plane1", 2, planepnt1, planenormal1);
  gregions.push_back(&plane1);
    
  // Create a geometric model with these regions
    
  JaliGeometry::GeometricModel gm(dim, gregions);
    

  const Jali::MeshFramework_t frameworks[] = {Jali::MSTK, Jali::Simple};
  const char *framework_names[] = {"MSTK", "Simple"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::MeshFramework_t);
  Jali::MeshFramework_t the_framework;
  for (int i = 0; i < numframeworks; i++) {
    // Set the framework
    the_framework = frameworks[i];
    if (!Jali::framework_available(the_framework)) continue;
    
    bool parallel = (nproc > 1);
    if (!Jali::framework_generates(the_framework, parallel, dim))
      continue;
    
    std::cerr << "Testing mesh sets with " << framework_names[i] <<
        std::endl;
    
    // Create the mesh
    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;
    
    bool faces_requested = true;
    bool edges_requested = (the_framework == Jali::MSTK) ? true : false;
    bool sides_requested = false;
    bool wedges_requested = false;
    bool corners_requested = false;
    
    int ierr = 0;
    int aerr = 0;
    try {
      factory.framework(the_framework);

      std::vector<Jali::Entity_kind> entitylist;
      entitylist.push_back(Jali::Entity_kind::FACE);
      if (edges_requested) entitylist.push_back(Jali::Entity_kind::EDGE);
      if (sides_requested) entitylist.push_back(Jali::Entity_kind::SIDE);
      if (wedges_requested) entitylist.push_back(Jali::Entity_kind::WEDGE);
      if (corners_requested) entitylist.push_back(Jali::Entity_kind::CORNER);
      factory.included_entities(entitylist);

      factory.partitioner(Jali::Partitioner_type::BLOCK);
      
      factory.geometric_model(&gm);

      // Create a mesh
      
      mesh = factory(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 8, 8, 8);
      
    } catch (const Errors::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);

    // Now retrieve the sets and make sure we get what we expect

    // Check that the mesh has the right number of entities in the box
    // region - should be 4x4x4

    bool create_if_missing = true;
    std::shared_ptr<Jali::MeshSet> mset =
        mesh->find_meshset_from_region("box1", Jali::Entity_kind::CELL,
                                       create_if_missing);
    CHECK(mset != nullptr);

    int nboxcells = 0, nboxnodes = 0;
    if (parallel) {
      
      // Check local number of cells in box1 set

      int nboxcells_owned_local =
          mesh->get_set_size("box1", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_OWNED);
      if (nproc == 4)
        CHECK_EQUAL(16, nboxcells_owned_local);
      MPI_Allreduce(&nboxcells_owned_local, &nboxcells, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nboxcells_ghost_local =
          mesh->get_set_size("box1", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_GHOST);
      if (nproc == 4)
        CHECK_EQUAL(20, nboxcells_ghost_local);
          
      int nboxcells_all_local =
          mesh->get_set_size("box1", Jali::Entity_kind::CELL,
                             Jali::Entity_type::ALL);
      if (nproc == 4)
        CHECK_EQUAL(36, nboxcells_all_local);
      

      // Check local number of nodes in box1 set

      int nboxnodes_owned_local =
          mesh->get_set_size("box1", Jali::Entity_kind::NODE,
                             Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nboxnodes_owned_local, &nboxnodes, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nboxnodes_all_local =
          mesh->get_set_size("box1", Jali::Entity_kind::NODE,
                             Jali::Entity_type::ALL);
      if (nproc == 4)
        CHECK_EQUAL(80, nboxnodes_all_local);
    } else {
      nboxcells = mesh->get_set_size("box1", Jali::Entity_kind::CELL,
                                     Jali::Entity_type::PARALLEL_OWNED);
      nboxnodes = mesh->get_set_size("box1", Jali::Entity_kind::NODE,
                                     Jali::Entity_type::PARALLEL_OWNED);
    }

    // Check global number of cells and nodes in box1 set
    CHECK_EQUAL(64, nboxcells);
    CHECK_EQUAL(125, nboxnodes);

    // Make sure the centroid of each returned cell is in box1

    Jali::Entity_ID_List boxcells;
    mesh->get_set_entities("box1", Jali::Entity_kind::CELL,
                           Jali::Entity_type::ALL, &boxcells);
    bool cells_in_box = true;
    for (auto const& c : boxcells) {
      JaliGeometry::Point ccen = mesh->cell_centroid(c);
      for (int i = 0; i < 3; i++) {
        if (ccen[i] < box1lo[i] || ccen[i] > box1hi[i]) {
          cells_in_box = false;
          break;
        }
      }
      if (!cells_in_box) break;
    }
    CHECK(cells_in_box);

    // Make sure the nodes returned in set lie in the box1
    Jali::Entity_ID_List boxnodes;
    mesh->get_set_entities("box1", Jali::Entity_kind::NODE,
                           Jali::Entity_type::ALL, &boxnodes);
    bool nodes_in_box = true;
    for (auto const& n : boxnodes) {
      JaliGeometry::Point nxyz;
      mesh->node_get_coordinates(n, &nxyz);
      for (int i = 0; i < 3; i++) {
        if (nxyz[i] < box1lo[i] || nxyz[i] > box1hi[i]) {
          nodes_in_box = false;
          break;
        }
      }
      if (!nodes_in_box) break;
    }
    CHECK(nodes_in_box);


    // Check the number of cells in the union of box1, box2 and box3

    mset = mesh->find_meshset_from_region("ureg", Jali::Entity_kind::CELL,
                                          create_if_missing);
    CHECK(mset != nullptr);

    int nuregcells = 0;
    int nuregcells_local = 0;
    nuregcells_local = mesh->get_set_size("ureg", Jali::Entity_kind::CELL,
                                          Jali::Entity_type::PARALLEL_OWNED);
    if (parallel)
      MPI_Allreduce(&nuregcells_local, &nuregcells, 1, MPI_INT, MPI_SUM,
                     MPI_COMM_WORLD);
    else
      nuregcells = nuregcells_local;

    int nuregcells_expected_local = 0;
    for (auto const& c : mesh->cells<Jali::Entity_type::PARALLEL_OWNED>()) {
      JaliGeometry::Point ccen = mesh->cell_centroid(c);
      if (box1.inside(ccen) || box2.inside(ccen) || box3.inside(ccen))
        nuregcells_expected_local++;
    }
    CHECK_EQUAL(nuregcells_expected_local, nuregcells_local);

    if (parallel) {
      int nuregcells_expected = 0;
      MPI_Allreduce(&nuregcells_expected_local, &nuregcells_expected, 1,
                     MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      CHECK_EQUAL(nuregcells_expected, nuregcells);
    }

    
    // Check the number of cells in the subtraction of box2 from box1

    mset = mesh->find_meshset_from_region("sreg", Jali::Entity_kind::CELL,
                                          create_if_missing);
    CHECK(mset != nullptr);

    int nsregcells = 0;
    int nsregcells_local = 0;
      nsregcells_local = mesh->get_set_size("sreg", Jali::Entity_kind::CELL,
                                            Jali::Entity_type::PARALLEL_OWNED);
    if (parallel)
      MPI_Allreduce(&nsregcells_local, &nsregcells, 1, MPI_INT, MPI_SUM,
                     MPI_COMM_WORLD);
    else
      nsregcells = nsregcells_local;

    int nsregcells_expected_local = 0;
    for (auto const& c : mesh->cells<Jali::Entity_type::PARALLEL_OWNED>()) {
      JaliGeometry::Point ccen = mesh->cell_centroid(c);
      if (box1.inside(ccen) && !box2.inside(ccen))
        nsregcells_expected_local++;
    }
    CHECK_EQUAL(nsregcells_expected_local, nsregcells_local);

    if (parallel) {
      int nsregcells_expected = 0;
      MPI_Allreduce(&nsregcells_expected_local, &nsregcells_expected, 1,
                     MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      CHECK_EQUAL(nsregcells_expected, nsregcells);
    }

    

    // Check the number of cells in the intersection of box1 and box3

    mset = mesh->find_meshset_from_region("ireg", Jali::Entity_kind::CELL,
                                          create_if_missing);
    CHECK(mset != nullptr);

    int niregcells = 0;
    int niregcells_local = 0;
    niregcells_local = mesh->get_set_size("ireg", Jali::Entity_kind::CELL,
                                          Jali::Entity_type::PARALLEL_OWNED);
    if (parallel)
      MPI_Allreduce(&niregcells_local, &niregcells, 1, MPI_INT, MPI_SUM,
                     MPI_COMM_WORLD);
    else
      niregcells = niregcells_local;

    int niregcells_expected_local = 0;
    for (auto const& c : mesh->cells<Jali::Entity_type::PARALLEL_OWNED>()) {
      JaliGeometry::Point ccen = mesh->cell_centroid(c);
      if (box1.inside(ccen) && box3.inside(ccen))
        niregcells_expected_local++;
    }
    CHECK_EQUAL(niregcells_expected_local, niregcells_local);

    if (parallel) {
      int niregcells_expected = 0;
      MPI_Allreduce(&niregcells_expected_local, &niregcells_expected, 1,
                     MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      CHECK_EQUAL(niregcells_expected, niregcells);
    }


    // Check the number of cells in the complement of box1, box2 and box3

    mset = mesh->find_meshset_from_region("creg", Jali::Entity_kind::CELL,
                                          create_if_missing);
    CHECK(mset != nullptr);

    int ncregcells = 0;
    int ncregcells_local = 0;
    ncregcells_local = mesh->get_set_size("creg", Jali::Entity_kind::CELL,
                                          Jali::Entity_type::PARALLEL_OWNED);
    if (parallel)
      MPI_Allreduce(&ncregcells_local, &ncregcells, 1, MPI_INT, MPI_SUM,
                     MPI_COMM_WORLD);
    else
      ncregcells = ncregcells_local;

    int ncregcells_expected_local = 0;
    for (auto const& c : mesh->cells<Jali::Entity_type::PARALLEL_OWNED>()) {
      JaliGeometry::Point ccen = mesh->cell_centroid(c);
      if (!box1.inside(ccen) && !box2.inside(ccen) && !box3.inside(ccen))
        ncregcells_expected_local++;
    }
    CHECK_EQUAL(ncregcells_expected_local, ncregcells_local);

    if (parallel) {
      int ncregcells_expected = 0;
      MPI_Allreduce(&ncregcells_expected_local, &ncregcells_expected, 1,
                     MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      CHECK_EQUAL(ncregcells_expected, ncregcells);
    }


    // Check that the mesh returns the right number of faces for the
    // left plane

    mset = mesh->find_meshset_from_region("plane1", Jali::Entity_kind::FACE,
                                          create_if_missing);
    CHECK(mset != nullptr);

    Jali::Entity_ID_List planefaces;
    mesh->get_set_entities("plane1", Jali::Entity_kind::FACE,
                           Jali::Entity_type::PARALLEL_OWNED, &planefaces);

    // In the parallel case, with 4 procs we don't know if the face
    // will be cut zero, one or two times by the partitioner

    if (!parallel)
      CHECK_EQUAL(64, planefaces.size());
    else {
      int nfaces_local = planefaces.size();
      int nfaces = 0;
      MPI_Allreduce(&nfaces_local, &nfaces, 1, MPI_INT, MPI_SUM,
                     MPI_COMM_WORLD);
      CHECK_EQUAL(64, nfaces);
    }

    for (auto const& f : planefaces) {
      Jali::Entity_ID_List fcells;
      mesh->face_get_cells(f, Jali::Entity_type::ALL, &fcells);
      CHECK_EQUAL(1, fcells.size());

      // Get outward normal of face (normalized)
      
      JaliGeometry::Point fnormal = mesh->face_normal(f, false, fcells[0]);
      fnormal /= norm(fnormal);
      
      // Should match the plane normal
      
      for (int i = 0; i < dim; ++i)
        CHECK_CLOSE(planenormal1[i], fnormal[i], 1.0e-10);
    }
  }

}


//! Test the reading of labeled sets from an Exodus II file - Only
//! MSTK, not Simple mesh can read Exodus II files

TEST(MESH_SETS_LABELED_3D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  std::string filename = "test/hex_3x3x3_sets.exo";

  int dim = 3;

  // The mesh being read has 3 material blocks, 21 sidesets
  // (facesets - different geometry/material combinations) and 21
  // nodesets (again different geometry/material combinations) AND 3
  // element sets (different from material blocks). Test the
  // retrieval of at least some of them

  std::vector<JaliGeometry::RegionPtr> gregions;

  // Two material blocks

  JaliGeometry::LabeledSetRegion lsrgn1("mat1", 1, "CELL", filename,
                                        "Exodus II", "10000");
  gregions.push_back(&lsrgn1);

  JaliGeometry::LabeledSetRegion lsrgn2("mat3", 2, "CELL", filename,
                                        "Exodus II", "30000");
  gregions.push_back(&lsrgn2);

  // Two cell/element sets

  JaliGeometry::LabeledSetRegion lsrgn3("cellset2", 3, "CELL", filename,
                                        "Exodus II", "2");
  gregions.push_back(&lsrgn3);

  JaliGeometry::LabeledSetRegion lsrgn4("cellset3", 4, "CELL", filename,
                                        "Exodus II", "3");
  gregions.push_back(&lsrgn4);

    
  // Capture 2 sidesets

  JaliGeometry::LabeledSetRegion lsrgn5("face101", 5, "FACE", filename,
                                        "Exodus II", "101");
  gregions.push_back(&lsrgn5);

  JaliGeometry::LabeledSetRegion lsrgn6("face10005", 6, "FACE", filename,
                                        "Exodus II", "10005");
  gregions.push_back(&lsrgn6);


  // Capture one nodeset (nodes on a face of the middle region)

  JaliGeometry::LabeledSetRegion lsrgn7("nodeset20004", 7, "NODE", filename,
                                        "Exodus II", "20004");
  gregions.push_back(&lsrgn7);
    
  // Create a geometric model with these regions
    
  JaliGeometry::GeometricModel gm(dim, gregions);
    

  const Jali::MeshFramework_t frameworks[] = {Jali::MSTK};
  const char *framework_names[] = {"MSTK"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::MeshFramework_t);
  Jali::MeshFramework_t the_framework;
  for (int i = 0; i < numframeworks; i++) {
    // Set the framework
    the_framework = frameworks[i];
    if (!Jali::framework_available(the_framework)) continue;

    bool parallel = (nproc > 1);

    std::cerr << "Testing labeled mesh sets with " << framework_names[i] <<
        std::endl;

    // Create the mesh
    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;

    bool faces_requested = true;
    bool edges_requested = false;
    bool sides_requested = false;
    bool wedges_requested = false;
    bool corners_requested = false;

    int ierr = 0;
    int aerr = 0;
    try {
      factory.framework(the_framework);

      std::vector<Jali::Entity_kind> entitylist;
      entitylist.push_back(Jali::Entity_kind::FACE);
      if (edges_requested) entitylist.push_back(Jali::Entity_kind::EDGE);
      if (sides_requested) entitylist.push_back(Jali::Entity_kind::SIDE);
      if (wedges_requested) entitylist.push_back(Jali::Entity_kind::WEDGE);
      if (corners_requested) entitylist.push_back(Jali::Entity_kind::CORNER);
      factory.included_entities(entitylist);

      factory.partitioner(Jali::Partitioner_type::BLOCK);

      factory.geometric_model(&gm);

      // Create a mesh
    
      mesh = factory(filename);

    } catch (const Errors::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);


    // Now retrieve the sets and make sure we get what we expect
    bool create_if_missing = true;
    std::shared_ptr<Jali::MeshSet> mset =
        mesh->find_meshset_from_region("mat1", Jali::Entity_kind::CELL,
                                       create_if_missing);
    CHECK(mset != nullptr);

    int nsetcells = 0;
    if (parallel) {
      // Check local number of cells in box1 set
      
      int nsetcells_owned_local =
          mesh->get_set_size("mat1", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nsetcells_owned_local, &nsetcells, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nsetcells_ghost_local =
          mesh->get_set_size("mat1", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_GHOST);
          
      int nsetcells_all_local =
          mesh->get_set_size("mat1", Jali::Entity_kind::CELL,
                             Jali::Entity_type::ALL);
      CHECK_EQUAL(nsetcells_all_local,
                  nsetcells_owned_local + nsetcells_ghost_local);
    } else {
      nsetcells = mesh->get_set_size("mat1", Jali::Entity_kind::CELL,
                                     Jali::Entity_type::PARALLEL_OWNED);
    }

    // Check global number of cells and nodes in mat1
    CHECK_EQUAL(9, nsetcells);

    mset = mesh->find_meshset_from_region("mat3", Jali::Entity_kind::CELL,
                                          create_if_missing);
    CHECK(mset != nullptr);

    nsetcells = 0;
    if (parallel) {
      // Check local number of cells in mat3

      int nsetcells_owned_local =
          mesh->get_set_size("mat3", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nsetcells_owned_local, &nsetcells, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nsetcells_ghost_local =
          mesh->get_set_size("mat3", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_GHOST);
          
      int nsetcells_all_local =
          mesh->get_set_size("mat3", Jali::Entity_kind::CELL,
                             Jali::Entity_type::ALL);
      CHECK_EQUAL(nsetcells_all_local,
                  nsetcells_owned_local + nsetcells_ghost_local);
    } else {
      nsetcells = mesh->get_set_size("mat3", Jali::Entity_kind::CELL,
                                     Jali::Entity_type::PARALLEL_OWNED);
    }

    // Check global number of cells and nodes in mat3
    CHECK_EQUAL(9, nsetcells);

    mset = mesh->find_meshset_from_region("cellset2", Jali::Entity_kind::CELL,
                                          create_if_missing);
    CHECK(mset != nullptr);

    nsetcells = 0;
    if (parallel) {
      // Check local number of cells in cell set 2

      int nsetcells_owned_local =
          mesh->get_set_size("cellset2", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nsetcells_owned_local, &nsetcells, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nsetcells_ghost_local =
          mesh->get_set_size("cellset2", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_GHOST);
          
      int nsetcells_all_local =
          mesh->get_set_size("cellset2", Jali::Entity_kind::CELL,
                             Jali::Entity_type::ALL);
      CHECK_EQUAL(nsetcells_all_local,
                  nsetcells_owned_local + nsetcells_ghost_local);
    } else {
      nsetcells = mesh->get_set_size("cellset2", Jali::Entity_kind::CELL,
                                     Jali::Entity_type::PARALLEL_OWNED);
    }

    // Check global number of cells and nodes in cell set 2
    CHECK_EQUAL(9, nsetcells);


    mset = mesh->find_meshset_from_region("cellset3", Jali::Entity_kind::CELL,
                                          create_if_missing);
    CHECK(mset != nullptr);

    nsetcells = 0;
    if (parallel) {
      // Check local number of cells in cell set 3

      int nsetcells_owned_local =
          mesh->get_set_size("cellset3", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nsetcells_owned_local, &nsetcells, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nsetcells_ghost_local =
          mesh->get_set_size("cellset3", Jali::Entity_kind::CELL,
                             Jali::Entity_type::PARALLEL_GHOST);
          
      int nsetcells_all_local =
          mesh->get_set_size("cellset3", Jali::Entity_kind::CELL,
                             Jali::Entity_type::ALL);
      CHECK_EQUAL(nsetcells_all_local,
                  nsetcells_owned_local + nsetcells_ghost_local);
    } else {
      nsetcells = mesh->get_set_size("cellset3", Jali::Entity_kind::CELL,
                                     Jali::Entity_type::PARALLEL_OWNED);
    }

    // Check global number of cells and nodes in cell set 3
    CHECK_EQUAL(9, nsetcells);


    mset = mesh->find_meshset_from_region("face101", Jali::Entity_kind::FACE,
                                          create_if_missing);
    CHECK(mset != nullptr);

    int nsetfaces = 0;
    if (parallel) {
      // Check local number of cells in side set 101

      int nsetfaces_owned_local =
          mesh->get_set_size("face101", Jali::Entity_kind::FACE,
                             Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nsetfaces_owned_local, &nsetfaces, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nsetfaces_ghost_local =
          mesh->get_set_size("face101", Jali::Entity_kind::FACE,
                             Jali::Entity_type::PARALLEL_GHOST);
      
      int nsetfaces_all_local =
          mesh->get_set_size("face101", Jali::Entity_kind::FACE,
                             Jali::Entity_type::ALL);
      CHECK_EQUAL(nsetfaces_all_local,
                  nsetfaces_owned_local + nsetfaces_ghost_local);
    } else {
      nsetfaces = mesh->get_set_size("face101", Jali::Entity_kind::FACE,
                                     Jali::Entity_type::PARALLEL_OWNED);
    }

    // Check global number of cells and nodes in cell set 3
    CHECK_EQUAL(9, nsetfaces);


    mset = mesh->find_meshset_from_region("face10005", Jali::Entity_kind::FACE,
                                          create_if_missing);
    CHECK(mset != nullptr);

    nsetfaces = 0;
    if (parallel) {
      // Check local number of cells in side set 10005

      int nsetfaces_owned_local =
          mesh->get_set_size("face10005", Jali::Entity_kind::FACE,
                             Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nsetfaces_owned_local, &nsetfaces, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nsetfaces_ghost_local =
          mesh->get_set_size("face10005", Jali::Entity_kind::FACE,
                             Jali::Entity_type::PARALLEL_GHOST);
          
      int nsetfaces_all_local =
          mesh->get_set_size("face10005", Jali::Entity_kind::FACE,
                             Jali::Entity_type::ALL);
      CHECK_EQUAL(nsetfaces_all_local,
                  nsetfaces_owned_local + nsetfaces_ghost_local);
    } else {
      nsetfaces = mesh->get_set_size("face10005", Jali::Entity_kind::FACE,
                                     Jali::Entity_type::PARALLEL_OWNED);
    }

    // Check global number of cells and nodes in side set 10005
    CHECK_EQUAL(3, nsetfaces);

    
    mset = mesh->find_meshset_from_region("nodeset20004",
                                          Jali::Entity_kind::NODE,
                                          create_if_missing);
    CHECK(mset != nullptr);

    int nsetnodes = 0;
    if (parallel) {
      // Check local number of nodes in nodeset 20004

      int nsetnodes_owned_local =
          mesh->get_set_size("nodeset20004", Jali::Entity_kind::NODE,
                             Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nsetnodes_owned_local, &nsetnodes, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nsetnodes_ghost_local =
          mesh->get_set_size("nodeset20004", Jali::Entity_kind::NODE,
                             Jali::Entity_type::PARALLEL_GHOST);
          
      int nsetnodes_all_local =
          mesh->get_set_size("nodeset20004", Jali::Entity_kind::NODE,
                             Jali::Entity_type::ALL);
      CHECK_EQUAL(nsetnodes_all_local,
                  nsetnodes_owned_local + nsetnodes_ghost_local);
    } else {
      nsetnodes = mesh->get_set_size("nodeset20004", Jali::Entity_kind::NODE,
                                     Jali::Entity_type::PARALLEL_OWNED);
    }

    // Check global number of nodes in node set 20004
    CHECK_EQUAL(8, nsetnodes);
    
  }
}
