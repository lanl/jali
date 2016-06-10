//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//
/*!
 * @file   test_boundary_ghosts.cc
 * @author Rao V. Garimella
 * @date   Mon Jun 8, 2016
 *
 * @brief  Test functionality of boundary ghosts (degenerate cells at exterior boundaries)
 *
 * We will test cells and sides (and assume wedges and corners are good)
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "Geometry.hh"

TEST(MESH_BOUNDARY_GHOSTS_2D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::Framework frameworks[] = {Jali::MSTK};
  const char *framework_names[] = {"MSTK"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::Framework);
  for (int fr = 0; fr < numframeworks; fr++) {

    Jali::Framework the_framework = frameworks[fr];
    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing for boundary ghosts with " << framework_names[fr] <<
        std::endl;

    // Create the mesh

    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Jali::FrameworkPreference prefs(factory.preference());
      prefs.clear();
      prefs.push_back(the_framework);
      factory.preference(prefs);

      std::vector<Jali::Entity_kind> entitylist = {Jali::Entity_kind::EDGE,
                                                   Jali::Entity_kind::FACE,
                                                   Jali::Entity_kind::SIDE};
      factory.included_entities(entitylist);
      factory.boundary_ghosts_requested(true);
      factory.partitioner(Jali::Partitioner_type::BLOCK);

      mesh = factory(0.0, 0.0, 1.0, 1.0, 4, 4);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);

    // Check cell counts

    int ncells_owned = mesh->num_cells<Jali::Entity_type::PARALLEL_OWNED>();
    int ncells_ghost = mesh->num_cells<Jali::Entity_type::PARALLEL_GHOST>();
    int ncells_bndry_ghost =
        mesh->num_cells<Jali::Entity_type::BOUNDARY_GHOST>();
    int ncells_all = mesh->num_cells();
    
    CHECK_EQUAL(16/nproc, ncells_owned);
    if (nproc > 1) {
      CHECK(ncells_ghost > 0);
      if (nproc == 4)
        CHECK_EQUAL(5, ncells_ghost);
    } else
      CHECK(!ncells_ghost);
    
    // 16 exterior faces, so 16 boundary ghost cells distributed over
    // however many procs
    CHECK_EQUAL(16/nproc, ncells_bndry_ghost);
    CHECK_EQUAL(ncells_owned + ncells_ghost + ncells_bndry_ghost,
                ncells_all);

    // Check side counts

    int nsides_owned = mesh->num_sides<Jali::Entity_type::PARALLEL_OWNED>();
    int nsides_ghost = mesh->num_sides<Jali::Entity_type::PARALLEL_GHOST>();
    int nsides_bndry_ghost =
        mesh->num_sides<Jali::Entity_type::BOUNDARY_GHOST>();
    int nsides_all = mesh->num_sides();
    
    CHECK_EQUAL(64/nproc, nsides_owned);
    if (nproc > 1) {
      CHECK(nsides_ghost > 0);
      if (nproc == 4)
        CHECK_EQUAL(20, nsides_ghost);
    } else
      CHECK(!nsides_ghost);
    
    // 16 exterior faces, so 16 boundary ghost cells with one sides each
    // Expect 16 boundary ghost sides distributed over however many procs
    CHECK_EQUAL(16/nproc, nsides_bndry_ghost);
    CHECK_EQUAL(nsides_owned + nsides_ghost + nsides_bndry_ghost,
                nsides_all);


    // Every face (even ones on the exterior boundary) should have two
    // cells connected to it

    for (auto const& f : mesh->faces()) {
      if (mesh->entity_get_type(Jali::Entity_kind::FACE, f) ==
          Jali::Entity_type::PARALLEL_GHOST) continue;

      Jali::Entity_ID_List fcells;
      mesh->face_get_cells(f, Jali::Entity_type::PARALLEL_ALL, &fcells);
      CHECK_EQUAL(2, fcells.size());
    }

    // Nodes on the exterior boundary should return boundary ghost cells

    for (auto const& n : mesh->nodes()) {
      if (mesh->entity_get_type(Jali::Entity_kind::NODE, n) ==
          Jali::Entity_type::PARALLEL_GHOST) continue;

      JaliGeometry::Point coord(2);
      mesh->node_get_coordinates(n, &coord);
      if (fabs(coord[0]) < 1.0e-12 || fabs(coord[0]-1) < 1.0e-12 ||
          fabs(coord[1]) < 1.0e-12 || fabs(coord[1]-1) < 1.0e-12) {
        // node is on exterior boundary - check that it knows about
        // boundary ghosts

        Jali::Entity_ID_List ncells;
        mesh->node_get_cells(n, Jali::Entity_type::PARALLEL_ALL, &ncells);

        bool found_bndry_ghost = false;
        for (auto const& c : ncells) {
          if (mesh->entity_get_type(Jali::Entity_kind::CELL, c) ==
              Jali::Entity_type::BOUNDARY_GHOST) {
            found_bndry_ghost = true;
            break;
          }
        }
        CHECK(found_bndry_ghost);
      }
    }

    // All owned and boundary ghost sides should have an opposite side

    for (auto const & s : mesh->sides<Jali::Entity_type::PARALLEL_OWNED>()) {
      Jali::Entity_ID s2 = mesh->side_get_opposite_side(s);
      CHECK(s2 >= 0);
      CHECK(s2 < nsides_all);
    }
    for (auto const & s : mesh->sides<Jali::Entity_type::BOUNDARY_GHOST>()) {
      Jali::Entity_ID s2 = mesh->side_get_opposite_side(s);
      CHECK(s2 >= 0);
      CHECK(s2 < nsides_all);
    }
  }
}


TEST(MESH_BOUNDARY_GHOSTS_3D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::Framework frameworks[] = {Jali::MSTK};
  const char *framework_names[] = {"MSTK"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::Framework);
  for (int fr = 0; fr < numframeworks; fr++) {

    Jali::Framework the_framework = frameworks[fr];
    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing for boundary ghosts with " << framework_names[fr] <<
        std::endl;

    // Create the mesh

    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Jali::FrameworkPreference prefs(factory.preference());
      prefs.clear();
      prefs.push_back(the_framework);
      factory.preference(prefs);

      std::vector<Jali::Entity_kind> entitylist = {Jali::Entity_kind::EDGE,
                                                   Jali::Entity_kind::FACE,
                                                   Jali::Entity_kind::SIDE};
      factory.included_entities(entitylist);
      factory.boundary_ghosts_requested(true);
      factory.partitioner(Jali::Partitioner_type::BLOCK);

      mesh = factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);

    // Check cell counts

    int ncells_owned = mesh->num_cells<Jali::Entity_type::PARALLEL_OWNED>();
    int ncells_ghost = mesh->num_cells<Jali::Entity_type::PARALLEL_GHOST>();
    int ncells_bndry_ghost =
        mesh->num_cells<Jali::Entity_type::BOUNDARY_GHOST>();
    int ncells_all = mesh->num_cells();
    
    CHECK_EQUAL(64/nproc, ncells_owned);
    if (nproc > 1) {
      CHECK(ncells_ghost > 0);
      if (nproc == 4)
        CHECK_EQUAL(20, ncells_ghost);  // 4x4x4 split along two directions
    } else
      CHECK(!ncells_ghost);
    
    // 96 exterior faces, so 96 boundary ghost cells distributed over
    // however many procs
    CHECK_EQUAL(96/nproc, ncells_bndry_ghost);
    CHECK_EQUAL(ncells_owned + ncells_ghost + ncells_bndry_ghost,
          ncells_all);

    // Check side counts

    int nsides_owned = mesh->num_sides<Jali::Entity_type::PARALLEL_OWNED>();
    int nsides_ghost = mesh->num_sides<Jali::Entity_type::PARALLEL_GHOST>();
    int nsides_bndry_ghost =
        mesh->num_sides<Jali::Entity_type::BOUNDARY_GHOST>();
    int nsides_all = mesh->num_sides();
    
    CHECK_EQUAL((64*24)/nproc, nsides_owned);  // 24 sides per cell, 64 cells
    if (nproc > 1) {
      CHECK(nsides_ghost > 0);
      if (nproc == 4)
        CHECK_EQUAL((20*24), nsides_ghost);
    } else
      CHECK(!nsides_ghost);
    
    // 96 exterior quad faces, so 96 boundary ghost cells with 4 sides each
    CHECK_EQUAL((96*4)/nproc, nsides_bndry_ghost);
    CHECK_EQUAL(nsides_owned + nsides_ghost + nsides_bndry_ghost,
                nsides_all);


    // Every face (even ones on the exterior boundary) should have two
    // cells connected to it

    for (auto const& f : mesh->faces()) {
      if (mesh->entity_get_type(Jali::Entity_kind::FACE, f) ==
          Jali::Entity_type::PARALLEL_GHOST) continue;

      Jali::Entity_ID_List fcells;
      mesh->face_get_cells(f, Jali::Entity_type::PARALLEL_ALL, &fcells);
      CHECK_EQUAL(2, fcells.size());
    }

    // Nodes on the exterior boundary should return boundary ghost cells

    for (auto const& n : mesh->nodes()) {
      if (mesh->entity_get_type(Jali::Entity_kind::NODE, n) ==
          Jali::Entity_type::PARALLEL_GHOST) continue;

      JaliGeometry::Point coord(2);
      mesh->node_get_coordinates(n, &coord);
      if (fabs(coord[0]) < 1.0e-12 || fabs(coord[0]-1) < 1.0e-12 ||
          fabs(coord[1]) < 1.0e-12 || fabs(coord[1]-1) < 1.0e-12 ||
          fabs(coord[2]) < 1.0e-12 || fabs(coord[2]-1) < 1.0e-12) {
        // node is on exterior boundary - check that it knows about
        // boundary ghosts

        Jali::Entity_ID_List ncells;
        mesh->node_get_cells(n, Jali::Entity_type::PARALLEL_ALL, &ncells);

        bool found_bndry_ghost = false;
        for (auto const& c : ncells) {
          if (mesh->entity_get_type(Jali::Entity_kind::CELL, c) ==
              Jali::Entity_type::BOUNDARY_GHOST) {
            found_bndry_ghost = true;
            break;
          }
        }
        CHECK(found_bndry_ghost);
      }
    }

    // All owned and boundary ghost sides should have an opposite side

    for (auto const & s : mesh->sides<Jali::Entity_type::PARALLEL_OWNED>()) {
      Jali::Entity_ID s2 = mesh->side_get_opposite_side(s);
      CHECK(s2 >= 0);
      CHECK(s2 < nsides_all);
    }
    for (auto const & s : mesh->sides<Jali::Entity_type::BOUNDARY_GHOST>()) {
      Jali::Entity_ID s2 = mesh->side_get_opposite_side(s);
      CHECK(s2 >= 0);
      CHECK(s2 < nsides_all);
    }
  }
}
