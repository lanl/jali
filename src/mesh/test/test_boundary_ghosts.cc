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
#include "Geometry.hh"

TEST(MESH_BOUNDARY_GHOSTS_2D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::MeshFramework_t frameworks[] = {Jali::MSTK};
  const char *framework_names[] = {"MSTK"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::MeshFramework_t);
  for (int fr = 0; fr < numframeworks; fr++) {

    Jali::MeshFramework_t the_framework = frameworks[fr];
    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing for boundary ghosts with " << framework_names[fr] <<
        std::endl;

    // Create the mesh

    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      factory.framework(the_framework);

      std::vector<Jali::Entity_kind> entitylist = {Jali::Entity_kind::EDGE,
                                                   Jali::Entity_kind::FACE,
                                                   Jali::Entity_kind::SIDE};
      factory.included_entities(entitylist);
      factory.boundary_ghosts_requested(true);
      factory.partitioner(Jali::Partitioner_type::BLOCK);

      mesh = factory(0.0, 0.0, 1.0, 1.0, 4, 4);

    } catch (const Errors::Message& e) {
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
      mesh->face_get_cells(f, Jali::Entity_type::ALL, &fcells);
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
        mesh->node_get_cells(n, Jali::Entity_type::ALL, &ncells);

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

  const Jali::MeshFramework_t frameworks[] = {Jali::MSTK};
  const char *framework_names[] = {"MSTK"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::MeshFramework_t);
  for (int fr = 0; fr < numframeworks; fr++) {

    Jali::MeshFramework_t the_framework = frameworks[fr];
    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing for boundary ghosts with " << framework_names[fr] <<
        std::endl;

    // Create the mesh

    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      factory.framework(the_framework);

      std::vector<Jali::Entity_kind> entitylist = {Jali::Entity_kind::EDGE,
                                                   Jali::Entity_kind::FACE,
                                                   Jali::Entity_kind::SIDE};
      factory.included_entities(entitylist);
      factory.boundary_ghosts_requested(true);
      factory.partitioner(Jali::Partitioner_type::BLOCK);

      mesh = factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);

    } catch (const Errors::Message& e) {
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
      mesh->face_get_cells(f, Jali::Entity_type::ALL, &fcells);
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
        mesh->node_get_cells(n, Jali::Entity_type::ALL, &ncells);

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
