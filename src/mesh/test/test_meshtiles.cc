//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//
// -------------------------------------------------------------
/**
 * @file   test_tiles.cc
 * @author Rao V. Garimella
 * @date   Mar 4, 2016
 *
 * @brief  Unit tests for iterating through tiles of a mesh
 *
 *
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include "Mesh.hh"
#include "MeshTile.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "Point.hh"

TEST(MESH_TILES_MPI) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::Framework frameworks[] = {Jali::MSTK, Jali::Simple};
  const char *framework_names[] = {"MSTK", "Simple"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::Framework);
  Jali::Framework the_framework;
  for (int i = 0; i < numframeworks; i++) {
    // Set the framework
    the_framework = frameworks[i];
    if (!Jali::framework_available(the_framework)) continue;

    int dim = 3;
    if (nproc > 1) {
      if (!Jali::framework_generates(the_framework, true, dim))
        continue;
    } else {
      if (!Jali::framework_generates(the_framework, false, dim))
        continue;
    }
    std::cerr << "Testing mesh tile with " << framework_names[i] <<
        std::endl;

    // Create the mesh
    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;

    bool faces_requested = true;
    bool edges_requested = (the_framework == Jali::MSTK) ? true : false;
    bool wedges_requested = (the_framework == Jali::MSTK) ? true : false;
    bool corners_requested = (the_framework == Jali::MSTK) ? true : false;

    int ierr = 0;
    int aerr = 0;
    int num_tiles = 8;  // number of tiles in mesh on each compute node
    try {
      Jali::FrameworkPreference prefs(factory.preference());
      prefs.clear();
      prefs.push_back(the_framework);
      factory.preference(prefs);

      // Create a mesh with tiles

      mesh = factory(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 8, 8, 8, NULL,
                     faces_requested, edges_requested,
                     wedges_requested, corners_requested, num_tiles);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);

    std::vector<Jali::Entity_ID> allcells;

    int ncells_on_proc = mesh->num_cells<Jali::Parallel_type::OWNED>();
    int ncells_per_tile = ncells_on_proc/num_tiles;


    // Iterate through the tiles and retrieve info about them

    int ntiles = 0;
    for (auto const & t : mesh->tiles()) {
      ntiles++;

      // Make sure we have the right number of cells (typically mesh
      // tiles do not have any ghost cells)

      int ncells_owned = t->num_cells<Jali::Parallel_type::OWNED>();
      CHECK_EQUAL(ncells_per_tile, ncells_owned);
      int ncells_ghost = t->num_cells<Jali::Parallel_type::GHOST>();
      CHECK_EQUAL(0, ncells_ghost);

      // Put them into a list making sure no duplicates were encountered

      int ncells = 0;
      for (auto c : t->cells<Jali::Parallel_type::OWNED>()) {
        CHECK(std::find(allcells.begin(), allcells.end(), c) == allcells.end());
        allcells.push_back(c);
        ncells++;
      }

      // check if we have the right number of nodes (each tile should
      // have 4 cells in a square pattern which gives us 9 nodes -
      // this will fail if the assumption of the tile structure is not
      // true)

      int nnodes_owned = t->num_nodes<Jali::Parallel_type::OWNED>();
      int nnodes_ghost = t->num_nodes<Jali::Parallel_type::GHOST>();
//      CHECK_EQUAL(9, nnodes_owned + nnodes_ghost);

      // check that we can iterate through the nodes correctly

      int nnodes = 0;
      for (auto n : t->nodes<Jali::Parallel_type::OWNED>())
        nnodes++;
      for (auto n : t->nodes<Jali::Parallel_type::GHOST>())
        nnodes++;
//    CHECK_EQUAL(9, nnodes);
      CHECK_EQUAL(nnodes_owned+nnodes_ghost, nnodes);

      if (edges_requested) {
        // Assuming the same cell pattern as before, we should have 12 edges
        int nedges_owned = t->num_edges<Jali::Parallel_type::OWNED>();
        int nedges_ghost = t->num_edges<Jali::Parallel_type::GHOST>();
//        CHECK_EQUAL(12, nedges_owned + nedges_ghost);

        // check that we can iterate through the edges correctly
        int nedges = 0;
        for (auto e : t->edges<Jali::Parallel_type::OWNED>())
          nedges++;
        for (auto e : t->edges<Jali::Parallel_type::GHOST>())
          nedges++;
//       CHECK_EQUAL(12, nedges);
      }

      if (faces_requested) {
        // Assuming the same cell pattern as before, we should have 12 faces
        int nfaces_owned = t->num_faces<Jali::Parallel_type::OWNED>();
        int nfaces_ghost = t->num_faces<Jali::Parallel_type::GHOST>();
//        CHECK_EQUAL(12, nfaces_owned + nfaces_ghost);

        // check that we can iterate through the faces correctly
        int nfaces = 0;
        for (auto f : t->faces<Jali::Parallel_type::OWNED>())
          nfaces++;
        for (auto f : t->faces<Jali::Parallel_type::GHOST>())
          nfaces++;
//        CHECK_EQUAL(12, nfaces);
        CHECK_EQUAL(nfaces_owned + nfaces_ghost, nfaces);
      }

      if (wedges_requested) {
        // We should have 48 Parallel_type::OWNED wedges per cell

        int nwedges = 0;
        for (auto w : t->wedges<Jali::Parallel_type::OWNED>())
          nwedges++;

        CHECK_EQUAL(48*ncells_per_tile, nwedges);
        int nwedges_owned = t->num_wedges<Jali::Parallel_type::OWNED>();
        CHECK_EQUAL(nwedges, nwedges_owned);
        int nwedges_ghost = t->num_wedges<Jali::Parallel_type::GHOST>();
        CHECK_EQUAL(0, nwedges_ghost);
      }

      if (corners_requested) {
        // We should have 8 corners per cell

        int ncorners = 0;
        for (auto cn : t->corners<Jali::Parallel_type::OWNED>())
          ncorners++;

        CHECK_EQUAL(8*ncells_per_tile, ncorners);
        int ncorners_owned = t->num_corners<Jali::Parallel_type::OWNED>();
        CHECK_EQUAL(ncorners, ncorners_owned);
        int ncorners_ghost = t->num_corners<Jali::Parallel_type::GHOST>();
        CHECK_EQUAL(0, ncorners_ghost);
      }
    }

    CHECK_EQUAL(ncells_on_proc, allcells.size());
  }
}
