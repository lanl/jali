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
 * @brief  Unit tests for mesh tile functionality
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

TEST(ONE_MESH_TILE) {

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

    std::cerr << "Testing single mesh tile with " << framework_names[i] <<
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
    try {
      Jali::FrameworkPreference prefs(factory.preference());
      prefs.clear();
      prefs.push_back(the_framework);
      factory.preference(prefs);

      // Create a mesh WITHOUT tiles
      mesh = factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5, NULL,
                     faces_requested, edges_requested, wedges_requested,
                     corners_requested);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }
    
    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);


    // Create a list of cells in the center of the mesh (inside the box
    // {{0.2,0.2,0.2}, {0.8,0.8,0.8}} which should result in 27 cells)

    std::vector<int> cell_list;
    for (auto const & c : mesh->cells()) {
      JaliGeometry::Point ccen = mesh->cell_centroid(c);
      if (ccen[0] > 0.2 && ccen[0] < 0.8 &&
          ccen[1] > 0.2 && ccen[1] < 0.8 &&
          ccen[2] > 0.2 && ccen[2] < 0.8)
        cell_list.emplace_back(c);
    }

    CHECK_EQUAL(27, cell_list.size());

    // Create the tile using the cell list

    std::shared_ptr<Jali::MeshTile> tile =
        make_meshtile(*mesh, cell_list, faces_requested, edges_requested,
                      wedges_requested, corners_requested);

    // Try to query the tile and see if it matches expectations


    // Does the tile think it belongs to the mesh it was created from

    //    CHECK(*mesh == tile->mesh());   // This comparison does not compile

    // Are the expected cells in the tile?

    std::vector<int> expected_cells = {31, 32, 33, 36, 37, 38, 41, 42, 43,
                                       56, 57, 58, 61, 62, 63, 66, 67, 68,
                                       81, 82, 83, 86, 87, 88, 91, 92, 93};

    // CHECK_EQUAL macro not able to parse the templated call correctly
    // so do it indirectly

    int nent = tile->num_cells();
    CHECK_EQUAL(expected_cells.size(), nent);

    int j = 0;
    for (auto const & cellid : tile->cells())
      CHECK_EQUAL(expected_cells[j++], cellid);

    // Are the expected nodes in the tile?

    std::vector<int> expected_nodes = {43, 44, 45, 46, 49, 50, 51, 52,
                                       55, 56, 57, 58, 61, 62, 63, 64,
                                       79, 80, 81, 82, 85, 86, 87, 88,
                                       91, 92, 93, 94, 97, 98, 99, 100,
                                       115, 116, 117, 118, 121, 122, 123, 124,
                                       127, 128, 129, 130, 133, 134, 135, 136,
                                       151, 152, 153, 154, 157, 158, 159, 160,
                                       163, 164, 165, 166, 169, 170, 171, 172};

    nent = tile->num_nodes();
    CHECK_EQUAL(expected_nodes.size(), nent);

    for (auto const & nodeid : tile->nodes())
      CHECK(std::find(expected_nodes.begin(), expected_nodes.end(), nodeid)
            != expected_nodes.end());


    // Check faces of the tile (should have 108 faces)

    nent = tile->num_faces();
    CHECK_EQUAL(108, nent);

    for (auto const & faceid : tile->faces()) {
      JaliGeometry::Point fcen = mesh->face_centroid(faceid);
      CHECK(fcen[0] > 0.19999 && fcen[0] < 0.80001 &&
            fcen[1] > 0.19999 && fcen[1] < 0.80001 &&
            fcen[2] > 0.19999 && fcen[2] < 0.80001);
    }


    // Check edges of the tile (should have 144)
    if (edges_requested) {
      nent = tile->num_edges();
      CHECK_EQUAL(144, nent);
    }

    if (wedges_requested) {
      // Check wedges of the tile (should have 27*48 = 1296)

      nent = tile->num_wedges();
      CHECK_EQUAL(1296, nent);

      // Check that the wedge centroids lie inside the box used to
      // filter the cells

      for (auto const & wedgeid : tile->wedges()) {
        std::vector<JaliGeometry::Point> wpnts;
        mesh->wedge_get_coordinates(wedgeid, &wpnts);
        JaliGeometry::Point wcen =
            (wpnts[0] + wpnts[1] + wpnts[2] + wpnts[3])/4.0;
        CHECK(wcen[0] > 0.2 && wcen[0] < 0.8 &&
              wcen[1] > 0.2 && wcen[1] < 0.8 &&
              wcen[2] > 0.2 && wcen[2] < 0.8);
      }
    }

    if (corners_requested) {
      // Check corners of the tile (should have 27*8 = 216)

      nent = tile->num_corners();
      CHECK_EQUAL(216, nent);

      // Check that the corner centroids are inside the box used to
      // filter the cells

      for (auto const & cornerid : tile->corners()) {
        std::vector<JaliGeometry::Point> cnpnts;
        mesh->corner_get_coordinates(cornerid, &cnpnts);
        JaliGeometry::Point cncen =
            (cnpnts[0] + cnpnts[1] + cnpnts[2] + cnpnts[3] +
             cnpnts[4] + cnpnts[5] + cnpnts[6] + cnpnts[7])/8.0;
        CHECK(cncen[0] > 0.2 && cncen[0] < 0.8 &&
              cncen[1] > 0.2 && cncen[1] < 0.8 &&
              cncen[2] > 0.2 && cncen[2] < 0.8);
      }
    }
  }

}

