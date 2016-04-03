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
    int num_tiles_requested = 7;  // number of tiles in mesh on each processor
    try {
      Jali::FrameworkPreference prefs(factory.preference());
      prefs.clear();
      prefs.push_back(the_framework);
      factory.preference(prefs);

      std::vector<Jali::Entity_kind> entitylist;
      entitylist.push_back(Jali::Entity_kind::FACE);
      if (edges_requested) entitylist.push_back(Jali::Entity_kind::EDGE);
      if (wedges_requested) entitylist.push_back(Jali::Entity_kind::WEDGE);
      if (corners_requested) entitylist.push_back(Jali::Entity_kind::CORNER);
      factory.included_entities(entitylist);

      // Create a mesh with tiles

      factory.num_tiles(num_tiles_requested);
      factory.num_ghost_layers_tile(1);

      mesh = factory(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 4, 4, 4);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);


    auto const& meshtiles = mesh->tiles();
    CHECK_EQUAL(num_tiles_requested, meshtiles.size());

    // Check that each cell has one and only one owner

    std::vector<int> cell_owner(mesh->num_cells(), -1);
    std::vector<int> cell_num_owners(mesh->num_cells(), 0);
    for (auto const& t : meshtiles) {
      int tileID = t->ID();
      for (auto const& c : t->cells<Jali::Parallel_type::OWNED>()) {
        cell_owner[c] = tileID;
        cell_num_owners[c]++;
      }
    }
    for (auto const& c : mesh->cells()) {
      CHECK(cell_owner[c] >= 0 && cell_owner[c] < num_tiles_requested);
      CHECK_EQUAL(1, cell_num_owners[c]);
    }

    // Check that the master tile of a ghost cell really owns the cell
      
    for (auto const& t : meshtiles) {
      int tileID = t->ID();
      for (auto const& c : t->cells<Jali::Parallel_type::GHOST>())
        CHECK_EQUAL(cell_owner[c], mesh->master_tile_ID_of_cell(c));
    }

    std::vector<int> node_owner(mesh->num_nodes(), -1);
    std::vector<int> node_num_owners(mesh->num_nodes(), 0);
    for (auto const& t : meshtiles) {
      int tileID = t->ID();
      for (auto const& n : t->nodes<Jali::Parallel_type::OWNED>()) {
        node_owner[n] = tileID;
        node_num_owners[n]++;
      }
    }
    for (auto const& n : mesh->nodes()) {
      CHECK(node_owner[n] >= 0 && node_owner[n] < num_tiles_requested);
      CHECK_EQUAL(1, node_num_owners[n]);
    }
    for (auto const& t : meshtiles) {
      int tileID = t->ID();
      for (auto const& n : t->nodes<Jali::Parallel_type::GHOST>())
        CHECK_EQUAL(node_owner[n], mesh->master_tile_ID_of_node(n));
    }

    // Check that each face has one and only one owner

    std::vector<int> face_owner(mesh->num_faces(), -1);
    std::vector<int> face_num_owners(mesh->num_faces(), 0);
    for (auto const& t : meshtiles) {
      int tileID = t->ID();
      for (auto const& f : t->faces<Jali::Parallel_type::OWNED>()) {
        face_owner[f] = tileID;
        face_num_owners[f]++;
      }
    }
    for (auto const& f : mesh->faces()) {
      CHECK(face_owner[f] >= 0 && face_owner[f] < num_tiles_requested);
      CHECK_EQUAL(1, face_num_owners[f]);
    }

    // Check that the master tile of a ghost face really owns the face
      
    for (auto const& t : meshtiles) {
      int tileID = t->ID();
      for (auto const& f : t->faces<Jali::Parallel_type::GHOST>())
        CHECK_EQUAL(face_owner[f], mesh->master_tile_ID_of_face(f));
    }


    // Check that each edge has one and only one owner

    if (edges_requested) {
      std::vector<int> edge_owner(mesh->num_edges(), -1);
      std::vector<int> edge_num_owners(mesh->num_edges(), 0);
      for (auto const& t : meshtiles) {
        int tileID = t->ID();
        for (auto const& e : t->edges<Jali::Parallel_type::OWNED>()) {
          edge_owner[e] = tileID;
          edge_num_owners[e]++;
        }
      }
      for (auto const& e : mesh->edges()) {
        CHECK(edge_owner[e] >= 0 && edge_owner[e] < num_tiles_requested);
        CHECK_EQUAL(1, edge_num_owners[e]);
      }

      // Check that the master tile of a ghost edge really owns the edge
      
      for (auto const& t : meshtiles) {
        int tileID = t->ID();
        for (auto const& e : t->edges<Jali::Parallel_type::GHOST>())
          CHECK_EQUAL(edge_owner[e], mesh->master_tile_ID_of_edge(e));
      }
    }


    // Check that each wedge has one and only one owner

    if (wedges_requested) {
      std::vector<int> wedge_owner(mesh->num_wedges(), -1);
      std::vector<int> wedge_num_owners(mesh->num_wedges(), 0);
      for (auto const& t : meshtiles) {
        int tileID = t->ID();
        for (auto const& w : t->wedges<Jali::Parallel_type::OWNED>()) {
          wedge_owner[w] = tileID;
          wedge_num_owners[w]++;
        }
      }
      for (auto const& w : mesh->wedges()) {
        CHECK(wedge_owner[w] >= 0 && wedge_owner[w] < num_tiles_requested);
        CHECK_EQUAL(1, wedge_num_owners[w]);
      }

      // Check that the master tile of a ghost wedge really owns the wedge
      
      for (auto const& t : meshtiles) {
        int tileID = t->ID();
        for (auto const& w : t->wedges<Jali::Parallel_type::GHOST>())
          CHECK_EQUAL(wedge_owner[w], mesh->master_tile_ID_of_wedge(w));
      }
    }


    if (corners_requested) {

      // Check that each corner has one and only one owner

      std::vector<int> corner_owner(mesh->num_corners(), -1);
      std::vector<int> corner_num_owners(mesh->num_corners(), 0);
      for (auto const& t : meshtiles) {
        int tileID = t->ID();
        for (auto const& cn : t->corners<Jali::Parallel_type::OWNED>()) {
          corner_owner[cn] = tileID;
          corner_num_owners[cn]++;
        }
      }
      for (auto const& cn : mesh->corners()) {
        CHECK(corner_owner[cn] >= 0 && corner_owner[cn] < num_tiles_requested);
        CHECK_EQUAL(1, corner_num_owners[cn]);
      }

      // Check that the master tile of a ghost corner really owns the corner
      
      for (auto const& t : meshtiles) {
        int tileID = t->ID();
        for (auto const& cn : t->corners<Jali::Parallel_type::GHOST>())
          CHECK_EQUAL(corner_owner[cn], mesh->master_tile_ID_of_corner(cn));
      }
    }

  }
}
