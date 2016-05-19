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
#include "BoxRegion.hh"
#include "PlaneRegion.hh"
#include "GeometricModel.hh"

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
    bool parallel = (nproc > 1);
    if (!Jali::framework_generates(the_framework, parallel, dim))
      continue;

    std::cerr << "Testing mesh tile with " << framework_names[i] <<
        std::endl;

    // Create the mesh
    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;

    bool faces_requested = true;
    bool edges_requested = (the_framework == Jali::MSTK) ? true : false;
    bool sides_requested = (the_framework == Jali::MSTK) ? true : false;
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
      if (sides_requested) entitylist.push_back(Jali::Entity_kind::SIDE);
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


    // Check that each side has one and only one owner

    if (sides_requested) {
      std::vector<int> side_owner(mesh->num_sides(), -1);
      std::vector<int> side_num_owners(mesh->num_sides(), 0);
      for (auto const& t : meshtiles) {
        int tileID = t->ID();
        for (auto const& s : t->sides<Jali::Parallel_type::OWNED>()) {
          side_owner[s] = tileID;
          side_num_owners[s]++;
        }
      }
      for (auto const& s : mesh->sides()) {
        CHECK(side_owner[s] >= 0 && side_owner[s] < num_tiles_requested);
        CHECK_EQUAL(1, side_num_owners[s]);
      }

      // Check that the master tile of a ghost side really owns the side
      
      for (auto const& t : meshtiles) {
        int tileID = t->ID();
        for (auto const& s : t->sides<Jali::Parallel_type::GHOST>())
          CHECK_EQUAL(side_owner[s], mesh->master_tile_ID_of_side(s));
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


//! Test retrieval of set entities on tiles

TEST(MESH_TILES_SETS) {

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
    bool parallel = (nproc > 1);
    if (!Jali::framework_generates(the_framework, parallel, dim))
      continue;

    std::cerr << "Testing mesh tile with " << framework_names[i] <<
        std::endl;

    // Create the mesh
    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;

    bool faces_requested = true;

    int ierr = 0;
    int aerr = 0;
    int num_tiles_requested = 8;  // number of tiles in mesh on each processor

    Jali::FrameworkPreference prefs(factory.preference());
    prefs.clear();
    prefs.push_back(the_framework);
    factory.preference(prefs);
    
    factory.included_entities({Jali::Entity_kind::FACE});
    
    // Request a mesh with tiles
    
    factory.num_tiles(num_tiles_requested);
    factory.num_ghost_layers_tile(1);
    
    
    // Create a box region which does not quite cover the whole
    // domain - only a central cube that is a bit over half the size
    // of the domain in each direction
    
    std::vector<JaliGeometry::RegionPtr> gregions;
    JaliGeometry::Point boxlo(-0.51, -0.51, -0.51), boxhi(0.51, 0.51, 0.51);
    JaliGeometry::BoxRegion box1("box1", 1, boxlo, boxhi);
    gregions.push_back(&box1);
    
    // Create a plane region consisting of the left boundary
    
    JaliGeometry::Point planepnt1(-1.0, 0.0, 0.0);
    JaliGeometry::Point planenormal1(-1.0, 0.0, 0.0);
    JaliGeometry::PlaneRegion plane1("plane1", 2, planepnt1, planenormal1);
    gregions.push_back(&plane1);
    
    // Create a geometric model with these regions
    
    JaliGeometry::GeometricModel gm(3, gregions);
    
    factory.geometric_model(&gm);

    factory.partitioner(Jali::Partitioner_type::BLOCK);

    // Create a mesh - in the parallel case, no tile will cross an MPI
    // boundary
    
    mesh = factory(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 8, 8, 8);

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);

    // Sanity check - mesh has the right number of entities in the box
    // region - should be 4x4x4

    Jali::Entity_ID_List boxcells;
    mesh->get_set_entities("box1",  Jali::Entity_kind::CELL,
                          Jali::Parallel_type::OWNED, &boxcells);
    CHECK_EQUAL(64, boxcells.size());
   

    auto const& meshtiles = mesh->tiles();
    CHECK_EQUAL(num_tiles_requested, meshtiles.size());

    // Do all the tiles have the same number of cells?

    for (auto const& t : meshtiles) {
      int tileID = t->ID();
      CHECK_EQUAL(64, t->num_cells<Jali::Parallel_type::OWNED>());
      CHECK_EQUAL(125, t->num_cells());
    }

    // Check that we are able to retrieve set entities on the tiles
    // and that they are the right entities

    std::vector<int> cell_owner(mesh->num_cells(), -1);
    std::vector<int> cell_num_owners(mesh->num_cells(), 0);
    for (auto const& t : meshtiles) {
      int tileID = t->ID();

      // Get the tile's entities inside the box

      t->get_set_entities("box1", Jali::Entity_kind::CELL,
                          Jali::Parallel_type::OWNED, &boxcells);

      // Make sure we got as many we expected - each tile has 4x4x4
      // cells and 2x2x2 cells of each tile are inside the box region

      CHECK_EQUAL(8, boxcells.size());

      // Make sure each cell is really in the box

      for (auto const& c : boxcells) {
        std::vector<JaliGeometry::Point> cellpnts;
        mesh->cell_get_coordinates(c, &cellpnts);
        for (auto const& p : cellpnts)
          CHECK(box1.inside(p));
      }

      // Get the number of owned+ghost cells

      t->get_set_entities("box1", Jali::Entity_kind::CELL,
                          Jali::Parallel_type::ALL, &boxcells);

      // Make sure we got as many we expected - each 4x4x4 set of
      // cells of the tile within the box will touch other tiles along
      // 3 surfaces and 3 lines and 1 vertex. So we can expect 3(4) +
      // 3(2) + 1 = 19 ghost cells. Along with the 8 owned cells,
      // thats 27 cells

      CHECK_EQUAL(27, boxcells.size());


      // Bounding box of tile

      JaliGeometry::Point tilelo(99.0, 99.0, 99.0);
      JaliGeometry::Point tilehi(-99.0, -99.0, -99.0);

      for (auto const& c : t->cells<Jali::Parallel_type::OWNED>()) {
        std::vector<JaliGeometry::Point> cellpnts;
        mesh->cell_get_coordinates(c, &cellpnts);
        for (auto const& p : cellpnts) {
          for (int i = 0; i < 3; ++i) {
            if (p[i] < tilelo[i]) tilelo[i] = p[i];
            if (p[i] > tilehi[i]) tilehi[i] = p[i];
          }
        }
      }

      // If the X-coordinate of the bounding box of the tile is
      // coincident with the left boundary, the querying the tile for
      // faces of plane1 should return a non-null set

      Jali::Entity_ID_List planefaces;
      t->get_set_entities("plane1", Jali::Entity_kind::FACE,
                          Jali::Parallel_type::OWNED, &planefaces);

      if (fabs(tilelo[0]-(-1)) < 1.0e-12) {
        CHECK_EQUAL(16, planefaces.size());

        for (auto const& f : planefaces) {
          Jali::Entity_ID_List fcells;
          mesh->face_get_cells(f, Jali::Parallel_type::ALL, &fcells);
          CHECK_EQUAL(1, fcells.size());

          // Get outward normal of face (normalized)

          JaliGeometry::Point fnormal = mesh->face_normal(f, false, fcells[0]);
          fnormal /= norm(fnormal);

          // Should match the plane normal

          for (int i = 0; i < 3; ++i)
            CHECK_CLOSE(planenormal1[i], fnormal[i], 1.0e-10);
        }
      }
      else
        CHECK_EQUAL(0, planefaces.size());
    }
  }
}
