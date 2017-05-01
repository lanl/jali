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
#include "Point.hh"

TEST(ONE_MESH_TILE) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::MeshFramework_t frameworks[] = {Jali::MSTK, Jali::Simple};
  const char *framework_names[] = {"MSTK", "Simple"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::MeshFramework_t);
  Jali::MeshFramework_t the_framework;
  for (int i = 0; i < numframeworks; i++) {
    // Set the framework
    the_framework = frameworks[i];
    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing single mesh tile with " << framework_names[i] <<
        std::endl;

    for (int num_halo_layers = 0; num_halo_layers < 3; ++num_halo_layers) {

      std::cerr << "Testing tile with " << num_halo_layers << " halo layers\n";

      // Create the mesh
      Jali::MeshFactory factory(MPI_COMM_WORLD);
      std::shared_ptr<Jali::Mesh> mesh;

      bool faces_requested = true;
      bool edges_requested = (the_framework == Jali::MSTK) ? true : false;
      bool wedges_requested = (the_framework == Jali::MSTK) ? true : false;
      bool sides_requested = (the_framework == Jali::MSTK) ? true : false;
      bool corners_requested = (the_framework == Jali::MSTK) ? true : false;

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

        // Create a mesh WITHOUT tiles (factory default is 0)
        mesh = factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

      } catch (const Errors::Message& e) {
        std::cerr << ": mesh error: " << e.what() << std::endl;
        ierr++;
      } catch (const std::exception& e) {
        std::cerr << ": error: " << e.what() << std::endl;
        ierr++;
      }

      MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      CHECK_EQUAL(aerr, 0);



      // Create a list of cells in the center of the mesh (inside the box
      // {{0.3,0.3,0.3}, {0.7,0.7,0.7}} which should result in 64 cells)

      double delta = 0.1;
      double min_owned = 0.299;
      double max_owned = 0.701;
      double min_all = min_owned - num_halo_layers*delta;
      double max_all = max_owned + num_halo_layers*delta;

      std::vector<int> expected_owned_cells, expected_ghost_cells;
      std::vector<int> expected_all_cells;
      for (auto const & c : mesh->cells()) {
        auto ccen = mesh->cell_centroid(c);
        if (ccen[0] > min_all && ccen[0] < max_all &&
            ccen[1] > min_all && ccen[1] < max_all &&
            ccen[2] > min_all && ccen[2] < max_all) {

          expected_all_cells.push_back(c);

          if (ccen[0] > min_owned && ccen[0] < max_owned &&
              ccen[1] > min_owned && ccen[1] < max_owned &&
              ccen[2] > min_owned && ccen[2] < max_owned)
            expected_owned_cells.push_back(c);
          else
            expected_ghost_cells.push_back(c);
        }
      }



      // Create the tile using the cell list

      auto tile = make_meshtile(*mesh, expected_owned_cells, num_halo_layers,
                                faces_requested, edges_requested,
                                sides_requested, wedges_requested,
                                corners_requested);




      // Try to query the tile and see if it matches expectations


      // Are the expected owned cells in the tile? Constant number of
      // owned cells regardless of how many halo layers are requested
      // - so check only for num_halo_layers=0

      int nent, j;

      if (num_halo_layers == 0) {
        nent = tile->num_cells<Jali::Entity_type::PARALLEL_OWNED>();
        CHECK_EQUAL(expected_owned_cells.size(), nent);

        auto const& tile_owned_cells = tile->cells<Jali::Entity_type::PARALLEL_OWNED>();

        j = 0;
        for (auto const& c : expected_owned_cells) {
          if (std::find(tile_owned_cells.begin(), tile_owned_cells.end(), c) !=
              tile_owned_cells.end())
            j++;
          else
            std::cerr << "Tile does not contain owned cell " << c << "\n";
        }
        CHECK_EQUAL(j, nent);
      }

      // Are the expected ghost cells in the tile?

      nent = tile->num_cells<Jali::Entity_type::PARALLEL_GHOST>();
      CHECK_EQUAL(expected_ghost_cells.size(), nent);

      auto const& tile_ghost_cells = tile->cells<Jali::Entity_type::PARALLEL_GHOST>();

      j = 0;
      for (auto const& c : expected_ghost_cells) {
        if (std::find(tile_ghost_cells.begin(), tile_ghost_cells.end(), c) !=
            tile_ghost_cells.end())
          j++;
        else
          std::cerr << "Tile does not contain ghost cell " << c << "\n";
      }
      CHECK_EQUAL(j, nent);


      // Are the expected owned+ghost cells in the tile?

      nent = tile->num_cells();  // Default is to get all cells
      CHECK_EQUAL(expected_all_cells.size(), nent);

      auto const& tile_all_cells = tile->cells();

      j = 0;
      for (auto const & c : expected_all_cells) {
        if (std::find(tile_all_cells.begin(), tile_all_cells.end(), c) !=
            tile_all_cells.end())
          j++;
        else
          std::cerr << "Tile does not contain owned or ghost cell " << c <<
              "\n";
      }
      CHECK_EQUAL(j, nent);




      // Process the nodes of the tile

      std::vector<int> expected_owned_nodes;
      std::vector<int> expected_ghost_nodes;
      std::vector<int> expected_all_nodes;
      for (auto const & n : mesh->nodes()) {
        JaliGeometry::Point nodexyz;
        mesh->node_get_coordinates(n, &nodexyz);
        if (nodexyz[0] > min_all && nodexyz[0] < max_all &&
            nodexyz[1] > min_all && nodexyz[1] < max_all &&
            nodexyz[2] > min_all && nodexyz[2] < max_all) {

          expected_all_nodes.push_back(n);

          if (nodexyz[0] > min_owned && nodexyz[0] < max_owned &&
              nodexyz[1] > min_owned && nodexyz[1] < max_owned &&
              nodexyz[2] > min_owned && nodexyz[2] < max_owned)
            expected_owned_nodes.push_back(n);
          else
            expected_ghost_nodes.push_back(n);
        }
      }


      // Are the expected owned nodes in the tile? Since there is only
      // one tile in the mesh, this is all the nodes of the owned
      // cells. Owned nodes list will be constant regardless of how
      // many halo layers are requested so check only if num_halo_layers = 0

      if (num_halo_layers == 0) {
        nent = tile->num_nodes<Jali::Entity_type::PARALLEL_OWNED>();
        CHECK_EQUAL(expected_owned_nodes.size(), nent);

        auto const& tile_owned_nodes = tile->nodes<Jali::Entity_type::PARALLEL_OWNED>();
        j = 0;
        for (auto const& n : expected_owned_nodes) {
          if (std::find(tile_owned_nodes.begin(), tile_owned_nodes.end(), n) !=
              tile_owned_nodes.end())
            j++;
          else
            std::cerr << "Tile does not contain owned node " << n << "\n";
        }
        CHECK_EQUAL(j, nent);
      }

      // Are the expected ghost nodes in the tile?

      nent = tile->num_nodes<Jali::Entity_type::PARALLEL_GHOST>();
      CHECK_EQUAL(expected_ghost_nodes.size(), nent);

      auto const& tile_ghost_nodes = tile->nodes<Jali::Entity_type::PARALLEL_GHOST>();
      j = 0;
      for (auto const & n : expected_ghost_nodes) {
        if (std::find(tile_ghost_nodes.begin(), tile_ghost_nodes.end(), n) !=
            tile_ghost_nodes.end())
          j++;
        else
          std::cerr << "Tile does not contain ghost node " << n << "\n";
      }
      CHECK_EQUAL(j, nent);

      // Are the expected owned+ghost nodes in the tile?

      nent = tile->num_nodes();
      CHECK_EQUAL(expected_all_nodes.size(), nent);

      auto const& tile_all_nodes = tile->nodes();
      j = 0;
      for (auto const & n : expected_all_nodes) {
        if (std::find(tile_all_nodes.begin(), tile_all_nodes.end(), n) !=
            tile_all_nodes.end())
          j++;
        else
          std::cerr << "Tile does not contain owned or ghost node " << n <<
              "\n";
      }
      CHECK_EQUAL(j, nent);




      // Process the faces

      std::vector<int> expected_owned_faces;
      std::vector<int> expected_ghost_faces;
      std::vector<int> expected_all_faces;
      for (auto const& f : mesh->faces()) {
        auto fcen = mesh->face_centroid(f);
        if (fcen[0] > min_all && fcen[0] < max_all &&
            fcen[1] > min_all && fcen[1] < max_all &&
            fcen[2] > min_all && fcen[2] < max_all) {

          expected_all_faces.push_back(f);

          if (fcen[0] > min_owned && fcen[0] < max_owned &&
              fcen[1] > min_owned && fcen[1] < max_owned &&
              fcen[2] > min_owned && fcen[2] < max_owned)
            expected_owned_faces.push_back(f);
          else
            expected_ghost_faces.push_back(f);
        }
      }


      // Are the expected owned faces in the tile? Since there is only
      // one tile in the mesh, all the faces of owned cells will be
      // owned. This list will remain the same regardless of how many
      // halo layers are requested, so check only if num_halo_layers=0

      if (num_halo_layers == 0) {
        nent = tile->num_faces<Jali::Entity_type::PARALLEL_OWNED>();
        CHECK_EQUAL(expected_owned_faces.size(), nent);

        auto const& tile_owned_faces = tile->faces<Jali::Entity_type::PARALLEL_OWNED>();
        j = 0;
        for (auto const& f : expected_owned_faces) {
          if (std::find(tile_owned_faces.begin(), tile_owned_faces.end(), f) !=
              tile_owned_faces.end())
            j++;
          else
            std::cerr << "Tile does not contain owned face " << f << "\n";
        }
        CHECK_EQUAL(j, nent);
      }

      // Are the expected ghost faces in the tile?

      nent = tile->num_faces<Jali::Entity_type::PARALLEL_GHOST>();
      CHECK_EQUAL(expected_ghost_faces.size(), nent);

      auto const& tile_ghost_faces = tile->faces<Jali::Entity_type::PARALLEL_GHOST>();
      j = 0;
      for (auto const & f : expected_ghost_faces) {
        if (std::find(tile_ghost_faces.begin(), tile_ghost_faces.end(), f) !=
            tile_ghost_faces.end())
          j++;
        else
          std::cerr << "Tile does not contain ghost face " << f << "\n";
      }
      CHECK_EQUAL(j, nent);

      // Are the expected owned+ghost faces in the tile?

      nent = tile->num_faces();
      CHECK_EQUAL(expected_all_faces.size(), nent);

      auto const& tile_all_faces = tile->faces();
      j = 0;
      for (auto const & f : expected_all_faces) {
        if (std::find(tile_all_faces.begin(), tile_all_faces.end(), f) !=
            tile_all_faces.end())
          j++;
        else
          std::cerr << "Tile does not contain owned or ghost face " << f <<
              "\n";
      }
      CHECK_EQUAL(nent, j);





      if (edges_requested) {

        // Process the edges

        std::vector<int> expected_owned_edges;
        std::vector<int> expected_ghost_edges;
        std::vector<int> expected_all_edges;
        for (auto const & e : mesh->edges()) {
          auto ecen = mesh->edge_centroid(e);
          if (ecen[0] > min_all && ecen[0] < max_all &&
              ecen[1] > min_all && ecen[1] < max_all &&
              ecen[2] > min_all && ecen[2] < max_all) {

            expected_all_edges.push_back(e);

            if (ecen[0] > min_owned && ecen[0] < max_owned &&
                ecen[1] > min_owned && ecen[1] < max_owned &&
                ecen[2] > min_owned && ecen[2] < max_owned)
              expected_owned_edges.push_back(e);
            else
              expected_ghost_edges.push_back(e);
          }
        }


        // Are the expected owned edges in the tile? Since there is
        // only one tile in the mesh, all the edges of owned cells
        // will be owned. This list will remain the same regardless of
        // the number of halo layers requested, so check only if
        // num_halo_layers = 0

        if (num_halo_layers == 0) {
          nent = tile->num_edges<Jali::Entity_type::PARALLEL_OWNED>();
          CHECK_EQUAL(expected_owned_edges.size(), nent);

          auto const& tile_owned_edges = tile->edges<Jali::Entity_type::PARALLEL_OWNED>();
          j = 0;
          for (auto const & e : expected_owned_edges) {
            if (std::find(tile_owned_edges.begin(), tile_owned_edges.end(), e) !=
                tile_owned_edges.end())
              j++;
            else
              std::cerr << "Tile does not contain owned edge " << e << "\n";
          }
          CHECK_EQUAL(j, nent);
        }

        // Are the expected ghost nodes in the tile?

        nent = tile->num_edges<Jali::Entity_type::PARALLEL_GHOST>();
        CHECK_EQUAL(expected_ghost_edges.size(), nent);

        auto const& tile_ghost_edges = tile->edges();
        j = 0;
        for (auto const& e : expected_ghost_edges) {
          if (std::find(tile_ghost_edges.begin(), tile_ghost_edges.end(), e) !=
              tile_ghost_edges.end())
            j++;
          else
            std::cerr << "Tile does not contain ghost edge " << e << "\n";
        }
        CHECK_EQUAL(j, nent);

        // Are the expected owned+ghost edges in the tile?

        nent = tile->num_edges();
        CHECK_EQUAL(expected_all_edges.size(), nent);

        auto const& tile_all_edges = tile->edges();
        j = 0;
        for (auto const & e : expected_all_edges) {
          if (std::find(tile_all_edges.begin(), tile_all_edges.end(), e) !=
              tile_all_edges.end())
            j++;
          else
            std::cerr << "Tile does not contain owned or ghost edge " << e <<
                "\n";
        }
        CHECK_EQUAL(j, nent);
      }



      if (sides_requested) {

        std::vector<int> expected_owned_sides;
        std::vector<int> expected_ghost_sides;
        std::vector<int> expected_all_sides;
        for (auto const & s : mesh->sides()) {
          std::vector<JaliGeometry::Point> spnts;
          mesh->side_get_coordinates(s, &spnts);
          auto scen = (spnts[0] + spnts[1] + spnts[2] + spnts[3])/4.0;
          if (scen[0] > min_all && scen[0] < max_all &&
              scen[1] > min_all && scen[1] < max_all &&
              scen[2] > min_all && scen[2] < max_all) {

            expected_all_sides.push_back(s);

            if (scen[0] > min_owned && scen[0] < max_owned &&
                scen[1] > min_owned && scen[1] < max_owned &&
                scen[2] > min_owned && scen[2] < max_owned)
              expected_owned_sides.push_back(s);
            else
              expected_ghost_sides.push_back(s);
          }
        }


        // Are the expected owned sides in the tile? All the sides
        // of owned cells will be owned. This list will remain the
        // same regardless of the number of halo layers requested, so
        // check only only if num_halo_layers == 0

        if (num_halo_layers == 0) {
          nent = tile->num_sides<Jali::Entity_type::PARALLEL_OWNED>();
          CHECK_EQUAL(expected_owned_sides.size(), nent);

          auto const& tile_owned_sides = tile->sides<Jali::Entity_type::PARALLEL_OWNED>();
          j = 0;
          for (auto const& s : expected_owned_sides) {
            if (std::find(tile_owned_sides.begin(), tile_owned_sides.end(), s) !=
                tile_owned_sides.end())
              j++;
            else
              std::cerr << "Tile does not contain owned side " << s << "\n";
          }
          CHECK_EQUAL(j, nent);
        }

        // Are the expected ghost sides in the tile?

        nent = tile->num_sides<Jali::Entity_type::PARALLEL_GHOST>();
        CHECK_EQUAL(expected_ghost_sides.size(), nent);

        auto const& tile_ghost_sides = tile->sides<Jali::Entity_type::PARALLEL_GHOST>();
        j = 0;
        for (auto const& s : expected_ghost_sides) {
          if (std::find(tile_ghost_sides.begin(), tile_ghost_sides.end(), s)
              != tile_ghost_sides.end())
            j++;
          else
            std::cerr << "Tile does not contain ghost side " << s << "\n";
        }
        CHECK_EQUAL(j, nent);

        // Are the expected owned+ghost sides in the tile?

        nent = tile->num_sides();
        CHECK_EQUAL(expected_all_sides.size(), nent);

        auto const& tile_all_sides = tile->sides();
        j = 0;
        for (auto const& s : expected_all_sides) {
          if (std::find(tile_all_sides.begin(), tile_all_sides.end(), s) !=
              tile_all_sides.end())
            j++;
          else
            std::cerr << "Tile does not contain owned or ghost side " << s <<
                "\n";
        }
        CHECK_EQUAL(j, nent);
      }


      if (wedges_requested) {

        std::vector<int> expected_owned_wedges;
        std::vector<int> expected_ghost_wedges;
        std::vector<int> expected_all_wedges;
        for (auto const & w : mesh->wedges()) {
          std::vector<JaliGeometry::Point> wpnts;
          mesh->wedge_get_coordinates(w, &wpnts);
          auto wcen = (wpnts[0] + wpnts[1] + wpnts[2] + wpnts[3])/4.0;
          if (wcen[0] > min_all && wcen[0] < max_all &&
              wcen[1] > min_all && wcen[1] < max_all &&
              wcen[2] > min_all && wcen[2] < max_all) {

            expected_all_wedges.push_back(w);

            if (wcen[0] > min_owned && wcen[0] < max_owned &&
                wcen[1] > min_owned && wcen[1] < max_owned &&
                wcen[2] > min_owned && wcen[2] < max_owned)
              expected_owned_wedges.push_back(w);
            else
              expected_ghost_wedges.push_back(w);
          }
        }


        // Are the expected owned wedges in the tile? All the wedges
        // of owned cells will be owned. This list will remain the
        // same regardless of the number of halo layers requested, so
        // check only only if num_halo_layers == 0

        if (num_halo_layers == 0) {
          nent = tile->num_wedges<Jali::Entity_type::PARALLEL_OWNED>();
          CHECK_EQUAL(expected_owned_wedges.size(), nent);

          auto const& tile_owned_wedges = tile->wedges<Jali::Entity_type::PARALLEL_OWNED>();
          j = 0;
          for (auto const& w : expected_owned_wedges) {
            if (std::find(tile_owned_wedges.begin(), tile_owned_wedges.end(), w) !=
                tile_owned_wedges.end())
              j++;
            else
              std::cerr << "Tile does not contain owned wedge " << w << "\n";
          }
          CHECK_EQUAL(j, nent);
        }

        // Are the expected ghost wedges in the tile?

        nent = tile->num_wedges<Jali::Entity_type::PARALLEL_GHOST>();
        CHECK_EQUAL(expected_ghost_wedges.size(), nent);

        auto const& tile_ghost_wedges = tile->wedges<Jali::Entity_type::PARALLEL_GHOST>();
        j = 0;
        for (auto const& w : expected_ghost_wedges) {
          if (std::find(tile_ghost_wedges.begin(), tile_ghost_wedges.end(), w)
              != tile_ghost_wedges.end())
            j++;
          else
            std::cerr << "Tile does not contain ghost wedge " << w << "\n";
        }
        CHECK_EQUAL(j, nent);

        // Are the expected owned+ghost wedges in the tile?

        nent = tile->num_wedges();
        CHECK_EQUAL(expected_all_wedges.size(), nent);

        auto const& tile_all_wedges = tile->wedges();
        j = 0;
        for (auto const& w : expected_all_wedges) {
          if (std::find(tile_all_wedges.begin(), tile_all_wedges.end(), w) !=
              tile_all_wedges.end())
            j++;
          else
            std::cerr << "Tile does not contain owned or ghost wedge " << w <<
                "\n";
        }
        CHECK_EQUAL(j, nent);
      }





      if (corners_requested) {
        std::vector<int> expected_owned_corners;
        std::vector<int> expected_ghost_corners;
        std::vector<int> expected_all_corners;
        for (auto const & cn : mesh->corners()) {
          std::vector<JaliGeometry::Point> cnpnts;
          mesh->corner_get_coordinates(cn, &cnpnts);
          JaliGeometry::Point cncen =
              (cnpnts[0] + cnpnts[1] + cnpnts[2] + cnpnts[3] +
               cnpnts[4] + cnpnts[5] + cnpnts[6] + cnpnts[7])/8.0;
          if (cncen[0] > min_all && cncen[0] < max_all &&
              cncen[1] > min_all && cncen[1] < max_all &&
              cncen[2] > min_all && cncen[2] < max_all) {

            expected_all_corners.push_back(cn);

            if (cncen[0] > min_owned && cncen[0] < max_owned &&
                cncen[1] > min_owned && cncen[1] < max_owned &&
                cncen[2] > min_owned && cncen[2] < max_owned)
              expected_owned_corners.push_back(cn);
            else
              expected_ghost_corners.push_back(cn);
          }
        }


        // Are the expected owned corners in the tile? All the corners
        // of owned cells will be owned. This list will remain the
        // same regardless of the number of halo layers requested, so
        // check only if num_halo_layers = 0

        if (num_halo_layers == 0) {
          nent = tile->num_corners<Jali::Entity_type::PARALLEL_OWNED>();
          CHECK_EQUAL(expected_owned_corners.size(), nent);

          auto const& tile_owned_corners = tile->corners<Jali::Entity_type::PARALLEL_OWNED>();
          j = 0;
          for (auto const & cn : expected_owned_corners) {
            if (std::find(tile_owned_corners.begin(), tile_owned_corners.end(),
                          cn) != tile_owned_corners.end())
              j++;
            else
              std::cerr << "Tile does not contain owned corner " << cn << "\n";
          }
          CHECK_EQUAL(nent, j);
        }

        // Are the expected ghost corners in the tile?

        nent = tile->num_corners<Jali::Entity_type::PARALLEL_GHOST>();
        CHECK_EQUAL(expected_ghost_corners.size(), nent);

        auto const& tile_ghost_corners = tile->corners<Jali::Entity_type::PARALLEL_GHOST>();
        j = 0;
        for (auto const & cn : expected_ghost_corners) {
          if (std::find(tile_ghost_corners.begin(), tile_ghost_corners.end(),
                        cn) != tile_ghost_corners.end())
            j++;
          else
            std::cerr << "Tile does not contain ghost corner " << cn << "\n";
        }
        CHECK_EQUAL(j, nent);

        // Are the expected owned+ghost corners in the tile?

        nent = tile->num_corners();
        CHECK_EQUAL(expected_all_corners.size(), nent);

        auto const& tile_all_corners = tile->corners();
        int j = 0;
        for (auto const& cn : expected_all_corners) {
          if (std::find(tile_all_corners.begin(), tile_all_corners.end(), cn)
              != tile_all_corners.end())
            j++;
          else
            std::cerr << "Tile does not contain owned or ghost corner " <<
                cn << "\n";
        }
        CHECK_EQUAL(j, nent);
      }
    }  // for (num_halo_layers = 0, 2)

  }  // for (i = 0, numframeworks)

}

