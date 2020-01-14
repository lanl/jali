/*
 Copyright (c) 2019, Triad National Security, LLC
 All rights reserved.

 Copyright 2019. Triad National Security, LLC. This software was
 produced under U.S. Government contract 89233218CNA000001 for Los
 Alamos National Laboratory (LANL), which is operated by Triad
 National Security, LLC for the U.S. Department of Energy. 
 All rights in the program are reserved by Triad National Security,
 LLC, and the U.S. Department of Energy/National Nuclear Security
 Administration. The Government is granted for itself and others acting
 on its behalf a nonexclusive, paid-up, irrevocable worldwide license
 in this material to reproduce, prepare derivative works, distribute
 copies to the public, perform publicly and display publicly, and to
 permit others to do so

 
 This is open source software distributed under the 3-clause BSD license.
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:
 
 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
 3. Neither the name of Triad National Security, LLC, Los Alamos
    National Laboratory, LANL, the U.S. Government, nor the names of its
    contributors may be used to endorse or promote products derived from this
    software without specific prior written permission.

 
 THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND
 CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*!
 * @file   test_wedges.cc
 * @author Rao V. Garimella
 * @date   Mon Oct 5, 2015
 *
 * @brief  Test functionality of wedges (simplices forming a decomposition of the cell)
 *
 *
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Geometry.hh"

TEST(MESH_WEDGES_2D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::MeshFramework_t frameworks[] = {Jali::MSTK};
  const char *framework_names[] = {"MSTK"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::MeshFramework_t);
  Jali::MeshFramework_t the_framework;
  for (int fr = 0; fr < numframeworks; fr++) {

    the_framework = frameworks[fr];
    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing wedge operators with " << framework_names[fr] <<
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
                                                   Jali::Entity_kind::WEDGE};
      factory.included_entities(entitylist);

      mesh = factory(0.0, 0.0, 1.0, 1.0, 2, 2);

    } catch (const Errors::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);

    int nwedges_owned = mesh->num_entities(Jali::Entity_kind::WEDGE,
                                           Jali::Entity_type::PARALLEL_OWNED);
    int nwedges_ghost = mesh->num_entities(Jali::Entity_kind::WEDGE,
                                           Jali::Entity_type::PARALLEL_GHOST);
    CHECK(nwedges_owned > 0);
    if (nproc > 1)
      CHECK(nwedges_ghost);
    else
      CHECK(!nwedges_ghost);

    nwedges_owned = mesh->num_wedges<Jali::Entity_type::PARALLEL_OWNED>();
    nwedges_ghost = mesh->num_wedges<Jali::Entity_type::PARALLEL_GHOST>();
    CHECK(nwedges_owned > 0);
    if (nproc > 1)
      CHECK(nwedges_ghost);
    else
      CHECK(!nwedges_ghost);

    double dp;

    for (auto const & c : mesh->cells()) {
      std::vector<Jali::Entity_ID> cwedges;
      mesh->cell_get_wedges(c, &cwedges);

      // Quad elements have 8 wedges

      CHECK_EQUAL(cwedges.size(), 8);

      double cellvol = mesh->cell_volume(c);

      for (auto const & w : cwedges) {

        // Since the wedge came from the cell 'c' the cell of of the
        // wedge must be 'c'

        Jali::Entity_ID wc = mesh->wedge_get_cell(w);
        CHECK_EQUAL(c, wc);

        JaliGeometry::Point ccen = mesh->cell_centroid(c);

        // Make sure the wedge knows which edge its associated with

        Jali::Entity_ID edg = mesh->wedge_get_edge(w);
        CHECK(edg >= 0);

        // Make sure the wedge knows which face its associated with

        Jali::Entity_ID f = mesh->wedge_get_face(w);
        CHECK(f >= 0);
        JaliGeometry::Point fcen = mesh->face_centroid(f);

        // Make sure the wedge knows which node its associated with

        Jali::Entity_ID n = mesh->wedge_get_node(w);
        CHECK(n >= 0);
        JaliGeometry::Point npnt;
        mesh->node_get_coordinates(n, &npnt);

        // Get the normal to facet0 of w (facet0 lies on the face)

        JaliGeometry::Point normal0 = mesh->wedge_facet_normal(w, 0);
        double norm0 = JaliGeometry::norm(normal0);

        // Make sure it points out of the cell by comparing with
        // the outward facing face normal

        int fdir;
        JaliGeometry::Point fnormal = mesh->face_normal(f, false, c, &fdir);
        dp = normal0*fnormal;
        CHECK(dp > 0);

        // Get the normal to facet1 of w (facet1 is perpendicular to the edge)

        JaliGeometry::Point normal1 = mesh->wedge_facet_normal(w, 1);
        double norm1 = JaliGeometry::norm(normal1);

        // Make sure that it points away from the node of the wedge by
        // comparing with the edge vector going from node n to the
        // opposite node

        int edir;
        JaliGeometry::Point evec = mesh->edge_vector(f, false, n, &edir);
        dp = normal1*evec;
        CHECK(dp > 0);

        // Make sure the link between the wedge and the opposite wedge
        // are correct

        Jali::Entity_ID w2 = mesh->wedge_get_opposite_wedge(w);
        CHECK(w2 >= -1);
        if (w2 != -1) {
          CHECK_EQUAL(w, mesh->wedge_get_opposite_wedge(w2));
          CHECK_EQUAL(f, mesh->wedge_get_face(w2));
          CHECK_EQUAL(edg, mesh->wedge_get_edge(w2));
          CHECK_EQUAL(n, mesh->wedge_get_node(w2));
          CHECK(c != mesh->wedge_get_cell(w2));

          // Also make sure the two wedges have equal and opposing
          // normals for their facet 0

          JaliGeometry::Point normal20 = mesh->wedge_facet_normal(w2, 0);
          double norm20 = JaliGeometry::norm(normal20);
          CHECK_CLOSE(norm0, norm20, 1.0e-6);
          dp = (normal0*normal20)/(norm0*norm20);
          CHECK_CLOSE(-1.0, dp, 1.0e-6);

          // Also make sure the two wedges have equal and coincident
          // normals to facet 1

          JaliGeometry::Point normal21 = mesh->wedge_facet_normal(w2, 1);
          double norm21 = JaliGeometry::norm(normal21);
          CHECK_CLOSE(norm1, norm21, 1.0e-6);
          dp = (normal1*normal21)/(norm1*norm21);
          CHECK_CLOSE(1.0, dp, 1.0e-6);
        }

        // Make sure the link between the edge and the adjacent wedge
        // are correct

        Jali::Entity_ID w3 = mesh->wedge_get_adjacent_wedge(w);
        CHECK_EQUAL(w, mesh->wedge_get_adjacent_wedge(w3));
        CHECK_EQUAL(f, mesh->wedge_get_face(w3));
        CHECK_EQUAL(edg, mesh->wedge_get_edge(w3));
        CHECK_EQUAL(c, mesh->wedge_get_cell(w3));

        // Also make sure the two wedges have equal and coincident
        // normals for their facet 0

        JaliGeometry::Point normal30 = mesh->wedge_facet_normal(w3, 0);
        double norm30 = JaliGeometry::norm(normal30);
        CHECK_CLOSE(norm0, norm30, 1.0e-6);
        dp = (normal0*normal30)/(norm0*norm30);
        CHECK_CLOSE(1.0, dp, 1.0e-6);

        // Also make sure the two wedges have equal and opposite
        // normals to facet 1

        JaliGeometry::Point normal31 = mesh->wedge_facet_normal(w3, 1);
        double norm31 = JaliGeometry::norm(normal31);
        CHECK_CLOSE(norm1, norm31, 1.0e-6);
        dp = (normal1*normal31)/(norm1*norm31);
        CHECK_CLOSE(-1.0, dp, 1.0e-6);


        // Get the opposite wedge (w4) of the adjacent wedge (w3) of w
        // and make sure that it is adjacent to the opposite wedge
        // (w2) of w.

        Jali::Entity_ID w4 = mesh->wedge_get_opposite_wedge(w3);
        CHECK(w4 >= -1);
        if (w4 != -1 && w2 != -1)
          CHECK_EQUAL(w4, mesh->wedge_get_adjacent_wedge(w2));

        // Get wedge coordinates and make sure they match up with the
        // expected coordinates (node, edge center, face center, cell
        // center)

        std::vector<JaliGeometry::Point> wcoords;
        mesh->wedge_get_coordinates(w, &wcoords);

        for (int i = 0; i < 2; ++i) {
          CHECK_EQUAL(wcoords[0][i], npnt[i]);
          CHECK_EQUAL(wcoords[1][i], fcen[i]);
          CHECK_EQUAL(wcoords[2][i], ccen[i]);
        }

        // Since there are 8 wedges in a quad element, its volume should
        // be 1/8th that of the cell

        double volume = mesh->wedge_volume(w);
        CHECK_CLOSE(volume, cellvol/cwedges.size(), 1.0e-06);

        // Now get the wedge coordinate in the positive volume order
        std::vector<JaliGeometry::Point> wcoords2;
        mesh->wedge_get_coordinates(w, &wcoords2, true);

        if (wcoords[1][0] == wcoords2[2][0] &&
            wcoords[1][1] == wcoords2[2][1] &&
            wcoords[2][0] == wcoords2[1][0] &&
            wcoords[2][1] == wcoords2[1][1]) {

          // coordinates are flipped - verify that the coordinate flipping
          // is warranted by computing the wedge volume using the natural
          // ordering and checking that its the opposite sign of the volume
          // returned by wedge_volume operator

          JaliGeometry::Point vec0 = wcoords[1]-wcoords[0];
          JaliGeometry::Point vec1 = wcoords[2]-wcoords[0];
          JaliGeometry::Point cpvec = vec0^vec1;
          double altvolume = cpvec[0];
          CHECK(altvolume*volume < 0.0);

          vec0 = wcoords2[1]-wcoords2[0];
          vec1 = wcoords2[2]-wcoords2[0];
          cpvec = vec0^vec1;
          altvolume = cpvec[0];
          CHECK(altvolume*volume > 0.0);
        }
      }  // for (w : cwedges)
    }
  }

}

TEST(MESH_WEDGES_3D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::MeshFramework_t frameworks[] = {Jali::MSTK};
  const char *framework_names[] = {"MSTK"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::MeshFramework_t);

  Jali::MeshFramework_t the_framework;
  for (int fr = 0; fr < numframeworks; fr++) {

    // Set the framework
    the_framework = frameworks[fr];
    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing wedge operators with " << framework_names[fr] <<
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
                                                   Jali::Entity_kind::WEDGE};
      factory.included_entities(entitylist);


      mesh = factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

    } catch (const Errors::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);

    double dp;
    for (auto const & c : mesh->cells()) {
      std::vector<Jali::Entity_ID> cwedges;
      mesh->cell_get_wedges(c, &cwedges);

      // Hex elements have 48 wedges

      CHECK_EQUAL(cwedges.size(), 48);

      double cellvol = mesh->cell_volume(c);

      for (auto const & w : cwedges) {

        // Since the wedge came from the cell 'c' the cell of of the
        // wedge must be 'c'

        Jali::Entity_ID wc = mesh->wedge_get_cell(w);
        CHECK_EQUAL(c, wc);

        JaliGeometry::Point ccen = mesh->cell_centroid(c);

        // Make sure the wedge knows which edge its associated with

        Jali::Entity_ID edg = mesh->wedge_get_edge(w);
        CHECK(e >= 0);

        JaliGeometry::Point ecen = mesh->edge_centroid(e);

        // Make sure the wedge knows which face its associated with

        Jali::Entity_ID f = mesh->wedge_get_face(w);
        CHECK(f >= 0);
        JaliGeometry::Point fcen = mesh->face_centroid(f);

        // Make sure the wedge knows which node its associated with

        Jali::Entity_ID n = mesh->wedge_get_node(w);
        CHECK(n >= 0);
        JaliGeometry::Point npnt;
        mesh->node_get_coordinates(n, &npnt);

        // Get the normal to facet0 of w (facet0 lies on the face)

        JaliGeometry::Point normal0 = mesh->wedge_facet_normal(w, 0);
        double norm0 = JaliGeometry::norm(normal0);

        // Make sure it points out of the cell by comparing with
        // the outward facing face normal

        int fdir;
        JaliGeometry::Point fnormal = mesh->face_normal(f, false, c, &fdir);
        dp = normal0*fnormal;
        CHECK(dp > 0);

        // Get the normal to facet1 of w (facet1 is perpendicular to the edge)

        JaliGeometry::Point normal1 = mesh->wedge_facet_normal(w, 1);
        double norm1 = JaliGeometry::norm(normal1);

        // Make sure that it points away from the node of the wedge by
        // comparing with the edge vector going from node n to the
        // opposite node

        int edir;
        JaliGeometry::Point evec = mesh->edge_vector(edg, false, n, &edir);
        dp = normal1*evec;
        CHECK(dp > 0);

        // Make sure the link between the wedge and the opposite wedge
        // are correct

        Jali::Entity_ID w2 = mesh->wedge_get_opposite_wedge(w);
        CHECK(w2 >= -1);
        if (w2 != -1) {
          CHECK_EQUAL(w, mesh->wedge_get_opposite_wedge(w2));
          CHECK_EQUAL(f, mesh->wedge_get_face(w2));
          CHECK_EQUAL(edg, mesh->wedge_get_edge(w2));
          CHECK_EQUAL(n, mesh->wedge_get_node(w2));
          CHECK(c != mesh->wedge_get_cell(w2));

          // Also make sure the two wedges have equal and opposing
          // normals for their facet 0

          JaliGeometry::Point normal20 = mesh->wedge_facet_normal(w2, 0);
          double norm20 = JaliGeometry::norm(normal20);
          CHECK_CLOSE(norm0, norm20, 1.0e-6);
          dp = (normal0*normal20)/(norm0*norm20);
          CHECK_CLOSE(-1.0, dp, 1.0e-6);

          // Also make sure the two wedges have equal and coincident
          // normals to facet 1

          JaliGeometry::Point normal21 = mesh->wedge_facet_normal(w2, 1);
          double norm21 = JaliGeometry::norm(normal21);
          CHECK_CLOSE(norm1, norm21, 1.0e-6);
          dp = (normal1*normal21)/(norm1*norm21);
          CHECK_CLOSE(1.0, dp, 1.0e-6);
        }

        // Make sure the link between the edge and the adjacent wedge
        // are correct

        Jali::Entity_ID w3 = mesh->wedge_get_adjacent_wedge(w);
        CHECK_EQUAL(w, mesh->wedge_get_adjacent_wedge(w3));
        CHECK_EQUAL(f, mesh->wedge_get_face(w3));
        CHECK_EQUAL(edg, mesh->wedge_get_edge(w3));
        CHECK_EQUAL(c, mesh->wedge_get_cell(w3));

        // Also make sure the two wedges have equal and coincident
        // normals for their facet 0

        JaliGeometry::Point normal30 = mesh->wedge_facet_normal(w3, 0);
        double norm30 = JaliGeometry::norm(normal30);
        CHECK_CLOSE(norm0, norm30, 1.0e-6);
        dp = (normal0*normal30)/(norm0*norm30);
        CHECK_CLOSE(1.0, dp, 1.0e-6);

        // Also make sure the two wedges have equal and opposite
        // normals to facet 1

        JaliGeometry::Point normal31 = mesh->wedge_facet_normal(w3, 1);
        double norm31 = JaliGeometry::norm(normal31);
        CHECK_CLOSE(norm1, norm31, 1.0e-6);
        dp = (normal1*normal31)/(norm1*norm31);
        CHECK_CLOSE(-1.0, dp, 1.0e-6);


        // Get the opposite wedge (w4) of the adjacent wedge (w3) of w
        // and make sure that it is adjacent to the opposite wedge
        // (w2) of w.

        Jali::Entity_ID w4 = mesh->wedge_get_opposite_wedge(w3);
        CHECK(w4 >= -1);
        if (w4 != -1 && w2 != -1)
          CHECK_EQUAL(w4, mesh->wedge_get_adjacent_wedge(w2));

        // Get wedge coordinates in the natural order and make sure
        // they match up with the expected coordinates (node, edge
        // center, face center, cell center)

        std::vector<JaliGeometry::Point> wcoords;
        mesh->wedge_get_coordinates(w, &wcoords);

        for (int i = 0; i < 3; ++i) {
          CHECK_EQUAL(wcoords[0][i], npnt[i]);
          CHECK_EQUAL(wcoords[1][i], ecen[i]);
          CHECK_EQUAL(wcoords[2][i], fcen[i]);
          CHECK_EQUAL(wcoords[3][i], ccen[i]);
        }

        // Since there are 48 wedges in a hex element, its volume should
        // be 1/48th that of the cell

        double volume = mesh->wedge_volume(w);
        CHECK_CLOSE(volume, cellvol/cwedges.size(), 1.0e-06);

        // Now get the wedge coordinate in the positive volume order
        std::vector<JaliGeometry::Point> wcoords2;
        mesh->wedge_get_coordinates(w, &wcoords2, true);

        if (wcoords[1][0] == wcoords2[2][0] &&
            wcoords[1][1] == wcoords2[2][1] &&
            wcoords[1][2] == wcoords2[2][2] &&
            wcoords[2][0] == wcoords2[1][0] &&
            wcoords[2][1] == wcoords2[1][1] &&
            wcoords[2][2] == wcoords2[1][2]) {
          // coordinates are flipped - verify that the coordinate flipping
          // is warranted by computing the wedge volume using the natural
          // ordering and checking that its the opposite sign of the volume
          // returned by wedge_volume operator

          JaliGeometry::Point vec0 = wcoords[1]-wcoords[0];
          JaliGeometry::Point vec1 = wcoords[2]-wcoords[0];
          JaliGeometry::Point vec2 = wcoords[3]-wcoords[0];
          JaliGeometry::Point cpvec = vec0^vec1;
          double altvolume = cpvec*vec2;
          CHECK(altvolume < 0.0);

          // Check that we get a +ve volume using the corrected ordering

          vec0 = wcoords2[1]-wcoords2[0];
          vec1 = wcoords2[2]-wcoords2[0];
          vec2 = wcoords2[3]-wcoords2[0];
          cpvec = vec0^vec1;
          altvolume = cpvec*vec2;
          CHECK(altvolume > 0.0);
        }
      }
    }
  }

}

