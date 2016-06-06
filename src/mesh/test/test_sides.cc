//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//
/*!
 * @file   test_sides.cc
 * @author Rao V. Garimella
 * @date   Wed May 10, 2015
 *
 * @brief  Test functionality of sides (simplices forming a decomposition of the cell)
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
#include "FrameworkTraits.hh"
#include "Geometry.hh"

TEST(MESH_SIDES_2D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::Framework frameworks[] = {Jali::MSTK};
  const char *framework_names[] = {"MSTK"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::Framework);
  Jali::Framework the_framework;
  for (int fr = 0; fr < numframeworks; fr++) {

    the_framework = frameworks[fr];
    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing side operators with " << framework_names[fr] <<
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

      mesh = factory(0.0, 0.0, 1.0, 1.0, 2, 2);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);

    int nsides_owned = mesh->num_entities(Jali::Entity_kind::SIDE,
                                          Jali::Entity_type::PARALLEL_OWNED);
    int nsides_ghost = mesh->num_entities(Jali::Entity_kind::SIDE,
                                          Jali::Entity_type::PARALLEL_GHOST);
    CHECK(nsides_owned > 0);
    if (nproc > 1)
      CHECK(nsides_ghost);
    else
      CHECK(!nsides_ghost);

    nsides_owned = mesh->num_sides<Jali::Entity_type::PARALLEL_OWNED>();
    nsides_ghost = mesh->num_sides<Jali::Entity_type::PARALLEL_GHOST>();
    CHECK(nsides_owned > 0);
    if (nproc > 1)
      CHECK(nsides_ghost);
    else
      CHECK(!nsides_ghost);


    double dp;

    for (auto const & c : mesh->cells()) {
      std::vector<Jali::Entity_ID> csides;
      mesh->cell_get_sides(c, &csides);

      // Quad elements have 4 sides

      CHECK_EQUAL(csides.size(), 4);

      double cellvol = mesh->cell_volume(c);

      for (auto const & s : csides) {

        // Since the side came from the cell 'c' the cell of the
        // side must be 'c'

        Jali::Entity_ID wc = mesh->side_get_cell(s);
        CHECK_EQUAL(c, wc);

        JaliGeometry::Point ccen = mesh->cell_centroid(c);

        // Make sure the side knows which edge its associated with

        Jali::Entity_ID e = mesh->side_get_edge(s);
        CHECK(e >= 0);

        // Make sure the side knows which face its associated with

        Jali::Entity_ID f = mesh->side_get_face(s);
        CHECK(f >= 0);
        JaliGeometry::Point fcen = mesh->face_centroid(f);

        // Make sure the side knows which nodes its associated with

        Jali::Entity_ID n0 = mesh->side_get_node(s, 0);
        CHECK(n0 >= 0);
        JaliGeometry::Point npnt0;
        mesh->node_get_coordinates(n0, &npnt0);

        // Make sure the side knows which nodes its associated with

        Jali::Entity_ID n1 = mesh->side_get_node(s, 1);
        CHECK(n1 >= 0);
        JaliGeometry::Point npnt1;
        mesh->node_get_coordinates(n1, &npnt1);

        // Get the normal to facet of s (facet lies on the face)

        JaliGeometry::Point normal = mesh->side_facet_normal(s);
        double norm = JaliGeometry::norm(normal);

        // Make sure it points out of the cell by comparing with
        // the outward facing face normal

        int fdir;
        JaliGeometry::Point fnormal = mesh->face_normal(f, false, c, &fdir);
        dp = normal*fnormal;
        CHECK(dp > 0);

        // Make sure the link between the side and the opposite side
        // are correct

        Jali::Entity_ID s2 = mesh->side_get_opposite_side(s);
        CHECK(s2 >= -1);
        if (s2 != -1) {
          CHECK_EQUAL(s, mesh->side_get_opposite_side(s2));
          CHECK_EQUAL(f, mesh->side_get_face(s2));
          CHECK_EQUAL(e, mesh->side_get_edge(s2));
          CHECK_EQUAL(n0, mesh->side_get_node(s2, 1));
          CHECK_EQUAL(n1, mesh->side_get_node(s2, 0));
          CHECK(c != mesh->side_get_cell(s2));

          // Also make sure the two sides have equal and opposing
          // normals for their facet

          JaliGeometry::Point normal2 = mesh->side_facet_normal(s2);
          double norm2 = JaliGeometry::norm(normal2);
          CHECK_CLOSE(norm, norm2, 1.0e-6);
          dp = (normal*normal2)/(norm*norm2);
          CHECK_CLOSE(-1.0, dp, 1.0e-6);
        }

        // Get side coordinates and make sure they match up with the
        // expected coordinates (node, edge center, face center, cell
        // center)

        std::vector<JaliGeometry::Point> scoords;
        mesh->side_get_coordinates(s, &scoords);

        for (int i = 0; i < 2; ++i) {
          CHECK_EQUAL(scoords[0][i], npnt0[i]);
          CHECK_EQUAL(scoords[1][i], npnt1[i]);
          CHECK_EQUAL(scoords[2][i], ccen[i]);
        }

        // Since there are 4 sides in a quad element, its volume should
        // be 1/4th that of the cell

        double volume = mesh->side_volume(s);
        CHECK_CLOSE(volume, cellvol/csides.size(), 1.0e-06);

        // Now get the side coordinate in the positive volume order
        std::vector<JaliGeometry::Point> scoords2;
        mesh->side_get_coordinates(s, &scoords2, true);

        if (scoords[0][0] == scoords2[1][0] &&
            scoords[0][1] == scoords2[1][1] &&
            scoords[1][0] == scoords2[0][0] &&
            scoords[1][1] == scoords2[0][1]) {

          // coordinates are flipped - verify that the coordinate flipping
          // is warranted by computing the side volume using the natural
          // ordering and checking that its the opposite sign of the volume
          // returned by side_volume operator

          JaliGeometry::Point vec0 = scoords[1]-scoords[0];
          JaliGeometry::Point vec1 = scoords[2]-scoords[0];
          JaliGeometry::Point cpvec = vec0^vec1;
          double altvolume = cpvec[0];
          CHECK(altvolume*volume < 0.0);

          vec0 = scoords2[1]-scoords2[0];
          vec1 = scoords2[2]-scoords2[0];
          cpvec = vec0^vec1;
          altvolume = cpvec[0];
          CHECK(altvolume*volume > 0.0);
        }
      }  // for (w : csides)
    }
  }

}

TEST(MESH_SIDES_3D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::Framework frameworks[] = {Jali::MSTK};
  const char *framework_names[] = {"MSTK"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::Framework);

  Jali::Framework the_framework;
  for (int fr = 0; fr < numframeworks; fr++) {

    // Set the framework
    the_framework = frameworks[fr];
    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing side operators with " << framework_names[fr] <<
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


      mesh = factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

    } catch (const Jali::Message& e) {
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
      std::vector<Jali::Entity_ID> csides;
      mesh->cell_get_sides(c, &csides);

      // Hex elements have 24 sides

      CHECK_EQUAL(24, csides.size());

      double cellvol = mesh->cell_volume(c);

      for (auto const & s : csides) {

        // Since the side came from the cell 'c' the cell of of the
        // side must be 'c'

        Jali::Entity_ID sc = mesh->side_get_cell(s);
        CHECK_EQUAL(c, sc);

        JaliGeometry::Point ccen = mesh->cell_centroid(c);

        // Make sure the side knows which edge its associated with

        Jali::Entity_ID e = mesh->side_get_edge(s);
        CHECK(e >= 0);

        JaliGeometry::Point ecen = mesh->edge_centroid(e);

        // Make sure the side knows which face its associated with

        Jali::Entity_ID f = mesh->side_get_face(s);
        CHECK(f >= 0);
        JaliGeometry::Point fcen = mesh->face_centroid(f);

        // Make sure the side knows which nodes its associated with

        Jali::Entity_ID n0 = mesh->side_get_node(s, 0);
        CHECK(n0 >= 0);
        JaliGeometry::Point npnt0;
        mesh->node_get_coordinates(n0, &npnt0);

        Jali::Entity_ID n1 = mesh->side_get_node(s, 1);
        CHECK(n1 >= 0);
        JaliGeometry::Point npnt1;
        mesh->node_get_coordinates(n1, &npnt1);

        // Get the normal to facet of s

        JaliGeometry::Point normal = mesh->side_facet_normal(s);
        double norm = JaliGeometry::norm(normal);

        // Make sure it points out of the cell by comparing with
        // the outward facing face normal

        int fdir;
        JaliGeometry::Point fnormal = mesh->face_normal(f, false, c, &fdir);
        dp = normal*fnormal;
        CHECK(dp > 0);

        // Make sure the link between the side and the opposite side
        // is correct

        Jali::Entity_ID s2 = mesh->side_get_opposite_side(s);
        CHECK(s2 >= -1);
        if (s2 != -1) {
          CHECK_EQUAL(s, mesh->side_get_opposite_side(s2));
          CHECK_EQUAL(f, mesh->side_get_face(s2));
          CHECK_EQUAL(e, mesh->side_get_edge(s2));
          CHECK_EQUAL(n0, mesh->side_get_node(s2, 1));
          CHECK_EQUAL(n1, mesh->side_get_node(s2, 0));
          CHECK(c != mesh->side_get_cell(s2));

          // Also make sure the two sides have equal and opposing
          // normals for their facet

          JaliGeometry::Point normal2 = mesh->side_facet_normal(s2, 0);
          double norm2 = JaliGeometry::norm(normal2);
          CHECK_CLOSE(norm, norm2, 1.0e-6);
          dp = (normal*normal2)/(norm*norm2);
          CHECK_CLOSE(-1.0, dp, 1.0e-6);
        }

        // Get side coordinates in the natural order and make sure
        // they match up with the expected coordinates (node 0, node
        // 1, face center, cell center)

        std::vector<JaliGeometry::Point> scoords;
        mesh->side_get_coordinates(s, &scoords);

        for (int i = 0; i < 3; ++i) {
          CHECK_EQUAL(scoords[0][i], npnt0[i]);
          CHECK_EQUAL(scoords[1][i], npnt1[i]);
          CHECK_EQUAL(scoords[2][i], fcen[i]);
          CHECK_EQUAL(scoords[3][i], ccen[i]);
        }

        // Since there are 24 sides in a hex element, and the hex is
        // cubical, each side volume should be 1/24th that of the cell

        double volume = mesh->side_volume(s);
        CHECK_CLOSE(volume, cellvol/csides.size(), 1.0e-06);

        // Now get the side coordinate in the positive volume order
        std::vector<JaliGeometry::Point> scoords2;
        mesh->side_get_coordinates(s, &scoords2, true);

        if (scoords[0][0] == scoords2[1][0] &&
            scoords[0][1] == scoords2[1][1] &&
            scoords[0][2] == scoords2[1][2] &&
            scoords[1][0] == scoords2[0][0] &&
            scoords[1][1] == scoords2[0][1] &&
            scoords[1][2] == scoords2[0][2]) {
          // coordinates are flipped - verify that the coordinate flipping
          // is warranted by computing the side volume using the natural
          // ordering and checking that its the opposite sign of the volume
          // returned by side_volume operator

          JaliGeometry::Point vec0 = scoords[1]-scoords[0];
          JaliGeometry::Point vec1 = scoords[2]-scoords[0];
          JaliGeometry::Point vec2 = scoords[3]-scoords[0];
          JaliGeometry::Point cpvec = vec0^vec1;
          double altvolume = cpvec*vec2;
          CHECK(altvolume < 0.0);

          // Check that we get a +ve volume using the corrected ordering

          vec0 = scoords2[1]-scoords2[0];
          vec1 = scoords2[2]-scoords2[0];
          vec2 = scoords2[3]-scoords2[0];
          cpvec = vec0^vec1;
          altvolume = cpvec*vec2;
          CHECK(altvolume > 0.0);
        }
      }
    }
  }

}

