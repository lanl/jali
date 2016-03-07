/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   test_entity_iterators.cc
 * @author Rao V. Garimella
 * @date   Tue Sep 15, 2015
 *
 * @brief
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

TEST(ENTITY_ITERATORS)
{

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::Framework frameworks[] = {
    Jali::MSTK, Jali::Simple
  };
  const char *framework_names[] = {
    "MSTK", "Simple"
  };

  const int numframeworks = sizeof(frameworks)/sizeof(Jali::Framework);


  Jali::Framework the_framework;
  for (int i = 0; i < numframeworks; i++) {


    // Set the framework

    the_framework = frameworks[i];

    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing entity iterators with " << framework_names[i] << std::endl;


    // Create the mesh

    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::unique_ptr<Jali::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Jali::FrameworkPreference prefs(factory.preference());
      prefs.clear();
      prefs.push_back(the_framework);

      factory.preference(prefs);

      mesh = factory(0.0,0.0,0.0,1.0,1.0,1.0,2,2,2,NULL,true,true,true,true);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr,&aerr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    CHECK_EQUAL(aerr,0);


    int id1, id2, id3, nent;

    id1 = 0;
    for (auto const & nodeid : mesh->nodes<Jali::OWNED>()) {
      CHECK_EQUAL(id1, nodeid);
      ++id1;
    }

    // Old style call
    CHECK_EQUAL(id1, mesh->num_entities(Jali::NODE, Jali::OWNED));

    // CHECK_EQUAL macro not able to parse the templated call correctly
    CHECK_EQUAL(id1, mesh->num_nodes<Jali::OWNED>());

    id2 = 0;
    for (auto const & nodeid : mesh->nodes<Jali::GHOST>()) {
      CHECK_EQUAL(id2, nodeid);
      ++id2;
    }
    CHECK_EQUAL(id2, mesh->num_entities(Jali::NODE, Jali::GHOST));
    CHECK_EQUAL(id2, mesh->num_nodes<Jali::GHOST>());


    id3 = 0;
    for (auto const & nodeid : mesh->nodes<Jali::ALL>()) {
      CHECK_EQUAL(id3, nodeid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::NODE, Jali::ALL));
    CHECK_EQUAL(id1+id2, mesh->num_nodes<Jali::ALL>());
    CHECK_EQUAL(id3, mesh->num_nodes<Jali::ALL>());



    id1 = 0;
    for (auto const & edgeid : mesh->edges<Jali::OWNED>()) {
      CHECK_EQUAL(id1, edgeid);
      ++id1;
    }
    CHECK_EQUAL(id1, mesh->num_entities(Jali::EDGE, Jali::OWNED));
    nent = mesh->num_edges<Jali::OWNED>();
    CHECK_EQUAL(id1, nent);

    id2 = 0;
    for (auto const & edgeid : mesh->edges<Jali::GHOST>()) {
      CHECK_EQUAL(id2, edgeid);
      ++id2;
    }
    CHECK_EQUAL(id2, mesh->num_entities(Jali::EDGE, Jali::GHOST));
    CHECK_EQUAL(id2, mesh->num_edges<Jali::GHOST>());

    id3 = 0;
    for (auto const & edgeid : mesh->edges<Jali::ALL>()) {
      CHECK_EQUAL(id3, edgeid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::EDGE, Jali::ALL));
    CHECK_EQUAL(id1+id2, mesh->num_edges<Jali::ALL>());
    CHECK_EQUAL(id3, mesh->num_edges<Jali::ALL>());




    id1 = 0;
    for (auto const & faceid : mesh->faces<Jali::OWNED>()) {
      CHECK_EQUAL(id1, faceid);
      ++id1;
    }
    CHECK_EQUAL(id1, mesh->num_entities(Jali::FACE,Jali::OWNED));
    CHECK_EQUAL(id1, mesh->num_faces<Jali::OWNED>());

    id2 = 0;
    for (auto const & faceid : mesh->faces<Jali::GHOST>()) {
      CHECK_EQUAL(id2, faceid);
      ++id2;
    }
    CHECK_EQUAL(id2, mesh->num_entities(Jali::FACE, Jali::GHOST));
    CHECK_EQUAL(id2, mesh->num_faces<Jali::GHOST>());

    id3 = 0;
    for (auto const & faceid : mesh->faces<Jali::ALL>()) {
      CHECK_EQUAL(id3, faceid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::FACE, Jali::ALL));
    CHECK_EQUAL(id1+id2, mesh->num_faces<Jali::ALL>());
    CHECK_EQUAL(id3, mesh->num_faces<Jali::ALL>());



    id1 = 0;
    for (auto const & wedgeid : mesh->wedges<Jali::OWNED>()) {
      CHECK_EQUAL(id1, wedgeid);
      ++id1;
    }
    CHECK_EQUAL(id1, mesh->num_entities(Jali::WEDGE, Jali::OWNED));
    CHECK_EQUAL(id1, mesh->num_wedges<Jali::OWNED>());

    id2 = 0;
    for (auto const & wedgeid : mesh->wedges<Jali::GHOST>()) {
      CHECK_EQUAL(id2, wedgeid);
      ++id2;
    }
    CHECK_EQUAL(id2, mesh->num_entities(Jali::WEDGE, Jali::GHOST));
    CHECK_EQUAL(id2, mesh->num_wedges<Jali::GHOST>());

    id3 = 0;
    for (auto const & wedgeid : mesh->wedges<Jali::ALL>()) {
      CHECK_EQUAL(id3, wedgeid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::WEDGE, Jali::ALL));
    CHECK_EQUAL(id1+id2, mesh->num_wedges<Jali::ALL>());
    CHECK_EQUAL(id3, mesh->num_wedges<Jali::ALL>());




    id1 = 0;
    for (auto const & cornerid : mesh->corners<Jali::OWNED>()) {
      CHECK_EQUAL(id1, cornerid);
      ++id1;
    }
    CHECK_EQUAL(id1, mesh->num_entities(Jali::CORNER, Jali::OWNED));
    CHECK_EQUAL(id1, mesh->num_corners<Jali::OWNED>());

    id2 = 0;
    for (auto const & cornerid : mesh->corners<Jali::GHOST>()) {
      CHECK_EQUAL(id2, cornerid);
      ++id2;
    }
    CHECK_EQUAL(id2, mesh->num_entities(Jali::CORNER, Jali::GHOST));
    CHECK_EQUAL(id2, mesh->num_corners<Jali::GHOST>());

    id3 = 0;
    for (auto const & cornerid : mesh->corners<Jali::ALL>()) {
      CHECK_EQUAL(id3, cornerid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::CORNER, Jali::ALL));
    CHECK_EQUAL(id1 + id2, mesh->num_corners<Jali::ALL>());
    CHECK_EQUAL(id3, mesh->num_corners<Jali::ALL>());



    id1 = 0;
    for (auto const & cellid : mesh->cells<Jali::OWNED>()) {
      CHECK_EQUAL(id1, cellid);
      ++id1;
    }
    CHECK_EQUAL(id1, mesh->num_entities(Jali::CELL, Jali::OWNED));
    CHECK_EQUAL(id1, mesh->num_cells<Jali::OWNED>());

    id2 = 0;
    for (auto const & cellid : mesh->cells<Jali::GHOST>()) {
      CHECK_EQUAL(id2, cellid);
      ++id2;
    }
    CHECK_EQUAL(id2, mesh->num_entities(Jali::CELL, Jali::GHOST));
    CHECK_EQUAL(id2, mesh->num_cells<Jali::GHOST>());

    id3 = 0;
    for (auto const & cellid : mesh->cells<Jali::ALL>()) {
      CHECK_EQUAL(id3, cellid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::CELL, Jali::ALL));
    CHECK_EQUAL(id1 + id2, mesh->num_cells<Jali::ALL>());
    CHECK_EQUAL(id3, mesh->num_cells<Jali::ALL>());

  } // for each framework i

}

