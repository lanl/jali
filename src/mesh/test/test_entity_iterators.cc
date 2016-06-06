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

TEST(ENTITY_ITERATORS) {
  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::Framework frameworks[] = {Jali::MSTK, Jali::Simple};
  const char *framework_names[] = {"MSTK", "Simple"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::Framework);
  Jali::Framework the_framework;
  for (int i = 0; i < numframeworks; i++) {
    the_framework = frameworks[i];
    if (!Jali::framework_available(the_framework)) continue;
    std::cerr << "Testing entity iterators with " << framework_names[i] <<
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

      if (the_framework == Jali::MSTK) {
        std::vector<Jali::Entity_kind> entitylist = {Jali::Entity_kind::EDGE,
                                                     Jali::Entity_kind::FACE,
                                                     Jali::Entity_kind::SIDE,
                                                     Jali::Entity_kind::WEDGE,
                                                     Jali::Entity_kind::CORNER};
        factory.included_entities(entitylist);
      } else {
        std::vector<Jali::Entity_kind> entitylist = {Jali::Entity_kind::FACE};
        factory.included_entities(entitylist);
      }

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


    int id1, id2, id3, nent;

    id1 = 0;
    for (auto const & nodeid : mesh->nodes<Jali::Entity_type::PARALLEL_OWNED>()) {
      CHECK_EQUAL(id1, nodeid);
      ++id1;
    }

    // Old style call
    CHECK_EQUAL(id1, mesh->num_entities(Jali::Entity_kind::NODE,
                                        Jali::Entity_type::PARALLEL_OWNED));

    // CHECK_EQUAL macro not able to parse the templated call correctly
    CHECK_EQUAL(id1, mesh->num_nodes<Jali::Entity_type::PARALLEL_OWNED>());

    id2 = 0;
    for (auto const & nodeid : mesh->nodes<Jali::Entity_type::PARALLEL_GHOST>()) {
      CHECK_EQUAL(id2, nodeid);
      ++id2;
    }
    CHECK_EQUAL(id2, mesh->num_entities(Jali::Entity_kind::NODE,
                                        Jali::Entity_type::PARALLEL_GHOST));
    CHECK_EQUAL(id2, mesh->num_nodes<Jali::Entity_type::PARALLEL_GHOST>());


    id3 = 0;
    for (auto const & nodeid : mesh->nodes<Jali::Entity_type::PARALLEL_ALL>()) {
      CHECK_EQUAL(id3, nodeid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::NODE,
                                              Jali::Entity_type::PARALLEL_ALL));
    CHECK_EQUAL(id1+id2, mesh->num_nodes<Jali::Entity_type::PARALLEL_ALL>());
    CHECK_EQUAL(id3, mesh->num_nodes<Jali::Entity_type::PARALLEL_ALL>());



    if (the_framework == Jali::MSTK) {
      id1 = 0;
      for (auto const & edgeid : mesh->edges<Jali::Entity_type::PARALLEL_OWNED>()) {
        CHECK_EQUAL(id1, edgeid);
        ++id1;
      }
      CHECK_EQUAL(id1, mesh->num_entities(Jali::Entity_kind::EDGE,
                                          Jali::Entity_type::PARALLEL_OWNED));
      nent = mesh->num_edges<Jali::Entity_type::PARALLEL_OWNED>();
      CHECK_EQUAL(id1, nent);
      
      id2 = 0;
      for (auto const & edgeid : mesh->edges<Jali::Entity_type::PARALLEL_GHOST>()) {
        CHECK_EQUAL(id2, edgeid);
        ++id2;
      }
      CHECK_EQUAL(id2, mesh->num_entities(Jali::Entity_kind::EDGE,
                                          Jali::Entity_type::PARALLEL_GHOST));
      CHECK_EQUAL(id2, mesh->num_edges<Jali::Entity_type::PARALLEL_GHOST>());
      
      id3 = 0;
      for (auto const & edgeid : mesh->edges<Jali::Entity_type::PARALLEL_ALL>()) {
        CHECK_EQUAL(id3, edgeid);
        ++id3;
      }
      CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::EDGE,
                                                Jali::Entity_type::PARALLEL_ALL));
      CHECK_EQUAL(id1+id2, mesh->num_edges<Jali::Entity_type::PARALLEL_ALL>());
      CHECK_EQUAL(id3, mesh->num_edges<Jali::Entity_type::PARALLEL_ALL>());
    }



    id1 = 0;
    for (auto const & faceid : mesh->faces<Jali::Entity_type::PARALLEL_OWNED>()) {
      CHECK_EQUAL(id1, faceid);
      ++id1;
    }
    CHECK_EQUAL(id1, mesh->num_entities(Jali::Entity_kind::FACE,
                                        Jali::Entity_type::PARALLEL_OWNED));
    CHECK_EQUAL(id1, mesh->num_faces<Jali::Entity_type::PARALLEL_OWNED>());

    id2 = 0;
    for (auto const & faceid : mesh->faces<Jali::Entity_type::PARALLEL_GHOST>()) {
      CHECK_EQUAL(id2, faceid);
      ++id2;
    }
    CHECK_EQUAL(id2, mesh->num_entities(Jali::Entity_kind::FACE,
                                        Jali::Entity_type::PARALLEL_GHOST));
    CHECK_EQUAL(id2, mesh->num_faces<Jali::Entity_type::PARALLEL_GHOST>());

    id3 = 0;
    for (auto const & faceid : mesh->faces<Jali::Entity_type::PARALLEL_ALL>()) {
      CHECK_EQUAL(id3, faceid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::FACE,
                                              Jali::Entity_type::PARALLEL_ALL));
    CHECK_EQUAL(id1+id2, mesh->num_faces<Jali::Entity_type::PARALLEL_ALL>());
    CHECK_EQUAL(id3, mesh->num_faces<Jali::Entity_type::PARALLEL_ALL>());



    if (the_framework == Jali::MSTK) {
      id1 = 0;
      for (auto const & sideid : mesh->sides<Jali::Entity_type::PARALLEL_OWNED>()) {
        CHECK_EQUAL(id1, sideid);
        ++id1;
      }
      CHECK_EQUAL(id1, mesh->num_entities(Jali::Entity_kind::SIDE,
                                          Jali::Entity_type::PARALLEL_OWNED));
      CHECK_EQUAL(id1, mesh->num_sides<Jali::Entity_type::PARALLEL_OWNED>());
      
      id2 = 0;
      for (auto const & sideid : mesh->sides<Jali::Entity_type::PARALLEL_GHOST>()) {
        CHECK_EQUAL(id2, sideid);
        ++id2;
      }
      CHECK_EQUAL(id2, mesh->num_entities(Jali::Entity_kind::SIDE,
                                          Jali::Entity_type::PARALLEL_GHOST));
      CHECK_EQUAL(id2, mesh->num_sides<Jali::Entity_type::PARALLEL_GHOST>());
      
      id3 = 0;
      for (auto const & sideid : mesh->sides<Jali::Entity_type::PARALLEL_ALL>()) {
        CHECK_EQUAL(id3, sideid);
        ++id3;
      }
      CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::SIDE,
                                                Jali::Entity_type::PARALLEL_ALL));
      CHECK_EQUAL(id1+id2, mesh->num_sides<Jali::Entity_type::PARALLEL_ALL>());
      CHECK_EQUAL(id3, mesh->num_sides<Jali::Entity_type::PARALLEL_ALL>());

      
      
      id1 = 0;
      for (auto const & wedgeid : mesh->wedges<Jali::Entity_type::PARALLEL_OWNED>()) {
        CHECK_EQUAL(id1, wedgeid);
      ++id1;
      }
      CHECK_EQUAL(id1, mesh->num_entities(Jali::Entity_kind::WEDGE,
                                          Jali::Entity_type::PARALLEL_OWNED));
      CHECK_EQUAL(id1, mesh->num_wedges<Jali::Entity_type::PARALLEL_OWNED>());
      
      id2 = 0;
      for (auto const & wedgeid : mesh->wedges<Jali::Entity_type::PARALLEL_GHOST>()) {
        CHECK_EQUAL(id2, wedgeid);
        ++id2;
      }
      CHECK_EQUAL(id2, mesh->num_entities(Jali::Entity_kind::WEDGE,
                                          Jali::Entity_type::PARALLEL_GHOST));
      CHECK_EQUAL(id2, mesh->num_wedges<Jali::Entity_type::PARALLEL_GHOST>());
      
      id3 = 0;
      for (auto const & wedgeid : mesh->wedges<Jali::Entity_type::PARALLEL_ALL>()) {
        CHECK_EQUAL(id3, wedgeid);
        ++id3;
      }
      CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::WEDGE,
                                                Jali::Entity_type::PARALLEL_ALL));
      CHECK_EQUAL(id1+id2, mesh->num_wedges<Jali::Entity_type::PARALLEL_ALL>());
      CHECK_EQUAL(id3, mesh->num_wedges<Jali::Entity_type::PARALLEL_ALL>());
      
      

      id1 = 0;
      for (auto const & cornerid : mesh->corners<Jali::Entity_type::PARALLEL_OWNED>()) {
        CHECK_EQUAL(id1, cornerid);
        ++id1;
      }
      CHECK_EQUAL(id1, mesh->num_entities(Jali::Entity_kind::CORNER,
                                          Jali::Entity_type::PARALLEL_OWNED));
      CHECK_EQUAL(id1, mesh->num_corners<Jali::Entity_type::PARALLEL_OWNED>());
      
      id2 = 0;
      for (auto const & cornerid : mesh->corners<Jali::Entity_type::PARALLEL_GHOST>()) {
        CHECK_EQUAL(id2, cornerid);
        ++id2;
      }
      CHECK_EQUAL(id2, mesh->num_entities(Jali::Entity_kind::CORNER,
                                          Jali::Entity_type::PARALLEL_GHOST));
      CHECK_EQUAL(id2, mesh->num_corners<Jali::Entity_type::PARALLEL_GHOST>());
      
      id3 = 0;
      for (auto const & cornerid : mesh->corners<Jali::Entity_type::PARALLEL_ALL>()) {
        CHECK_EQUAL(id3, cornerid);
        ++id3;
      }
      CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::CORNER,
                                                Jali::Entity_type::PARALLEL_ALL));
      CHECK_EQUAL(id1 + id2, mesh->num_corners<Jali::Entity_type::PARALLEL_ALL>());
      CHECK_EQUAL(id3, mesh->num_corners<Jali::Entity_type::PARALLEL_ALL>());
    }


    id1 = 0;
    for (auto const & cellid : mesh->cells<Jali::Entity_type::PARALLEL_OWNED>()) {
      CHECK_EQUAL(id1, cellid);
      ++id1;
    }
    CHECK_EQUAL(id1, mesh->num_entities(Jali::Entity_kind::CELL,
                                        Jali::Entity_type::PARALLEL_OWNED));
    CHECK_EQUAL(id1, mesh->num_cells<Jali::Entity_type::PARALLEL_OWNED>());

    id2 = 0;
    for (auto const & cellid : mesh->cells<Jali::Entity_type::PARALLEL_GHOST>()) {
      CHECK_EQUAL(id2, cellid);
      ++id2;
    }
    CHECK_EQUAL(id2, mesh->num_entities(Jali::Entity_kind::CELL,
                                        Jali::Entity_type::PARALLEL_GHOST));
    CHECK_EQUAL(id2, mesh->num_cells<Jali::Entity_type::PARALLEL_GHOST>());

    id3 = 0;
    for (auto const & cellid : mesh->cells<Jali::Entity_type::PARALLEL_ALL>()) {
      CHECK_EQUAL(id3, cellid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::CELL,
                                              Jali::Entity_type::PARALLEL_ALL));
    CHECK_EQUAL(id1 + id2, mesh->num_cells<Jali::Entity_type::PARALLEL_ALL>());
    CHECK_EQUAL(id3, mesh->num_cells<Jali::Entity_type::PARALLEL_ALL>());
  }  // for each framework i

}

