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
#include "errors.hh"

TEST(ENTITY_ITERATORS) {
  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  const Jali::MeshFramework_t frameworks[] = {Jali::MSTK, Jali::Simple};
  const char *framework_names[] = {"MSTK", "Simple"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::MeshFramework_t);
  Jali::MeshFramework_t the_framework;
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
      factory.framework(the_framework);

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

    } catch (const Errors::Message& e) {
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
    for (auto const & nodeid : mesh->nodes<Jali::Entity_type::ALL>()) {
      CHECK_EQUAL(id3, nodeid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::NODE,
                                              Jali::Entity_type::ALL));
    CHECK_EQUAL(id1+id2, mesh->num_nodes<Jali::Entity_type::ALL>());
    CHECK_EQUAL(id3, mesh->num_nodes<Jali::Entity_type::ALL>());



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
      for (auto const & edgeid : mesh->edges<Jali::Entity_type::ALL>()) {
        CHECK_EQUAL(id3, edgeid);
        ++id3;
      }
      CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::EDGE,
                                                Jali::Entity_type::ALL));
      CHECK_EQUAL(id1+id2, mesh->num_edges<Jali::Entity_type::ALL>());
      CHECK_EQUAL(id3, mesh->num_edges<Jali::Entity_type::ALL>());
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
    for (auto const & faceid : mesh->faces<Jali::Entity_type::ALL>()) {
      CHECK_EQUAL(id3, faceid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::FACE,
                                              Jali::Entity_type::ALL));
    CHECK_EQUAL(id1+id2, mesh->num_faces<Jali::Entity_type::ALL>());
    CHECK_EQUAL(id3, mesh->num_faces<Jali::Entity_type::ALL>());



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
      for (auto const & sideid : mesh->sides<Jali::Entity_type::ALL>()) {
        CHECK_EQUAL(id3, sideid);
        ++id3;
      }
      CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::SIDE,
                                                Jali::Entity_type::ALL));
      CHECK_EQUAL(id1+id2, mesh->num_sides<Jali::Entity_type::ALL>());
      CHECK_EQUAL(id3, mesh->num_sides<Jali::Entity_type::ALL>());

      
      
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
      for (auto const & wedgeid : mesh->wedges<Jali::Entity_type::ALL>()) {
        CHECK_EQUAL(id3, wedgeid);
        ++id3;
      }
      CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::WEDGE,
                                                Jali::Entity_type::ALL));
      CHECK_EQUAL(id1+id2, mesh->num_wedges<Jali::Entity_type::ALL>());
      CHECK_EQUAL(id3, mesh->num_wedges<Jali::Entity_type::ALL>());
      
      

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
      for (auto const & cornerid : mesh->corners<Jali::Entity_type::ALL>()) {
        CHECK_EQUAL(id3, cornerid);
        ++id3;
      }
      CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::CORNER,
                                                Jali::Entity_type::ALL));
      CHECK_EQUAL(id1 + id2, mesh->num_corners<Jali::Entity_type::ALL>());
      CHECK_EQUAL(id3, mesh->num_corners<Jali::Entity_type::ALL>());
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
    for (auto const & cellid : mesh->cells<Jali::Entity_type::ALL>()) {
      CHECK_EQUAL(id3, cellid);
      ++id3;
    }
    CHECK_EQUAL(id1 + id2, mesh->num_entities(Jali::Entity_kind::CELL,
                                              Jali::Entity_type::ALL));
    CHECK_EQUAL(id1 + id2, mesh->num_cells<Jali::Entity_type::ALL>());
    CHECK_EQUAL(id3, mesh->num_cells<Jali::Entity_type::ALL>());
  }  // for each framework i

}

