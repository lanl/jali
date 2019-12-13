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
 * @file   test_meshsets.cc
 * @author Rao V. Garimella
 * @date   Nov 2, 2016
 *
 * @brief  Unit tests for retrieving entities in a mesh
 *
 * We are testing retrieval of cells, faces, edges (only for MSTK) and
 * nodes in the mesh. The serial test is trivial, the parallel test
 * exercises the underlying machinery in MSTK
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Point.hh"
#include "BoxRegion.hh"
#include "PlaneRegion.hh"
#include "LogicalRegion.hh"
#include "LabeledSetRegion.hh"
#include "GeometricModel.hh"

TEST(MESH_SETS_3D) {
  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  int dim = 3;

  const Jali::MeshFramework_t frameworks[] = {Jali::MSTK, Jali::Simple};
  const char *framework_names[] = {"MSTK", "Simple"};
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::MeshFramework_t);
  Jali::MeshFramework_t the_framework;
  for (int i = 0; i < numframeworks; i++) {
    // Set the framework
    the_framework = frameworks[i];
    if (!Jali::framework_available(the_framework)) continue;
    
    bool parallel = (nproc > 1);
    if (!Jali::framework_generates(the_framework, parallel, dim))
      continue;
    
    std::cerr << "Testing mesh sets with " << framework_names[i] <<
        std::endl;
    
    // Create the mesh
    Jali::MeshFactory factory(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh;
    
    int ierr = 0;
    int aerr = 0;
    try {
      factory.framework(the_framework);

      std::vector<Jali::Entity_kind> entitylist;
      entitylist.push_back(Jali::Entity_kind::FACE);
      if (the_framework == Jali::MSTK)
        entitylist.push_back(Jali::Entity_kind::EDGE);
      entitylist.push_back(Jali::Entity_kind::NODE);

      factory.included_entities(entitylist);

      factory.partitioner(Jali::Partitioner_type::BLOCK);
      
      // Create a structured mesh
      
      mesh = factory(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 2, 2, 1);

    } catch (const Errors::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK_EQUAL(aerr, 0);

    // Now retrieve the entity counds and make sure we get what we expect

    int ncells = 0, nnodes = 0, nedges = 0, nfaces = 0;
    if (parallel) {
      
      // Check local number of cells in mesh

      int ncells_owned_local = mesh->num_entities(Jali::Entity_kind::CELL,
                                                  Jali::Entity_type::PARALLEL_OWNED);
      if (nproc == 4)
        CHECK_EQUAL(1, ncells_owned_local);
      MPI_Allreduce(&ncells_owned_local, &ncells, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int ncells_ghost_local = mesh->num_entities(Jali::Entity_kind::CELL,
                                                  Jali::Entity_type::PARALLEL_GHOST);
      if (nproc == 4)
        CHECK_EQUAL(3, ncells_ghost_local);
          
      int ncells_all_local = mesh->num_entities(Jali::Entity_kind::CELL,
                                                Jali::Entity_type::ALL);
      if (nproc == 4)
        CHECK_EQUAL(4, ncells_all_local);
      
      // Check local number of faces in mesh

      int nfaces_owned_local = mesh->num_entities(Jali::Entity_kind::FACE,
                                                  Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nfaces_owned_local, &nfaces, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      // Check local number of edges in mesh

      int nedges_owned_local = mesh->num_entities(Jali::Entity_kind::EDGE,
                                                  Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nedges_owned_local, &nedges, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      // Check local number of nodes in mesh

      int nnodes_owned_local = mesh->num_entities(Jali::Entity_kind::NODE,
                                                  Jali::Entity_type::PARALLEL_OWNED);
      MPI_Allreduce(&nnodes_owned_local, &nnodes, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);

      int nnodes_all_local = mesh->num_entities(Jali::Entity_kind::NODE,
                                                   Jali::Entity_type::ALL);
      if (nproc == 4)
        CHECK_EQUAL(18, nnodes_all_local);
    } else {
      ncells = mesh->num_entities(Jali::Entity_kind::CELL,
                                     Jali::Entity_type::PARALLEL_OWNED);
      nfaces = mesh->num_entities(Jali::Entity_kind::FACE,
                                     Jali::Entity_type::PARALLEL_OWNED);
      nedges = mesh->num_entities(Jali::Entity_kind::EDGE,
                                     Jali::Entity_type::PARALLEL_OWNED);
      nnodes = mesh->num_entities(Jali::Entity_kind::NODE,
                                     Jali::Entity_type::PARALLEL_OWNED);
    }

    // Check global number of cells and nodes in mesh
    CHECK_EQUAL(4, ncells);
    CHECK_EQUAL(20, nfaces);
    if (the_framework == Jali::MSTK) CHECK_EQUAL(33, nedges);
    CHECK_EQUAL(18, nnodes);
  }  // for each framework
}

