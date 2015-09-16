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
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

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
    Jali::Mesh *mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Jali::FrameworkPreference prefs(factory.preference());
      prefs.clear(); 
      prefs.push_back(the_framework);

      factory.preference(prefs);

      mesh = factory(0.0,0.0,0.0,1.0,1.0,1.0,2,2,2);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr,&aerr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    CHECK_EQUAL(aerr,0);


    int id1, id2;

    id1 = 0;
    Jali::Mesh::node_iterator nit = mesh->begin_owned_nodes();
    while (nit != mesh->end_owned_nodes()) {
      int nodeid = *nit;
      CHECK_EQUAL(id1,nodeid);
      ++nit;
      ++id1;
    }
    CHECK_EQUAL(id1,mesh->num_entities(Jali::NODE,Jali::OWNED));

    id2 = 0;
    nit = mesh->begin_ghost_nodes();
    while (nit != mesh->end_ghost_nodes()) {
      int nodeid = *nit;
      CHECK_EQUAL(id2,nodeid);
      ++nit;
      ++id2;
    }
    CHECK_EQUAL(id2,mesh->num_entities(Jali::NODE,Jali::GHOST));

    CHECK_EQUAL(id1+id2,mesh->num_entities(Jali::NODE,Jali::ALL));
             


    id1 = 0;
    Jali::Mesh::face_iterator fit = mesh->begin_owned_faces();
    while (fit != mesh->end_owned_faces()) {
      int faceid = *fit;
      CHECK_EQUAL(id1,faceid);
      ++fit;
      ++id1;
    }
    CHECK_EQUAL(id1,mesh->num_entities(Jali::FACE,Jali::OWNED));

    id2 = 0;
    fit = mesh->begin_ghost_faces();
    while (fit != mesh->end_ghost_faces()) {
      int faceid = *fit;
      CHECK_EQUAL(id2,faceid);
      ++fit;
      ++id2;
    }
    CHECK_EQUAL(id2,mesh->num_entities(Jali::FACE,Jali::GHOST));
    
    CHECK_EQUAL(id1+id2,mesh->num_entities(Jali::FACE,Jali::ALL));
                

    id1 = 0;
    Jali::Mesh::cell_iterator cit = mesh->begin_owned_cells();
    while (cit != mesh->end_owned_cells()) {
      int cellid = *cit;
      CHECK_EQUAL(id1,cellid);
      ++cit;
      ++id1;
    }
    CHECK_EQUAL(id1,mesh->num_entities(Jali::CELL,Jali::OWNED));

    id2 = 0;
    cit = mesh->begin_ghost_cells();
    while (cit != mesh->end_ghost_cells()) {
      int cellid = *cit;
      CHECK_EQUAL(id2,cellid);
      ++cit;
      ++id2;
    }
    CHECK_EQUAL(id2,mesh->num_entities(Jali::CELL,Jali::GHOST));

    CHECK_EQUAL(id1+id2,mesh->num_entities(Jali::CELL,Jali::ALL));

  } // for each framework i

}

