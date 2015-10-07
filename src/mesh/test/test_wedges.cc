// -------------------------------------------------------------
/**
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
#include "FrameworkTraits.hh"
#include "Geometry.hh"

TEST(MESH_WEDGES_2D)
{

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  const Jali::Framework frameworks[] = {  
    Jali::MSTK
  };
  const char *framework_names[] = {
    "MSTK"
  };
  
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::Framework);


  Jali::Framework the_framework;
  for (int fr = 0; fr < numframeworks; fr++) {


    // Set the framework

    the_framework = frameworks[fr];

    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing wedge operators with " << framework_names[fr] << std::endl;


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

      mesh = factory(0.0,0.0,1.0,1.0,2,2,NULL,true,true,true,false);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr,&aerr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    CHECK_EQUAL(aerr,0);

    int ncells = mesh->num_entities(Jali::CELL,Jali::ALL);
    for (Jali::Entity_ID c = 0; c < ncells; ++c) {
      std::vector<Jali::Entity_ID> cwedges;
      mesh->cell_get_wedges(c,&cwedges);

      // Quad elements have 8 wedges

      CHECK_EQUAL(cwedges.size(),8);

      double cellvol = mesh->cell_volume(c);
            
      Jali::Entity_ID_List::iterator itw = cwedges.begin();
      while (itw != cwedges.end()) {
        Jali::Entity_ID w = *itw;

        // Since the wedge came from the cell 'c' the cell of of the
        // wedge must be 'c'

        Jali::Entity_ID wc = mesh->wedge_get_cell(w);        
        CHECK_EQUAL(c,wc);

        JaliGeometry::Point ccen = mesh->cell_centroid(c);

        // Make sure the wedge knows which edge its associated with

        Jali::Entity_ID e = mesh->wedge_get_edge(w);
        CHECK(e != -1);

        // Make sure the wedge knows which face its associated with

        Jali::Entity_ID f = mesh->wedge_get_face(w);
        CHECK(f != -1);
        JaliGeometry::Point fcen = mesh->face_centroid(f);

        // Make sure the wedge knows which node its associated with

        Jali::Entity_ID n = mesh->wedge_get_node(w);
        CHECK(n != -1);
        JaliGeometry::Point npnt;
        mesh->node_get_coordinates(n,&npnt);

        // Get the normals to facet0 and facet1 of w;

        JaliGeometry::Point normal0 = mesh->wedge_facet_normal(w,0);
        double norm0 = norm(normal0);
        JaliGeometry::Point normal1 = mesh->wedge_facet_normal(w,1);
        double norm1 = norm(normal1);

        double dp;

        // Make sure the link between the edge and the opposite wedge
        // are correct

        Jali::Entity_ID w2 = mesh->wedge_get_opposite_wedge(w);
        if (w2 != -1) {
          CHECK_EQUAL(w,mesh->wedge_get_opposite_wedge(w2));
          CHECK_EQUAL(f,mesh->wedge_get_face(w2));
          CHECK_EQUAL(e,mesh->wedge_get_edge(w2));
          CHECK_EQUAL(n,mesh->wedge_get_node(w2));
          CHECK(c != mesh->wedge_get_cell(w2));

          // Also make sure the two wedges have equal and opposing
          // normals for their facet 0

          JaliGeometry::Point normal20 = mesh->wedge_facet_normal(w2,0);
          double norm20 = norm(normal20);
          CHECK_CLOSE(norm0,norm20,1.0e-6);
          dp = (normal0*normal20)/(norm0*norm20);
          CHECK_CLOSE(-1.0,dp,1.0e-6);

          // Also make sure the two wedges have equal and coincident
          // normals to facet 1
          
          JaliGeometry::Point normal21 = mesh->wedge_facet_normal(w2,1);
          double norm21 = norm(normal21);          
          CHECK_CLOSE(norm1,norm21,1.0e-6);
          dp = (normal1*normal21)/(norm1*norm21);
          CHECK_CLOSE(1.0,dp,1.0e-6);
        }

        // Make sure the link between the edge and the adjacent wedge
        // are correct

        Jali::Entity_ID w3 = mesh->wedge_get_adjacent_wedge(w);
        CHECK_EQUAL(w,mesh->wedge_get_adjacent_wedge(w3));
        CHECK_EQUAL(f,mesh->wedge_get_face(w3));
        CHECK_EQUAL(e,mesh->wedge_get_edge(w3));
        CHECK_EQUAL(c,mesh->wedge_get_cell(w3));

        // Also make sure the two wedges have equal and coincident
        // normals for their facet 0

        JaliGeometry::Point normal30 = mesh->wedge_facet_normal(w3,0);
        double norm30 = norm(normal30);
        CHECK_CLOSE(norm0,norm30,1.0e-6);
        dp = (normal0*normal30)/(norm0*norm30);
        CHECK_CLOSE(1.0,dp,1.0e-6);

        // Also make sure the two wedges have equal and opposite
        // normals to facet 1
          
        JaliGeometry::Point normal31 = mesh->wedge_facet_normal(w3,1);
        double norm31 = norm(normal31);          
        CHECK_CLOSE(norm1,norm31,1.0e-6);
        dp = (normal1*normal31)/(norm1*norm31);
        CHECK_CLOSE(-1.0,dp,1.0e-6);


        // Get the opposite wedge (w4) of the adjacent wedge (w3) of w
        // and make sure that it is adjacent to the opposite wedge
        // (w2) of w.

        Jali::Entity_ID w4 = mesh->wedge_get_opposite_wedge(w3);
        if (w4 != -1 && w2 != -1) 
          CHECK_EQUAL(w4,mesh->wedge_get_adjacent_wedge(w2));
        
        // Get wedge coordinates and make sure they match up with the
        // expected coordinates (node, edge center, face center, cell
        // center)

        std::vector<JaliGeometry::Point> wcoords;
        mesh->wedge_get_coordinates(w,&wcoords);

        for (int i = 0; i < 2; ++i) {
          CHECK_EQUAL(wcoords[0][i],npnt[i]);
          CHECK_EQUAL(wcoords[1][i],fcen[i]);
          CHECK_EQUAL(wcoords[2][i],ccen[i]);
        }

        // Since there are 48 wedges in a hex element, its volume should
        // be 1/48th that of the cell

        double volume = mesh->wedge_volume(w);
        CHECK_CLOSE(volume,cellvol/cwedges.size(),1.0e-06);

        ++itw;
      } // while (itw != cwedges.end())
    }

    // Check node to wedge links

    int nnodes = mesh->num_entities(Jali::NODE,Jali::ALL);
    
  }

}

TEST(MESH_WEDGES_3D)
{

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  const Jali::Framework frameworks[] = {  
    Jali::MSTK
  };
  const char *framework_names[] = {
    "MSTK"
  };
  
  const int numframeworks = sizeof(frameworks)/sizeof(Jali::Framework);


  Jali::Framework the_framework;
  for (int fr = 0; fr < numframeworks; fr++) {


    // Set the framework

    the_framework = frameworks[fr];

    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing wedge operators with " << framework_names[fr] << std::endl;


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

      mesh = factory(0.0,0.0,0.0,1.0,1.0,1.0,2,2,2,NULL,true,true,true,false);

    } catch (const Jali::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    MPI_Allreduce(&ierr,&aerr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    CHECK_EQUAL(aerr,0);

    int ncells = mesh->num_entities(Jali::CELL,Jali::ALL);
    for (Jali::Entity_ID c = 0; c < ncells; ++c) {
      std::vector<Jali::Entity_ID> cwedges;
      mesh->cell_get_wedges(c,&cwedges);

      // Hex elements have 48 wedges

      CHECK_EQUAL(cwedges.size(),48);

      double cellvol = mesh->cell_volume(c);
            
      Jali::Entity_ID_List::iterator itw = cwedges.begin();
      while (itw != cwedges.end()) {
        Jali::Entity_ID w = *itw;

        // Since the wedge came from the cell 'c' the cell of of the
        // wedge must be 'c'

        Jali::Entity_ID wc = mesh->wedge_get_cell(w);        
        CHECK_EQUAL(c,wc);

        JaliGeometry::Point ccen = mesh->cell_centroid(c);

        // Make sure the wedge knows which edge its associated with

        Jali::Entity_ID e = mesh->wedge_get_edge(w);
        CHECK(e != -1);

        JaliGeometry::Point ecen = mesh->edge_centroid(e);

        // Make sure the wedge knows which face its associated with

        Jali::Entity_ID f = mesh->wedge_get_face(w);
        CHECK(f != -1);
        JaliGeometry::Point fcen = mesh->face_centroid(f);

        // Make sure the wedge knows which node its associated with

        Jali::Entity_ID n = mesh->wedge_get_node(w);
        CHECK(n != -1);
        JaliGeometry::Point npnt;
        mesh->node_get_coordinates(n,&npnt);

        // Get the normals to facet0 and facet1 of w;

        JaliGeometry::Point normal0 = mesh->wedge_facet_normal(w,0);
        double norm0 = norm(normal0);
        JaliGeometry::Point normal1 = mesh->wedge_facet_normal(w,1);
        double norm1 = norm(normal1);

        double dp;

        // Make sure the link between the edge and the opposite wedge
        // are correct

        Jali::Entity_ID w2 = mesh->wedge_get_opposite_wedge(w);
        if (w2 != -1) {
          CHECK_EQUAL(w,mesh->wedge_get_opposite_wedge(w2));
          CHECK_EQUAL(f,mesh->wedge_get_face(w2));
          CHECK_EQUAL(e,mesh->wedge_get_edge(w2));
          CHECK_EQUAL(n,mesh->wedge_get_node(w2));
          CHECK(c != mesh->wedge_get_cell(w2));

          // Also make sure the two wedges have equal and opposing
          // normals for their facet 0

          JaliGeometry::Point normal20 = mesh->wedge_facet_normal(w2,0);
          double norm20 = norm(normal20);
          CHECK_CLOSE(norm0,norm20,1.0e-6);
          dp = (normal0*normal20)/(norm0*norm20);
          CHECK_CLOSE(-1.0,dp,1.0e-6);

          // Also make sure the two wedges have equal and coincident
          // normals to facet 1
          
          JaliGeometry::Point normal21 = mesh->wedge_facet_normal(w2,1);
          double norm21 = norm(normal21);          
          CHECK_CLOSE(norm1,norm21,1.0e-6);
          dp = (normal1*normal21)/(norm1*norm21);
          CHECK_CLOSE(1.0,dp,1.0e-6);
        }

        // Make sure the link between the edge and the adjacent wedge
        // are correct

        Jali::Entity_ID w3 = mesh->wedge_get_adjacent_wedge(w);
        CHECK_EQUAL(w,mesh->wedge_get_adjacent_wedge(w3));
        CHECK_EQUAL(f,mesh->wedge_get_face(w3));
        CHECK_EQUAL(e,mesh->wedge_get_edge(w3));
        CHECK_EQUAL(c,mesh->wedge_get_cell(w3));

        // Also make sure the two wedges have equal and coincident
        // normals for their facet 0

        JaliGeometry::Point normal30 = mesh->wedge_facet_normal(w3,0);
        double norm30 = norm(normal30);
        CHECK_CLOSE(norm0,norm30,1.0e-6);
        dp = (normal0*normal30)/(norm0*norm30);
        CHECK_CLOSE(1.0,dp,1.0e-6);

        // Also make sure the two wedges have equal and opposite
        // normals to facet 1
          
        JaliGeometry::Point normal31 = mesh->wedge_facet_normal(w3,1);
        double norm31 = norm(normal31);          
        CHECK_CLOSE(norm1,norm31,1.0e-6);
        dp = (normal1*normal31)/(norm1*norm31);
        CHECK_CLOSE(-1.0,dp,1.0e-6);


        // Get the opposite wedge (w4) of the adjacent wedge (w3) of w
        // and make sure that it is adjacent to the opposite wedge
        // (w2) of w.

        Jali::Entity_ID w4 = mesh->wedge_get_opposite_wedge(w3);
        if (w4 != -1 && w2 != -1) 
          CHECK_EQUAL(w4,mesh->wedge_get_adjacent_wedge(w2));
        
        // Get wedge coordinates and make sure they match up with the
        // expected coordinates (node, edge center, face center, cell
        // center)

        std::vector<JaliGeometry::Point> wcoords;
        mesh->wedge_get_coordinates(w,&wcoords);

        for (int i = 0; i < 3; ++i) {
          CHECK_EQUAL(wcoords[0][i],npnt[i]);
          CHECK_EQUAL(wcoords[1][i],ecen[i]);
          CHECK_EQUAL(wcoords[2][i],fcen[i]);
          CHECK_EQUAL(wcoords[3][i],ccen[i]);
        }

        // Since there are 48 wedges in a hex element, its volume should
        // be 1/48th that of the cell

        double volume = mesh->wedge_volume(w);
        CHECK_CLOSE(volume,cellvol/cwedges.size(),1.0e-06);

        ++itw;
      }
    }
 }

}

