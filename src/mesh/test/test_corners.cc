// -------------------------------------------------------------
/**
 * @File   test_wedges.cc
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

TEST(MESH_CORNERS_2D) {

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
  for (int i = 0; i < numframeworks; i++) {


    // Set the framework

    the_framework = frameworks[i];

    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing wedge operators with " << framework_names[i] << std::endl;


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

      mesh = factory(0.0,0.0,1.0,1.0,2,2,NULL,true,true,true,true);

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
      std::vector<Jali::Entity_ID> ccorners;
      mesh->cell_get_corners(c,&ccorners);

      // Quad elements have 4 corners

      CHECK_EQUAL(4,ccorners.size());

      double cellvol = mesh->cell_volume(c);
            
      Jali::Entity_ID_List::iterator itc = ccorners.begin();
      while (itc != ccorners.end()) {
        Jali::Entity_ID cn = *itc;

        // Since the corner came from the cell 'c' the cell of of the
        // corner must be 'c'

        Jali::Entity_ID wc = mesh->corner_get_cell(cn);        
        CHECK_EQUAL(c,wc);

        JaliGeometry::Point ccen = mesh->cell_centroid(c);

        // Make sure the corner knows which node its associated with

        Jali::Entity_ID n = mesh->corner_get_node(cn);
        CHECK(n != -1);
        JaliGeometry::Point npnt;
        mesh->node_get_coordinates(n,&npnt);

        // Get the wedges of the corner

        Jali::Entity_ID_List cnwedges;
        
        mesh->corner_get_wedges(cn,&cnwedges);
        CHECK_EQUAL(2,cnwedges.size());   // corner in a 2D cell has 2 wedges

        double volwedges = 0;

        Jali::Entity_ID_List::iterator itw = cnwedges.begin();
        while (itw != cnwedges.end()) {
          Jali::Entity_ID w = *itw;

          // Make sure that the node of the wedge is the same as the
          // node of the corner

          CHECK_EQUAL(n,mesh->wedge_get_node(w));
          
          // Add up the volume of the wedges

          volwedges += mesh->wedge_volume(w);

          ++itw;
        }


        double cnvolume = mesh->corner_volume(cn);

        CHECK_CLOSE(cellvol/4.0,cnvolume,1.0e-06);
        CHECK_CLOSE(volwedges,cnvolume,1.0e-06);


        // Get the facets of the corner and calculate the "volume" a
        // different way (by getting the facets, connecting them to
        // the geometric center of the corner polyhedron to form tris
        // and adding up those tri "volumes")

        std::vector<JaliGeometry::Point> cncoords;

        mesh->corner_get_coordinates(cn,&cncoords);

        CHECK_EQUAL(4,cncoords.size());

        // Find geometric center of corner

        JaliGeometry::Point cncen(mesh->space_dimension());
        cncen.set(0.0);
        
        std::vector<JaliGeometry::Point>::iterator itp = cncoords.begin();
        while (itp != cncoords.end()) {
          JaliGeometry::Point xyz = *itp;
          cncen += xyz;
          ++itp;
        }
        cncen /= 4;

        // Now form a tet with each of the corner facets and the
        // center of the corner and compute its volume; Add this volume
        // to the alternate form of the corner volume

        double cnvolume2 = 0.0;
        for (int i = 0; i < 4; i++) {
          JaliGeometry::Point vec0 = cncoords[(i+1)%4] - cncoords[i];
          JaliGeometry::Point vec1 = cncen - cncoords[i];
          JaliGeometry::Point cpvec = vec0^vec1;
          double trivol = cpvec[0];
          cnvolume2 += trivol;
        }

        CHECK_CLOSE(cnvolume,cnvolume2,1.0e-06);

        ++itc;
      } // while (itc != ccorners.end())
    } // for c = 0, ncells
  }

}



TEST(MESH_CORNERS_3D) {

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
  for (int i = 0; i < numframeworks; i++) {


    // Set the framework

    the_framework = frameworks[i];

    if (!Jali::framework_available(the_framework)) continue;

    std::cerr << "Testing wedge operators with " << framework_names[i] << std::endl;


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

    int ncells = mesh->num_entities(Jali::CELL,Jali::ALL);
    for (Jali::Entity_ID c = 0; c < ncells; ++c) {
      std::vector<Jali::Entity_ID> ccorners;
      mesh->cell_get_corners(c,&ccorners);

      // Hex elements have 8 corners

      CHECK_EQUAL(8,ccorners.size());

      double cellvol = mesh->cell_volume(c);
            
      Jali::Entity_ID_List::iterator itc = ccorners.begin();
      while (itc != ccorners.end()) {
        Jali::Entity_ID cn = *itc;

        // Since the corner came from the cell 'c' the cell of of the
        // corner must be 'c'

        Jali::Entity_ID wc = mesh->corner_get_cell(cn);        
        CHECK_EQUAL(c,wc);

        JaliGeometry::Point ccen = mesh->cell_centroid(c);

        // Make sure the corner knows which node its associated with

        Jali::Entity_ID n = mesh->corner_get_node(cn);
        CHECK(n != -1);
        JaliGeometry::Point npnt;
        mesh->node_get_coordinates(n,&npnt);

        // Get the wedges of the corner

        Jali::Entity_ID_List cnwedges;
        
        mesh->corner_get_wedges(cn,&cnwedges);
        CHECK_EQUAL(6,cnwedges.size());   // 6 wedges in corner at trivalent node of 3D cell

        double volwedges = 0;

        Jali::Entity_ID_List::iterator itw = cnwedges.begin();
        while (itw != cnwedges.end()) {
          Jali::Entity_ID w = *itw;

          // Make sure that the node of the wedge is the same as the
          // node of the corner

          CHECK_EQUAL(n,mesh->wedge_get_node(w));
          
          // Add up the volume of the wedges

          volwedges += mesh->wedge_volume(w);

          ++itw;
        }


        double cnvolume = mesh->corner_volume(cn);

        CHECK_CLOSE(cellvol/8.0,cnvolume,1.0e-06);
        CHECK_CLOSE(volwedges,cnvolume,1.0e-06);


        // Get the facets of the corner and calculate the "volume" a
        // different way (by getting the facets, connecting them to
        // the geometric center of the corner polyhedron to form tris
        // and adding up those tri "volumes")

        std::vector<JaliGeometry::Point> cncoords;
        std::vector< std::array<int,3> > facets;

        mesh->corner_get_facetization(cn,&cncoords,&facets);

        CHECK_EQUAL(8,cncoords.size()); // 8 nodes in the corner
        CHECK_EQUAL(12,facets.size());  // 6 "faces" with 2 facets each

        // Find geometric center of corner

        JaliGeometry::Point cncen(mesh->space_dimension());
        cncen.set(0.0);
        
        std::vector<JaliGeometry::Point>::iterator itp = cncoords.begin();
        while (itp != cncoords.end()) {
          JaliGeometry::Point xyz = *itp;
          cncen += xyz;
          ++itp;
        }
        cncen /= cncoords.size();

        // Now form a tet with each of the corner facets and the
        // center of the corner and compute its volume; Add this volume
        // to the alternate form of the corner volume

        double cnvolume2 = 0.0;
        std::vector< std::array<int,3> >::iterator itf = facets.begin();
        while (itf != facets.end()) {
          std::array<int,3> fctpnts = *itf;
          JaliGeometry::Point vec0 = cncoords[fctpnts[1]] - cncoords[fctpnts[0]];
          JaliGeometry::Point vec1 = cncoords[fctpnts[2]] - cncoords[fctpnts[0]];
          JaliGeometry::Point vec2 = cncen - cncoords[fctpnts[0]];
          JaliGeometry::Point cpvec = vec0^vec1; // outward facing normal
          double tetvol = -(cpvec*vec2)/6.0;
          cnvolume2 += tetvol;
          ++itf;
        }

        CHECK_CLOSE(cnvolume,cnvolume2,1.0e-06);

        ++itc;
      } // while (itc != ccorners.end())
    } // for c = 0, ncells
  }

} // MESH_CORNERS_3D

