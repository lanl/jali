// -------------------------------------------------------------
// file: test_mesh_factory.cc
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 18, 2011 by William A. Perkins
// Last Change: Tue Aug  2 10:58:57 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <iostream>
#include <UnitTest++.h>

#include "dbc.hh"
#include "../MeshFactory.hh"
#include "../FrameworkTraits.hh"

// Check to see if we have some files to read

#ifndef BOGUS_TEST_FILE
#define BOGUS_TEST_FILE "bogus.file"
#endif

#ifndef EXODUS_TEST_FILE
#error EXODUS_TEST_FILE must be defined
#endif

#ifndef NEMESIS_TEST_FILE
#error EXODUS_TEST_FILE must be defined
#endif

#ifndef MOAB_TEST_FILE
#error EXODUS_TEST_FILE must be defined
#endif

// -------------------------------------------------------------
// check_preference
// -------------------------------------------------------------

static void
check_preference(Jali::MeshFactory& mesh_factory,
                 const Jali::Framework& f)
{
  Jali::FrameworkPreference pref;
  pref.push_back(f);
  if (Jali::framework_available(f)) {
    mesh_factory.preference(pref);
    CHECK(mesh_factory.preference().front() == f);
  } else {
    CHECK_THROW(mesh_factory.preference(pref), Jali::Message);
  }
}


SUITE (MeshFramework)
{

  // This tests setting the Mesh Factory framework preferences.  If
  // only one framework is preferred, and it is not available, an
  // exception should be thrown, while setting preferences
  TEST (PreferenceThrow)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);
    Jali::FrameworkPreference pref(mesh_factory.preference());

    // The Simple framework should always be there
    check_preference(mesh_factory, Jali::Simple);
    check_preference(mesh_factory, Jali::MOAB);
    check_preference(mesh_factory, Jali::STKMESH);
    check_preference(mesh_factory, Jali::MSTK);
  }

  TEST (Generate)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    Jali::FrameworkPreference pref;
    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    double x0(0.0), y0(0.0), z0(0.0);
    double x1(10.0), y1(10.0), z1(10.0);
    int nx(10), ny(10), nz(10);

    // The Simple framework is always available, but will only
    // generate in serial

    pref.clear(); pref.push_back(Jali::Simple);
    mesh_factory.preference(pref);

    if (parallel) {
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  Jali::Message);
      mesh.reset();
    } else {
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(mesh);
      mesh.reset();
    }

    // The STK, if available, framework will always generate

    if (framework_available(Jali::STKMESH)) {
      pref.clear(); pref.push_back(Jali::STKMESH);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(mesh);
      mesh.reset();
    }

    // The MSTK framework, if available, will always generate

    if (framework_available(Jali::MSTK)) {
      pref.clear(); pref.push_back(Jali::MSTK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(mesh);
      mesh.reset();
    }

    // The MOAB framework cannot generate


    if (framework_available(Jali::MOAB)) {
      pref.clear(); pref.push_back(Jali::MOAB);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  Jali::Message);
      mesh.reset();
    }
  }

  TEST (Generate2D)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    Jali::FrameworkPreference pref;
    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    double x0(0.0), y0(0.0);
    double x1(10.0), y1(10.0);
    int nx(10), ny(10);

    // The MSTK framework, if available, will always generate

    if (framework_available(Jali::MSTK)) {
      pref.clear(); pref.push_back(Jali::MSTK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0,
                          x1, y1,
                          nx, ny);
      CHECK(mesh);
      mesh.reset();
    }

    // The Simple framework is always available, but
    // cannot generate 2D meshes

    pref.clear(); pref.push_back(Jali::Simple);
    mesh_factory.preference(pref);

    CHECK_THROW(mesh = mesh_factory(x0, y0,
                                    x1, y1,
                                    nx, ny),
                Jali::Message);
    mesh.reset();

    // The STK, even if available cannot generate 2D meshes

    if (framework_available(Jali::STKMESH)) {
      pref.clear(); pref.push_back(Jali::STKMESH);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0,
                                      x1, y1,
                                      nx, ny),
                  Jali::Message);
      mesh.reset();
    }

    // The MOAB framework cannot generate


    if (framework_available(Jali::MOAB)) {
      pref.clear(); pref.push_back(Jali::MOAB);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0,
                                      x1, y1,
                                      nx, ny),
                  Jali::Message);
      mesh.reset();
    }
  }


  // The Simple framework cannot read anything, even if it exists
  TEST (ReadSimple) {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    Jali::FrameworkPreference pref;
    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    pref.clear(); pref.push_back(Jali::Simple);
    mesh_factory.preference(pref);
    CHECK_THROW(mesh = mesh_factory(BOGUS_TEST_FILE),
                Jali::Message);
    CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                Jali::Message);
    CHECK_THROW(mesh = mesh_factory(EXODUS_TEST_FILE),
                Jali::Message);
    CHECK_THROW(mesh = mesh_factory(NEMESIS_TEST_FILE),
                Jali::Message);
  }

  // Try to read a MOAB HDF5 file, which can only be read by the MOAB
  // framework
  TEST (ReadMOABHDF5)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    if (Jali::framework_available(Jali::MOAB)) {
      mesh = mesh_factory(MOAB_TEST_FILE);
      CHECK(mesh);
      mesh.reset();
    } else {
      CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                  Jali::Message);
    }

    // Try it with another framework just for grins

    if (Jali::framework_available(Jali::STKMESH)) {
      Jali::FrameworkPreference pref;
      pref.push_back(Jali::STKMESH);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                  Jali::Message);
    }
  }

  TEST (ReadExodus)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    bool available =
        (Jali::framework_available(Jali::STKMESH) && !parallel) ||
        Jali::framework_available(Jali::MSTK);

    if (available) {
      mesh = mesh_factory(EXODUS_TEST_FILE);
      CHECK(mesh);
    } else {
      CHECK_THROW(mesh = mesh_factory(EXODUS_TEST_FILE),
                  Jali::Message);
    }
  }

  TEST (ReadNemesis)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);
    if ((Jali::framework_available(Jali::STKMESH) ||
         Jali::framework_available(Jali::MSTK)) &&
        parallel) {
      mesh = mesh_factory(NEMESIS_TEST_FILE);
      CHECK(mesh);
    } else {
      CHECK_THROW(mesh = mesh_factory(NEMESIS_TEST_FILE),
                  Jali::Message);
    }
  }


  // Make sure mesh factory options can be set and transmitted
  // correctly to the mesh constructors

  TEST (Options)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    // Check which entities are included by default in mesh generation
    // Should be NODES, FACES and CELLS

    std::vector<Jali::Entity_kind> entity_list;
    entity_list = mesh_factory.included_entities();

    CHECK(std::find(entity_list.begin(), entity_list.end(),
                    Jali::Entity_kind::NODE) != entity_list.end());
    CHECK(std::find(entity_list.begin(), entity_list.end(),
                    Jali::Entity_kind::FACE) != entity_list.end());
    CHECK(std::find(entity_list.begin(), entity_list.end(),
                    Jali::Entity_kind::CELL) != entity_list.end());

    // Turn edges on 
    mesh_factory.included_entities({Jali::Entity_kind::EDGE,
            Jali::Entity_kind::FACE, Jali::Entity_kind::CELL});

    mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
    CHECK(mesh);
    CHECK(mesh->num_edges());
    mesh.reset();
    mesh_factory.reset_options();

    // Turn wedges on (faces and edges are turned on when wedges are turned on)
    mesh_factory.included_entities({Jali::Entity_kind::WEDGE});

    mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
    CHECK(mesh);
    CHECK(mesh->num_faces());
    CHECK(mesh->num_edges());
    CHECK(mesh->num_wedges());
    mesh.reset();
    mesh_factory.reset_options();

    // Turn corners on (wedges, edges and faces get turned on as well)
    mesh_factory.included_entities({Jali::Entity_kind::CORNER});

    mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
    CHECK(mesh);
    CHECK(mesh->num_faces());
    CHECK(mesh->num_edges());  
    CHECK(mesh->num_wedges());
    CHECK(mesh->num_corners());
    mesh.reset();
    mesh_factory.reset_options();

    
    // Generate 1D meshes with CARTESIAN and record the volume of the
    // cell farthest away from the origin (since it is a 5 cell mesh,
    // this would be cell 4)
    std::vector<double> x = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    mesh = mesh_factory(x);
    CHECK(mesh);
    int c = 4;
    double cart_vol = mesh->cell_volume(c);
    mesh.reset();
    mesh_factory.reset_options();

    mesh_factory.mesh_geometry(JaliGeometry::Geom_type::SPHERICAL);
    mesh = mesh_factory(x);
    CHECK(mesh);
    double sph_vol = mesh->cell_volume(c);

    CHECK(sph_vol != cart_vol);

    if (parallel) {
      Jali::FrameworkPreference pref;
      if (Jali::framework_available(Jali::MSTK) ||
          Jali::framework_available(Jali::STKMESH)) {

        if (Jali::framework_available(Jali::MSTK))
          pref.push_back(Jali::MSTK);
        if (Jali::framework_available(Jali::STKMESH))
          pref.push_back(Jali::STKMESH);

        // Make sure we have ghost cells if we asked for them
        mesh_factory.num_ghost_layers_distmesh(1);
        mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
        CHECK(mesh);
        
        CHECK(mesh->num_cells<Jali::Parallel_type::GHOST>());
        mesh.reset();
        
        // Make sure we don't have ghost cells if we didn't ask for them
        mesh_factory.num_ghost_layers_distmesh(0);
        mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
        CHECK(mesh);
        
        CHECK(mesh->num_cells<Jali::Parallel_type::GHOST>() == 0);
        mesh.reset();
        mesh_factory.reset_options();
      
        // Make sure we can choose partitioners if we want to.  The
        // default partitioner is METIS. Instead we choose ZOLTAN's
        // Recursive Coordinate Bisection algorithm to partition a
        // mesh that is thin along the Z-direction.  The algorithm
        // should naturally partition only along the X-Y directions
        // and not Z

        mesh_factory.partitioner(Jali::Partitioner_type::ZOLTAN_RCB); 
        mesh = mesh_factory(0.0, 0.0, 0.0, 100.0, 100.0, 1.0, 4, 4, 10);
        CHECK(mesh);

        // Make sure any horizontal face that has only one cell
        // connected to it is only on the top or bottom boundary of
        // the domain. If not, it would indicate that the partitioning
        // took pace along the Z-direction

        for (auto const& f : mesh->faces()) {
          JaliGeometry::Point fnormal = mesh->face_normal(f);
          if (fnormal[2] == 0.0) continue;  // z-component of normal is 0
          
          Jali::Entity_ID_List fregs;
          mesh->face_get_cells(f, Jali::Parallel_type::OWNED, &fregs);

          if (fregs.size() == 2) continue;  // Not a face on a mesh boundary

          JaliGeometry::Point fcen = mesh->face_centroid(f);
          CHECK(fabs(fcen[2]) < 1.0e-06 || fabs(fcen[2]-1.0) < 1.0e-06);
        }
        mesh.reset();
        mesh_factory.reset_options();
      }
    }
  }
}
