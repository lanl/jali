// Copyright Los Alamos National Security, LLC 2016-

#include <iostream>
#include <UnitTest++.h>

#include "dbc.hh"
#include "../MeshFactory.hh"

// Check to see if we have some files to read

SUITE (MeshFramework)
{
  

  TEST (Generate)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    double x0(0.0), y0(0.0), z0(0.0);
    double x1(10.0), y1(10.0), z1(10.0);
    int nx(10), ny(10), nz(10);

    // The Simple framework is always available, but will only
    // generate in serial

    mesh_factory.framework(Jali::Simple);

    if (parallel) {
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0, x1, y1, z1, nx, ny, nz),
                  Errors::Message);
      mesh.reset();
    } else {
      mesh = mesh_factory(x0, y0, z0, x1, y1, z1, nx, ny, nz);
      CHECK(mesh);
      mesh.reset();
    }

    // The MSTK framework will always generate

    mesh_factory.framework(Jali::MSTK);
    mesh = mesh_factory(x0, y0, z0, x1, y1, z1, nx, ny, nz);
    CHECK(mesh);
    mesh.reset();

    // STKmesh can generate in serial and parallel (but its not available)

    CHECK_THROW(mesh_factory.framework(Jali::STKMESH), Errors::Message);
    mesh.reset();
    
    //      mesh = mesh_factory(x0, y0, z0, x1, y1, z1, nx, ny, nz);
    //      CHECK(mesh);
    //      mesh.reset();
    //    }

    // The MOAB framework is not available 


    CHECK_THROW(mesh_factory.framework(Jali::MOAB), Errors::Message);
    mesh.reset();
    
    //      CHECK_THROW(mesh = mesh_factory(x0, y0, z0, x1, y1, z1, nx, ny, nz),
    //                  Errors::Message);
    //      mesh.reset();
  }

  TEST (Generate2D)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    double x0(0.0), y0(0.0);
    double x1(10.0), y1(10.0);
    int nx(10), ny(10);

    // The MSTK framework, if available, will always generate

    if (Jali::framework_available(Jali::MSTK)) {
      mesh_factory.framework(Jali::MSTK);
      mesh = mesh_factory(x0, y0, x1, y1, nx, ny);
      CHECK(mesh);
      mesh.reset();
    }

    // The Simple framework is always available, but cannot generate
    // 2D meshes

    mesh_factory.framework(Jali::Simple);

    CHECK_THROW(mesh = mesh_factory(x0, y0, x1, y1, nx, ny), Errors::Message);
    mesh.reset();

    // The STKmesh is not available

    CHECK(!Jali::framework_available(Jali::STKMESH));
    CHECK_THROW(mesh_factory.framework(Jali::STKMESH), Errors::Message);
    //  CHECK_THROW(mesh = mesh_factory(x0, y0, x1, y1, nx, ny), Errors::Message);
    //  mesh.reset();

    // The MOAB framework is not available
    
    CHECK(!Jali::framework_available(Jali::MOAB));
    CHECK_THROW(mesh_factory.framework(Jali::MOAB), Errors::Message);
    //  CHECK_THROW(mesh = mesh_factory(x0, y0, x1, y1, nx, ny), Errors::Message);
    //  mesh.reset();

  }


  // The Simple framework cannot read anything, even if it exists
  TEST (ReadSimple) {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    mesh_factory.framework(Jali::Simple);
    CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                Errors::Message);
    CHECK_THROW(mesh = mesh_factory(EXODUS_TEST_FILE),
                Errors::Message);
    CHECK_THROW(mesh = mesh_factory(NEMESIS_TEST_FILE),
                Errors::Message);
  }

  // Try to read a MOAB HDF5 file, which can only be read by the MOAB
  // framework - but the MOAB framework is not available in Jali
  TEST (ReadMOABHDF5) {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    if (Jali::framework_available(Jali::MOAB)) {
      mesh_factory.framework(Jali::MOAB);
      mesh = mesh_factory(MOAB_TEST_FILE);
      CHECK(mesh);
      mesh.reset();
    } else {
      CHECK_THROW(mesh_factory.framework(Jali::MOAB), Errors::Message);
    }
  }

  TEST (ReadExodus) {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    if (Jali::framework_available(Jali::MSTK)) {
      mesh_factory.framework(Jali::MSTK);
      CHECK(mesh = mesh_factory(EXODUS_TEST_FILE));
    } else if (Jali::framework_available(Jali::MSTK)) {
      mesh_factory.framework(Jali::STKMESH);
      CHECK(mesh = mesh_factory(EXODUS_TEST_FILE));
    }
  }

  // Only MSTK can read parallel Exodus II (Nemesis I) files as implemented
  TEST (ReadNemesis) {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    bool parallel(nproc > 1);

    std::shared_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    if (Jali::framework_available(Jali::MSTK) && parallel) {
      mesh_factory.framework(Jali::MSTK);
      CHECK(mesh = mesh_factory(NEMESIS_TEST_FILE));
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

    // Turn faces, edges on 
    mesh_factory.included_entities({Jali::Entity_kind::EDGE,
            Jali::Entity_kind::FACE, Jali::Entity_kind::CELL});

    mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
    CHECK(mesh);
    CHECK(mesh->num_edges());
    CHECK(mesh->num_faces());
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


    // Turn on all entities using Entity_kind::ALL_KIND and check that
    // they are present

    mesh_factory.included_entities(Jali::Entity_kind::ALL_KIND);
    
    mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
    CHECK(mesh);    
    CHECK(mesh->num_faces());
    CHECK(mesh->num_edges());
    CHECK(mesh->num_wedges());
    CHECK(mesh->num_corners());
    CHECK(mesh->num_cells());
    mesh.reset();
    mesh_factory.reset_options();


    // Generate 1D meshes with CARTESIAN and record the volume of the
    // cell farthest away from the origin (since it is a 5 cell mesh,
    // this would be cell 4) - ONLY IN Simple FRAMEWORK

    int c = 4;
    double cart_vol = 0.0, sph_vol = 0.0;
    std::vector<double> x = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    mesh_factory.framework(Jali::Simple);
    
    if (parallel)
      CHECK_THROW(mesh_factory(x), Errors::Message);
    else {
      mesh = mesh_factory(x);
      CHECK(mesh);
      cart_vol = mesh->cell_volume(c);      
      mesh.reset();
    }
    mesh_factory.reset_options();

    mesh_factory.mesh_geometry(JaliGeometry::Geom_type::SPHERICAL);
    mesh_factory.framework(Jali::Simple);
    if (parallel)
      CHECK_THROW(mesh_factory(x), Errors::Message);
    else {
      mesh = mesh_factory(x);
      CHECK(mesh);
      sph_vol = mesh->cell_volume(c);
      mesh.reset();
    }
    mesh_factory.reset_options();

    if (!parallel)
      CHECK(sph_vol != cart_vol);

      
    if (parallel) {
      if (Jali::framework_available(Jali::MSTK) ||
          Jali::framework_available(Jali::STKMESH)) {

        if (Jali::framework_available(Jali::MSTK))
          mesh_factory.framework(Jali::MSTK);
        if (Jali::framework_available(Jali::STKMESH))
          mesh_factory.framework(Jali::STKMESH);

        // Make sure we have ghost cells if we asked for them
        mesh_factory.num_ghost_layers_distmesh(1);
        mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
        CHECK(mesh);
        
        CHECK(mesh->num_cells<Jali::Entity_type::PARALLEL_GHOST>());
        mesh.reset();
        
        // Make sure we don't have ghost cells if we didn't ask for them
        mesh_factory.num_ghost_layers_distmesh(0);
        mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
        CHECK(mesh);
        
        CHECK(mesh->num_cells<Jali::Entity_type::PARALLEL_GHOST>() == 0);
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
        // took place along the Z-direction

        for (auto const& f : mesh->faces()) {
          JaliGeometry::Point fnormal = mesh->face_normal(f);
          if (fnormal[2] == 0.0) continue;  // z-component of normal is 0
          
          Jali::Entity_ID_List fregs;
          mesh->face_get_cells(f, Jali::Entity_type::PARALLEL_OWNED, &fregs);

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
