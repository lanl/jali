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
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
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
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    bool parallel(nproc > 1);
        
    Jali::FrameworkPreference pref;
    std::unique_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    double x0( 0.0), y0( 0.0), z0( 0.0);
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
      mesh = NULL;
    } else {
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(mesh.get());
      mesh = NULL;
    }

    // The STK, if available, framework will always generate

    if (framework_available(Jali::STKMESH)) {
      pref.clear(); pref.push_back(Jali::STKMESH);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(mesh.get());
      mesh = NULL;
    }

    // The MSTK framework, if available, will always generate

    if (framework_available(Jali::MSTK)) {
      pref.clear(); pref.push_back(Jali::MSTK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(mesh.get());
      mesh = NULL;
    }

    // The MOAB framework cannot generate


    if (framework_available(Jali::MOAB)) {
      pref.clear(); pref.push_back(Jali::MOAB);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  Jali::Message);
      mesh = NULL;
    }
  }

  TEST (Generate2D)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    bool parallel(nproc > 1);
    
    Jali::FrameworkPreference pref;
    std::unique_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    double x0( 0.0), y0( 0.0);
    double x1(10.0), y1(10.0);
    int nx(10), ny(10);

    // The MSTK framework, if available, will always generate

    if (framework_available(Jali::MSTK)) {
      pref.clear(); pref.push_back(Jali::MSTK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0,
                          x1, y1,
                          nx, ny);
      CHECK(mesh.get());
      mesh = NULL;
    }

    // The Simple framework is always available, but 
    // cannot generate 2D meshes

    pref.clear(); pref.push_back(Jali::Simple);
    mesh_factory.preference(pref);

    CHECK_THROW(mesh = mesh_factory(x0, y0,
                                    x1, y1,
                                    nx, ny),
                Jali::Message);
    mesh = NULL;

    // The STK, even if available cannot generate 2D meshes

    if (framework_available(Jali::STKMESH)) {
      pref.clear(); pref.push_back(Jali::STKMESH);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0,
                                      x1, y1,
                                      nx, ny),
                  Jali::Message);
      mesh = NULL;
    }

    // The MOAB framework cannot generate


    if (framework_available(Jali::MOAB)) {
      pref.clear(); pref.push_back(Jali::MOAB);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0,
                                      x1, y1,
                                      nx, ny),
                  Jali::Message);
      mesh = NULL;
    }
  }


  // The Simple framework cannot read anything, even if it exists
  TEST (ReadSimple) {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    bool parallel(nproc > 1);
    
    Jali::FrameworkPreference pref;
    std::unique_ptr<Jali::Mesh> mesh;
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
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    bool parallel(nproc > 1);

    std::unique_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    if (Jali::framework_available(Jali::MOAB)) {
      mesh = mesh_factory(MOAB_TEST_FILE);
      CHECK(mesh.get());
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
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    bool parallel(nproc > 1);
    
    std::unique_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

    bool available =
        (Jali::framework_available(Jali::STKMESH) && !parallel) ||
        Jali::framework_available(Jali::MSTK);

    if (available) {
      mesh = mesh_factory(EXODUS_TEST_FILE);
      CHECK(mesh.get());
    } else {
      CHECK_THROW(mesh = mesh_factory(EXODUS_TEST_FILE),
                  Jali::Message);
    }
  }

  TEST (ReadNemesis) 
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    bool parallel(nproc > 1);
    
    std::unique_ptr<Jali::Mesh> mesh;
    Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);
    if ((Jali::framework_available(Jali::STKMESH) ||
         Jali::framework_available(Jali::MSTK)) && 
        parallel) {
      mesh = mesh_factory(NEMESIS_TEST_FILE);
      CHECK(mesh.get());
    } else {
      CHECK_THROW(mesh = mesh_factory(NEMESIS_TEST_FILE),
                  Jali::Message);
    }
  }      

}
