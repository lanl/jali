// -------------------------------------------------------------
/**
 * @file   test_mesh_file.cc
 * @author William A. Perkins
 * @date Mon May 16 14:03:23 2011
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 14, 2011 by William A. Perkins
// Last Change: Mon May 16 14:03:23 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <iostream>
#include <UnitTest++.h>

#include "dbc.hh"
#include "../MeshFileType.hh"
#include "../MeshException.hh"

SUITE (MeshFileType)
{
  TEST (ExodusII)
  {
    // EXODUS_TEST_FILE is macro defined by cmake
    std::string fname(EXODUS_TEST_FILE);

    Jali::Format f;
    try {
      f = Jali::file_format(MPI_COMM_WORLD, fname);
    } catch (const Jali::Message& e) {
      throw e;
    }

    CHECK(f == Jali::ExodusII);
  }

  TEST (Nemesis)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

    // NEMESIS_TEST_FILE is macro defined by cmake
    std::string fname(NEMESIS_TEST_FILE);

    Jali::Format f;
    if (nproc > 1 && nproc <= 4) {
      int ierr[1];
      ierr[0] = 0;
      try {
        f = Jali::file_format(MPI_COMM_WORLD, fname);
      } catch (const Jali::Message& e) {
        throw e;
      }

      CHECK(f == Jali::Nemesis);
    } else {
      CHECK_THROW(f = Jali::file_format(MPI_COMM_WORLD, fname),
                  Jali::FileMessage);
    }
  }

  TEST (MOABHD5)
  {
    // MOAB_TEST_FILE is macro defined by cmake
    std::string fname(MOAB_TEST_FILE);

    Jali::Format f;
    try {
      f = Jali::file_format(MPI_COMM_WORLD, fname);
    } catch (const Jali::Message& e) {
      throw e;
    }

    CHECK(f == Jali::MOABHDF5);
  }

  TEST (PathFailure)
  {
    std::string fname("/some/bogus/path.exo");

    Jali::Format f;

    CHECK_THROW(Jali::file_format(MPI_COMM_WORLD, fname),
                Jali::FileMessage);
  }

  TEST (MagicNumberFailure)
  {
    std::string fname(BOGUS_TEST_FILE);

    Jali::Format f;

    CHECK_THROW(Jali::file_format(MPI_COMM_WORLD, fname),
                Jali::FileMessage);
  }

}

