// -------------------------------------------------------------
/**
 * @file   test_mesh_framework.cc
 * @author William A. Perkins
 * @date Tue Oct  4 06:15:36 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 14, 2011 by William A. Perkins
// Last Change: Tue Oct  4 06:15:36 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <algorithm>
#include <UnitTest++.h>

#include "dbc.hh"
#include "../MeshFramework.hh"
#include "../FrameworkTraits.hh"

SUITE (Framework) 
{
  TEST (DefaultPreference) {
    
    Jali::FrameworkPreference pref(Jali::default_preference());

    CHECK(std::find(pref.begin(), pref.end(), Jali::Simple) != pref.end());
    
#ifdef HAVE_MOAB_MESH
    CHECK(std::find(pref.begin(), pref.end(), Jali::MOAB) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), Jali::MOAB) == pref.end());
#endif

#ifdef HAVE_STK_MESH
    CHECK(std::find(pref.begin(), pref.end(), Jali::STKMESH) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), Jali::STKMESH) == pref.end());
#endif

#ifdef HAVE_MSTK_MESH
    CHECK(std::find(pref.begin(), pref.end(), Jali::MSTK) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), Jali::MSTK) == pref.end());
#endif

  }

  TEST (AvailablePreference) 
  {
    Jali::FrameworkPreference pref;
    
    pref.clear();
    pref.push_back(Jali::MOAB);
    pref = Jali::available_preference(pref);
#ifdef HAVE_MOAB_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
    pref.clear();
    pref.push_back(Jali::STKMESH);
    pref = Jali::available_preference(pref);
#ifdef HAVE_STK_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
    pref.clear();
    pref.push_back(Jali::MSTK);
    pref = Jali::available_preference(pref);
#ifdef HAVE_MSTK_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
  }

  TEST (Readability)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    bool parallel(nproc > 1);
    
    CHECK(!Jali::framework_reads(Jali::Simple, Jali::ExodusII, parallel));
    CHECK(!Jali::framework_reads(Jali::Simple, Jali::Nemesis, parallel));
    CHECK(!Jali::framework_reads(Jali::Simple, Jali::MOABHDF5, parallel));

    CHECK(Jali::framework_reads(Jali::MOAB, Jali::ExodusII, parallel));
    CHECK(!Jali::framework_reads(Jali::MOAB, Jali::Nemesis, parallel));
    CHECK(Jali::framework_reads(Jali::MOAB, Jali::MOABHDF5, parallel));

    if (parallel) {
      CHECK(!Jali::framework_reads(Jali::STKMESH, Jali::ExodusII, parallel));
      CHECK(Jali::framework_reads(Jali::STKMESH, Jali::Nemesis, parallel));
    } else {
      CHECK(Jali::framework_reads(Jali::STKMESH, Jali::ExodusII, parallel));
      CHECK(!Jali::framework_reads(Jali::STKMESH, Jali::Nemesis, parallel));
    }
    CHECK(!Jali::framework_reads(Jali::STKMESH, Jali::MOABHDF5, parallel));

  }

  TEST (Generatability)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    bool parallel(nproc > 1);
    
    CHECK(!Jali::framework_generates(Jali::MOAB, parallel,3));
    CHECK(Jali::framework_generates(Jali::MSTK, parallel,3));
    CHECK(Jali::framework_generates(Jali::STKMESH, parallel,3));
    if (parallel) {
      CHECK(!Jali::framework_generates(Jali::Simple, parallel,3));
    } 


    //    CHECK(!Jali::framework_generates(Jali::MOAB, parallel,2));
    //    CHECK(Jali::framework_generates(Jali::MSTK, parallel,2));
    //    CHECK(!Jali::framework_generates(Jali::STKMESH, parallel,2));
    //    if (parallel) {
    //      CHECK(!Jali::framework_generates(Jali::Simple, parallel,2));
    //    } 

  }

}
