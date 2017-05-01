/*
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.

Copyright 2017. Los Alamos National Security, LLC. This software was
produced under U.S. Government contract DE-AC52-06NA25396 for Los
Alamos National Laboratory (LANL), which is operated by Los Alamos
National Security, LLC for the U.S. Department of Energy. The
U.S. Government has rights to use, reproduce, and distribute this
software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so
as not to confuse it with the version available from LANL.
 
Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

1.  Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3.  Neither the name of Los Alamos National Security, LLC, Los Alamos
National Laboratory, LANL, the U.S. Government, nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   Domain.cc
 * @author Rao V. Garimella
 * @date Mon Aug  1 10:05:25 2011
 *
 * @brief
 *
 *
 */

#include "Domain.hh"

namespace JaliGeometry {

// -------------------------------------------------------------
//  class Domain
// -------------------------------------------------------------

// -------------------------------------------------------------
// Domain:: constructors / destructor
// -------------------------------------------------------------

// Constructor

Domain::Domain(const unsigned int dim): spatial_dimension_(dim)
{
  if (dim < 2u || dim > 3u) {
    std::cerr << "Only 2D and 3D domains are supported" << std::endl;
    throw std::exception();
  }

  GeometricModels.clear();
  FreeRegions.clear();
}

// Copy constructor

Domain::Domain(const Domain& old)
{
  int i, ng, nr;

  spatial_dimension_ = old.spatial_dimension();

  ng = old.Num_Geometric_Models();
  for (i = 0; i < ng; i++) {
    GeometricModelPtr g = old.Geometric_Model_i(i);
    GeometricModels.push_back(g);
  }

  nr = old.Num_Free_Regions();
  for (i = 0; i < nr; i++) {
    RegionPtr r = old.Free_Region_i(i);
    FreeRegions.push_back(r);
  }
}

// Destructor

Domain::~Domain(void)
{
  GeometricModels.clear();
  FreeRegions.clear();
}



// Constructor with lists of geometric models and free regions

Domain::Domain(const unsigned int dim,
               const std::vector<GeometricModelPtr>& in_geometric_models,
               const std::vector<RegionPtr>& in_Regions) :
  spatial_dimension_(dim), GeometricModels(in_geometric_models),
  FreeRegions(in_Regions)
{
  if (dim != 2 || dim != 3) {
    std::cerr << "Only 2D and 3D domains are supported" << std::endl;
    throw std::exception();
  }
}



// Add a geometric model

void Domain::Add_Geometric_Model(const GeometricModelPtr& gm)
{
  // Make sure spatial dimension of domain and geometric model are the same

  if (spatial_dimension_ != gm->dimension()) {
    std::cerr << "Spatial dimension of domain and geometric model mismatch" << std::endl;
    throw std::exception();
  }

  GeometricModels.push_back(gm);
}


// Add a Free Region

void Domain::Add_Free_Region(const RegionPtr& regptr)
{
  if (spatial_dimension_ < regptr->dimension()) {
    std::cerr << "Spatial dimension of domain is less than that of the free region" << std::endl;
    throw std::exception();
  }

  FreeRegions.push_back(regptr);
}

// Number of geometric models

int Domain::Num_Geometric_Models(void) const
{
  return GeometricModels.size();
}

// Get the i'th Decomposition

GeometricModelPtr Domain::Geometric_Model_i(const int i) const
{
  return GeometricModels[i];
}

// Number of Free Regions

int Domain::Num_Free_Regions(void) const
{
  return FreeRegions.size();
}

// Get the i'th Free Region

RegionPtr Domain::Free_Region_i(const int i) const
{
  return FreeRegions[i];
}

} // namespace JaliGeometry
