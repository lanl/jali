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


/**
 * @file   Domain.hh
 * @author Rao Garimella
 * @date   Mon Aug  1 09:57:42 2011
 *
 * @brief  Declaration of the Domain class
 *
 *
 */

#ifndef _Domain_hh_
#define _Domain_hh_

#include <vector>

#include "Region.hh"
#include "GeometricModel.hh"

namespace JaliGeometry {

// -------------------------------------------------------------
//  class Domain
// -------------------------------------------------------------
/// A class that represent the entire geometric domain
/**
 * The geometric domain contains one or more objects describing the specifics
 * of how the domain is decomposed into subdomains. In its simplest form,
 * the geometric domain is a list of geometric regions (see
 * JaliGeometry::Region) that may or may not tile the domain with or
 * without overlap. In its most sophisticated form, the domain has a list
 * of geometric models, each of which is a particular non-overlapping, tiling
 * decomposition of the domain into geometric regions. Other free floating
 * regions may also exist in the domain for purposes like post-processing
 */

class Domain {
public:

  // Default constructor.

  Domain(const unsigned int dim);

  // Copy constructor

  Domain(const Domain& old);


  // Constructor with Geometric Model List

  Domain(const unsigned int dim,
         const std::vector<GeometricModelPtr>& in_geometric_models,
         const std::vector<RegionPtr>& in_Regions);

  // Destructor

  virtual ~Domain(void);


  inline
  unsigned int spatial_dimension() const
  {
    return spatial_dimension_;
  }


  // Add a Geometric Model

  void Add_Geometric_Model(const GeometricModelPtr& gm);


  // Add a Free Region

  void Add_Free_Region(const RegionPtr& regptr);


  // Number of Geometric Models

  int Num_Geometric_Models(void) const;


  // Get the i'th Geometric Model

  GeometricModelPtr Geometric_Model_i(const int i) const;


  // Number of Free Regions

  int Num_Free_Regions(void) const;


  // Get the i'th Free Region

  RegionPtr Free_Region_i(const int i) const;

private:

  // Dimension of domain

  unsigned int spatial_dimension_;


  // List of geometric models in domain

  std::vector<GeometricModelPtr> GeometricModels;


  // List of Free Region pointers (regions that are not part of any
  // geometric model)

  std::vector<RegionPtr> FreeRegions;

};

} // namespace JaliGeometry

#endif

