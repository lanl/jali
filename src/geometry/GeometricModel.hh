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
 * @file   GeometricModel.hh
 * @author Rao Garimella
 * @date   Sep 15, 2011
 *
 * @brief  Declaration of the GeometricModel class
 *
 *
 */

#ifndef _GeometricModel_hh_
#define _GeometricModel_hh_

#include <vector>

#include "Region.hh"
//#include "RegionFactory.hh"

namespace JaliGeometry {

// -------------------------------------------------------------
//  class GeometricModel
// -------------------------------------------------------------
// A class that represent a geometric model or more specifically,
// a particular decomposition of the domain into subdomains
/**
 * The geometric model is an object that contains a list of
 * geometric regions that tile the domain (no gaps, no overlaps)
 **/

class GeometricModel {
public:

  // constructor.

  GeometricModel(const unsigned int dim);

  // Copy constructor

  GeometricModel(const GeometricModel& old);


  // Constructor from parameter list

//  GeometricModel(const unsigned int dim, Teuchos::ParameterList gm_param_list,
//                 const Epetra_MpiComm *comm, const VerboseObject *verbobj=NULL);


  // Constructor from a list of regions

  GeometricModel(const unsigned int dim,
                 const std::vector<RegionPtr>& in_Regions);


  // Destructor

  ~GeometricModel(void);


  // Topological Dimension of geometric model

  inline
  unsigned int dimension() const {
    return topo_dimension_;
  }


  // Add a Region to a GeometricModel

  void Add_Region(const RegionPtr& r);


  // Number of Regions

  int Num_Regions(void) const;


  // Get the i'th region of the model

  RegionPtr Region_i(const int i) const;


  // Get a region by its ID
  RegionPtr FindRegion(const int id) const;


  // Get a region by its ID
  RegionPtr FindRegion(const std::string name) const;


  // Check if regions cover the domain extents.
  // This will work perfectly for domains with rectangular regions
  // but not so for other types of regions

  bool Rough_Check_Tiling(void) const;

private:

  // Topological dimension of the model

  unsigned int topo_dimension_;

  // List of regions in this geometric model

  std::vector<RegionPtr> Regions;

};


  typedef GeometricModel *GeometricModelPtr;

} // namespace JaliGeometry

#endif

