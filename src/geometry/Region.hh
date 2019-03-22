/*
 Copyright (c) 2019, Triad National Security, LLC
 All rights reserved.

 Copyright 2019. Triad National Security, LLC. This software was
 produced under U.S. Government contract 89233218CNA000001 for Los
 Alamos National Laboratory (LANL), which is operated by Triad
 National Security, LLC for the U.S. Department of Energy. 
 All rights in the program are reserved by Triad National Security,
 LLC, and the U.S. Department of Energy/National Nuclear Security
 Administration. The Government is granted for itself and others acting
 on its behalf a nonexclusive, paid-up, irrevocable worldwide license
 in this material to reproduce, prepare derivative works, distribute
 copies to the public, perform publicly and display publicly, and to
 permit others to do so

 
 This is open source software distributed under the 3-clause BSD license.
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:
 
 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
 3. Neither the name of Triad National Security, LLC, Los Alamos
    National Laboratory, LANL, the U.S. Government, nor the names of its
    contributors may be used to endorse or promote products derived from this
    software without specific prior written permission.

 
 THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND
 CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/**
 * @file   Region.hh
 * @author William A. Perkins
 * @date Mon Aug  1 09:57:42 2011
 *
 * @brief  Declaration of the abstract Region class
 *
 *
 */

#ifndef _Region_hh_
#define _Region_hh_

#include <vector>

#include "Point.hh"

namespace JaliGeometry {


enum class Region_type {
    BOX = 0,
    PLANE,
    LABELEDSET,
    LAYER,
    SURFACE,
    POINT,
    COLORFUNCTION,
    LOGICAL,
    POLYGON
};

enum class LifeCycle_type {PERMANENT = 0, TEMPORARY};

enum class Bool_type {
    NOBOOLEAN = -1,
    COMPLEMENT,
    UNION,
    INTERSECT,
    SUBTRACT
};

// -------------------------------------------------------------
//  class Region
// -------------------------------------------------------------
/// An class that represent a geometric region
/**
 * A Region is just some arbitrary subset of space, that can be
 * specified in a myriad of ways.  At a minimum, there is a need to be
 * able to determine if a point is inside that space.  Other needs to
 * be added later.
 *
 * The region class does not use a constructor based on the XML parameter
 * list because it has to create derived region classes based on the shape
 * parameter of the region specification.
 *
 *
 */

class Region {
public:

  /// Default constructor.
  Region(void);

  /// Constructor with name and ID
  Region(const std::string name, const unsigned int id,
         const unsigned int dim = 3,
         const LifeCycle_type lifecycle = LifeCycle_type::PERMANENT);

  /// Copy constructor
  Region(const Region& old);

  /// Destructor
  virtual ~Region(void);


  /// Set the dimension of the region
  inline
  void set_dimension(const unsigned int dim)
  {
    topo_dimension_ = dim;
  }

  /// Name of the region
  inline
  std::string name(void) const
  {
    return name_;
  }

  /// Integer identifier of the region
  inline
  unsigned int id(void) const
  {
    return id_;
  }

  // Topological dimension of region (0 - point, 1 - curve, 2 - surface, 3 - volume)
  inline
  unsigned int dimension(void) const
  {
    return topo_dimension_;
  }

  // Get the Lifecycle of this region - Do mesh entity sets derived from
  // it have to be kept around or are they temporary and can be destroyed
  // as soon as they are used?

  inline
  LifeCycle_type lifecycle(void) const
  {
    return lifecycle_;
  }

  // Type of the region
  virtual Region_type type() const = 0;

  /// Is the specified point inside the Region
  /// Does being on the boundary count as inside or not?
  virtual bool inside(const Point& p) const = 0;


  /// Get the extents of the Region
  /// void extents(Point *pmin, Point *pmax) const;

private:

  // Lifecycle (Temporary or Permanent)

  LifeCycle_type lifecycle_;

  // Topological dimension of region (0, 1, 2, 3)

  unsigned int topo_dimension_;

  // Name of identifier

  std::string name_;

  // Integer identifier of region

  unsigned int id_;

  // lower/left/back corner or min xyz

  // Point min_pnt;

  // Upper/right/front corner or max xyz

  // Point max_pnt;

};

typedef Region *RegionPtr;

// A thing to hold some Region

typedef std::vector< RegionPtr > RegionVector;

} // namespace JaliGeometry

#endif

