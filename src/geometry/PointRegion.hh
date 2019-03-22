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
 * @file   PointRegion.hh
 * @author Rao Garimella
 * @date
 *
 * @brief  Declaration of PointRegion class
 *
 *
 */

#ifndef _PointRegion_hh_
#define _PointRegion_hh_

#include <string>

#include "Region.hh"


namespace JaliGeometry {

// -------------------------------------------------------------
//  class PointRegion
// -------------------------------------------------------------
/// A point in space

class PointRegion : public Region {
public:

  PointRegion(const std::string name, const unsigned int id, const Point& p,
              const LifeCycle_type lifecycle = LifeCycle_type::PERMANENT);

  /// Protected copy constructor to avoid unwanted copies.
  PointRegion(const PointRegion& old);

  /// Destructor
  ~PointRegion(void);

  // Type of the region
  inline Region_type type() const { return Region_type::POINT; }

  /// Get the point defining the region
  const Point& point(void) const { return p_; }

  /// Is the specified point inside this region - in this case it
  /// means coincident with the region point

  bool inside(const Point& p) const;

protected:

  const Point p_;              /* point */

};

/// A smart pointer to PointRegion instances
// typedef Teuchos::RCP<PointRegion> PointRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef PointRegion *PointRegionPtr;

} // namespace JaliGeometry


#endif
