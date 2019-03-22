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
 * @file   BoxRegion.hh
 * @author Rao Garimella, William A. Perkins
 * @date Wed Sep 28 08:54:19 2011
 *
 * @brief  Declaration of BoxRegion class (adapted from RectangularRegion)
 *
 *
 */

#ifndef _BoxRegion_hh_
#define _BoxRegion_hh_

#include "Region.hh"
#include "dbc.hh"
#include "errors.hh"

namespace JaliGeometry {

// -------------------------------------------------------------
//  class BoxRegion
// -------------------------------------------------------------
/// A rectangular region in space, defined by two points

class BoxRegion : public Region {
public:

  /// Default constructor uses two corner points (order not important).

  BoxRegion(const std::string name, const unsigned int id, const Point& p0,
            const Point& p1,
            const LifeCycle_type lifecycle = LifeCycle_type::PERMANENT);

  /// Protected copy constructor to avoid unwanted copies.
  BoxRegion(const BoxRegion& old);

  /// Destructor
  ~BoxRegion(void);


  // Type of the region
  inline
  Region_type type() const { return Region_type::BOX; }

  /// Get the first point defining the region
  inline
  const Point& point0(void) const { return p0_; }

  /// Get the second point defining the region
  inline
  const Point& point1(void) const { return p1_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

  /// corners
  inline
  void corners(Point *lo_corner, Point *hi_corner) const
  {
    if (lo_corner == NULL) {
      Errors::Message mesg("lo_corner not specified");
      Exceptions::Jali_throw(mesg);
    }
    if (hi_corner == NULL) {
      Errors::Message mesg("hi_corner not specified");
      Exceptions::Jali_throw(mesg);
    }
    //    lo_corner->init(p0_.dim());
    lo_corner->set(p0_);
    //    hi_corner->init(p1_.dim());
    hi_corner->set(p1_);
  }

  // Is the box degenerate - zero length in one or more directions and
  // if so in how many directions?
  bool is_degenerate(int *ndeg) const;

protected:

  const Point p0_;              /**< one corner of the region  */
  const Point p1_;              /**< the other corner of the region */

  /// Is the specified value between the two values (inclusive, order not important)
  static bool between_(const double& x, const double& x0, const double& x1);

};

/// A smart pointer to BoxRegion instances
//
// typedef Teuchos::RCP<BoxRegion> BoxRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef BoxRegion *BoxRegionPtr;

} // namespace JaliGeometry


#endif
