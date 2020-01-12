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
 * @file   BoxRegion.cc
 * @author Rao Garimella, William A. Perkins
 * @date Fri Jul 29 12:28:10 2011
 *
 * @brief  Implementation of BoxRegion class (Adapted from RectangularRegion.cc)
 *
 *
 */

#include "BoxRegion.hh"
#include "errors.hh"

namespace JaliGeometry {

// -------------------------------------------------------------
//  class BoxRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// BoxRegion:: constructors / destructor
// -------------------------------------------------------------
BoxRegion::BoxRegion(const std::string name, const unsigned int id,
                     const Point& p0, const Point& p1,
                     const LifeCycle_type lifecycle)
  : Region(name,id,p0.dim(),lifecycle), p0_(p0), p1_(p1)
{

  if (p0_.dim() != p1_.dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in dimensions of corner points of BoxRegion \"" << Region::name() << "\"\n Perhaps the region is improperly defined?\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::Jali_throw(mesg);
  }

  // Check if this is a reduced dimensionality box (e.g. even though
  // it is in 3D space it is a 2D box)

  int dim = p0.dim();
  for (int i = 0; i < p0.dim(); i++)
    if (p0[i] == p1[i]) dim--;

  if (dim < p0.dim()) Region::set_dimension(dim);
}

BoxRegion::BoxRegion(const BoxRegion& old)
  : Region(old), p0_(old.p0_), p1_(old.p1_)
{
  // empty
}

BoxRegion::~BoxRegion(void)
{

}

// -------------------------------------------------------------
// BoxRegion::between_
// -------------------------------------------------------------
bool
BoxRegion::between_(const double& x, const double& x0, const double& x1)
{
  double tol = 1.0e-08;
  return  (x0+tol >= x && x >= x1-tol) || (x1+tol >= x && x >= x0-tol);
}

// -------------------------------------------------------------
// BoxRegion::inside
// -------------------------------------------------------------
bool
BoxRegion::inside(const Point& p) const
{

  if (p.dim() != p0_.dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in corner dimension of BoxRegion \"" << Region::name() << "\" and query point.\n Perhaps the region is improperly defined?\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::Jali_throw(mesg);
  }

  bool result(true);
  for (int i = 0; i < p.dim(); ++i) {
    result = result && between_(p[i], p0_[i], p1_[i]);
  }
  return result;
}

// -------------------------------------------------------------
// BoxRegion::is_degenerate (also indicate in how many dimensions)
// -------------------------------------------------------------
bool
BoxRegion::is_degenerate(int *ndeg) const
{
  *ndeg = 0;
  for (int i = 0; i < p0_.dim(); ++i) {
    if (p0_[i] == p1_[i]) (*ndeg)++;
  }
  if (*ndeg) return true;

  return false;
}

} // namespace JaliGeometry
