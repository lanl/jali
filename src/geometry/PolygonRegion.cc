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
 * @file   PolygonRegion.cc
 * @author Rao Garimella
 * @date Fri Jul 29 12:28:10 2011
 *
 * @brief  Implementation of PolygonRegion class
 *
 *
 */

#include "PolygonRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace JaliGeometry {

// -------------------------------------------------------------
//  class PolygonRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// Polygon:: constructors / destructor
// -------------------------------------------------------------
PolygonRegion::PolygonRegion(const std::string name, const unsigned int id,
                             const unsigned int num_points,
                             const std::vector<Point>& points,
                             const LifeCycle_type lifecycle)
  : Region(name,id,points[0].dim(),lifecycle), num_points_(num_points),
    points_(points),normal_(points[0].dim()),elim_dir_(0)
{
  init();
}

PolygonRegion::PolygonRegion(const PolygonRegion& old)
  : Region(old), num_points_(old.num_points_), points_(old.points_),
    normal_(old.normal_), elim_dir_(old.elim_dir_)
{
  // empty
}

PolygonRegion::~PolygonRegion(void)
{

}

void PolygonRegion::init() {

  if (num_points_ < dimension()) {
    std::stringstream tempstr;
    tempstr << "\nDimension " << dimension() <<
      " regions need to be specified by at least " << dimension() <<
      " points\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::Jali_throw(mesg);
  }

//  if (dimension() == 2 && num_points_ > 2) {
//      (verbobj->os()) << "\nDimension " << dimension() <<
//        " regions specified by more points (" << num_points_ << ") " <<
//        "than needed\n" << "Using only the first two\n";
//    }
//  }

  if (dimension() == 2) {
    Point vec = points_[1] - points_[0];
    vec /= norm(vec);
    normal_.set(vec[1],-vec[0]);

    elim_dir_ = (vec[0] < vec[1]) ? 0 : 1;
  }
  else if (dimension() == 3) {
    Point vec0 = points_[2]-points_[1];
    Point vec1 = points_[0]-points_[1];
    normal_ = vec0^vec1;
    normal_ /= norm(normal_);

#ifdef ENABLE_DBC
    for (int i = 3; i < num_points_; i++) {
      vec0 = points_[(i+1)%num_points_]-points_[i];
      vec1 = points_[(i-1+num_points_)%num_points_]-points_[i];
      Point nrml = vec0^vec1;
      nrml /= norm(nrml);

      double dp = nrml*normal_;
      if (fabs(dp-1.0) > 1.0e-06) {
        Errors::Message mesg("Polygon region is not exactly planar");
        Exceptions::Jali_throw(mesg);
      }
    }
#endif

    /* Determine which direction to eliminate while projecting to one
       of the coordinate planes - to do this we have to find the
       direction in which the polygon is the smallest or in other
       words, the direction in which the normal to the polygon is the
       largest */

    int dmax = -1;
    double maxlen = -1;
    for (int i = 0; i < 3; i++)
      if (normal_[i] > maxlen) {
        maxlen = normal_[i];
        dmax = i;
      }

    elim_dir_ = dmax;
  }
  else {
    std::stringstream tempstr;
    tempstr << "Cannot handle regions of dimension " << dimension() << "\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::Jali_throw(mesg);
  }
}

// -------------------------------------------------------------
// PolygonRegion::inside
// -------------------------------------------------------------
bool
PolygonRegion::inside(const Point& p) const
{

#ifdef ENABLE_DBC
  if (p.dim() != points_[0].dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in corner dimension of Polygon \"" << Region::name() << "\" and query point.\n Perhaps the region is improperly defined?\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::Jali_throw(mesg);
  }
#endif

  /* First check if the point is on the infinite line/plane */

  double d(0.0), res(0.0);

  for (int i = 0; i < p.dim(); ++i) {
    res += normal_[i]*p[i];
    d += normal_[i]*points_[0][i];
  }
  res -= d;

  if (fabs(res) > 1.0e-12)
    return false;


  bool result(false);
  if (dimension() == 2) {
    // Now check if it lies in the line segment

    // vector from start of segment to point
    Point vec0 = p-points_[0];

    // segment vector
    Point vec1 = points_[1]-points_[0];

    // Normalize
    double slen = norm(vec1);
    vec1 /= slen;

    double dp = vec0*vec1;

    // projection of vec0 along segment lies inside the segment
    if (dp >= 0 && dp <= slen) {

      // projected point along line segment
      Point p1 = points_[0] + dp*vec1;

      // vector between point and its projection
      Point dvec = p - p1;

      // distance between point and its projection
      double d_sqr = L22(dvec);

      // Is the distance 0? Point is inside segment
      if (d_sqr < 1.0e-16)
        result = true;
    }
  }
  else {
    /* Now check if the point is in the polygon */

    /*
       Find the indices of coordinates on the projection plane

       if elim_dir_ is 0, then d0 = 1, d1 = 2 (YZ plane)
       if elim_dir_ is 1, then d0 = 2, d1 = 0 (XZ plane)
       if elim_dir_ is 2, then d0 = 0, d1 = 1 (XY plane)
    */

    double d0 = (elim_dir_+1)%3;
    double d1 = (elim_dir_+2)%3;

    /* Now apply the Jordan curve theorem to do the in/out test */
    /* odd number of crossings - point is inside                */

    double u, v;
    u = p[d0]; v = p[d1];

    for (int i = 0; i < num_points_; i++) {
      int iplus1 = (i+1)%num_points_;
      double u_i = points_[i][d0];
      double v_i = points_[i][d1];
      double u_iplus1 = points_[iplus1][d0];
      double v_iplus1 = points_[iplus1][d1];

      // don't compute - v_iplus1-v_i could be zero
      //      double slope = (u_iplus1-u_i)/(v_iplus1-v_i);

      if (((v_i > v && v_iplus1 <= v) || (v_iplus1 > v && v_i <= v)) &&
          (u <= (u_i + (v-v_i)*(u_iplus1-u_i)/(v_iplus1-v_i))))
        result = !result;
    }

    // The above check is not guaranteed to give an +ve result if the point is on
    // the boundary. So do an additional check

    if (!result) {

      for (int i = 0; i < num_points_; i++) {

        int iplus1 = (i+1)%num_points_;

        Point p_i(2);
        p_i.set(points_[i][d0],points_[i][d1]);

        Point p_iplus1(2);
        p_iplus1.set(points_[iplus1][d0],points_[iplus1][d1]);

        // vector from first point of segment to query point
        Point vec0(2);
        vec0.set(u-p_i[0],v-p_i[1]);

        // line segment vector
        Point vec1(2);
        vec1 = p_iplus1 - p_i;

        // unit vector along segment
        double slen = norm(vec1);
        vec1 /= slen;

        double dp = vec0*vec1;

        // projection of vec0 along segment lies outside the segment
        if (dp < 0 || dp > slen)
          continue;

        // projected point along line segment
        Point p1(2);
        p1 = p_i + dp*vec1;

        // vector between projected point and query point
        Point dvec(2);
        dvec.set(u-p1[0],v-p1[1]);

        // distance between point and its projection
        double d_sqr = L22(dvec);

        // Is the distance 0? Point is inside segment
        if (d_sqr < 1.0e-16) {
          result = true;
          break;
        }
      }
    }
  }

  return result;
}


} // namespace JaliGeometry
