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
/*
This is the geometry component of the Jali code.
Authors: Rao Garimella
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef   JALI_GEOMETRY_POINT_HH_
#define   JALI_GEOMETRY_POINT_HH_

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

namespace JaliGeometry {

class Point {
 public:
  Point() {
    d = 0;
    xyz[0] = xyz[1] = xyz[2] = 0.0;
  }
  Point(const Point& p) {
    d = p.d;
    std::copy(p.xyz, p.xyz+d, xyz);
  }
  explicit Point(const int N) {
    d = N;
    xyz[0] = xyz[1] = xyz[2] = 0.0;
  }
  Point(const double& x, const double& y) {
    d = 2;
    xyz[0] = x;
    xyz[1] = y;
  }
  Point(const double& x, const double& y, const double& z) {
    d = 3;
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }
  ~Point() {}

  // main members

  // Not necessary - the constructor and set functions do the
  // equivalent things
  //
  //  void init(const int N) {
  //    d = N;
  //    xyz[0] = xyz[1] = xyz[2] = 0.0;
  //  }

  void set(const double& val) {
    assert(this);
    assert(d > 0);
    for (int i = 0; i < d; i++) xyz[i] = val;
  }
  void set(const Point& p) {
    d = p.d;
    std::copy(p.xyz, p.xyz+d, xyz);
  }
  void set(const double *val) {
    assert(val);
    assert(d > 0);
    std::copy(val, val+d, xyz);
  }
  void set(const int N, const double *val) {
    assert(val);
    d = N;
    std::copy(val, val+d, xyz);
  }
  void set(const double& x, const double& y) {
    d = 2;
    xyz[0] = x;
    xyz[1] = y;
  }
  void set(const double& x, const double& y, const double& z) {
    d = 3;
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }

  int is_valid() { return (d == 1 || d == 2 || d == 3) ? 1 : 0; }

  // access members
  double& operator[] (const int i) { return xyz[i]; }
  const double& operator[] (const int i) const { return xyz[i]; }

  double x() const { return xyz[0]; }
  double y() const { return (d > 1) ? xyz[1] : 0.0; }
  double z() const { return (d > 2) ? xyz[2] : 0.0; }

  int dim() const { return d; }

  // operators
  Point& operator=(const Point& p) {
    d = p.d;
    std::copy(p.xyz, p.xyz+d, xyz);
    return *this;
  }

  Point& operator+=(const Point& p) {
    for (int i = 0; i < d; i++) xyz[i] += p[i];
    return *this;
  }
  Point& operator-=(const Point& p) {
    for (int i = 0; i < d; i++) xyz[i] -= p[i];
    return *this;
  }
  Point& operator*=(const double& c) {
    for (int i = 0; i < d; i++) xyz[i] *= c;
    return *this;
  }
  Point& operator/=(const double& c) {
    for (int i = 0; i < d; i++) xyz[i] /= c;
    return *this;
  }

  friend Point operator*(const double& r, const Point& p) {
    if (p.d == 1) {
      Point newp = Point(p.d);
      newp.set(r*p[0]);
      return newp;
    }
    else if (p.d == 2)
      return Point(r*p[0], r*p[1]);
    else
      return Point(r*p[0], r*p[1], r*p[2]);
  }
  friend Point operator*(const Point& p, const double& r) { return r*p; }
  friend double operator*(const Point& p, const Point& q) {
    double s = 0.0;
    for (int i = 0; i < p.d; i++ ) s += p[i]*q[i];
    return s;
  }

  friend Point operator/(const Point& p, const double& r) { return p * (1.0/r); }

  friend Point operator+(const Point& p, const Point& q) {
    if (p.d == 1) {
      Point newp = Point(p.d);
      newp.set(p[0]+q[0]);
      return newp;
    }
    else if (p.d == 2)
      return Point(p[0]+q[0], p[1]+q[1]);
    else
      return Point(p[0]+q[0], p[1]+q[1], p[2]+q[2]);
  }
  friend Point operator-(const Point& p, const Point& q) {
    if (p.d == 1) {
      Point newp = Point(p.d);
      newp.set(p[0]-q[0]);
      return newp;
    }
    else if (p.d == 2)
      return Point(p[0]-q[0], p[1]-q[1]);
    else
      return Point(p[0]-q[0], p[1]-q[1], p[2]-q[2]);
  }
  friend Point operator-(const Point& p) {
    if (p.d == 1) {
      Point newp = Point(p.d);
      newp.set(-p[0]);
      return newp;
    }
    else if (p.d == 2)
      return Point(-p[0], -p[1]);
    else
      return Point(-p[0], -p[1], -p[2]);
  }

  friend Point operator^(const Point& p, const Point& q) {
    Point pq(p.d);
    if (p.d == 2) {
      pq[0] = p[0] * q[1] - q[0] * p[1];
    } else if (p.d == 3) {
      pq[0] = p[1] * q[2] - p[2] * q[1];
      pq[1] = p[2] * q[0] - p[0] * q[2];
      pq[2] = p[0] * q[1] - p[1] * q[0];
    }
    return pq;
  }

  friend std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << p.x();
    if (p.d > 1) os << " " << p.y();
    if (p.d > 2) os << " " << p.z();
    return os;
  }

 private:
  int d;
  double xyz[3];
};  // end class Point

/* miscellaneous */
inline double L22(const Point& p) { return p*p; }
inline double norm(const Point& p) { return sqrt(p*p); }

typedef std::vector<Point> Point_List;

}  // namespace JaliGeometry

#endif

