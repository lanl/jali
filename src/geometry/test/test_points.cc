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
#include <UnitTest++.h>

#include <iostream>


#include "../Point.hh"

#include "mpi.h"


TEST(Point)
{
  double x = 0.4, y = 0.6, z = 0.8;
  double xyz[3] = {1, 2, 3.5};

  // Create different types of points

  JaliGeometry::Point p0(2);
  p0.set(x,y);
  CHECK_EQUAL(x,p0.x());
  CHECK_EQUAL(y,p0.y());

  JaliGeometry::Point p1(2);
  p1.set(2,xyz);
  CHECK_EQUAL(xyz[0],p1.x());
  CHECK_EQUAL(xyz[1],p1.y());


  x = x + 0.1; y = y + 0.1;
  JaliGeometry::Point p2(x,y);
  CHECK_EQUAL(x,p2.x());
  CHECK_EQUAL(y,p2.y());

  JaliGeometry::Point p3(p2);
  CHECK_EQUAL(x,p3.x());
  CHECK_EQUAL(y,p3.y());

  JaliGeometry::Point p4 = p2;
  CHECK_EQUAL(x,p4.x());
  CHECK_EQUAL(y,p4.y());

  JaliGeometry::Point p5(3);
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;
  p5.set(3,xyz);
  CHECK_EQUAL(x,p5.x());
  CHECK_EQUAL(y,p5.y());
  CHECK_EQUAL(z,p5.z());


  // Create a bunch of 3D points and delete some of them to trigger
  // different patterns of deletion in the underlying coordinate
  // buffer

  int num3 = 25;
  JaliGeometry::Point *points3[25];
  for (int i = 0; i < num3; ++i) {
    points3[i] = new JaliGeometry::Point(x,y,z);
  }

  // Now delete some points in 3D - point 3, 4, 5, 8, 7, 16, 18, 17
  // Deletion of 3 is a stand alone deletion
  // Deletion of 4 followed by point 5 will cause merging of two holes
  // Deletion of 8 followed by point 7 will cause merging of two holes
  // Deletion of 16, 18, 17 will cause bridging of three holes

  // no obvious way of checking the underlying coordinate buffer
  // machinery except through a debugger

  delete points3[2];

  delete points3[4];
  delete points3[5];

  delete points3[8];
  delete points3[7];

  delete points3[16];
  delete points3[18];
  delete points3[17];



  // Create a bunch of mixed dim points

  int num = 25;
  std::vector<JaliGeometry::Point *> points;
  points.resize(25);

  for (int i = 0; i < num; ++i) {
    if (i%2)
      points[i] = new JaliGeometry::Point(x,y);
    else {
      points[i] = new JaliGeometry::Point(3);
      points[i]->set(x,y,z);
    }
    x = x + 1.0;
    y = y + 0.5;
    z = z + 2.0;
  }


  // Don't bother deleting other points. They will get deleted when exiting


  // we should test the mathematical operations as well here

}

