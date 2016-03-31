#ifndef _JALI_GEOM_H
#define _JALI_GEOM_H


// using namespace std;
#include <vector>

#include "Point.hh"


namespace JaliGeometry {

// Mesh Geometry Type

enum class Geom_type {
  CARTESIAN,
  CYLINDRICAL,
  SPHERICAL
};

// Return the volume and centroid of a general polyhedron
//
// ccoords  - vertices of the polyhedron (in no particular order)
// nf       - number of faces of polyhedron
// nfnodes  - number of nodes for each face
// fcoords  - linear array of face coordinates in in ccw manner
//            assuming normal of face is pointing out (
//
// So if the polyhedron has 5 faces with 5,3,3,3 and 3 nodes each
// then entries 1-5 of fcoords describes face 1, entries 6-9
// describes face 2 and so on
//
// So much common work has to be done for computing the centroid
// and volume calculations that they have been combined into one
//
// The volume of all polyhedra except tets is computed as a sum of
// volumes of tets created by connecting the polyhedron center to
// a face center and an edge of the face

void polyhed_get_vol_centroid(const std::vector<Point> ccoords,
                              const unsigned int nf,
                              const std::vector<unsigned int> nfnodes,
                              const std::vector<Point> fcoords,
                              double *volume,
                              Point *centroid);

// Is point in polyhed

bool point_in_polyhed(const Point testpnt,
                      const std::vector<Point> ccoords,
                      const unsigned int nf,
                      const std::vector<unsigned int> nfnodes,
                      const std::vector<Point> fcoords);

// Compute area, centroid and normal of polygon

// In 2D, the area is computed by a contour integral around the
// perimeter. In 3D, the area is computed by connecting a
// "center" point of the polygon to the edges of the polygon and
// summing the areas of the resulting triangles
//
// The normal of a 3D polygon is computed as the sum of the area
// weighted normals of the triangular facets

void polygon_get_area_centroid_normal(const std::vector<Point> coords,
                                      double *area, Point *centroid,
                                      Point *normal);

// Get area weighted normal of polygon
// In 2D, the normal is unambiguous - the normal is evaluated at one corner
// In 3D, the procedure evaluates the normal at each corner and averages it

//  Point polygon_get_normal(const std::vector<Point> coords);


// Is point in polygon

bool point_in_polygon(const Point testpnt,
                      const std::vector<Point> coords);

// Compute volume and centroid of 1d segment, accounting for geometry
void segment_get_vol_centroid(const std::vector<Point> ccoords,
                              Geom_type my_geom_type,
                              double *volume, Point* centroid);

// Compute the face area in a 1d mesh
void face1d_get_area(const std::vector<Point> fcoords,
                     Geom_type my_geom_type,
                     double *area);

}  // namespace JaliGeometry


#endif
