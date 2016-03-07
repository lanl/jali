/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   MeshFactory.hh
 * @author William A. Perkins
 * @date Wed Sep 28 09:10:15 2011
 *
 * @brief  declaration of the MeshFactory class
 *
 *
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 10, 2011 by William A. Perkins
// Last Change: Wed Sep 28 09:10:15 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _MeshFactory_hh_
#define _MeshFactory_hh_

#include <string>
#include <vector>

#include <mpi.h>

#include "MeshException.hh"
#include "MeshFramework.hh"
#include "Mesh.hh"

#include "GeometricModel.hh"

namespace Jali {

// -------------------------------------------------------------
//  class MeshFactory
// -------------------------------------------------------------
class MeshFactory {
 protected:

  /// The parallel environment
  const MPI_Comm my_comm;

  /// A list of preferred mesh frameworks to consider
  FrameworkPreference my_preference;

 private:

  /// private, undefined copy constructor to avoid unwanted copies
  MeshFactory(MeshFactory& old);

  /// Create a mesh by reading the specified file (or set of files)
  Mesh *create(const std::string& filename,
               const JaliGeometry::GeometricModelPtr &gm =
               (JaliGeometry::GeometricModelPtr) NULL,
               const bool request_faces = true,
               const bool request_edges = false,
               const bool request_wedges = false,
               const bool request_corners = false,
               const int num_tiles = 0);


  /// Create a hexahedral mesh of the specified dimensions
  Mesh *create(double x0, double y0, double z0,
               double x1, double y1, double z1,
               int nx, int ny, int nz,
               const JaliGeometry::GeometricModelPtr &gm =
               (JaliGeometry::GeometricModelPtr) NULL,
               const bool request_faces = true,
               const bool request_edges = false,
               const bool request_wedges = false,
               const bool request_corners = false,
               const int num_tiles = 0);


  /// Create a quadrilateral mesh of the specified dimensions
  Mesh *create(double x0, double y0,
               double x1, double y1,
               int nx, int ny,
               const JaliGeometry::GeometricModelPtr &gm =
               (JaliGeometry::GeometricModelPtr) NULL,
               const bool request_faces = true,
               const bool request_edges = false,
               const bool request_wedges = false,
               const bool request_corners = false,
               const int num_tiles = 0);


  /// Create a mesh by extract subsets of entities from an existing mesh
  Mesh *create(const Mesh *inmesh,
               const std::vector<std::string> setnames,
               const Entity_kind setkind,
               const bool flatten = false,
               const bool extrude = false,
               const bool request_faces = true,
               const bool request_edges = false,
               const bool request_wedges = false,
               const bool request_corners = false,
               const int num_tiles = 0);

 public:

  /// Default constructor.
  explicit MeshFactory(const MPI_Comm& communicator);

  /// Destructor
  ~MeshFactory(void);

  /// Get the framework preference
  const FrameworkPreference& preference(void) const
  { return my_preference; }

  /// Set the framework preference
  void preference(const FrameworkPreference& pref);

  /// Create a mesh by reading the specified file (or set of files) -- operator
  std::unique_ptr<Mesh> operator() (const std::string& filename,
                                    const JaliGeometry::GeometricModelPtr &gm =
                                    (JaliGeometry::GeometricModelPtr) NULL,
                                    const bool request_faces = true,
                                    const bool request_edges = false,
                                    const bool request_wedges = false,
                                    const bool request_corners = false,
                                    const int num_tiles = 0) {

    return std::unique_ptr<Mesh>(create(filename, gm, request_faces,
                                        request_edges, request_wedges,
                                        request_corners, num_tiles));
  }

  /// Create a hexahedral mesh of the specified dimensions -- operator
  std::unique_ptr<Mesh> operator() (double x0, double y0, double z0,
                                    double x1, double y1, double z1,
                                    int nx, int ny, int nz,
                                    const JaliGeometry::GeometricModelPtr &gm =
                                    (JaliGeometry::GeometricModelPtr) NULL,
                                    const bool request_faces = true,
                                    const bool request_edges = false,
                                    const bool request_wedges = false,
                                    const bool request_corners = false,
                                    const int num_tiles = 0) {

    return std::unique_ptr<Mesh>(create(x0, y0, z0, x1, y1, z1, nx, ny, nz, gm,
                                        request_faces, request_edges,
                                        request_wedges, request_corners,
                                        num_tiles));
  }

  /// Create a quadrilateral mesh of the specified dimensions -- operator
  std::unique_ptr<Mesh> operator() (double x0, double y0,
                                    double x1, double y1,
                                    int nx, int ny,
                                    const JaliGeometry::GeometricModelPtr &gm =
                                    (JaliGeometry::GeometricModelPtr) NULL,
                                    const bool request_faces = true,
                                    const bool request_edges = false,
                                    const bool request_wedges = false,
                                    const bool request_corners = false,
                                    const int num_tiles = 0)  {

    return std::unique_ptr<Mesh>(create(x0, y0, x1, y1, nx, ny, gm,
                                        request_faces, request_edges,
                                        request_wedges, request_corners,
                                        num_tiles));
  }

  /// Create a mesh by extract subsets of entities from an existing mesh
  std::unique_ptr<Mesh> operator() (const Mesh *inmesh,
                                    const std::vector<std::string> setnames,
                                    const Entity_kind setkind,
                                    const bool flatten = false,
                                    const bool extrude = false,
                                    const bool request_faces = true,
                                    const bool request_edges = false,
                                    const bool request_wedges = false,
                                    const bool request_corners = false,
                                    const int num_tiles = 0) {

    return std::unique_ptr<Mesh>(create(inmesh, setnames, setkind, flatten,
                                        extrude, request_faces, request_edges,
                                        request_wedges, request_corners,
                                        num_tiles));
  }

};

}  // namespace Jali

#endif
