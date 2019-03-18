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


/*!
 * @file   MeshFactory.hh
 * @author Rao Garimella
 *
 * @brief  declaration of the MeshFactory class
 *
 */

#ifndef _MESHFACTORY_HH_
#define _MESHFACTORY_HH_

#include <mpi.h>

#include <string>
#include <vector>
#include <memory>
#include <utility>

#include "MeshDefs.hh"
#include "Mesh.hh"
#include "GeometricModel.hh"
#include "Geometry.hh"
#include "errors.hh"
#include "exceptions.hh"

namespace Jali {

/// A type to identify available mesh frameworks
enum MeshFramework_t {
  Simple = 1,
  MSTK,
  MOAB,
  STKMESH
};

/// A type to identify mesh file formats
enum MeshFormat_t {
  ExodusII = 1,
  MOABHDF5,
  FLAGX3D
};

/// Get a name for a given framework
std::string framework_name(MeshFramework_t const & f);

/// Check if a framework is available for use
bool framework_available(MeshFramework_t const & f);

/// Check if a framework can generate a mesh
bool framework_generates(MeshFramework_t const & f, bool const parallel,
                         int const dim);

bool framework_reads(MeshFramework_t const& f, bool const parallel,
                     MeshFormat_t const& format);

bool framework_extracts(MeshFramework_t const& f, bool const parallel,
                        int const dim);


/* @brief Factory class to construct meshes using chosen mesh framework */

class MeshFactory {
 public:

  /// Default constructor.
  explicit MeshFactory(MPI_Comm const& communicator) : comm_(communicator) {}

  /// Delete the copy constructor
  MeshFactory(MeshFactory& factory) = delete;

  /// Destructor
  ~MeshFactory(void) {}

  /// Reset options to default values;
  void reset_options(void);

  // Get/Set methods for mesh initialization parameters

  /// Get the framework being used
  MeshFramework_t framework(void) const {
    return framework_;
  }

  /// Set the framework to use
  void framework(MeshFramework_t const& framework) {
    if (framework == Simple || framework == MSTK) {
      framework_ = framework;
    } else {
      std::stringstream mesgstrm;
      mesgstrm << framework_name(framework) << " not available\n";
      Errors::Message errmsg(mesgstrm.str());
      Exceptions::Jali_throw(errmsg);
    }
  }

  /// Get the partitioner for meshes to be created (default METIS)
  Partitioner_type partitioner(void) const {
    return partitioner_;
  }

  /// Set the partitioner
  void partitioner(Partitioner_type partitioner) {
    partitioner_ = partitioner;
  }

  /// Get the geometry type for the meshes to be created (default CARTESIAN)
  JaliGeometry::Geom_type mesh_geometry(void) const {
    return geom_type_;
  }

  /// Set the geometry type for the meshes to be created
  void mesh_geometry(JaliGeometry::Geom_type geom_type) {
    geom_type_ = geom_type;
  }

  /// Get the number of ghost layers around distributed mesh
  /// partitions in the meshes to be created (default 1)
  int num_ghost_layers_distmesh(void) const {
    return num_ghost_layers_distmesh_;
  }

  /// Set the number of ghost layers around distributed mesh
  /// partitions in the meshes to be created
  void num_ghost_layers_distmesh(int num_layers) {
    num_ghost_layers_distmesh_ = num_layers;
  }

  /// Are ghost/virtual elements outside external boundaries of the
  /// meshes to be created (default false)
  bool boundary_ghosts_requested(void) const {
    return request_boundary_ghosts_;
  }

  /// Set the number of ghost/virtual element layers outside
  /// external boundarieis in the meshes to be created
  void boundary_ghosts_requested(bool requested_or_not) {
    request_boundary_ghosts_ = requested_or_not;
  }

  /// Get the number of tiles to be created (default 0)

  int num_tiles(void) const {
    return num_tiles_;
  }

  /// Set the number of tiles to be created
  
  void num_tiles(int n) {
    num_tiles_ = n;
  }

  /// Get the number of ghost layers around on-node mesh tiles in the
  /// meshes to be created (default 0)

  int num_ghost_layers_tile(void) const {
    return num_ghost_layers_tile_;
  }

  /// Set the number of ghost layers around on-node mesh tiles in
  /// the meshes to be created
  
  void num_ghost_layers_tile(int num_layers) {
    num_ghost_layers_tile_ = num_layers;
  }

  /// Request that the GIDs be made contiguous
  void contiguous_gids(bool make_contiguous) {
    contiguous_gids_ = make_contiguous;
  }

  /// @brief Get explicitly represented entity kinds 
  ///
  /// Get the types of entities that are explicitly requested in the
  /// meshes to be created (Nodes are always included. Cells are usually
  /// included unless explicitly excluded for particle meshes - NOT
  /// IMPLEMENTED YET)
  ///
  
  std::vector<Entity_kind> included_entities(void) const {
    std::vector<Entity_kind> list;
    list.push_back(Entity_kind::NODE);  // always present
    if (request_edges_) list.push_back(Entity_kind::EDGE);
    if (request_faces_) list.push_back(Entity_kind::FACE);
    if (request_sides_) list.push_back(Entity_kind::SIDE);
    if (request_wedges_) list.push_back(Entity_kind::WEDGE);
    if (request_corners_) list.push_back(Entity_kind::CORNER);
    if (request_cells_) list.push_back(Entity_kind::CELL);
    return list;
  }

  /// Set the type of entities that are explicitly INCLUDED/REQUESTED in the
  /// meshes to be created
  
  void included_entities(std::vector<Entity_kind> const& list) {
    for (auto const& e : list)
      included_entities(e);
  }

  /// Set the type of entities that are explicitly INCLUDED/REQUESTED in the
  /// meshes to be created
  void included_entities(Entity_kind const e) {
    if (e == Entity_kind::ALL_KIND) {
      request_edges_ = true;
      request_faces_ = true;
      request_sides_ = true;
      request_wedges_ = true;
      request_corners_ = true;
      request_cells_ = true;
    } else {
      switch (e) {
        case Entity_kind::NODE: break;  // included by default - nothing to do
        case Entity_kind::EDGE: request_edges_ = true; break;
        case Entity_kind::FACE: request_faces_ = true; break;
        case Entity_kind::SIDE: request_sides_ = true; break;
        case Entity_kind::WEDGE: request_wedges_ = true; break;
        case Entity_kind::CORNER: request_corners_ = true; break;
        case Entity_kind::CELL: request_cells_ = true; break;
        default:
          std::cerr << "Invalid entity kind asked to be explicitly included\n";
      }
    }
  }

  /// Set the type of entities that are explicitly EXCLUDED in the
  /// meshes to be created
  void excluded_entities(std::vector<Entity_kind> const& list) {
    for (auto const& e : list) {
      switch (e) {
        case Entity_kind::NODE:
          std::cerr << "MeshFactory::exclude_entities: " <<
              "Cannot turn off node creation in meshes\n";
          break;
        case Entity_kind::EDGE: request_edges_ = false; break;
        case Entity_kind::SIDE: request_sides_ = false; break;
        case Entity_kind::WEDGE: request_wedges_ = false; break;
        case Entity_kind::CORNER: request_corners_ = false; break;
        case Entity_kind::FACE:
          std::cerr << "MeshFactory::excluded_entities: " <<
              "Cannot turn off face creation in meshes CURRENTLY - " <<
              "Needed for cell geometric quantity computations\n"; break;
        case Entity_kind::CELL:
          std::cerr << "MeshFactory::excluded_entities: " <<
              "Cannot turn off cell creation in meshes CURRENTLY\n";
          break;
        default:
          std::cerr << "Invalid entity kind asked to be explicitly excluded\n";
      }
    }
  }

  /// Get a pointer to the geometric model underlying the mesh to be created
  JaliGeometry::GeometricModelPtr geometric_model(void) const {
    return geometric_model_;
  }

  /// Set the pointer to the geometric model undelying the mesh to be created
  void geometric_model(JaliGeometry::GeometricModelPtr gm) {
    geometric_model_ = gm;
  }

  /// Create a mesh by reading the specified file (or set of files) -- operator
  std::shared_ptr<Mesh> operator() (std::string const& filename) {
    return create(filename);
  }

  /// Create a hexahedral mesh of the specified dimensions -- operator
  std::shared_ptr<Mesh> operator() (double const x0, double const y0,
                                    double const z0,
                                    double const x1, double const y1,
                                    double const z1,
                                    int const nx, int const ny, int const nz) {
    return create(x0, y0, z0, x1, y1, z1, nx, ny, nz);
  }

  /// Create a quadrilateral mesh of the specified dimensions -- operator
  std::shared_ptr<Mesh> operator() (double const x0, double const y0,
                                    double const x1, double const y1,
                                    int const nx, int const ny) {
    return create(x0, y0, x1, y1, nx, ny);
  }

  /// Create a 1d mesh -- operator
  std::shared_ptr<Mesh> operator() (std::vector<double> const& x) {
    return create(x);
  }

  /// Create a 1d mesh -- operator
  std::shared_ptr<Mesh> operator() (double const x0, double const x1,
                                    int const nx) {
    double dX = (x1-x0)/((double)nx);
    double myX = x0;

    std::vector<double> x(nx);
    for (auto it = x.begin(); it != x.end(); it++) {
      *it = myX;
      myX += dX;
    }

    return create(x);
  }

  /// Create a mesh by extract subsets of entities from an existing mesh
  std::shared_ptr<Mesh> operator() (std::shared_ptr<Mesh> const inmesh,
                                    std::vector<std::string> const& setnames,
                                    Entity_kind const setkind,
                                    bool const flatten = false,
                                    bool const extrude = false) {
    return create(inmesh, setnames, setkind, flatten, extrude);
  }

 private:

  /// Create a mesh by reading the specified file (or set of files)
  std::shared_ptr<Mesh> create(std::string const& filename);

  /// Create a hexahedral mesh of the specified dimensions
  //  No need of geom_type argument as its always CARTESIAN

  std::shared_ptr<Mesh> create(double const x0, double const y0,
                               double const z0,
                               double const x1, double const y1,
                               double const z1,
                               int const nx, int const ny, int const nz);

  /// Create a quadrilateral mesh of the specified dimensions
  std::shared_ptr<Mesh> create(double const x0, double const y0,
                               double const x1, double const y1,
                               int const nx, int const ny);
  
  /// Create a 1d mesh
  std::shared_ptr<Mesh> create(std::vector<double> const& x);


  /// Create a mesh by extract subsets of entities from an existing mesh
  std::shared_ptr<Mesh> create(std::shared_ptr<Mesh> const inmesh,
                               std::vector<std::string> const& setnames,
                               Entity_kind const setkind,
                               bool const flatten = false,
                               bool const extrude = false);


  /// The parallel environment
  MPI_Comm const comm_;

  /// The framework to use for creating the mesh
  MeshFramework_t const default_framework_ = MSTK;
  MeshFramework_t framework_ = default_framework_;

  /// What type of entities to include/exclude in the meshes to be created
  /// (nodes must ALWAYS be present)
  bool const request_edges_default_ = false;
  bool request_edges_ = request_edges_default_;
  bool const request_faces_default_ = true;  // cell geometry computation needs faces
  bool request_faces_ = request_faces_default_;
  bool const request_cells_default_ = true;
  bool request_cells_ = request_cells_default_;
  bool const request_sides_default_ = false;
  bool request_sides_ = request_sides_default_;
  bool const request_wedges_default_ = false;
  bool request_wedges_ = request_wedges_default_;
  bool const request_corners_default_ = false;
  bool request_corners_ = request_corners_default_;

  /// Whether the mesh should contain boundary ghosts (dummy cells
  /// outside a boundary to ease some numerical computations

  bool request_boundary_ghosts_default_ = false;
  bool request_boundary_ghosts_ = request_boundary_ghosts_default_;
  
  /// Number of on-node mesh tiles
  int const num_tiles_default_ = 0;
  int num_tiles_ = num_tiles_default_;

  /// Number of ghost/halo layers at the tile level on compute node
  int const num_ghost_layers_tile_default_ = 0;
  int num_ghost_layers_tile_ = num_ghost_layers_tile_default_;

  /// Number of ghost/halo layers for mesh partitions across compute nodes
  int const num_ghost_layers_distmesh_default_ = 1;
  int num_ghost_layers_distmesh_ = num_ghost_layers_distmesh_default_;

  /// Partitioner type
  /// SHOULD CHANGE THIS TO Partitioner_type::BLOCK
  Partitioner_type const partitioner_default_ = Partitioner_type::INDEX;
  Partitioner_type partitioner_ = partitioner_default_;

  /// Geometry type
  JaliGeometry::Geom_type const geom_type_default_ =
      JaliGeometry::Geom_type::CARTESIAN;
  JaliGeometry::Geom_type geom_type_ = geom_type_default_;

  /// Geometric model  
  JaliGeometry::GeometricModel *geometric_model_ = nullptr;

  /// Should GIDs be made contiguous?
  bool contiguous_gids_ = false;
};

}  // namespace Jali

#endif
