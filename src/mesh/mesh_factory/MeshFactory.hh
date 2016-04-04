/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/*!
 * @file   MeshFactory.hh
 * @author William A. Perkins, Rao Garimella
 * @date Wed Sep 28 09:10:15 2011
 *
 * @brief  declaration of the MeshFactory class
 *
 *
 */

#ifndef _MeshFactory_hh_
#define _MeshFactory_hh_

#include <string>
#include <vector>

#include <mpi.h>
#include <memory>
#include <utility>

#include "MeshException.hh"
#include "MeshFramework.hh"
#include "Mesh.hh"

#include "GeometricModel.hh"
#include "Geometry.hh"

namespace Jali {

// -------------------------------------------------------------
//  class MeshFactory
// -------------------------------------------------------------
class MeshFactory {

 protected:

  /// The parallel environment
  const MPI_Comm my_comm_;

  /// A list of preferred mesh frameworks to consider
  FrameworkPreference my_preference_;

  /// What type of entities to include/exclude in the meshes to be created
  /// (nodes must ALWAYS be present)
  bool request_edges_;
  bool request_faces_;
  bool request_cells_;
  bool request_wedges_;
  bool request_corners_;

  /// Number of on-node mesh tiles
  int num_tiles_;

  /// Number of ghost/halo layers at the tile level on compute node
  int num_ghost_layers_tile_;

  /// Number of ghost/halo layers for mesh partitions across compute nodes
  int num_ghost_layers_distmesh_;

  /// Partitioner type
  Partitioner_type partitioner_;

  /// Geometry type
  JaliGeometry::Geom_type geom_type_;

  /// Geometric model
  JaliGeometry::GeometricModel *geometric_model_;

 private:

  /// private, undefined copy constructor to avoid unwanted copies
  MeshFactory(MeshFactory& old);

  /// Create a mesh by reading the specified file (or set of files)
  std::shared_ptr<Mesh> create(const std::string& filename);

  /// Create a hexahedral mesh of the specified dimensions
  //  No need of geom_type argument as its always CARTESIAN

  std::shared_ptr<Mesh> create(double x0, double y0, double z0,
                               double x1, double y1, double z1,
                               int nx, int ny, int nz);

  /// Create a quadrilateral mesh of the specified dimensions
  std::shared_ptr<Mesh> create(double x0, double y0,
                               double x1, double y1,
                               int nx, int ny);

  /// Create a 1d mesh
  std::shared_ptr<Mesh> create(const std::vector<double>& x);



  /// Create a mesh by extract subsets of entities from an existing mesh
  std::shared_ptr<Mesh> create(const std::shared_ptr<Mesh> inmesh,
                               const std::vector<std::string> setnames,
                               const Entity_kind setkind,
                               const bool flatten = false,
                               const bool extrude = false);

 public:

  /// Default constructor.
  explicit MeshFactory(const MPI_Comm& communicator);

  /// Destructor
  ~MeshFactory(void);

  /// Reset options to default values;
  void reset_options(void);

  // Get/Set methods for mesh initialization parameters

  /// Get the framework preference
  const FrameworkPreference& preference(void) const
  { return my_preference_; }

  /// Set the framework preference
  void preference(const FrameworkPreference& pref);

  /// Get the partitioner for meshes to be created (default METIS)
  Partitioner_type partitioner(void) const {return partitioner_;}

  /// Set the partitioner
  void partitioner(Partitioner_type partitioner) {partitioner_ = partitioner;}

  /// Get the geometry type for the meshes to be created (default CARTESIAN)
  JaliGeometry::Geom_type mesh_geometry(void) const {return geom_type_;}

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

  /// Get the number of tiles to be created (default 0)
  int num_tiles(void) const {return num_tiles_;}

  /// Set the number of tiles to be created
  void num_tiles(int n) {num_tiles_ = n;}

  /// Get the number of ghost layers around on-node mesh tiles in the
  /// meshes to be created (default 0)
  int num_ghost_layers_tile(void) const {return num_ghost_layers_tile_;}

  /// Set the number of ghost layers around distributed mesh
  /// partitions in the meshes to be created
  void num_ghost_layers_tile(int num_layers) {
    num_ghost_layers_tile_ = num_layers;
  }

  /*! @brief Get explicitly represented entity kinds 
   
   Get the types of entities that are explicitly requested in the
   meshes to be created (Nodes are always included. Cells are usually
   included unless explicitly excluded for particle meshes - NOT
   IMPLEMENTED YET)
  */
  std::vector<Entity_kind> included_entities(void) const {
    std::vector<Entity_kind> list;
    list.push_back(Entity_kind::NODE);  // always present
    if (request_edges_) list.push_back(Entity_kind::EDGE);
    if (request_faces_) list.push_back(Entity_kind::FACE);
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
      request_wedges_ = true;
      request_corners_ = true;
      request_cells_ = true;
    }
    else {
      switch (e) {
        case Entity_kind::NODE: break;  // included by default - nothing to do
        case Entity_kind::EDGE: request_edges_ = true; break;
        case Entity_kind::FACE: request_faces_ = true; break;
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
  std::shared_ptr<Mesh> operator() (const std::string& filename) {
    return create(filename);
  }

  /// Create a hexahedral mesh of the specified dimensions -- operator
  std::shared_ptr<Mesh> operator() (double x0, double y0, double z0,
                                    double x1, double y1, double z1,
                                    int nx, int ny, int nz) {
      return create(x0, y0, z0, x1, y1, z1, nx, ny, nz);
  }

  /// Create a quadrilateral mesh of the specified dimensions -- operator
  std::shared_ptr<Mesh> operator() (double x0, double y0,
                                    double x1, double y1,
                                    int nx, int ny) {
    return create(x0, y0, x1, y1, nx, ny);
  }

  /// Create a 1d mesh -- operator
  std::shared_ptr<Mesh> operator() (const std::vector<double>& x) {
    return create(x);
  }

  /// Create a 1d mesh -- operator
  std::shared_ptr<Mesh> operator() (double x0, double x1, int nx) {
    double dX = (x1-x0)/((double)nx);
    double myX = x0;

    std::vector<double> x(nx);
    for(auto it = x.begin(); it != x.end(); it++) {
      *it = myX;
      myX += dX;
    }

    return create(x);
  }

  /// Create a mesh by extract subsets of entities from an existing mesh
  std::shared_ptr<Mesh> operator() (const std::shared_ptr<Mesh> inmesh,
                                    const std::vector<std::string> setnames,
                                    const Entity_kind setkind,
                                    const bool flatten = false,
                                    const bool extrude = false) {
    return create(inmesh, setnames, setkind, flatten, extrude);
  }

};

}  // namespace Jali

#endif
