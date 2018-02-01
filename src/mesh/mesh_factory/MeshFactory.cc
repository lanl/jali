// -------------------------------------------------------------
// Copyright Los Alamos National Security, LLC
// -------------------------------------------------------------


#include "MeshFactory.hh"
#include "Geometry.hh"

#include "Mesh_simple.hh"

#ifdef HAVE_MSTK_MESH
#include "Mesh_MSTK.hh"
#endif

#ifdef HAVE_STK_MESH
#include "Mesh_STK.hh"
#endif

#ifdef HAVE_MOAB_MESH
#define USE_MPI
#include "Mesh_MOAB.hh"
#endif

namespace Jali {

/// Get a name for a given framework
std::string framework_name(MeshFramework_t const& f) {
  switch (f) {
    case (Simple):
      return "Simple";
      break;
    case (STKMESH):
      return "stk::mesh";
      break;
    case (MOAB):
      return "MOAB";
      break;
    case (MSTK):
      return "MSTK";
      break;
    default:
      Errors::Message mesg("Unknown framework");
      Exceptions::Jali_throw(mesg);
  }
}

/// Check if a framework is available for use
bool framework_available(MeshFramework_t const& f) {
  if (f == Simple || f == MSTK)
    return true;
  else
    return false;
}

/// Check if a framework can generate a mesh
bool framework_generates(MeshFramework_t const& f, bool const parallel,
                         int const dim) {
  switch (f) {
    case Jali::Simple:
      return (!parallel && (dim == 1 || dim == 3));
    case Jali::MSTK:
      return (dim == 2 || dim == 3);
    case Jali::STKMESH:
      return (dim == 3 && !parallel);
    default:
      return false;
  }
}

bool framework_reads(MeshFramework_t const& f, bool const parallel,
                     MeshFormat_t const& format) {
  switch (f) {
    case Jali::MSTK:
      return (format == Jali::ExodusII);
    case Jali::STKMESH:
      return (format == Jali::ExodusII && !parallel);
    case Jali::MOAB:
      return (format == Jali::MOABHDF5);
    default:
      return false;
  }
}

bool framework_extracts(MeshFramework_t const& f, bool const parallel,
                        int const dim) {
  switch (f) {
    case Jali::MSTK:
      return (dim == 2 || dim == 3);
    default:
      return false;
  }
}


// -------------------------------------------------------------
//  class MeshFactory
// -------------------------------------------------------------

// Reset mesh constructor options to their default values

void MeshFactory::reset_options(void) {
  /// The framework to use for creating the mesh
  MeshFramework_t framework_ = default_framework_;

  // What type of entities to include/exclude in the meshes to be created
  // (nodes must ALWAYS be present)
  request_edges_ = request_edges_default_;
  request_faces_ = request_faces_default_;
  request_cells_ = request_cells_default_;
  request_wedges_ = request_wedges_default_;
  request_corners_ = request_corners_default_;

  /// Number of on-node mesh tiles
  num_tiles_ = num_tiles_default_;

  /// Number of ghost/halo layers at the tile level on compute node
  num_ghost_layers_tile_ = num_ghost_layers_tile_default_;

  /// Number of ghost/halo layers for mesh partitions across compute nodes
  num_ghost_layers_distmesh_ = num_ghost_layers_distmesh_default_;

  /// Whether ghost/virtual elements outside external boundaries are requested
  request_boundary_ghosts_ = false;

  /// Partitioner type
  partitioner_ = partitioner_default_;

  /// Geometry type
  geom_type_ = geom_type_default_;

  /// Geometric model
  geometric_model_ = nullptr;
}

/**
 *
 * @brief Create a mesh by reading the specified file (or file set).
 * @param filename mesh file to read
 * @return mesh instance
 */

std::shared_ptr<Mesh>
MeshFactory::create(std::string const& filename) {
  Errors::Message errmsg("MeshFactory:: unable to create mesh");
  int ierr = 0, aerr = 0;

  std::shared_ptr<Mesh> result;
  try {
    switch (framework_) {
      case MSTK: {        
        result =
            std::make_shared<Mesh_MSTK>(filename, comm_, geometric_model_,
                                        request_faces_, request_edges_,
                                        request_sides_, request_wedges_,
                                        request_corners_,
                                        num_tiles_, num_ghost_layers_tile_,
                                        num_ghost_layers_distmesh_,
                                        request_boundary_ghosts_,
                                        partitioner_, geom_type_);
        if (geometric_model_ &&
            (geometric_model_->dimension() != result->space_dimension())) {
          errmsg.add_data("Geometric model and mesh dimension do not match");
          Exceptions::Jali_throw(errmsg);
        }
        return result;
        break;
      }
      default:
        errmsg.add_data("Chosen framework cannot import meshes");
        ierr = 1;
        break;
    }
  } catch (const Errors::Message& msg) {
    ierr = 1;
    errmsg.add_data(msg.what());
  } catch (const std::exception& stde) {
    ierr = 1;
    errmsg.add_data("internal error: ");
    errmsg.add_data(stde.what());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);
}

/**
 * @brief Create a mesh by generating a block of hexahedral cells.
 *
 * @param x0 minimum x-coordinate
 * @param y0 minimum y-coordinate
 * @param z0 minimum z-coordinate
 * @param x1 maximum x-coordinate
 * @param y1 maximum y-coordinate
 * @param z1 maximum z-coordinate
 * @param nx number of cells in the x-direction
 * @param ny number of cells in the y-direction
 * @param nz number of cells in the z-direction
 *
 * @return mesh instance
 */
std::shared_ptr<Mesh>
MeshFactory::create(double const x0, double const y0, double const z0,
                    double const x1, double const y1, double const z1,
                    int const nx, int const ny, int const nz) {
  std::stringstream mesgstr;
  std::shared_ptr<Mesh> result;
  Errors::Message errmsg("MeshFactory::create - Unable to create 3D mesh");
  int ierr = 0, aerr = 0;

  unsigned int dim = 3;

  if (geometric_model_ && (geometric_model_->dimension() != 3)) {
    Errors::Message mesg("Geometric model and mesh dimension do not match");
    Exceptions::Jali_throw(mesg);
  }

  if (nx <= 0 || ny <= 0 || nz <= 0) {
    ierr = 1;
    mesgstr << "invalid number of mesh cells requested: " << nx << "x" << ny
            << "x" << nz;
    errmsg.add_data(mesgstr.str());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);

  if (x1 - x0 <= 0.0 || y1 - y0 <= 0.0 || z1 - z0 <= 0.0) {
    ierr += 1;
    mesgstr << "invalid mesh dimensions requested: " << (x1 - x0) << "x" <<
        (y1 - y0) << "x" << (z1 - z0);
    errmsg.add_data(mesgstr.str());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);

  int numprocs;
  MPI_Comm_size(comm_, &numprocs);

  try {
    switch (framework_) {
      case Simple: {
        if (numprocs == 1) {
          result =
              std::make_shared<Mesh_simple>(x0, y0, z0,
                                            x1, y1, z1,
                                            nx, ny, nz,
                                            comm_, geometric_model_,
                                            request_faces_, request_edges_,
                                            request_sides_, request_wedges_,
                                            request_corners_,
                                            num_tiles_, num_ghost_layers_tile_,
                                            num_ghost_layers_distmesh_,
                                            request_boundary_ghosts_,
                                            partitioner_);
          return result;
        }
        else {
          ierr = 1;
          errmsg.add_data("Simple mesh cannot generate parallel meshes");
        }
        break;
      }
      case MSTK: {
        result =
            std::make_shared<Mesh_MSTK>(x0, y0, z0, x1, y1, z1, nx, ny, nz,
                                        comm_, geometric_model_,
                                        request_faces_, request_edges_,
                                        request_sides_, request_wedges_,
                                        request_corners_,
                                        num_tiles_, num_ghost_layers_tile_,
                                        num_ghost_layers_distmesh_,
                                        request_boundary_ghosts_,
                                        partitioner_);
        return result;
      }
      default:
        ierr = 1;
        errmsg.add_data("Chosen framework cannot generate meshes");
    }        
  } catch (const Errors::Message& msg) {
    ierr = 1;
    errmsg.add_data(msg.what());
  } catch (const std::exception& stde) {
    ierr = 1;
    errmsg.add_data("internal error: ");
    errmsg.add_data(stde.what());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);
}


/**
 * @brief Create a mesh by generating a block of quadrilateral cells.
 *
 * @param x0 origin x-coordinate
 * @param y0 origin y-coordinate
 * @param x1 maximum x-coordinate
 * @param y1 maximum y-coordinate
 * @param nx number of cells in the x-direction
 * @param ny number of cells in the y-direction
 *
 * @return mesh instance
 */

std::shared_ptr<Mesh>
MeshFactory::create(double const x0, double const y0,
                    double const x1, double const y1,
                    int const nx, int const ny) {
  std::shared_ptr<Mesh> result;
  Errors::Message errmsg("MeshFactory::create: error: ");
  int ierr = 0, aerr = 0;

  unsigned int dim = 2;

  if (geometric_model_ && (geometric_model_->dimension() != 2)) {
    Errors::Message mesg("Geometric model and mesh dimension do not match");
    Exceptions::Jali_throw(mesg);
  }

  if (nx <= 0 || ny <= 0) {
    ierr = 1;
    std::stringstream mesgstream;
    mesgstream << "invalid number of mesh cells requested: " << nx << " x " <<
        ny;
    errmsg.add_data(mesgstream.str());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);

  if (x1 - x0 <= 0.0 || y1 - y0 <= 0.0) {
    ierr = 1;
    std::stringstream mesgstream;
    mesgstream << "invalid mesh dimensions requested: " << (x1-x0) << " x " <<
        (y1-y0);
    errmsg.add_data(mesgstream.str());
  }

  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);

  int numprocs;
  MPI_Comm_size(comm_, &numprocs);

  try {
    switch (framework_) {
      case Simple: {
        if (numprocs == 1) {
          result =
              std::make_shared<Mesh_simple>(x0, y0, x1, y1, nx, ny,
                                            comm_, geometric_model_,
                                            request_faces_, request_edges_,
                                            request_sides_, request_wedges_,
                                            request_corners_,
                                            num_tiles_, num_ghost_layers_tile_,
                                            num_ghost_layers_distmesh_,
                                            request_boundary_ghosts_,
                                            partitioner_, geom_type_);
          return result;
        } else {
          ierr = 1;
          errmsg.add_data("Simple framework cannot generate mesh in parallel");
        }
        break;
      }
      case MSTK: {
        result =
            std::make_shared<Mesh_MSTK>(x0, y0, x1, y1, nx, ny,
                                        comm_, geometric_model_,
                                        request_faces_, request_edges_,
                                        request_sides_, request_wedges_,
                                        request_corners_,
                                        num_tiles_, num_ghost_layers_tile_,
                                        num_ghost_layers_distmesh_,
                                        request_boundary_ghosts_,
                                        partitioner_, geom_type_);
        return result;
      }
      default: {
        errmsg.add_data("Chosen framework cannnot generate meshes");        
      }
    }
  } catch (const Errors::Message& msg) {
    ierr = 1;
    errmsg.add_data(msg.what());
  } catch (const std::exception& stde) {
    ierr = 1;
    errmsg.add_data("internal error: ");
    errmsg.add_data(stde.what());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);
}


/**
 * @brief Create a mesh by generating a block of 1d cells.
 *
 * @param x vector of spatial coordinates of nodes
 * @return mesh instance
 */

std::shared_ptr<Mesh>
MeshFactory::create(std::vector<double> const& x) {
  std::shared_ptr<Mesh> result;
  Errors::Message errmsg("MeshFactory::create: error: ");
  int ierr = 0, aerr = 0;

  unsigned int dim = 1;

  if (geometric_model_ && (geometric_model_->dimension() != 1)) {
    ierr = 1;
    errmsg.add_data("Geometric model and mesh dimension do not match");
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);

  if (x.size() < 2) {
    ierr = 1;
    std::stringstream mesgstream;
    mesgstream << "invalid num nodes requested: " << x.size();
    errmsg.add_data(mesgstream.str());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);

  double delta = x.back() - x.front();
  if (delta <= 0.0) {
    ierr = 1;
    std::stringstream mesgstream;
    mesgstream << "Invalid mesh spacing requested (" << delta << ")";
    errmsg.add_data(mesgstream.str());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);

  int numprocs;
  MPI_Comm_size(comm_, &numprocs);

  try {
    switch (framework_) {
      case Simple: {
        if (numprocs == 1) {
          result =
              std::make_shared<Mesh_simple>(x, comm_, geometric_model_,
                                            request_faces_, request_edges_,
                                            request_sides_, request_wedges_,
                                            request_corners_,
                                            num_tiles_, num_ghost_layers_tile_,
                                            num_ghost_layers_distmesh_,
                                            request_boundary_ghosts_,
                                            partitioner_, geom_type_);
          return result;
        } else {
          ierr = 1;
          errmsg.add_data("Simple framework cannot generate parallel 1D mesh");
        }
      }
      default: {
        ierr = 1;
        errmsg.add_data("Chosen framework cannot generate 1D mesh");
      }
    }
  } catch (const Errors::Message& msg) {
    ierr = 1;
    errmsg.add_data(msg.what());
  } catch (const std::exception& stde) {
    ierr = 1;
    errmsg.add_data("internal error: ");
    errmsg.add_data(stde.what());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);
}

/**
 * @brief Create mesh by extracting subset of entities
 *
 * This creates a mesh by extracting subsets of entities from an existing
 * mesh possibly flattening it by removing the last dimension or (in the
 * future) extruding it when it makes sense
 *
 * @param inmesh  mesh from which to extract the new mesh
 * @param setnames  the names of entity sets to extract the mesh from
 * @param setkind   the kind of entities in the new mesh
 * @param flatten   whether the mesh should be dimensionally reduced (topological and spatial dimensionality of new mesh will be 1 less than the old mesh)
 * @param extrude  whether the extracted entities should be extruded in the z-direction
 *
 * @return
 */
std::shared_ptr<Mesh>
MeshFactory::create(std::shared_ptr<Mesh> const inmesh,
                    std::vector<std::string> const& setnames,
                    Entity_kind const setkind,
                    bool const flatten, bool const extrude) {
  std::shared_ptr<Mesh> result;
  Errors::Message errmsg("MeshFactory::create: error: ");
  int ierr = 0, aerr = 0;

  int dim = inmesh->manifold_dimension();
  int numprocs;
  MPI_Comm_size(comm_, &numprocs);

  try {
    switch (framework_) {
      case MSTK: {
        result =
            std::make_shared<Mesh_MSTK>(*inmesh,
                                        setnames, setkind,
                                        flatten, extrude,
                                        request_faces_, request_edges_,
                                        request_sides_, request_wedges_,
                                        request_corners_,
                                        num_tiles_, num_ghost_layers_tile_,
                                        num_ghost_layers_distmesh_,
                                        request_boundary_ghosts_,
                                        partitioner_, geom_type_);
        return result;
      }
      default: {
        ierr = 1;
        errmsg.add_data("Chosen framework cannot extract meshes");
      }
    }
  } catch (const Errors::Message& msg) {
    ierr = 1;
    errmsg.add_data(msg.what());
  } catch (const std::exception& stde) {
    ierr = 1;
    errmsg.add_data("internal error: ");
    errmsg.add_data(stde.what());
  }
  MPI_Allreduce(&ierr, &aerr, 1, MPI_INT, MPI_SUM, comm_);
  if (aerr > 0) Exceptions::Jali_throw(errmsg);
}

}  // namespace Jali
