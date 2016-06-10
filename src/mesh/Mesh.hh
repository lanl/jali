//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#ifndef _JALIMESH_H_
#define _JALIMESH_H_

#include <mpi.h>

#include <memory>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <cassert>
#include <typeinfo>

#include "MeshDefs.hh"
#include "Point.hh"
#include "GeometricModel.hh"
#include "Region.hh"
#include "Geometry.hh"

#include "MeshTile.hh"

#define JALI_CACHE_VARS 1  // Switch to 0 to turn caching off

//! \mainpage Jali
//!
//! Jali is a parallel unstructured mesh infrastructure library for
//! multiphysics applications. It simplifies the process of importing,
//! querying, manipulating and exporting unstructured mesh information
//! for physics and numerical methods developers. Jali is capable of
//! representing arbitrarily complex polytopal meshes in 2D and 3D. It
//! can answer topological queries about cells, faces, edges and nodes
//! as well subcell entities called corners and wedges
//! (iotas). Additionally, Jali can answer geometric queries about the
//! entities it supports. It can handle distributed meshes scalably
//! upto thousands of processors. Jali is built upon the open source
//! unstructured mesh infrastructure library, MSTK, developed at Los
//! Alamos National Laboratory since 2004. In addition, it is possible
//! to use Jali with STKmesh from Sandia National Labs and MOAB from
//! Argonne National Labs although support for these frameworks is
//! limited.
//!
//! In addition, to the mesh representation, the Jali library also
//! includes a rudimentary geometric model representation consisting of
//! simple geometric regions. Geometric regions are used to create mesh
//! sets (sets of cells, sets of faces, etc.) for specification of
//! material properties, boundary conditions and initial conditions.
//!
//! Finally, the Jali library contains a simple state manager for
//! storing and retrieving field data on mesh entities. Currently, the
//! state manager does not have any capability to synchronize data
//! across processors.
//!
//! The main classes of relevance to developers in the Jali package are:
//!  - Jali::Mesh                     in file       Mesh.hh
//!  - Jali::MeshFactory              in file       MeshFactory.hh
//!  - Jali::State                    in file       JaliState.h
//!  - JaliGeometry::GeometricModel   in file       GeometricModel.hh
//! .

//  RegionFactory will get included once we decide on a way to read XML
//  input files and generate parameter lists
//  - JaliGeometry::RegionFactory    in file       RegionFactory.hh

namespace Jali {

  //! \class Mesh.hh
  //! \brief Base mesh class
  //!
  //! Use the associated mesh factory to create an instance of a
  //! derived class based on a particular mesh framework (like MSTK,
  //! STKmesh etc.)
  //!
  //! **** IMPORTANT NOTE ABOUT CONSTANTNESS OF THIS CLASS ****
  //! Instantiating a const version of this class only guarantees that
  //! the underlying mesh topology and geometry does not change (the
  //! public interfaces conforms strictly to this definition). However,
  //! for purposes of memory savings we use lazy initialization and
  //! caching of face data, edge data, geometry quantities, columns
  //! etc., which means that these data may still change. We also
  //! cannot initialize the cached quantities in the constructor since
  //! they depend on initialization of data structures in the derived
  //! class - however, the base class gets constructed before the
  //! derived class gets constructed so it is not possible without more
  //! obscure acrobatics. This is why some of the caching data
  //! declarations are declared with the keyword 'mutable' and routines
  //! that modify the mutable data are declared with a constant
  //! qualifier.
  //!


class Mesh {
  
 public:
  
  //! \brief constructor
  //!
  //! constructor - cannot call directly. Code must set mesh framework
  //! preference to one of the available mesh frameworks (MSTK) and
  //! call the mesh_factory to make a mesh. If it is absolutely
  //! necessary, one can call the constructor of one of the available
  //! mesh frameworks directly
  
  Mesh(const bool request_faces = true,
       const bool request_edges = false,
       const bool request_sides = false,
       const bool request_wedges = false,
       const bool request_corners = false,
       const int num_tiles_ini = 0,
       const int num_ghost_layers_tile = 0,
       const int num_ghost_layers_distmesh = 1,
       const bool request_boundary_ghosts = false,
       const Partitioner_type partitioner = Partitioner_type::METIS,
       const JaliGeometry::Geom_type geom_type =
       JaliGeometry::Geom_type::CARTESIAN,
       const MPI_Comm incomm = MPI_COMM_WORLD) :
    spacedim(3), celldim(3), mesh_type_(Mesh_type::GENERAL),
    cell_geometry_precomputed(false), face_geometry_precomputed(false),
    edge_geometry_precomputed(false), side_geometry_precomputed(false),
    corner_geometry_precomputed(false),
    faces_requested(request_faces), edges_requested(request_edges),
    sides_requested(request_sides), wedges_requested(request_wedges),
    corners_requested(request_corners),
    num_tiles_ini_(num_tiles_ini),
    num_ghost_layers_tile_(num_ghost_layers_tile),
    num_ghost_layers_distmesh_(num_ghost_layers_distmesh),
    boundary_ghosts_requested_(request_boundary_ghosts),
    partitioner_pref_(partitioner),
    cell2face_info_cached(false), face2cell_info_cached(false),
    cell2edge_info_cached(false), face2edge_info_cached(false),
    side_info_cached(false), wedge_info_cached(false),
    corner_info_cached(false), type_info_cached(false),
    geometric_model_(NULL), comm(incomm),
    geomtype(geom_type) {
    
    if (corners_requested)  // corners are defined in terms of wedges
      wedges_requested = true;
    if (wedges_requested)   // wedges are defined in terms of sides
      sides_requested = true;
    if (sides_requested) {
      faces_requested = true;
      edges_requested = true;
    }
  }

  //! destructor
  //!
  //! destructor - must be virtual to downcast base class to derived class
  //! (I don't understand why but the stackoverflow prophets say so)

  virtual ~Mesh() {}

  //! MPI communicator being used

  inline
  MPI_Comm get_comm() const {
    return comm;
  }

  // Geometric type for the mesh - CARTESIAN, CYLINDRICAL, or SPHERICAL
  inline
  JaliGeometry::Geom_type geom_type() const {
    return geomtype;
  }

  //! Set the spatial dimension of points in the mesh - typically
  //! invoked by the constructor of a derived mesh class

  inline
  void set_space_dimension(const unsigned int dim) {
    spacedim = dim;
  }

  //! Spatial dimension of points in the mesh

  inline
  unsigned int space_dimension() const {
    return spacedim;
  }

  //! Set the topological dimension of mesh cells - typically
  //! invoked by the constructor of a derived mesh class

  inline
  void set_cell_dimension(const unsigned int dim) {
    celldim = dim;   // 3 is solid mesh, 2 is surface mesh, 1 is wire mesh
  }


  //! Topological dimension of mesh cells (hexes, tets, polyhedra have
  //! a cell dimension of 3; quads, triangles and polygons have a cell
  //! dimension of 2; line elements - NOT SUPPORTED - have a dimension
  //! of 1; particles - NOT SUPPORTED - have a dimension of 0)

  inline
  unsigned int cell_dimension() const {
    return celldim;
  }

  //! Set the pointer to a geometric model underpinning the mesh
  //! Typically, set by the constructor of a derived mesh class

  inline
  void set_geometric_model(const JaliGeometry::GeometricModelPtr &gm) {
    geometric_model_ = gm;
  }

  //! Return a pointer to a geometric model underpinning the mesh The
  //! geometric model consists of regions that are used to define mesh
  //! sets

  inline
  JaliGeometry::GeometricModelPtr geometric_model() const {
    return geometric_model_;
  }


  //! Set mesh type - Mesh_type::RECTANGULAR or Mesh_type::GENERAL

  inline
  void set_mesh_type(const Mesh_type mesh_type) {
    mesh_type_ = mesh_type;
  }

  //! Get mesh type - Mesh_type::RECTANGULAR or Mesh_type::GENERAL

  inline
  Mesh_type mesh_type() const {
    return mesh_type_;
  }

  //! Get type of entity - PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST

  Entity_type entity_get_type(const Entity_kind kind,
                              const Entity_ID entid) const;


  //! Parent entity in the source mesh if mesh was derived from another mesh

  virtual
  Entity_ID entity_get_parent(const Entity_kind kind, const Entity_ID entid)
      const;


  //! Get cell type - UNKNOWN, TRI, QUAD, POLYGON, TET, PRISM, PYRAMID, HEX,
  //! POLYHED
  //! See MeshDefs.hh

  virtual
  Cell_type cell_get_type(const Entity_ID cellid) const = 0;

  //
  // General mesh information
  // -------------------------
  //

  //! Number of entities of any kind (cell, face, node) and in a
  //! particular category (PARALLEL_OWNED, PARALLEL_GHOST, ALL)

  unsigned int num_entities(const Entity_kind kind,
                            const Entity_type type) const;

  //! Number of nodes of type (PARALLEL_OWNED, PARALLEL_GHOST, ALL)

  template<Entity_type type = Entity_type::ALL>
  unsigned int num_nodes() const;

  //! Number of edges of type (PARALLEL_OWNED, PARALLEL_GHOST, ALL)

  template<Entity_type type = Entity_type::ALL>
  unsigned int num_edges() const;

  //! Number of faces of type (PARALLEL_OWNED, PARALLEL_GHOST, ALL)

  template<Entity_type type = Entity_type::ALL>
  unsigned int num_faces() const;

  //! Number of sides of type (PARALLEL_OWNED, PARALLEL_GHOST, ALL)

  template<Entity_type type = Entity_type::ALL>
  unsigned int num_sides() const;

  //! Number of wedges of type (PARALLEL_OWNED, PARALLEL_GHOST, ALL)

  template<Entity_type type = Entity_type::ALL>
  unsigned int num_wedges() const;

  //! Number of corners of type (PARALLEL_OWNED, PARALLEL_GHOST, ALL)

  template<Entity_type type = Entity_type::ALL>
  unsigned int num_corners() const;

  //! Number of cells of type (PARALLEL_OWNED, PARALLEL_GHOST, ALL)

  template<Entity_type type = Entity_type::ALL>
  unsigned int num_cells() const;


  //! Global ID of any entity

  virtual
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const = 0;


  //! List of references to mesh tiles (collections of mesh cells)
  // Don't want to make the vector contain const references to tiles
  // because the tiles may be asked to add or remove some entities

  const std::vector<std::shared_ptr<MeshTile>> & tiles() {
    return meshtiles;
  }

  //! Number of mesh tiles (on a compute node)

  int num_tiles() const {return meshtiles.size();}

  //! Nodes of mesh (of a particular parallel type OWNED, GHOST or ALL)

  template<Entity_type type = Entity_type::ALL>
  const std::vector<Entity_ID> & nodes() const;

  //! Edges of mesh (of a particular parallel type OWNED, GHOST or ALL)

  template<Entity_type type = Entity_type::ALL>
  const std::vector<Entity_ID> & edges() const;

  //! Faces of mesh (of a particular parallel type OWNED, GHOST or ALL)

  template<Entity_type type = Entity_type::ALL>
  const std::vector<Entity_ID> & faces() const;

  //! Sides of mesh (of a particular parallel type OWNED, GHOST or ALL)

  template<Entity_type type = Entity_type::ALL>
  const std::vector<Entity_ID> & sides() const;

  //! Wedges of mesh (of a particular parallel type OWNED, GHOST or ALL)

  template<Entity_type type = Entity_type::ALL>
  const std::vector<Entity_ID> & wedges() const;

  //! Corners of mesh (of a particular parallel type OWNED, GHOST or ALL)

  template<Entity_type type = Entity_type::ALL>
  const std::vector<Entity_ID> & corners() const;

  //! Cells of mesh (of a particular parallel type OWNED, GHOST or ALL)

  template<Entity_type type = Entity_type::ALL>
  const std::vector<Entity_ID> & cells() const;


  // Master tile ID for entities

  int master_tile_ID_of_node(Entity_ID const nodeid) const {
    return tiles_initialized_ ? node_master_tile_ID_[nodeid] : -1;
  }
  int master_tile_ID_of_edge(Entity_ID const edgeid) const {
    return tiles_initialized_ ? edge_master_tile_ID_[edgeid] : -1;
  }
  int master_tile_ID_of_face(Entity_ID const faceid) const {
    return tiles_initialized_ ? face_master_tile_ID_[faceid] : -1;
  }
  int master_tile_ID_of_cell(Entity_ID const cellid) const {
    return tiles_initialized_ ? cell_master_tile_ID_[cellid] : -1;
  }
  int master_tile_ID_of_side(Entity_ID const sideid) const {
    return (tiles_initialized_ ?
            cell_master_tile_ID_[side_get_cell(sideid)] : -1);
  }
  int master_tile_ID_of_wedge(Entity_ID const wedgeid) const {
    return (tiles_initialized_ ?
            cell_master_tile_ID_[wedge_get_cell(wedgeid)] : -1);
  }
  int master_tile_ID_of_corner(Entity_ID const cornerid) const {
    return (tiles_initialized_ ?
            cell_master_tile_ID_[corner_get_cell(cornerid)] : -1);
  }


  //
  // Mesh Entity Adjacencies
  //-------------------------


  // Downward Adjacencies
  //---------------------

  //! Get number of faces in a cell

  unsigned int cell_get_num_faces(const Entity_ID cellid) const;

  //! Get faces of a cell.
  //!
  //! The Google coding guidelines regarding function arguments is purposely
  //! violated here to allow for a default input argument
  //!
  //! On a distributed mesh, this will return all the faces of the
  //! cell, OWNED or GHOST. If ordered = true, the faces will be
  //! returned in a standard order according to Exodus II convention
  //! for standard cells; in all other situations (ordered = false or
  //! non-standard cells), the list of faces will be in arbitrary order

  void cell_get_faces(const Entity_ID cellid,
                      Entity_ID_List *faceids,
                      const bool ordered = false) const;

  //! Get faces of a cell and directions in which the cell uses the face
  //!
  //! The Google coding guidelines regarding function arguments is purposely
  //! violated here to allow for a default input argument
  //!
  //! On a distributed mesh, this will return all the faces of the
  //! cell, OWNED or GHOST. If ordered = true, the faces will be
  //! returned in a standard order according to Exodus II convention
  //! for standard cells; in all other situations (ordered = false or
  //! non-standard cells), the list of faces will be in arbitrary order
  //!
  //! In 3D, direction is 1 if face normal points out of cell
  //! and -1 if face normal points into cell
  //! In 2D, direction is 1 if face/edge is defined in the same
  //! direction as the cell polygon, and -1 otherwise

  void cell_get_faces_and_dirs(const Entity_ID cellid,
                               Entity_ID_List *faceids,
                               std::vector<dir_t> *facedirs,
                               const bool ordered = false) const;


  //! Get edges of a cell (in no particular order)

  void cell_get_edges(const Entity_ID cellid,
                      Entity_ID_List *edgeids) const;

  //! Get edges and dirs of a 2D cell. This is to make the code cleaner
  //! for integrating over the cell in 2D where faces and edges are
  //! identical but integrating over the cells using face information
  //! is more cumbersome (one would have to take the face normals,
  //! rotate them and then get a consistent edge vector)

  void cell_2D_get_edges_and_dirs(const Entity_ID cellid,
                                  Entity_ID_List *edgeids,
                                  std::vector<dir_t> *edge_dirs) const;

  //! Get nodes of a cell (in no particular order)

  virtual
  void cell_get_nodes(const Entity_ID cellid,
                      Entity_ID_List *nodeids) const = 0;


  //! Get edges of a face and directions in which the face uses the edges
  //!
  //! On a distributed mesh, this will return all the edges of the
  //! face, OWNED or GHOST. If ordered = true, the edges will be
  //! returned in a ccw order around the face as it is naturally defined.
  //!
  //! IMPORTANT NOTE IN 2D CELLS: In meshes where the cells are two
  //! dimensional, faces and edges are identical. For such cells, this
  //! operator will return a single edge and a direction of 1. However,
  //! this direction cannot be relied upon to compute, say, a contour
  //! integral around the 2D cell.

  void face_get_edges_and_dirs(const Entity_ID faceid,
                               Entity_ID_List *edgeids,
                               std::vector<dir_t> *edgedirs,
                               const bool ordered = false) const;


  //! Get the local index of a face edge in a cell edge list
  //! Example:
  //!
  //! face_get_edges(face=5) --> {20, 21, 35, 9, 10}
  //! cell_get_edges(cell=18) --> {1, 2, 3, 5, 8, 9, 10, 13, 21, 35, 20, 37, 40}
  //! face_to_cell_edge_map(face=5,cell=18) --> {10, 8, 9, 5, 6}


  void face_to_cell_edge_map(const Entity_ID faceid,
                             const Entity_ID cellid,
                             std::vector<int> *map) const;

  //! Get nodes of face
  //! On a distributed mesh, all nodes (OWNED or GHOST) of the face
  //! are returned
  //! In 3D, the nodes of the face are returned in ccw order consistent
  //! with the face normal
  //! In 2D, nfnodes is 2

  virtual
  void face_get_nodes(const Entity_ID faceid,
                      Entity_ID_List *nodeids) const = 0;


  //! Get nodes of edge

  void edge_get_nodes(const Entity_ID edgeid,
                      Entity_ID *nodeid0, Entity_ID *nodeid1) const;

  
  //! Get sides of a cell (in no particular order)

  void cell_get_sides(const Entity_ID cellid,
                      Entity_ID_List *sideids) const;

  //! Get wedges of cell (in no particular order)

  void cell_get_wedges(const Entity_ID cellid,
                       Entity_ID_List *wedgeids) const;

  //! Get corners of cell (in no particular order)

  void cell_get_corners(const Entity_ID cellid,
                        Entity_ID_List *cornerids) const;

  //! Get corner at cell and node combination

  Entity_ID cell_get_corner_at_node(const Entity_ID cellid,
                                    const Entity_ID nodeid) const;


  //! Face of a side

  Entity_ID side_get_cell(const Entity_ID sideid) const;

  //! Face of a side

  Entity_ID side_get_face(const Entity_ID sideid) const;

  //! Edge of a side

  Entity_ID side_get_edge(const Entity_ID sideid) const;

  //! Sense in which side is using its edge (i.e. do p0, p1 of side,
  //! edge match or not)

  int side_get_edge_use(const Entity_ID sideid) const;

  //! Node of a side (each edge of a side has two nodes - inode (0,1)
  //! indicates which one to return)

  Entity_ID side_get_node(const Entity_ID sideid, int inode) const;

  //! Wedge of a side 
  //! Each side points to two wedges - iwedge (0, 1)
  //! indicates which one to return; the wedge returned will be
  //! consistent with the node returned by side_get_node.
  //! So, side_get_node(s,i) = wedge_get_node(side_get_wedge(s,i))

  Entity_ID side_get_wedge(const Entity_ID sideid, int iwedge) const;

  //! Face of a wedge

  Entity_ID wedge_get_face(const Entity_ID wedgeid) const;

  //! Edge of a wedge

  Entity_ID wedge_get_edge(const Entity_ID wedgeid) const;

  //! Node of a wedge

  Entity_ID wedge_get_node(const Entity_ID wedgeid) const;

  //! Node of a corner

  Entity_ID corner_get_node(const Entity_ID cornerid) const;

  //! Wedges of a corner

  void corner_get_wedges(const Entity_ID cornerid,
                         Entity_ID_List *wedgeids) const;

  //! Face get facets (or should we return a vector of standard pairs
  //! containing the wedge and a facet index?)

  void face_get_facets(const Entity_ID faceid,
                       Entity_ID_List *facetids) const;

  // Upward adjacencies
  //-------------------

  //! Cells of type 'type' connected to a node - The order of cells
  //! is not guaranteed to be the same for corresponding nodes on
  //! different processors

  virtual
  void node_get_cells(const Entity_ID nodeid,
                      const Entity_type type,
                      Entity_ID_List *cellids) const = 0;


  //! Faces of type 'type' connected to a node - The order of faces
  //! is not guaranteed to be the same for corresponding nodes on
  //! different processors

  virtual
  void node_get_faces(const Entity_ID nodeid,
                      const Entity_type type,
                      Entity_ID_List *faceids) const = 0;

  //! Wedges connected to a node - The wedges are returned in no
  //! particular order. Also, the order of nodes is not guaranteed to
  //! be the same for corresponding nodes on different processors

  void node_get_wedges(const Entity_ID nodeid,
                       const Entity_type type,
                       Entity_ID_List *wedgeids) const;

  //! Corners connected to a node - The corners are returned in no
  //! particular order. Also, the order of corners is not guaranteed to
  //! be the same for corresponding nodes on different processors

  void node_get_corners(const Entity_ID nodeid,
                        const Entity_type type,
                        Entity_ID_List *cornerids) const;

  //! Get faces of type of a particular cell that are connected to the
  //! given node - The order of faces is not guarnateed to be the same
  //! for corresponding nodes on different processors

  virtual
  void node_get_cell_faces(const Entity_ID nodeid,
                           const Entity_ID cellid,
                           const Entity_type type,
                            Entity_ID_List *faceids) const = 0;

  //! Cells connected to a face - The cells are returned in no
  //! particular order. Also, the order of cells is not guaranteed to
  //! be the same for corresponding faces on different processors

  void face_get_cells(const Entity_ID faceid,
                      const Entity_type type,
                      Entity_ID_List *cellids) const;

  //! Cell of a wedge

  Entity_ID wedge_get_cell(const Entity_ID wedgeid) const;

  //! Side of a wedge

  Entity_ID wedge_get_side(const Entity_ID wedgeid) const;

  //! Corner of a wedge

  Entity_ID wedge_get_corner(const Entity_ID wedgeid) const;

  //! wedges of a facet

  // void wedges_of_a_facet (const Entity_ID facetid, Entity_ID_List *wedgeids)
  //    const;

  //! Cell of a corner

  Entity_ID corner_get_cell(const Entity_ID cornerid) const;


  // Same level adjacencies
  //-----------------------

  //! Face connected neighboring cells of given cell of a particular type
  //! (e.g. a hex has 6 face neighbors)
  //!
  //! The order in which the cellids are returned cannot be
  //! guaranteed in general except when type = ALL, in which case
  //! the cellids will correcpond to cells across the respective
  //! faces given by cell_get_faces

  virtual
  void cell_get_face_adj_cells(const Entity_ID cellid,
                               const Entity_type type,
                               Entity_ID_List *fadj_cellids) const = 0;

  //! Node connected neighboring cells of given cell
  //! (a hex in a structured mesh has 26 node connected neighbors)
  //! The cells are returned in no particular order

  virtual
  void cell_get_node_adj_cells(const Entity_ID cellid,
                               const Entity_type type,
                               Entity_ID_List *cellids) const = 0;


  //! Opposite side in neighboring cell of a side. The two sides share
  //! facet 0 of wedge comprised of nodes 0,1 of the common edge and
  //! center point of the common face in 3D, and nodes 0,1 of the
  //! common edge in 2D. At boundaries, this routine returns -1

  Entity_ID side_get_opposite_side(const Entity_ID wedgeid) const;


  //! Opposite wedge in neighboring cell of a wedge. The two wedges
  //! share facet 0 of wedge comprised of the node, center point of
  //! the common edge and center point of the common face in 3D, and
  //! node and edge center in 2D. At boundaries, this routine returns
  //! -1

  Entity_ID wedge_get_opposite_wedge(const Entity_ID wedgeid) const;

  //! adjacent wedge along edge in the same cell. The two wedges share
  //! facet 1 of wedge comprised of edge center, face center and zone center
  //! in 3D, and node and zone center in 2D


  Entity_ID wedge_get_adjacent_wedge(const Entity_ID wedgeid) const;


  //
  // Mesh entity geometry
  //--------------
  //

  //! Node coordinates

  // Preferred operator
  virtual
  void node_get_coordinates(const Entity_ID nodeid,
                            JaliGeometry::Point *ncoord) const = 0;

  virtual
  void node_get_coordinates(const Entity_ID nodeid,
                            std::array<double, 3> *ncoord) const;
  virtual
  void node_get_coordinates(const Entity_ID nodeid,
                            std::array<double, 2> *ncoord) const;
  virtual
  void node_get_coordinates(const Entity_ID nodeid, double *ncoord) const;

  //! Face coordinates - conventions same as face_to_nodes call
  //! Number of nodes is the vector size divided by number of spatial dimensions

  virtual
  void face_get_coordinates(const Entity_ID faceid,
                            std::vector<JaliGeometry::Point> *fcoords)
      const = 0;

  //! Coordinates of cells in standard order (Exodus II convention)
  //!
  //! STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  //! For a general polyhedron this will return the node coordinates in
  //! arbitrary order
  //! Number of nodes is vector size divided by number of spatial dimensions

  virtual
  void cell_get_coordinates(const Entity_ID cellid,
                            std::vector<JaliGeometry::Point> *ccoords)
      const = 0;

  //! Coordinates of side
  //!
  //! The coordinates will be returned in a fixed
  //! order - If posvol_order is true, the node coordinates will be
  //! ordered such that they will result in a +ve volume calculation;
  //! otherwise, they will be returned in the natural order in which
  //! they are defined for the side. For sides, the natural order
  //! automatically gives positive volume for 2D and 3D - its only in
  //! 1D that some side coordinates have to be reordered (see wedges
  //! below for which one wedge of each side will give a -ve volume if
  //! the natural coordinate order is used)

  void side_get_coordinates(const Entity_ID sideid,
                            std::vector<JaliGeometry::Point> *scoords,
                            bool posvol_order = false) const;

  //! Coordinates of wedge
  //!
  //! If posvol_order = true, then the coordinates will be returned
  //! in an order that will result in a positive volume (in 3D this assumes
  //! that the computation for volume is done as (V01 x V02).V03 where V0i
  //! is a vector from coordinate 0 to coordinate i of the tet). If posvol_order
  //! is false, the coordinates will be returned in a fixed order - in 2D,
  //! this is node point, edge/face center, cell center and in 3D, this is
  //! node point, edge center, face center, cell center
  //!
  //! By default the coordinates are returned in the natural order
  //! (posvol_order = false)


  void wedge_get_coordinates(const Entity_ID wedgeid,
                             std::vector<JaliGeometry::Point> *wcoords,
                             bool posvol_order = false) const;


  //! Coordinates of corner points. In 2D, these are ordered in a ccw
  //! manner. In 3D, they are not ordered in any particular way and
  //! this routine may not be too useful since the topology of the
  //! corner is not guaranteed to be standard like in 2D. Its better
  //! to work with the wedges of the corner

  void corner_get_coordinates(const Entity_ID cornerid,
                              std::vector<JaliGeometry::Point> *cncoords) const;

  //! Get a facetized description of corner geometry in 3D. The facet
  //! points index into the pointcoords vector. Each facet is
  //! guaranteed to have its points listed such that its normal points
  //! out of the corner

  void
  corner_get_facetization(const Entity_ID cornerid,
                          std::vector<JaliGeometry::Point> *pointcoords,
                          std::vector<std::array<Entity_ID, 3>> *facetpoints)
      const;

  //! "facets" (line segments) describing a corner in 2D. The facet points are
  //! (0,1) (1,2) (2,3) and (3,4) referring to the point coordinates. They are
  //! guaranteed to be in ccw order around the quadrilateral corner

  void
  corner_get_facetization(const Entity_ID cornerid,
                          std::vector<JaliGeometry::Point> *pointcoords,
                          std::vector<std::array<Entity_ID, 2>>
                          *facetpoints) const;

// "facets" (points - node and cell center) describing a corner in 1D :)

  void
  corner_get_facetization(const Entity_ID cornerid,
                          std::vector<JaliGeometry::Point> *pointcoords,
                          std::vector<std::array<Entity_ID, 1>>
                          *facetpoints) const;

  // Mesh entity geometry
  //--------------
  //


  //! Volume/Area of cell

  double
  cell_volume(const Entity_ID cellid, const bool recompute = false) const;

  //! Area/length of face

  double
  face_area(const Entity_ID faceid, const bool recompute = false) const;

  //! Length of edge

  double
  edge_length(const Entity_ID edgeid, const bool recompute = false) const;

  //! Volume of side

  double
  side_volume(const Entity_ID sideid, const bool recompute = false) const;

  //! Volume of wedge

  double
  wedge_volume(const Entity_ID wedgeid, const bool recompute = false) const;

  //! Volume of a corner

  double
  corner_volume(const Entity_ID cornerid, const bool recompute = false) const;

  //! Centroid of cell

  JaliGeometry::Point
  cell_centroid(const Entity_ID cellid, const bool recompute = false) const;

  //! Centroid of face

  JaliGeometry::Point
  face_centroid(const Entity_ID faceid, const bool recompute = false) const;

  //! Centroid/center of edge(never cached)

  JaliGeometry::Point edge_centroid(const Entity_ID edgeid) const;

  //! Normal to face
  //! The vector is normalized and then weighted by the area of the face
  //!
  //! If recompute is TRUE, then the normal is recalculated using current
  //! face coordinates but not stored. (If the recomputed normal must be
  //! stored, then call recompute_geometric_quantities).
  //!
  //! If cellid is not specified, the normal is the natural normal of
  //! the face. This means that at boundaries, the normal may point in
  //! or out of the domain depending on how the face is defined. On the
  //! other hand, if cellid is specified, the normal is the outward
  //! normal with respect to the cell. In planar and solid meshes, the
  //! normal with respect to the cell on one side of the face is just
  //! the negative of the normal with respect to the cell on the other
  //! side. In general surfaces meshes, this will not be true at C1
  //! discontinuities

  //! if cellid is specified, then orientation returns the direction of
  //! the natural normal of the face with respect to the cell (1 is
  //! pointing out of the cell and -1 pointing in)


  JaliGeometry::Point face_normal(const Entity_ID faceid,
                                     const bool recompute = false,
                                     const Entity_ID cellid = -1,
                                     int *orientation = NULL) const;


  //! Edge vector - not normalized (or normalized and weighted by length
  //! of the edge)
  //!
  //! If recompute is TRUE, then the vector is recalculated using current
  //! edge coordinates but not stored. (If the recomputed vector must be
  //! stored, then call recompute_geometric_quantities).
  //!
  //! If pointid is specified, the vector is the natural direction of
  //! the edge (from point0 to point1).  On the other hand, if pointid
  //! is specified (has to be a point of the face), the vector is from
  //! specified point to opposite point of edge.
  //!
  //! if pointid is specified, then orientation returns the direction of
  //! the natural direction of the edge with respect to the point (1 is
  //! away from the point and -1 is towards)


  JaliGeometry::Point edge_vector(const Entity_ID edgeid,
                                     const bool recompute = false,
                                     const Entity_ID pointid = -1,
                                     int *orientation = NULL) const;

  //! Point in cell?

  bool point_in_cell(const JaliGeometry::Point &p,
                      const Entity_ID cellid) const;


  //! Outward normal to facet of side that is shared with side from
  //! neighboring cell.
  //!
  //! The vector is normalized and then weighted by the area of the
  //! face. If recompute is TRUE, then the normal is recalculated
  //! using current wedge coordinates but not stored. (If the
  //! recomputed normal must be stored, then call
  //! recompute_geometric_quantities).
  //!
  //! Normals of other facets are typically not used in discretizations

  JaliGeometry::Point side_facet_normal(const Entity_ID sideid,
                                        const bool recompute = false) const;

  //! Outward normal to facet of wedge.
  //!
  //! The vector is normalized and then weighted by the area of the
  //! face Typically, one would ask for only facet 0 (shared with
  //! wedge from neighboring cell) and facet 1 (shared with adjacent
  //! wedge in the same side). If recompute is TRUE, then the normal
  //! is recalculated using current wedge coordinates but not
  //! stored. (If the recomputed normal must be stored, then call
  //! recompute_geometric_quantities).
  //!
  //! Normals of other facets are typically not used in discretizations

  JaliGeometry::Point wedge_facet_normal(const Entity_ID wedgeid,
                                          const unsigned int which_facet,
                                          const bool recompute = false) const;


  //
  // Mesh modification
  //-------------------

  //! Set coordinates of node

  virtual
  void node_set_coordinates(const Entity_ID nodeid,
                             const JaliGeometry::Point ncoord) = 0;


  virtual
  void node_set_coordinates(const Entity_ID nodeid,
                             const double *ncoord) = 0;


  //! Update geometric quantities (volumes, normals, centroids, etc.)
  //! and cache them - called for initial caching or for update after
  //! mesh modification

  void update_geometric_quantities();

  //
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------
  //

  //! Is this is a valid name of a set containing entities of 'kind'

  bool valid_set_name(const std::string setname,
                       const Entity_kind kind) const;


  //! Get number of entities of type 'category' in set

  virtual
  unsigned int get_set_size(const Set_Name setname,
                             const Entity_kind kind,
                             const Entity_type type) const = 0;

  virtual
  unsigned int get_set_size(const char *setname,
                             const Entity_kind kind,
                             const Entity_type type) const = 0;


  //! Get list of entities of type 'category' in set

  virtual
  void get_set_entities(const Set_Name setname,
                         const Entity_kind kind,
                         const Entity_type type,
                         Entity_ID_List *entids) const = 0;

  virtual
  void get_set_entities(const char *setname,
                         const Entity_kind kind,
                         const Entity_type type,
                         Entity_ID_List *entids) const = 0;



  //! \brief Export to Exodus II file
  //! Export mesh to Exodus II file. If with_fields is true, the fields in
  //! JaliState are also exported out.

  virtual
  void write_to_exodus_file(const std::string exodusfilename,
                            const bool with_fields = true) const {}

  //! \brief Export to GMV file
  //! Export mesh to GMV file. If with_fields is true, the fields in
  //! JaliState are also exported out.

  virtual
  void write_to_gmv_file(const std::string gmvfilename,
                         const bool with_fields = true) const {}

  //! \brief Precompute and cache corners, wedges, edges, cells
  // WHY IS THIS VIRTUAL?

  virtual
  void cache_extra_variables();

 protected:

  int compute_cell_geometric_quantities() const;
  int compute_face_geometric_quantities() const;
  int compute_edge_geometric_quantities() const;
  int compute_side_geometric_quantities() const;
  int compute_corner_geometric_quantities() const;


  // get faces of a cell and directions in which it is used - this function
  // is implemented in each mesh framework. The results are cached in
  // the base class

  virtual
  void cell_get_faces_and_dirs_internal(const Entity_ID cellid,
                                        Entity_ID_List *faceids,
                                        std::vector<dir_t> *face_dirs,
                                        const bool ordered = false) const = 0;

  // Cells connected to a face - this function is implemented in each
  // mesh framework. The results are cached in the base class

  virtual
  void face_get_cells_internal(const Entity_ID faceid,
                                const Entity_type type,
                                Entity_ID_List *cellids) const = 0;


  // edges of a face - this function is implemented in each mesh
  // framework. The results are cached in the base class

  virtual
  void face_get_edges_and_dirs_internal(const Entity_ID faceid,
                                        Entity_ID_List *edgeids,
                                        std::vector<dir_t> *edge_dirs,
                                        const bool ordered = true) const = 0;

  // edges of a cell - this function is implemented in each mesh
  // framework. The results are cached in the base class.

  virtual
  void cell_get_edges_internal(const Entity_ID cellid,
                                Entity_ID_List *edgeids) const = 0;

  // edges and directions of a 2D cell - this function is implemented
  // in each mesh framework. The results are cached in the base class.

  virtual
  void
  cell_2D_get_edges_and_dirs_internal(const Entity_ID cellid,
                                      Entity_ID_List *edgeids,
                                      std::vector<dir_t> *edge_dirs) const = 0;


  // get nodes of an edge - virtual function that will be implemented
  // in each mesh framework. The results are cached in the base class

  virtual
  void edge_get_nodes_internal(const Entity_ID edgeid,
                               Entity_ID *enode0, Entity_ID *enode1) const = 0;

  //! \brief Get info about mesh fields on a particular type of entity

  //! Get info about the number of fields, their names and their types
  //! on a particular type of entity on the mesh - DESIGNED TO BE
  //! CALLED ONLY BY THE JALI STATE MANAGER FOR INITIALIZATION OF MESH
  //! STATE FROM THE MESH FILE

  virtual
  void get_field_info(Entity_kind on_what, int *num,
                      std::vector<std::string> *varnames,
                      std::vector<std::string> *vartypes) const {*num = 0;}

  //! \brief Retrieve a field on the mesh - cannot template virtual funcs

  //! Retrieve a field on the mesh. If the return value is false, it
  //! could be that (1) the field does not exist (2) it exists but is
  //! associated with a different type of entity (3) the variable type
  //! sent in was the wrong type (int instead of double or double
  //! instead of std::array<double,2> or std::array<double,2> instead
  //! of std::array<double,3> etc - DESIGNED TO BE CALLED ONLY BY THE
  //! JALI STATE MANAGER FOR INITIALIZATION OF MESH STATE FROM THE
  //! MESH FILE

  virtual
  bool get_field(std::string field_name, Entity_kind on_what, int *data) const
  {return false;}
  virtual
  bool get_field(std::string field_name, Entity_kind on_what,
                 double *data) const {return false;}
  virtual
  bool
  get_field(std::string field_name, Entity_kind on_what,
            std::array<double, (std::size_t)2> *data) const {return false;}
  virtual
  bool get_field(std::string field_name, Entity_kind on_what,
                 std::array<double, (std::size_t)3> *data) const {return false;}
  virtual
  bool get_field(std::string field_name, Entity_kind on_what,
                 std::array<double, (std::size_t)6> *data) const {return false;}

  //! \brief Store a field on the mesh - cannot template as its virtual

  //! Store a field on the mesh. If the return value is false, it
  //! means that the mesh already has a field of that name but its of
  //! a different type or its on a different type of entity - DESIGNED
  //! TO BE CALLED ONLY BY THE JALI STATE MANAGER FOR INITIALIZATION
  //! OF MESH STATE FROM THE MESH FILE

  virtual
  bool store_field(std::string field_name, Entity_kind on_what, int *data)
  {return false;}
  virtual
  bool store_field(std::string field_name, Entity_kind on_what, double *data)
  {return false;}
  virtual
  bool store_field(std::string field_name, Entity_kind on_what,
                   std::array<double, (std::size_t)2> *data) {return false;}
  virtual
  bool store_field(std::string field_name, Entity_kind on_what,
                   std::array<double, (std::size_t)3> *data) {return false;}
  virtual
  bool store_field(std::string field_name, Entity_kind on_what,
                   std::array<double, (std::size_t)6> *data) {return false;}

  // The following methods are declared const since they do not modify the
  // mesh but just modify cached variables declared as mutable

  int compute_cell_geometry(const Entity_ID cellid,
                            double *volume,
                            JaliGeometry::Point *centroid) const;
  int compute_face_geometry(const Entity_ID faceid,
                            double *area,
                            JaliGeometry::Point *centroid,
                            JaliGeometry::Point *normal0,
                            JaliGeometry::Point *normal1) const;
  int compute_edge_geometry(const Entity_ID edgeid,
                            double *length,
                            JaliGeometry::Point *edge_vector,
                            JaliGeometry::Point *centroid) const;

  // The outward_facet_normal is the area-weighted normal of the side
  // that lies on the boundary of the cell. This normal points out of
  // the cell. The mid_facet_normal is the normal of the common facet
  // between the two wedges of the side. This normal points out of
  // wedge 0 of the side and into wedge 1

  void compute_side_geometry(const Entity_ID sideid,
                             double *volume,
                             JaliGeometry::Point *outward_facet_normal,
                             JaliGeometry::Point *mid_facet_normal) const;

  void compute_corner_geometry(const Entity_ID cornerid,
                              double *volume) const;

  void cache_type_info() const;
  void cache_cell2face_info() const;
  void cache_face2cell_info() const;
  void cache_cell2edge_info() const;
  void cache_face2edge_info() const;
  void cache_edge2node_info() const;
  void cache_side_info() const;
  void cache_wedge_info() const;
  void cache_corner_info() const;

  void build_tiles();
  void add_tile(std::shared_ptr<MeshTile> tile2add);
  void init_tiles();
  int get_new_tile_ID() const { return meshtiles.size(); }

  // Set master tile ID for entities

  void set_master_tile_ID_of_node(Entity_ID const nodeid,
                                 int const tileid) {
    node_master_tile_ID_[nodeid] = tileid;
  }
  void set_master_tile_ID_of_edge(Entity_ID const edgeid,
                                 int const tileid) {
    edge_master_tile_ID_[edgeid] = tileid;
  }
  void set_master_tile_ID_of_face(Entity_ID const faceid,
                                 int const tileid) {
    face_master_tile_ID_[faceid] = tileid;
  }
  void set_master_tile_ID_of_cell(Entity_ID const cellid,
                                 int const tileid) {
    cell_master_tile_ID_[cellid] = tileid;
  }
  void set_master_tile_ID_of_wedge(Entity_ID const wedgeid,
                                   int const tileid) {}
  void set_master_tile_ID_of_side(Entity_ID const sideid,
                                   int const tileid) {}
  void set_master_tile_ID_of_corner(Entity_ID const cornerid,
                                    int const tileid) {}

  //! @brief Get the partitioning of a regular mesh such that each
  //! partition is a rectangular block
  //!
  //! @param dim Dimension of problem - 1, 2 or 3
  //! @param domain 2*dim values for min/max of domain
  //!  (xmin, xmax, ymin, ymax, zmin, zmax)
  //! @param num_cells_in_dir  number of cells in each direction
  //! @param num_blocks_requested number of blocks requested
  //! @param blocklimits min/max limits for each block
  //! @param blocknumcells num cells in each direction for blocks
  //!
  //! Returns 1 if successful, 0 otherwise
  
  int
  block_partition_regular_mesh(int const dim,
                               double const * const domain,
                               int const * const num_cells_in_dir,
                               int const num_blocks_requested,
                               std::vector<std::array<double, 6>> *blocklimits,
                               std::vector<std::array<int, 3>> *blocknumcells);


  // Data

  unsigned int celldim, spacedim;
  JaliGeometry::Geom_type geomtype = JaliGeometry::Geom_type::CARTESIAN;

  MPI_Comm comm;

  // MeshTile data (A meshtile is a list of cell indices that will be
  // processed together)

  const int num_tiles_ini_;
  const int num_ghost_layers_tile_;
  const int num_ghost_layers_distmesh_;
  const bool boundary_ghosts_requested_;
  const Partitioner_type partitioner_pref_;
  bool tiles_initialized_ = false;
  std::vector<std::shared_ptr<MeshTile>> meshtiles;
  std::vector<int> node_master_tile_ID_, edge_master_tile_ID_;
  std::vector<int> face_master_tile_ID_, cell_master_tile_ID_;

  // Some geometric quantities

  mutable std::vector<double> cell_volumes, face_areas, edge_lengths,
    side_volumes, corner_volumes;
  mutable std::vector<JaliGeometry::Point> cell_centroids,
    face_centroids, face_normal0, face_normal1, edge_vectors, edge_centroids;

  // outward facing normal from side to side in adjacent cell
  mutable std::vector<JaliGeometry::Point> side_outward_facet_normal;
  // Normal of the common facet of the two wedges - normal points out
  // of wedge 0 of side into wedge 1
  mutable std::vector<JaliGeometry::Point> side_mid_facet_normal;

  // Entity lists

  mutable std::vector<int> nodeids_owned_, nodeids_ghost_, nodeids_all_;
  mutable std::vector<int> edgeids_owned_, edgeids_ghost_, edgeids_all_;
  mutable std::vector<int> faceids_owned_, faceids_ghost_, faceids_all_;
  mutable std::vector<int> sideids_owned_, sideids_ghost_,
    sideids_boundary_ghost_, sideids_all_;
  mutable std::vector<int> wedgeids_owned_, wedgeids_ghost_,
    wedgeids_boundary_ghost_, wedgeids_all_;
  mutable std::vector<int> cornerids_owned_, cornerids_ghost_,
    cornerids_boundary_ghost_, cornerids_all_;
  mutable std::vector<int> cellids_owned_, cellids_ghost_,
    cellids_boundary_ghost_, cellids_all_;
  std::vector<int> dummy_list_;  // for unspecialized cases

  // Type info for essential entities - sides, wedges and corners will
  // get their type from their owning cell

  mutable std::vector<Entity_type> cell_type;
  mutable std::vector<Entity_type> face_type;  // if faces requested
  mutable std::vector<Entity_type> edge_type;  // if edges requested
  mutable std::vector<Entity_type> node_type;

  // Some standard topological relationships that are cached. The rest
  // are computed on the fly or obtained from the derived class

  mutable std::vector<Entity_ID_List> cell_face_ids;
  mutable std::vector<std::vector<dir_t>> cell_face_dirs;
  mutable std::vector<Entity_ID_List> face_cell_ids;
  mutable std::vector<Entity_ID_List> cell_edge_ids;
  mutable std::vector<Entity_ID_List> face_edge_ids;
  mutable std::vector<std::vector<dir_t>> face_edge_dirs;
  mutable std::vector<std::array<Entity_ID, 2>> edge_node_ids;

  // cell_2D_edge_dirs is an unusual topological relationship
  // requested by MHD discretization - It has no equivalent in 3D

  mutable std::vector<std::vector<dir_t>> cell_2D_edge_dirs;


  // Topological relationships involving standard and non-standard
  // entities (sides, corners and wedges). The non-standard entities
  // may be required for polyhedral elements and more accurate
  // discretizations.
  //
  // 1D:
  // A side is a line segment from a node to the cell. Wedges and
  // corners are the same as sides.
  //
  // 2D:
  // A side is a triangle formed by the two nodes of an edge/face and
  // the cell center. A wedge is half of a side formed by one node of
  // the edge, the edge center and the cell center. A corner is a
  // quadrilateral formed by the two wedges in a cell at a node
  //
  // 3D:
  // A side is a tet formed by the two nodes of an edge, a face center
  // and a cell center. A wedge is half a side, formed by a node of
  // the edge, the edge center, the face center and the cell center. A
  // corner is formed by all the wedges of a cell at a node.

  // Sides
  mutable std::vector<Entity_ID> side_cell_id;
  mutable std::vector<Entity_ID> side_face_id;
  mutable std::vector<Entity_ID> side_edge_id;
  mutable std::vector<bool> side_edge_use;  // true: side, edge - p0, p1 match
  mutable std::vector<std::array<Entity_ID, 2>> side_node_ids;
  mutable std::vector<Entity_ID> side_opp_side_id;

  // Wedges - most wedge info is derived from sides
  mutable std::vector<Entity_ID> wedge_corner_id;

  // some other one-many adjacencies
  mutable std::vector<std::vector<Entity_ID>> cell_side_ids;
  mutable std::vector<std::vector<Entity_ID>> cell_corner_ids;
  //  mutable std::vector<std::vector<Entity_ID>> edge_side_ids;
  mutable std::vector<std::vector<Entity_ID>> node_corner_ids;
  mutable std::vector<std::vector<Entity_ID>> corner_wedge_ids;

  // Rectangular or general
  mutable Mesh_type mesh_type_;

  // flags to indicate what data is current

  mutable bool faces_requested, edges_requested, sides_requested,
    wedges_requested, corners_requested;
  mutable bool type_info_cached;
  mutable bool cell2face_info_cached, face2cell_info_cached;
  mutable bool cell2edge_info_cached, face2edge_info_cached;
  mutable bool edge2node_info_cached;
  mutable bool side_info_cached, wedge_info_cached, corner_info_cached;
  mutable bool cell_geometry_precomputed, face_geometry_precomputed,
    edge_geometry_precomputed, side_geometry_precomputed,
    corner_geometry_precomputed;

  // Pointer to geometric model that contains descriptions of
  // geometric regions - These geometric regions are used to define
  // entity sets for properties, boundary conditions etc.

  JaliGeometry::GeometricModelPtr geometric_model_;


  //! Make the State class a friend so that it can access protected
  //! methods for retrieving and storing mesh fields

  friend class State;

  //! Make the MeshTile class a friend so it can access protected functions
  //! for getting a new tile ID

  friend class MeshTile;

  // Make the make_meshtile function a friend so that it can access
  // the protected functions init_tiles and add_tile

  friend
  std::shared_ptr<MeshTile> make_meshtile(Mesh& parent_mesh,
                                          std::vector<Entity_ID> const& cells,
                                          int const num_ghost_layers_tile,
                                          bool const request_faces,
                                          bool const request_edges,
                                          bool const request_sides,
                                          bool const request_wedges,
                                          bool const request_corners);

 private:

  /// Method to get partitioning of a mesh into num parts

  void get_partitioning(int const num_parts,
                        Partitioner_type const parttype,
                        std::vector<std::vector<int>> *partitions);

  /// Method to get crude partitioning by chopping up the index space

  void get_partitioning_by_index_space(int const num_parts,
                                       std::vector<std::vector<int>> *partitions);

  /// Method to get crude partitioning by subdivision into rectangular blocks

  void get_partitioning_by_blocks(int const num_parts,
                                       std::vector<std::vector<int>> *partitions);

  /// Method to get partitioning of a mesh into num parts using METIS

#ifdef Jali_HAVE_METIS
  void get_partitioning_with_metis(int const num_parts,
                                   std::vector<std::vector<int>> *partitions);
#endif

  /// Method to get partitioning of a mesh into num parts using ZOLTAN

#ifdef Jali_HAVE_ZOLTAN
  void get_partitioning_with_zoltan(int const num_parts,
                                    std::vector<std::vector<int>> *partitions);
#endif

};  // End class Mesh



// Templated version of num_entities with specializations on enum

template<Entity_type type> inline
unsigned int Mesh::num_nodes() const {
  std::cerr << "num_nodes: Not defined for parallel type " << type << "\n";
  return 0;
}
template<> inline
unsigned int
Mesh::num_nodes<Entity_type::PARALLEL_OWNED>() const {
  return nodeids_owned_.size();
}
template<> inline
unsigned int
Mesh::num_nodes<Entity_type::PARALLEL_GHOST>() const {
  return nodeids_ghost_.size();
}
template<> inline
unsigned int Mesh::num_nodes<Entity_type::ALL>() const {
  return num_nodes<Entity_type::PARALLEL_OWNED>() +
      num_nodes<Entity_type::PARALLEL_GHOST>();
}

template<Entity_type type> inline
unsigned int Mesh::num_edges() const {
  std::cerr << "num_edges: Not defined for parallel type " << type << "\n";
  return 0;
}
template<> inline
unsigned int
Mesh::num_edges<Entity_type::PARALLEL_OWNED>() const {
  return edgeids_owned_.size();
}
template<> inline
unsigned int
Mesh::num_edges<Entity_type::PARALLEL_GHOST>() const {
  return edgeids_ghost_.size();
}
template<> inline
unsigned int Mesh::num_edges<Entity_type::ALL>() const {
  return num_edges<Entity_type::PARALLEL_OWNED>() +
      num_edges<Entity_type::PARALLEL_GHOST>();
}

template<Entity_type type> inline
unsigned int Mesh::num_faces() const {
  std::cerr << "num_faces: Not defined for parallel type " << type << "\n";
  return 0;
}
template<> inline
unsigned int
Mesh::num_faces<Entity_type::PARALLEL_OWNED>() const {
  return faceids_owned_.size();
}
template<> inline
unsigned int
Mesh::num_faces<Entity_type::PARALLEL_GHOST>() const {
  return faceids_ghost_.size();
}
template<> inline
unsigned int Mesh::num_faces<Entity_type::ALL>() const {
  return (num_faces<Entity_type::PARALLEL_OWNED>() +
          num_faces<Entity_type::PARALLEL_GHOST>());
}

template<Entity_type type> inline
unsigned int Mesh::num_sides() const {
  std::cerr << "num_sides: Not defined for parallel type " << type << "\n";
  return 0;
}
template<> inline
unsigned int
Mesh::num_sides<Entity_type::PARALLEL_OWNED>() const {
  return sideids_owned_.size();
}
template<> inline
unsigned int
Mesh::num_sides<Entity_type::PARALLEL_GHOST>() const {
  return sideids_ghost_.size();
}
template<> inline
unsigned int
Mesh::num_sides<Entity_type::BOUNDARY_GHOST>() const {
  return sideids_boundary_ghost_.size();
}
template<> inline
unsigned int Mesh::num_sides<Entity_type::ALL>() const {
  return (num_sides<Entity_type::PARALLEL_OWNED>() +
          num_sides<Entity_type::PARALLEL_GHOST>() +
          num_sides<Entity_type::BOUNDARY_GHOST>());
}

template<Entity_type type> inline
unsigned int Mesh::num_wedges() const {
  std::cerr << "num_wedges: Not defined for parallel type " << type << "\n";
  return 0;
}
template<> inline
unsigned int
Mesh::num_wedges<Entity_type::PARALLEL_OWNED>() const {
  return wedgeids_owned_.size();
}
template<> inline
unsigned int
Mesh::num_wedges<Entity_type::PARALLEL_GHOST>() const {
return wedgeids_ghost_.size();
}
template<> inline
unsigned int
Mesh::num_wedges<Entity_type::BOUNDARY_GHOST>() const {
return wedgeids_boundary_ghost_.size();
}
template<> inline
unsigned int Mesh::num_wedges<Entity_type::ALL>() const {
  return (num_wedges<Entity_type::PARALLEL_OWNED>() +
          num_wedges<Entity_type::PARALLEL_GHOST>() +
          num_wedges<Entity_type::BOUNDARY_GHOST>());
}

template<Entity_type type> inline
unsigned int Mesh::num_corners() const {
  std::cerr << "num_corners: Not defined for type " << type << "\n";
  return 0;
}
template<> inline
unsigned int
Mesh::num_corners<Entity_type::PARALLEL_OWNED>() const {
  return cornerids_owned_.size();
}
template<> inline
unsigned int
Mesh::num_corners<Entity_type::PARALLEL_GHOST>() const {
  return cornerids_ghost_.size();
}
template<> inline
unsigned int Mesh::num_corners<Entity_type::BOUNDARY_GHOST>() const {
  return cornerids_boundary_ghost_.size();
}
template<> inline
unsigned int Mesh::num_corners<Entity_type::ALL>() const {
  return (num_corners<Entity_type::PARALLEL_OWNED>() +
          num_corners<Entity_type::PARALLEL_GHOST>() +
          num_corners<Entity_type::BOUNDARY_GHOST>());
}

template<Entity_type type> inline
unsigned int Mesh::num_cells() const {
  std::cerr << "num_cells: Not defined for parallel type " << type << "\n";
  return 0;
}
template<> inline
unsigned int
Mesh::num_cells<Entity_type::PARALLEL_OWNED>() const {
  return cellids_owned_.size();
}
template<> inline
unsigned int
Mesh::num_cells<Entity_type::PARALLEL_GHOST>() const {
  return cellids_ghost_.size();
}
template<> inline
unsigned int
Mesh::num_cells<Entity_type::BOUNDARY_GHOST>() const {
  return cellids_boundary_ghost_.size();
}
template<> inline
unsigned int Mesh::num_cells<Entity_type::ALL>() const {
  return (num_cells<Entity_type::PARALLEL_OWNED>() +
          num_cells<Entity_type::PARALLEL_GHOST>() +
          num_cells<Entity_type::BOUNDARY_GHOST>());
}


inline
unsigned int Mesh::num_entities(const Entity_kind kind,
                                const Entity_type type) const {
  switch (kind) {
    case Entity_kind::NODE:
      switch (type) {
        case Entity_type::PARALLEL_OWNED:
          return num_nodes<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_nodes<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_nodes<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::EDGE:
      switch (type) {
        case Entity_type::PARALLEL_OWNED:
          return num_edges<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_edges<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_edges<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::FACE:
      switch (type) {
        case Entity_type::PARALLEL_OWNED:
          return num_faces<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_faces<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_faces<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::SIDE:
      switch (type) {
        case Entity_type::PARALLEL_OWNED:
          return num_sides<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_sides<Entity_type::PARALLEL_GHOST>();
        case Entity_type::BOUNDARY_GHOST:
          return num_sides<Entity_type::BOUNDARY_GHOST>();
        case Entity_type::ALL:
          return num_sides<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::WEDGE:
      switch (type) {
        case Entity_type::PARALLEL_OWNED:
          return num_wedges<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_wedges<Entity_type::PARALLEL_GHOST>();
        case Entity_type::BOUNDARY_GHOST:
          return num_wedges<Entity_type::BOUNDARY_GHOST>();
        case Entity_type::ALL:
          return num_wedges<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::CORNER:
      switch (type) {
        case Entity_type::PARALLEL_OWNED:
          return num_corners<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_corners<Entity_type::PARALLEL_GHOST>();
        case Entity_type::BOUNDARY_GHOST:
          return num_corners<Entity_type::BOUNDARY_GHOST>();
        case Entity_type::ALL:
          return num_corners<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::CELL:
      switch (type) {
        case Entity_type::PARALLEL_OWNED:
          return num_cells<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_cells<Entity_type::PARALLEL_GHOST>();
        case Entity_type::BOUNDARY_GHOST:
          return num_cells<Entity_type::BOUNDARY_GHOST>();
        case Entity_type::ALL:
          return num_cells<Entity_type::ALL>();
        default: return 0;
      }
    default:
      return 0;
  }
}


// templated version of functions returning entity lists (default
// implementation prints error message - meaningful values returned
// through template specialization)

template<Entity_type type> inline
const std::vector<Entity_ID> & Mesh::nodes() const {
  std::cerr << "Mesh::nodes() - " <<
      "Meaningless to query for list of nodes of parallel type " <<
      type << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::nodes<Entity_type::PARALLEL_OWNED>() const {
  return nodeids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::nodes<Entity_type::PARALLEL_GHOST>() const {
  return nodeids_ghost_;
}
template<> inline
const std::vector<Entity_ID> & Mesh::nodes<Entity_type::ALL>() const {
  return nodeids_all_;
}

template<Entity_type type> inline
const std::vector<Entity_ID> & Mesh::edges() const {
  std::cerr << "Mesh::edges() - " <<
      "Meaningless to query for list of edges of parallel type " <<
      type << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::edges<Entity_type::PARALLEL_OWNED>() const {
  return edgeids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::edges<Entity_type::PARALLEL_GHOST>() const {
  return edgeids_ghost_;
}
template<> inline
const std::vector<Entity_ID> & Mesh::edges<Entity_type::ALL>() const {
  return edgeids_all_;
}

template<Entity_type type> inline
const std::vector<Entity_ID> & Mesh::faces() const {
  std::cerr << "Mesh::faces() - " <<
      "Meaningless to query for list of faces of parallel type " <<
      type << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::faces<Entity_type::PARALLEL_OWNED>() const {
  return faceids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::faces<Entity_type::PARALLEL_GHOST>() const {
  return faceids_ghost_;
}
template<> inline
const std::vector<Entity_ID> & Mesh::faces<Entity_type::ALL>() const {
  return faceids_all_;
}

template<Entity_type type> inline
const std::vector<Entity_ID> & Mesh::sides() const {
  std::cerr << "Mesh::sides() - " <<
      "Meaningless to query for list of sides of parallel type " <<
      type << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::sides<Entity_type::PARALLEL_OWNED>() const {
  return sideids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::sides<Entity_type::PARALLEL_GHOST>() const {
  return sideids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::sides<Entity_type::BOUNDARY_GHOST>() const {
  return sideids_boundary_ghost_;
}
template<> inline
const std::vector<Entity_ID> & Mesh::sides<Entity_type::ALL>() const {
  return sideids_all_;
}

template<Entity_type type> inline
const std::vector<Entity_ID> & Mesh::wedges() const {
  std::cerr << "Mesh::wedges() - " <<
      "Meaningless to query for list of wedges of parallel type " <<
      type << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::wedges<Entity_type::PARALLEL_OWNED>() const {
  return wedgeids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::wedges<Entity_type::PARALLEL_GHOST>() const {
  return wedgeids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::wedges<Entity_type::BOUNDARY_GHOST>() const {
  return wedgeids_boundary_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::wedges<Entity_type::ALL>() const {
  return wedgeids_all_;
}

template<Entity_type type> inline
const std::vector<Entity_ID> & Mesh::corners() const {
  std::cerr << "Mesh::corners() - " <<
      "Meaningless to query for list of corners of parallel type " <<
      type << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::corners<Entity_type::PARALLEL_OWNED>() const {
  return cornerids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::corners<Entity_type::PARALLEL_GHOST>() const {
  return cornerids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::corners<Entity_type::BOUNDARY_GHOST>() const {
  return cornerids_boundary_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::corners<Entity_type::ALL>() const {
  return cornerids_all_;
}

template<Entity_type type> inline
const std::vector<Entity_ID> & Mesh::cells() const {
  std::cerr << "Mesh::cells() - " <<
      "Meaningless to query for list of cells of parallel type " <<
      type << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::cells<Entity_type::PARALLEL_OWNED>() const {
  return cellids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::cells<Entity_type::PARALLEL_GHOST>() const {
  return cellids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::cells<Entity_type::BOUNDARY_GHOST>() const {
  return cellids_boundary_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
Mesh::cells<Entity_type::ALL>() const {
  return cellids_all_;
}



// Inline functions of the Mesh class

inline
void Mesh::cell_get_faces(const Entity_ID cellid, Entity_ID_List *faceids,
                          const bool ordered) const {
  cell_get_faces_and_dirs(cellid, faceids, NULL, ordered);
}

inline
void Mesh::edge_get_nodes(const Entity_ID edgeid, Entity_ID *nodeid0,
                          Entity_ID *nodeid1) const {
#ifdef JALI_CACHE_VARS
  *nodeid0 = edge_node_ids[edgeid][0];
  *nodeid1 = edge_node_ids[edgeid][1];
#else
  edge_get_nodes_internal(edgeid, nodeid0, nodeid1);
#endif
}

inline
Entity_ID Mesh::side_get_wedge(const Entity_ID sideid, int iwedge) const {
  assert(wedges_requested);
  return (iwedge ? 2*sideid + 1 : 2*sideid);
}

inline
Entity_ID Mesh::side_get_face(const Entity_ID sideid) const {
  assert(sides_requested);
  assert(side_info_cached);
  return side_face_id[sideid];
}

inline
Entity_ID Mesh::side_get_edge(const Entity_ID sideid) const {
  assert(sides_requested);
  assert(side_info_cached);
  return side_edge_id[sideid];
}

inline
int Mesh::side_get_edge_use(const Entity_ID sideid) const {
  assert(sides_requested);
  assert(side_info_cached);
  return static_cast<int>(side_edge_use[sideid]);
}

inline
Entity_ID Mesh::side_get_cell(const Entity_ID sideid) const {
  assert(sides_requested);
  assert(side_info_cached);
  return side_cell_id[sideid];
}

inline
Entity_ID Mesh::side_get_node(const Entity_ID sideid, const int inode) const {
  assert(sides_requested);
  assert(side_info_cached && edge2node_info_cached);
  assert(inode == 0 || inode == 1);

  Entity_ID edgeid = side_edge_id[sideid];
  Entity_ID enodes[2];
  edge_get_nodes(edgeid, &enodes[0], &enodes[1]);
  bool use = side_edge_use[sideid];
  return (use ? enodes[inode] : enodes[!inode]);
}

inline
Entity_ID Mesh::side_get_opposite_side(const Entity_ID sideid) const {
  assert(sides_requested);
  assert(side_info_cached);
  return side_opp_side_id[sideid];
}

inline
Entity_ID Mesh::wedge_get_cell(const Entity_ID wedgeid) const {
  assert(sides_requested && wedges_requested);
  assert(side_info_cached);
  int sideid = wedgeid/2;  // which side does wedge belong to
  return side_get_cell(sideid);
}

inline
Entity_ID Mesh::wedge_get_face(const Entity_ID wedgeid) const {
  assert(sides_requested && wedges_requested);
  assert(side_info_cached);
  Entity_ID sideid = static_cast<Entity_ID>(wedgeid/2);
  return side_get_face(sideid);
}

inline
Entity_ID Mesh::wedge_get_edge(const Entity_ID wedgeid) const {
  assert(sides_requested && wedges_requested);
  assert(side_info_cached);
  Entity_ID sideid = static_cast<Entity_ID>(wedgeid/2);
  return side_get_edge(sideid);
}

inline
Entity_ID Mesh::wedge_get_node(const Entity_ID wedgeid) const {
  assert(sides_requested && wedges_requested);
  assert(side_info_cached);
  Entity_ID sideid = static_cast<Entity_ID>(wedgeid/2);
  int iwedge = wedgeid%2;  // Is it wedge 0 or wedge 1 of side
  return side_get_node(sideid, iwedge);
}

inline
Entity_ID Mesh::wedge_get_corner(const Entity_ID wedgeid) const {
  assert(sides_requested && wedges_requested);
  assert(side_info_cached && wedge_info_cached);
  return wedge_corner_id[wedgeid];
}

inline
Entity_ID Mesh::wedge_get_adjacent_wedge(Entity_ID const wedgeid) const {
  assert(wedges_requested);
  
  // Wedges come in pairs; their IDs are (2*sideid) and (2*sideid+1)
  // If the wedge ID is an odd number, then the adjacent wedge ID is
  // wedge ID minus one; If it is an even number, the adjacent wedge
  // ID is wedge ID plus one

  return (wedgeid%2 ? wedgeid - 1 : wedgeid + 1);
}

inline
Entity_ID Mesh::wedge_get_opposite_wedge(Entity_ID const wedgeid) const {
  assert(wedges_requested);
  assert(side_info_cached);

  Entity_ID sideid = static_cast<Entity_ID>(wedgeid/2);
  int iwedge = wedgeid%2;  // Is it wedge 0 or wedge 1 of side

  Entity_ID oppsideid = side_opp_side_id[sideid];
  if (oppsideid == -1)
    return -1;
  else {
    // if wedge is wedge 0 of side, the opposite wedge will be wedge 1
    // of the opposite side
    Entity_ID adjwedgeid = iwedge ? 2*oppsideid : 2*oppsideid + 1;
    return adjwedgeid;
  }
}

inline
void Mesh::corner_get_wedges(const Entity_ID cornerid,
                             Entity_ID_List *cwedges) const {
  assert(corners_requested);
  assert(corner_info_cached);

  int nwedges = corner_wedge_ids[cornerid].size();
  (*cwedges).resize(nwedges);
  std::copy(corner_wedge_ids[cornerid].begin(),
            corner_wedge_ids[cornerid].end(),
            cwedges->begin());
}

inline
Entity_ID Mesh::corner_get_node(const Entity_ID cornerid) const {
  assert(corners_requested);
  assert(corner_info_cached && side_info_cached);
  assert(corner_wedge_ids[cornerid].size());

  // Instead of calling corner_get_wedges which involves a list copy,
  // we will directly access the first element of the corner_wedge_ids
  // array
  Entity_ID w0 = corner_wedge_ids[cornerid][0];
  return wedge_get_node(w0);
}

inline
Entity_ID Mesh::corner_get_cell(const Entity_ID cornerid) const {
  assert(corners_requested);
  assert(corner_info_cached && side_info_cached);
  assert(corner_wedge_ids[cornerid].size());

  // Instead of calling corner_get_wedges which involves a list copy,
  // we will directly access the first element of the corner_wedge_ids
  // array
  Entity_ID w0 = corner_wedge_ids[cornerid][0];
  return wedge_get_cell(w0);
}

// Inefficient fallback implementation - hopefully the derived class
// has a more direct implementation

inline
void Mesh::node_get_coordinates(const Entity_ID nodeid,
                                std::array<double, 3> *ncoord) const {
  assert(spacedim == 3);
  JaliGeometry::Point p;
  node_get_coordinates(nodeid, &p);
  (*ncoord)[0] = p[0];
  (*ncoord)[1] = p[1];
  (*ncoord)[2] = p[2];
}

// Inefficient fallback implementation - hopefully the derived class
// has a more direct implementation

inline
void Mesh::node_get_coordinates(const Entity_ID nodeid,
                                std::array<double, 2> *ncoord) const {
  assert(spacedim == 2);
  JaliGeometry::Point p;
  node_get_coordinates(nodeid, &p);
  (*ncoord)[0] = p[0];
  (*ncoord)[1] = p[1];
}
  
// Inefficient fallback implementation - hopefully the derived class
// has a more direct implementation

inline
void Mesh::node_get_coordinates(const Entity_ID nodeid, double *ncoord) const {
  assert(spacedim == 1);
  JaliGeometry::Point p;
  node_get_coordinates(nodeid, &p);
  *ncoord = p[0];
}


}  // end namespace Jali





#endif /* _JALI_MESH_H_ */
