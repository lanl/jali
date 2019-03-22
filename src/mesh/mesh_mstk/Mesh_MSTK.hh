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
#ifndef _JALI_MESH_MSTK_H_
#define _JALI_MESH_MSTK_H_

#include "MSTK.h"

#include <cstdint>
#include <memory>
#include <vector>
#include <sstream>
#include <typeinfo>
#include <algorithm>

#include "Mesh.hh"
#include "Point.hh"
#include "GeometricModel.hh"
#include "LabeledSetRegion.hh"
#include "PointRegion.hh"
#include "LogicalRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Jali {

// Mesh class based on the MSTK framework.
//
// Instantiating a const version of this class only guarantees that
// the underlying mesh topology and geometry does not change. For
// purposes of memory savings we use lazy initialization of face and
// edge lists (they are already present in the underlying MSTK mesh
// data structures), which means that the data structures holding the
// mesh information may still change.


class Mesh_MSTK : public Mesh {
 public:

  // Constructors that read the mesh from a file

  // The request_faces and request_edges arguments have to be at the
  // end and not in the middle because if we omit them and specify a
  // pointer argument like gm or verbosity_obj, then there is implicit
  // conversion of the pointer to bool, thereby defeating the intent
  // of the call and making the pointer argument seem NULL. In C++11,
  // we could "delete" the illegal version of the call effectively
  // blocking the implicit conversion.


  Mesh_MSTK(const std::string filename, const MPI_Comm& incomm,
            const JaliGeometry::GeometricModelPtr& gm =
            (JaliGeometry::GeometricModelPtr) NULL,
            const bool request_faces = true,
            const bool request_edges = false,
            const bool request_sides = false,
            const bool request_wedges = false,
            const bool request_corners = false,
            const int num_tiles = 0,
            const int num_ghost_layers_tile = 0,
            const int num_ghost_layers_distmesh = 1,
            const bool request_boundary_ghosts = false,
            const Partitioner_type partitioner = Partitioner_type::METIS,
            const bool contiguous_gids = false,
            const JaliGeometry::Geom_type geom_type =
            JaliGeometry::Geom_type::CARTESIAN);
  
  // Constructors that generate a mesh internally (regular hexahedral mesh only)

  // 3D
  Mesh_MSTK(const double x0, const double y0, const double z0,
            const double x1, const double y1, const double z1,
            const unsigned int nx, const unsigned int ny,
            const unsigned int nz,
            const MPI_Comm& incomm,
            const JaliGeometry::GeometricModelPtr& gm =
            (JaliGeometry::GeometricModelPtr) NULL,
            const bool request_faces = true,
            const bool request_edges = false,
            const bool request_sides = false,
            const bool request_wedges = false,
            const bool request_corners = false,
            const int num_tiles = 0,
            const int num_ghost_layers_tile = 0,
            const int num_ghost_layers_distmesh = 1,
            const bool request_boundary_ghosts = false,
            const Partitioner_type partitioner = Partitioner_type::METIS,
            const bool contiguous_gids = false);


  // 2D
  Mesh_MSTK(const double x0, const double y0,
            const double x1, const double y1,
            const int nx, const int ny,
            const MPI_Comm& comm,
            const JaliGeometry::GeometricModelPtr& gm =
            (JaliGeometry::GeometricModelPtr) NULL,
            const bool request_faces = true,
            const bool request_edges = false,
            const bool request_sides = false,
            const bool request_wedges = false,
            const bool request_corners = false,
            const int num_tiles = 0,
            const int num_ghost_layers_tile = 0,
            const int num_ghost_layers_distmesh = 1,
            const bool request_boundary_ghosts = false,
            const Partitioner_type partitioner = Partitioner_type::METIS,
            const bool contiguous_gids = false,
            const JaliGeometry::Geom_type geom_type =
            JaliGeometry::Geom_type::CARTESIAN);

  // Construct a mesh by extracting a subset of entities from another
  // mesh. The subset may be specified by a setname or a list of
  // entities. In some cases like extracting a surface mesh from a
  // volume mesh, constructor can be asked to "flatten" the mesh to a
  // lower dimensional space or to extrude the mesh to give higher
  // dimensional cells

  Mesh_MSTK(const std::shared_ptr<Mesh> inmesh,
            const std::vector<std::string>& setnames,
            const Entity_kind entity_kind,
            const bool flatten = false,
            const bool extrude = false,
            const bool request_faces = true,
            const bool request_edges = false,
            const bool request_sides = false,
            const bool request_wedges = false,
            const bool request_corners = false,
            const int num_tiles = 0,
            const int num_ghost_layers_tile = 0,
            const int num_ghost_layers_distmesh = 1,
            const bool request_boundary_ghosts = false,
            const Partitioner_type partitioner = Partitioner_type::METIS,
            const bool contiguous_gids = false,
            const JaliGeometry::Geom_type geom_type =
            JaliGeometry::Geom_type::CARTESIAN);

  Mesh_MSTK(const Mesh& inmesh,
            const std::vector<std::string>& setnames,
            const Entity_kind entity_kind,
            const bool flatten = false,
            const bool extrude = false,
            const bool request_faces = true,
            const bool request_edges = false,
            const bool request_sides = false,
            const bool request_wedges = false,
            const bool request_corners = false,
            const int num_tiles = 0,
            const int num_ghost_layers_tile = 0,
            const int num_ghost_layers_distmesh = 1,
            const bool request_boundary_ghosts = false,
            const Partitioner_type partitioner = Partitioner_type::METIS,
            const bool contiguous_gids = false,
            const JaliGeometry::Geom_type geom_type =
            JaliGeometry::Geom_type::CARTESIAN);

  Mesh_MSTK(const Mesh& inmesh,
            const std::vector<int>& entity_list,
            const Entity_kind entity_kind,
            const bool flatten = false,
            const bool extrude = false,
            const bool request_faces = true,
            const bool request_edges = false,
            const bool request_sides = false,
            const bool request_wedges = false,
            const bool request_corners = false,
            const int num_tiles = 0,
            const int num_ghost_layers_tile = 0,
            const int num_ghost_layers_distmesh = 1,
            const bool request_boundary_ghosts = false,
            const Partitioner_type partitioner = Partitioner_type::METIS,
            const bool contiguous_gids = false,
            const JaliGeometry::Geom_type geom_type =
            JaliGeometry::Geom_type::CARTESIAN);


  ~Mesh_MSTK();


  // Get cell type

  Cell_type cell_get_type(const Entity_ID cellid) const;


  // Parent entity in the source mesh if mesh was derived from another mesh

  Entity_ID entity_get_parent(const Entity_kind kind,
                              const Entity_ID entid) const;




  //
  // General mesh information
  // -------------------------
  //

  // Global ID of any entity

  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const;



  //
  // Mesh Entity Adjacencies
  //-------------------------

  // Get nodes of cell
  // On a distributed mesh, all nodes (OWNED or GHOST) of the cell
  // are returned
  // Nodes are returned in a standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD Entity_kind::CELL TYPES in 3D
  // For a general polyhedron this will return the nodes in
  // arbitrary order
  // In 2D, the nodes of the polygon will be returned in ccw order
  // consistent with the face normal

  void cell_get_nodes(const Entity_ID cellid,
                      Entity_ID_List *nodeids) const;


  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2

  void face_get_nodes(const Entity_ID faceid,
                      Entity_ID_List *nodeids) const;


  // Get nodes of edge On a distributed mesh all nodes (Entity_type::PARALLEL_OWNED or
  // Entity_type::PARALLEL_GHOST) of the face are returned

  void edge_get_nodes_internal(const Entity_ID edgeid, Entity_ID *point0,
                               Entity_ID *point1) const;

  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node

  void node_get_cells(const Entity_ID nodeid,
                      const Entity_type ptype,
                      Entity_ID_List *cellids) const;

  // Faces of type 'ptype' connected to a node

  void node_get_faces(const Entity_ID nodeid,
                      const Entity_type ptype,
                      Entity_ID_List *faceids) const;

  // Get faces of ptype of a particular cell that are connected to the
  // given node

  void node_get_cell_faces(const Entity_ID nodeid,
                           const Entity_ID cellid,
                           const Entity_type ptype,
                           Entity_ID_List *faceids) const;


  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = USED, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces

  void cell_get_face_adj_cells(const Entity_ID cellid,
                               const Entity_type ptype,
                               Entity_ID_List *fadj_cellids) const;

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order

  void cell_get_node_adj_cells(const Entity_ID cellid,
                               const Entity_type ptype,
                               Entity_ID_List *nadj_cellids) const;


  //
  // Mesh entity geometry
  //--------------
  //

  // Node coordinates - 3 in 3D and 2 in 2D

  void node_get_coordinates(const Entity_ID nodeid,
                            JaliGeometry::Point *ncoord) const;
  void node_get_coordinates(const Entity_ID nodeid,
                            std::array<double, 3> *ncoord) const;
  void node_get_coordinates(const Entity_ID nodeid,
                            std::array<double, 2> *ncoord) const;



  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size divided by number of spatial dimensions

  void face_get_coordinates(const Entity_ID faceid,
                            std::vector<JaliGeometry::Point> *fcoords) const;

  // Coordinates of cells in standard order (Exodus  II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD Entity_kind::CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions

  void cell_get_coordinates(const Entity_ID cellid,
                            std::vector<JaliGeometry::Point> *ccoords) const;

  // Modify the coordinates of a node

  void node_set_coordinates(const Entity_ID nodeid,
                            const JaliGeometry::Point coords);

  void node_set_coordinates(const Entity_ID nodeid, const double *coords);



  //
  // Boundary Conditions or Sets
  //----------------------------
  //


  // Get number of entities of type 'category' in set

  unsigned int get_set_size(const Set_Name setname,
                            const Entity_kind kind,
                            const Entity_type ptype) const;


  unsigned int get_set_size(const char *setname,
                            const Entity_kind kind,
                            const Entity_type ptype) const;

  // Get list of entities of type 'category' in set

  void get_set_entities(const Set_Name setname,
                        const Entity_kind kind,
                        const Entity_type ptype,
                        std::vector<Entity_ID> *entids) const;


  void get_set_entities(const char *setname,
                        const Entity_kind kind,
                        const Entity_type ptype,
                        std::vector<Entity_ID> *entids) const;


  //! \brief Export to Exodus II file
  //! Export mesh to Exodus II file. If with_fields is true, the fields in
  //! JaliState are also exported out.

  void write_to_exodus_file(const std::string exodusfilename,
                            bool with_fields = true) const {
    if (with_fields)
      MESH_ExportToFile(mesh, exodusfilename.c_str(), "exodusii", 0, NULL,
                        NULL, mpicomm);
    else
      MESH_ExportToFile(mesh, exodusfilename.c_str(), "exodusii", -1, NULL,
                        NULL, mpicomm);
  }


  //! \brief Export to GMV file
  //! Export mesh to GMV file. If with_fields is true, the fields in
  //! JaliState are also exported out.

  void write_to_gmv_file(const std::string gmvfilename,
                         bool with_fields = true) const {
    if (with_fields)
      MESH_ExportToFile(mesh, gmvfilename.c_str(), "gmv", 0, NULL, NULL,
                        mpicomm);
    else
      MESH_ExportToFile(mesh, gmvfilename.c_str(), "gmv", -1, NULL, NULL,
                        mpicomm);
  }

  // Run MSTK's internal checks - meant for debugging only
  // Returns true if everything is ok, false otherwise

  bool run_internal_mstk_checks() const;

 protected:

  //! Get the number of fields on entities of a particular type along
  //! with their names and variable types - DESIGNED TO BE CALLED ONLY
  //! BY THE JALI STATE MANAGER FOR INITIALIZATION OF MESH STATE FROM
  //! THE MESH FILE

  void get_field_info(Entity_kind on_what, int *num,
                      std::vector<std::string> *varnames,
                      std::vector<std::string> *vartypes) const;

  //! Retrieve an integer field on the mesh. If the return value is
  //! false, it could be that (1) the field does not exist (2) it
  //! exists but is associated with a different type of entity (3) the
  //! variable type sent in was the wrong type (int instead of double
  //! or double instead of std::array<double,2> or
  //! std::array<double,2> instead of std::array<double,3> etc -
  //! DESIGNED TO BE CALLED ONLY BY THE JALI STATE MANAGER FOR
  //! INITIALIZATION OF MESH STATE FROM THE MESH FILE

  bool get_field(std::string field_name, Entity_kind on_what, int *data) const;

  //! Retrieve an integer field on the mesh. If the return value is
  //! false, it could be that (1) the field does not exist (2) it
  //! exists but is associated with a different type of entity (3) the
  //! variable type sent in was the wrong type (int instead of double
  //! or double instead of std::array<double,2> or
  //! std::array<double,2> instead of std::array<double,3> etc -
  //! DESIGNED TO BE CALLED ONLY BY THE JALI STATE MANAGER FOR
  //! INITIALIZATION OF MESH STATE FROM THE MESH FILE

  bool get_field(std::string field_name, Entity_kind on_what,
                 double *data) const;

  //! Retrieve an array field on the mesh. If the return value is
  //! false, it could be that (1) the field does not exist (2) it
  //! exists but is associated with a different type of entity (3) the
  //! variable type sent in was the wrong type (int instead of double
  //! or double instead of std::array<double,2> or
  //! std::array<double,2> instead of std::array<double,3> etc -
  //! DESIGNED TO BE CALLED ONLY BY THE JALI STATE MANAGER FOR
  //! INITIALIZATION OF MESH STATE FROM THE MESH FILE

  //! HAVE TO DO THIS WONKY THING BECAUSE COMPILER COMPLAINS THAT IT
  //! CANNOT FIND A get_field IMPLEMENTATION WITH std::array<double,2ul>

  template<std::size_t N>
  bool get_field_internal(std::string field_name, Entity_kind on_what,
                          std::array<double, N> *data) const;
  bool get_field(std::string field_name, Entity_kind on_what,
                 std::array<double, 2> *data) const {
    return get_field_internal(field_name, on_what,data);
  }
  bool get_field(std::string field_name, Entity_kind on_what,
                 std::array<double, 3> *data) const {
    return get_field_internal(field_name, on_what,data);
  }
  bool get_field(std::string field_name, Entity_kind on_what,
                 std::array<double, 6> *data) const {
    return get_field_internal(field_name, on_what,data);
  }

  //! Store an integer field on the mesh. If the return value is false, it
  //! means that the mesh already has a field of that name but its of
  //! a different type or its on a different type of entity - DESIGNED
  //! TO BE CALLED ONLY BY THE JALI STATE MANAGER FOR INITIALIZATION
  //! OF MESH STATE FROM THE MESH FILE

  bool store_field(std::string field_name, Entity_kind on_what, int *data);

  //! Store an real field on the mesh. If the return value is false, it
  //! means that the mesh already has a field of that name but its of
  //! a different type or its on a different type of entity - DESIGNED
  //! TO BE CALLED ONLY BY THE JALI STATE MANAGER FOR INITIALIZATION
  //! OF MESH STATE FROM THE MESH FILE

  bool store_field(std::string field_name, Entity_kind on_what, double *data);

  // Store a vector or tensor field data with the mesh

  //! HAVE TO DO THIS WONKY THING BECAUSE COMPILER COMPLAINS THAT IT
  //! CANNOT FIND A get_field IMPLEMENTATION WITH std::array<double,2ul>

  template<std::size_t N>
  bool store_field_internal(std::string field_name, Entity_kind on_what,
                            std::array<double, N> *data);
  bool store_field(std::string field_name, Entity_kind on_what,
                   std::array<double, 2> *data) {
    return store_field_internal(field_name, on_what,data);
  }
  bool store_field(std::string field_name, Entity_kind on_what,
                   std::array<double, 3> *data) {
    return store_field_internal(field_name, on_what,data);
  }
  bool store_field(std::string field_name, Entity_kind on_what,
                   std::array<double, 6> *data) {
    return store_field_internal(field_name, on_what, data);
  }

  void get_labeled_set_entities(const JaliGeometry::LabeledSetRegionPtr rgn,
                                const Entity_kind kind,
                                Entity_ID_List *owned_entities,
                                Entity_ID_List *ghost_entities) const;

 private:

  // Private methods
  // ----------------------------

  void clear_internals_();

  void pre_create_steps_(const int space_dimension, const MPI_Comm& incomm,
                         const JaliGeometry::GeometricModelPtr& gm);
  void post_create_steps_();

  void init_mesh_from_file_(const std::string filename,
                            const Partitioner_type partitioner =
                            PARTITIONER_DEFAULT);

  void collapse_degen_edges();
  Cell_type MFace_Celltype(MFace_ptr f);
  Cell_type MRegion_Celltype(MRegion_ptr r);
  void label_celltype();

  void init_pvert_lists();
  void init_pedge_lists();
  void init_pedge_dirs();
  void init_pface_lists();
  void init_pface_dirs();
  void init_pface_dirs_3();
  void init_pface_dirs_2();
  void init_pcell_lists();

  void init_vertex_id2handle_maps();
  void init_edge_id2handle_maps();
  void init_face_id2handle_maps();
  void init_cell_id2handle_maps();
  void init_global_ids();

  void init_cell_map();
  void init_face_map();
  void init_edge_map();
  void init_node_map();

  void init_nodes();
  void init_edges();
  void init_faces();
  void init_cells();

  void create_boundary_ghosts();

  void init_set_info();
  void inherit_labeled_sets(MAttrib_ptr copyatt, List_ptr src_entities);

  // internal name of sets (particularly labeled sets)
  std::string
  internal_name_of_set(const JaliGeometry::RegionPtr r,
                       const Entity_kind entity_kind) const;

  // internal name of sets (particularly labeled sets of cells)
  std::string
  other_internal_name_of_set(const JaliGeometry::RegionPtr r,
                             const Entity_kind entity_kind) const;

  int  generate_regular_mesh(Mesh_ptr mesh, double x0, double y0, double z0,
                             double x1, double y1, double z1, int nx,
                             int ny, int nz);
  int  generate_regular_mesh(Mesh_ptr mesh, double x0, double y0,
                             double x1, double y1, int nx, int ny);

  void extract_mstk_mesh(const Mesh_MSTK& inmesh,
                         const List_ptr entity_ids,
                         const MType entity_dim,
                         const bool flatten = false,
                         const bool extrude = false,
                         const bool request_faces = true,
                         const bool request_edges = false,
                         const int num_ghost_layers_distmesh = 1,
                         const bool request_boundary_ghosts = false,
                         const Partitioner_type partitioner =
                         Partitioner_type::METIS);


  // Downward Adjacencies
  //---------------------

  // Get faces of a cell and directions in which the cell uses the face

  // The Jali coding guidelines regarding function arguments is purposely
  // violated here to allow for a default input argument

  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order

  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise

  void cell_get_faces_and_dirs_internal(const Entity_ID cellid,
                                        Entity_ID_List *faceids,
                                        std::vector<dir_t> *face_dirs,
                                        const bool ordered = false) const;

  void cell_get_faces_and_dirs_ordered(const Entity_ID cellid,
                                       Entity_ID_List *faceids,
                                       std::vector<dir_t> *face_dirs) const;

  void cell_get_faces_and_dirs_unordered(const Entity_ID cellid,
                                         Entity_ID_List *faceids,
                                         std::vector<dir_t> *face_dirs) const;


  // Cells connected to a face

  void face_get_cells_internal(const Entity_ID faceid,
                               const Entity_type ptype,
                               Entity_ID_List *cellids) const;


  // Get edges of a cell

  void cell_get_edges_internal(const Entity_ID cellid,
                               Entity_ID_List *edgeids) const;

  // Get edges and directions of a 2D cell

  void cell_2D_get_edges_and_dirs_internal(const Entity_ID cellid,
                                           Entity_ID_List *edgeids,
                                           std::vector<dir_t> *edgedirs) const;

  // Edges and edge directions of a face

  void face_get_edges_and_dirs_internal(const Entity_ID cellid,
                                        Entity_ID_List *edgeids,
                                        std::vector<dir_t> *edgedirs,
                                        bool ordered = true) const;

  // Map from Jali's mesh entity kind to MSTK's mesh type.

  MType entity_kind_to_mtype(const Entity_kind kind) const {

    // The first index is cell dimension (0,1,2,3) and the second index
    // is the entity kind
    //
    // map order in each row is NODE, EDGE, FACE, CELL
    //
    // So, for a 1D mesh, nodes are MVERTEX type in MSTK, edges and faces
    // are also MVERTEX type, and cells are MEDGE type
    //
    // For a 2D mesh, nodes are MVERTEX type, edges and faces are MEDGE
    // type, and cells are MFACE type

    static MType const
      kind2mtype[4][4] = {{MVERTEX, MVERTEX, MVERTEX, MVERTEX},   // 0d meshes
                          {MVERTEX, MVERTEX, MVERTEX, MEDGE},     // 1d meshes
                          {MVERTEX, MEDGE,   MEDGE,   MFACE},     // 2d meshes
                          {MVERTEX, MEDGE,   MFACE,   MREGION}};  // 3d meshes

    return kind2mtype[manifold_dimension()][(int)kind];
  }


  // Private data

  MPI_Comm mpicomm;
  int myprocid, numprocs;

  Mesh_ptr mesh;

  int serial_run;


  // Local handles to entity lists (Vertices, "Faces", "Cells")

  // For a surface mesh, "Faces" refers to mesh edges and "Cells"
  // refers to mesh faces
  //
  // For a solid mesh, "Faces" refers to mesh faces and "Cells"
  // refers to mesh regions


  // These are MSTK's definitions of types of parallel mesh entities
  // These definitions are slightly different from what Jali has defined
  //
  // There are 2 types of entities relevant to this code - Owned and Ghost
  //
  // 1. Entity_type::PARALLEL_OWNED - owned by this processor
  //
  // 2. Entity_type::PARALLEL_GHOST - not owned by this processor
  //
  // Parallell_type::ALL = Entity_type::PARALLEL_OWNED + Entity_type::PARALLEL_GHOST 

  MSet_ptr OwnedVerts, NotOwnedVerts;

  mutable MSet_ptr OwnedEdges, NotOwnedEdges;

  mutable MSet_ptr OwnedFaces, NotOwnedFaces;

  MSet_ptr OwnedCells, GhostCells, BoundaryGhostCells;

  // Flags to indicate if face and edge info is initialized

  mutable bool faces_initialized, edges_initialized;

  // Marker to indicate if a cell is a boundary ghost cell

  MAttrib_ptr boundary_ghost_att = NULL;

  // Deleted entity lists if some pre-processing had to be done
  // to the mesh to eliminate degenerate entities

  bool entities_deleted;
  List_ptr deleted_vertices, deleted_edges, deleted_faces, deleted_regions;

  // Local ID to MSTK handle map

  std::vector<MEntity_ptr> vtx_id_to_handle;
  mutable std::vector<MEntity_ptr> edge_id_to_handle;
  mutable std::vector<MEntity_ptr> face_id_to_handle;
  std::vector<MEntity_ptr> cell_id_to_handle;


  // flag whether to flip a face dir or not when returning nodes of a
  // face (relevant only on partition boundaries)

  mutable bool *faceflip;

  // flag whether to flip an edge dir or not when returning nodes of an edge
  // (relevant only on partition boundaries)

  mutable bool *edgeflip;

  // Attribute to precompute and store celltype

  MAttrib_ptr celltype_att;

  // Parent entity attribute - populated if the mesh is derived from
  // another mesh

  MAttrib_ptr rparentatt, fparentatt, eparentatt, vparentatt;

  const Mesh_MSTK *parent_mesh;

  // variables needed for mesh deformation

  double *meshxyz;
  double *target_cell_volumes, *min_cell_volumes, *target_weights;

  // whether to make GIDs continguous or not
  bool contiguous_gids_;
};


// Retrieve field data from the mesh - special implementation for
// std::array<double,N>. Other implementations are in .cc file
template<std::size_t N>
inline
bool Mesh_MSTK::get_field_internal(std::string field_name, Entity_kind on_what,
                                   std::array<double, N> *data) const {
  MAttrib_ptr mattrib = MESH_AttribByName(mesh, field_name.c_str());
  if (!mattrib) return false;

  if (entity_kind_to_mtype(on_what) != (int) MAttrib_Get_EntDim(mattrib))
    return false;

  MAttType atttype = MAttrib_Get_Type(mattrib);
  if (atttype != VECTOR && atttype != TENSOR)
    return false;

  MType enttype = MAttrib_Get_EntDim(mattrib);
  int ncomp = MAttrib_Get_NumComps(mattrib);

  if (ncomp != N) return false;

  int nent;
  switch (enttype) {
    case MVERTEX: nent = MESH_Num_Vertices(mesh); break;
    case MEDGE: nent = MESH_Num_Edges(mesh); break;
    case MFACE: nent = MESH_Num_Faces(mesh); break;
    case MREGION: nent = MESH_Num_Regions(mesh); break;
    default: nent = 0; break;
  }

  MEntity_ptr ment;
  for (int i = 0; i < nent; i++) {
    switch (enttype) {
      case MVERTEX: ment = MESH_Vertex(mesh, i); break;
      case MEDGE: ment = MESH_Edge(mesh, i); break;
      case MFACE: ment = MESH_Face(mesh, i); break;
      case MREGION: ment = MESH_Region(mesh, i); break;
      default: ment = NULL; break;
    }
    if (ment == NULL) continue;  // not needed since nent=0 but here anyway

    int ival;
    double rval;
    void *pval;
    MEnt_Get_AttVal(ment, mattrib, &ival, &rval, &pval);

    for (int j = 0; j < N; j++)
      data[i][j] = ((double *)pval)[j];

  }  // for (int i...);

  return true;
}  // Mesh_MSTK::get_mesh_field




// Store field data with the mesh - special implementation for
// std::array<double,N>. Other implementations are in .cc file

template<std::size_t N>
inline
bool Mesh_MSTK::store_field_internal(std::string field_name,
                                     Entity_kind on_what,
                                     std::array<double,N> *data) {
  MType mtype = entity_kind_to_mtype(on_what);
  MAttType atttype;

  MAttrib_ptr mattrib = MESH_AttribByName(mesh, field_name.c_str());
  if (mattrib) {
    if (mtype != (int) MAttrib_Get_EntDim(mattrib))
      return false;

    atttype = MAttrib_Get_Type(mattrib);
    if (atttype != VECTOR || atttype != TENSOR) {
      std::cerr << "Mesh_MSTK::store_field -" <<
          " found attribute with same name but different type" << std::endl;
      return false;
    }
    int ncomp = MAttrib_Get_NumComps(mattrib);
    if (ncomp != N) {
      std::cerr << "Mesh_MSTK::store_field - found vector/tensor attribute " <<
          "but with different number of components" << std::endl;
      return false;
    }
  }
  else {
    if (N != 2 && N != 3 && N != 6) return false;
    atttype =  (N == 2 || (N == 3 && Mesh::space_dimension() == 3)) ? VECTOR :
        TENSOR;
    mattrib = MAttrib_New(mesh, field_name.c_str(), atttype, mtype, N);
  }

  int nent;
  switch (mtype) {
    case MVERTEX: nent = MESH_Num_Vertices(mesh); break;
    case MEDGE: nent = MESH_Num_Edges(mesh); break;
    case MFACE: nent = MESH_Num_Faces(mesh); break;
    case MREGION: nent = MESH_Num_Regions(mesh); break;
    default: nent = 0; break;
  }

  MEntity_ptr ment;
  for (int i = 0; i < nent; i++) {
    switch (mtype) {
      case MVERTEX: ment = MESH_Vertex(mesh, i); break;
      case MEDGE: ment = MESH_Edge(mesh, i); break;
      case MFACE: ment = MESH_Face(mesh, i); break;
      case MREGION: ment = MESH_Region(mesh, i); break;
      default: ment = NULL; break;
    }
    if (ment == NULL) continue;  // not needed since nent=0 but here anyway

    void *pval;
    double *vec = new double[N];
    std::copy(&(data[i][0]), &(data[i][0]) + N, vec);
    pval = vec;

    MEnt_Set_AttVal(ment, mattrib, 0, 0.0, pval);
  } // for

  return true;
}  // Mesh_MSTK::store_mesh_field

}  // End namespace Jali

#endif /* _JALI_MESH_MSTK_H_ */
