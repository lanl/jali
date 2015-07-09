#ifndef _JALIMESH_H_
#define _JALIMESH_H_

#include <memory>
#include "mpi.h"

#include "MeshDefs.hh"
#include "Point.hh"
#include "GeometricModel.hh"
#include "Region.hh"


#include <map>

namespace Jali
{

  // Base mesh class 
  //
  // Use the associated mesh factory to create an instance of a
  // derived class based on a particular mesh framework (like MSTK,
  // STKmesh etc.)
  //
  // **** IMPORTANT NOTE ABOUT CONSTANTNESS OF THIS CLASS ****
  // Instantiating a const version of this class only guarantees that
  // the underlying mesh topology and geometry does not change (the
  // public interfaces conforms strictly to this definition). However,
  // for purposes of memory savings we use lazy initialization and
  // caching of face data, edge data, geometry quantities, columns
  // etc., which means that these data may still change. We also
  // cannot initialize the cached quantities in the constructor since
  // they depend on initialization of data structures in the derived
  // class - however, the base class gets constructed before the
  // derived class gets constructed so it is not possible without more
  // obscure acrobatics. This is why some of the caching data
  // declarations are declared with the keyword 'mutable' and routines
  // that modify the mutable data are declared with a constant
  // qualifier.
  //


class Mesh
{

 public:

  // constructor

  Mesh(const bool request_faces=true,
       const bool request_edges=false,
       const MPI_Comm incomm=MPI_COMM_WORLD)
    : spacedim(3), celldim(3), mesh_type_(GENERAL), 
      cell_geometry_precomputed(false), face_geometry_precomputed(false),
      edge_geometry_precomputed(false), 
      faces_requested(request_faces), edges_requested(request_edges),
      cell2face_info_cached(false), face2cell_info_cached(false), 
      cell2edge_info_cached(false), face2edge_info_cached(false), 
      geometric_model_(NULL), comm(incomm)
  {
  }

  // destructor - must be virtual to downcast base class to derived class
  // (I don't understand why but the stackoverflow prophets say so)

  virtual ~Mesh() {}

  inline
  MPI_Comm get_comm() const {
    return comm;
  }

  inline
  void set_space_dimension(const unsigned int dim) {
    spacedim = dim;
  }

  inline
  unsigned int space_dimension() const
  {
    return spacedim;
  }

  inline
  void set_cell_dimension(const unsigned int dim) {
    celldim = dim;   // 3 is solid mesh, 2 is surface mesh, 1 is wire mesh
  }

  inline
  unsigned int cell_dimension() const
  {
    return celldim;
  }

  inline
  void set_geometric_model(const JaliGeometry::GeometricModelPtr &gm) {
    geometric_model_ = gm;
  }

  inline
  JaliGeometry::GeometricModelPtr geometric_model() const
  {
    return geometric_model_;
  }


  // Set/Get mesh type - RECTANGULAR, GENERAL (See MeshDefs.hh)

  inline
  void set_mesh_type(const Mesh_type mesh_type) {
    mesh_type_ = mesh_type;
  }

  inline
  Mesh_type mesh_type() const {
    return mesh_type_;
  }


  // Get parallel type of entity - OWNED, GHOST, ALL (See MeshDefs.hh)

  virtual
  Parallel_type entity_get_ptype(const Entity_kind kind,
                                 const Entity_ID entid) const = 0;


  // Parent entity in the source mesh if mesh was derived from another mesh

  virtual
  Entity_ID entity_get_parent(const Entity_kind kind, const Entity_ID entid) const;


  // Get cell type - UNKNOWN, TRI, QUAD, POLYGON, TET, PRISM, PYRAMID, HEX, POLYHED 
  // See MeshDefs.hh

  virtual
  Cell_type cell_get_type(const Entity_ID cellid) const = 0;

  // Cell type name

  std::string cell_type_to_name(const Cell_type type);

  //
  // General mesh information
  // -------------------------
  //

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)

  virtual
  unsigned int num_entities (const Entity_kind kind,
                             const Parallel_type ptype) const = 0;


  // Global ID of any entity

  virtual
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const = 0;



  //
  // Mesh Entity Adjacencies
  //-------------------------


  // Downward Adjacencies
  //---------------------

  // Get faces of a cell.

  // The .... coding guidelines regarding function arguments is purposely
  // violated here to allow for a default input argument

  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order

  unsigned int cell_get_num_faces(const Entity_ID cellid) const;

  void cell_get_faces (const Entity_ID cellid,
                       Entity_ID_List *faceids,
                       const bool ordered=false) const;

  // Get faces of a cell and directions in which the cell uses the face 

  // The ... coding guidelines regarding function arguments is purposely
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

  void cell_get_faces_and_dirs (const Entity_ID cellid,
                                Entity_ID_List *faceids,
                                std::vector<int> *facedirs,
				const bool ordered=false) const;


  // Get edges of a cell

  void cell_get_edges (const Entity_ID cellid, 
                       Entity_ID_List *edgeids) const;

  // Get edges and dirs of a 2D cell. This is to make the code cleaner
  // for integrating over the cell in 2D where faces and edges are
  // identical but integrating over the cells using face information
  // is more cumbersome (one would have to take the face normals,
  // rotate them and then get a consistent edge vector)

  void cell_2D_get_edges_and_dirs (const Entity_ID cellid, 
                                   Entity_ID_List *edgeids,
                                   std::vector<int> *edge_dirs) const;


  // Get nodes of a cell

  virtual
  void cell_get_nodes (const Entity_ID cellid, 
                       Entity_ID_List *nodeids) const=0;


  // Get edges of a face and directions in which the face uses the edges 

  // On a distributed mesh, this will return all the edges of the
  // face, OWNED or GHOST. If ordered = true, the edges will be
  // returned in a ccw order around the face as it is naturally defined.
 
  // IMPORTANT NOTE IN 2D CELLS: In meshes where the cells are two
  // dimensional, faces and edges are identical. For such cells, this
  // operator will return a single edge and a direction of 1. However,
  // this direction cannot be relied upon to compute, say, a contour
  // integral around the 2D cell. 

  void face_get_edges_and_dirs (const Entity_ID faceid,
                                Entity_ID_List *edgeids,
                                std::vector<int> *edgedirs,
				const bool ordered=false) const;


  // Get the local index of a face edge in a cell edge list
  // Example:
  //
  // face_get_edges(face=5) --> {20, 21, 35, 9, 10}
  // cell_get_edges(cell=18) --> {1, 2, 3, 5, 8, 9, 10, 13, 21, 35, 20, 37, 40}
  // face_to_cell_edge_map(face=5,cell=18) --> {10, 8, 9, 5, 6}


  void face_to_cell_edge_map(const Entity_ID faceid, 
                             const Entity_ID cellid,
                             std::vector<int> *map) const;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2

  virtual
  void face_get_nodes (const Entity_ID faceid,
                       Entity_ID_List *nodeids) const = 0;


  // Get nodes of edge
  
  virtual
  void edge_get_nodes (const Entity_ID edgeid, 
		       Entity_ID *nodeid0, Entity_ID *nodeid1) const = 0;

  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node - The order of cells is
  // not guaranteed to be the same for corresponding nodes on
  // different processors

  virtual
  void node_get_cells (const Entity_ID nodeid,
                       const Parallel_type ptype,
                       Entity_ID_List *cellids) const = 0;


  // Faces of type 'ptype' connected to a node - The order of faces is
  // not guarnateed to be the same for corresponding nodes on
  // different processors


  virtual
  void node_get_faces (const Entity_ID nodeid,
                       const Parallel_type ptype,
                       Entity_ID_List *faceids) const = 0;

  // Get faces of ptype of a particular cell that are connected to the
  // given node - The order of faces is not guarnateed to be the same
  // for corresponding nodes on different processors

  virtual
  void node_get_cell_faces (const Entity_ID nodeid,
                            const Entity_ID cellid,
                            const Parallel_type ptype,
                            Entity_ID_List *faceids) const = 0;

  // Cells connected to a face - The cells are returned in no
  // particular order. Also, the order of cells is not guaranteed to
  // be the same for corresponding faces on different processors

  void face_get_cells (const Entity_ID faceid,
                       const Parallel_type ptype,
                       Entity_ID_List *cellids) const;



  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = ALL, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces

  virtual
  void cell_get_face_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List *fadj_cellids) const = 0;

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order

  virtual
  void cell_get_node_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List *cellids) const = 0;


  //
  // Mesh entity geometry
  //--------------
  //

  // Node coordinates - 3 in 3D and 2 in 2D

  virtual
  void node_get_coordinates (const Entity_ID nodeid,
                             JaliGeometry::Point *ncoord) const = 0;


  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size divided by number of spatial dimensions

  virtual
  void face_get_coordinates (const Entity_ID faceid,
                             std::vector<JaliGeometry::Point> *fcoords) const = 0;

  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions

  virtual
  void cell_get_coordinates (const Entity_ID cellid,
                             std::vector<JaliGeometry::Point> *ccoords) const = 0;


  // Mesh entity geometry
  //--------------
  //


  // Volume/Area of cell

  double cell_volume (const Entity_ID cellid, const bool recompute=false) const;

  // Area/length of face

  double face_area(const Entity_ID faceid, const bool recompute=false) const;

  // Length of edge

  double edge_length(const Entity_ID edgeid, const bool recompute=false) const;

  // Centroid of cell

  JaliGeometry::Point cell_centroid (const Entity_ID cellid, const bool recompute=false) const;

  // Centroid of face

  JaliGeometry::Point face_centroid (const Entity_ID faceid, const bool recompute=false) const;

  // Normal to face
  // The vector is normalized and then weighted by the area of the face
  //
  // If recompute is TRUE, then the normal is recalculated using current
  // face coordinates but not stored. (If the recomputed normal must be
  // stored, then call recompute_geometric_quantities). 
  //
  // If cellid is not specified, the normal is the natural normal of
  // the face. This means that at boundaries, the normal may point in
  // or out of the domain depending on how the face is defined. On the
  // other hand, if cellid is specified, the normal is the outward
  // normal with respect to the cell. In planar and solid meshes, the
  // normal with respect to the cell on one side of the face is just
  // the negative of the normal with respect to the cell on the other
  // side. In general surfaces meshes, this will not be true at C1
  // discontinuities

  // if cellid is specified, then orientation returns the direction of
  // the natural normal of the face with respect to the cell (1 is
  // pointing out of the cell and -1 pointing in)


  JaliGeometry::Point face_normal (const Entity_ID faceid, 
				     const bool recompute=false, 
				     const Entity_ID cellid=-1, 
				     int *orientation=NULL) const;


  // Edge vector - not normalized (or normalized and weighted by length
  // of the edge)
  //
  // If recompute is TRUE, then the vector is recalculated using current
  // edge coordinates but not stored. (If the recomputed vector must be
  // stored, then call recompute_geometric_quantities). 
  //
  // If pointid is specified, the vector is the natural direction of
  // the edge (from point0 to point1).  On the other hand, if pointid
  // is specified (has to be a point of the face), the vector is from
  // specified point to opposite point of edge.
  //
  // if pointid is specified, then orientation returns the direction of
  // the natural direction of the edge with respect to the point (1 is
  // away from the point and -1 is towards)


  JaliGeometry::Point edge_vector (const Entity_ID edgeid, 
				     const bool recompute=false, 
				     const Entity_ID pointid=-1,
				     int *orientation=NULL) const;

  // Point in cell

  bool point_in_cell (const JaliGeometry::Point &p, 
                      const Entity_ID cellid) const;


  //
  // Mesh modification
  //-------------------

  // Set coordinates of node

  virtual
  void node_set_coordinates (const Entity_ID nodeid,
                             const JaliGeometry::Point ncoord) = 0;


  virtual
  void node_set_coordinates (const Entity_ID nodeid,
                             const double *ncoord) = 0;

  //
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------
  //

  // Is this is a valid ID of a set containing entities of 'kind'

  bool valid_set_name (const std::string setname,
                       const Entity_kind kind) const;


  // Get number of entities of type 'category' in set

  virtual
  unsigned int get_set_size (const Set_Name setname,
                             const Entity_kind kind,
                             const Parallel_type ptype) const = 0;

  virtual
  unsigned int get_set_size (const char *setname,
                             const Entity_kind kind,
                             const Parallel_type ptype) const = 0;


  // Get list of entities of type 'category' in set

  virtual
  void get_set_entities (const Set_Name setname,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         Entity_ID_List *entids) const = 0;

  virtual
  void get_set_entities (const char *setname,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         Entity_ID_List *entids) const = 0;

 protected:

  int compute_cell_geometric_quantities() const;
  int compute_face_geometric_quantities() const;
  int compute_edge_geometric_quantities() const;


  // get faces of a cell and directions in which it is used - this function
  // is implemented in each mesh framework. The results are cached in 
  // the base class

  virtual
  void cell_get_faces_and_dirs_internal (const Entity_ID cellid,
                                         Entity_ID_List *faceids,
                                         std::vector<int> *face_dirs,
                                         const bool ordered=false) const = 0;

  // Cells connected to a face - this function is implemented in each
  // mesh framework. The results are cached in the base class
  
  virtual
  void face_get_cells_internal (const Entity_ID faceid,
                                const Parallel_type ptype,
                                Entity_ID_List *cellids) const = 0;


  // edges of a face - this function is implemented in each mesh
  // framework. The results are cached in the base class

  virtual
  void face_get_edges_and_dirs_internal (const Entity_ID faceid,
					 Entity_ID_List *edgeids,
					 std::vector<int> *edge_dirs,
					 const bool ordered=true) const = 0;

  // edges of a cell - this function is implemented in each mesh
  // framework. The results are cached in the base class. 

  virtual
  void cell_get_edges_internal (const Entity_ID cellid,
                                Entity_ID_List *edgeids) const = 0;

  // edges and directions of a 2D cell - this function is implemented
  // in each mesh framework. The results are cached in the base class.

  virtual
  void cell_2D_get_edges_and_dirs_internal (const Entity_ID cellid,
                                            Entity_ID_List *edgeids,
                                            std::vector<int> *edge_dirs) const = 0;


  // Export mesh alone to different format files

  virtual
  void write_to_exodus_file(const std::string exodusfilename) const = 0;

  virtual
  void write_to_gmv_file(const std::string gmvfilename) const = 0;

 private:

  unsigned int celldim, spacedim;

  MPI_Comm comm;

  mutable std::vector<double> cell_volumes, face_areas, edge_lengths;
  mutable std::vector<JaliGeometry::Point> cell_centroids,
    face_centroids, face_normal0, face_normal1, edge_vectors;
  mutable std::vector<Entity_ID_List> cell_face_ids;
  mutable std::vector< std::vector<int> > cell_face_dirs;
  mutable std::vector<Entity_ID_List > face_cell_ids;
  mutable std::vector< std::vector<Parallel_type> > face_cell_ptype;
  mutable std::vector<Entity_ID_List> cell_edge_ids;
  mutable std::vector< std::vector<int> > cell_2D_edge_dirs;
  mutable std::vector<Entity_ID_List> face_edge_ids;
  mutable std::vector< std::vector<int> > face_edge_dirs;
  mutable Mesh_type mesh_type_;

  // flags to indicate what data is current

  mutable bool faces_requested, edges_requested;
  mutable bool cell2face_info_cached, face2cell_info_cached;
  mutable bool cell2edge_info_cached, face2edge_info_cached;
  mutable bool cell_geometry_precomputed, face_geometry_precomputed,
    edge_geometry_precomputed;

  JaliGeometry::GeometricModelPtr geometric_model_;


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
			    JaliGeometry::Point *edge_vector) const;


  void cache_cell2face_info() const; 
  void cache_face2cell_info() const;
  void cache_cell2edge_info() const;
  void cache_face2edge_info() const;

}; // End class Mesh


inline
void Mesh::cell_get_faces(const Entity_ID cellid, Entity_ID_List *faceids,
                          const bool ordered) const {
  cell_get_faces_and_dirs(cellid, faceids, NULL, ordered);
}




} // end namespace Jali





#endif /* _JALI_MESH_H_ */
