/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <algorithm>

#include "Mesh_simple.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Jali
{

Mesh_simple::Mesh_simple (double x0, double y0, double z0,
 			  double x1, double y1, double z1,
			  int nx, int ny, int nz,
                          const MPI_Comm& comm,
			  const JaliGeometry::GeometricModelPtr gm,
                          const bool request_faces,
                          const bool request_edges,
                          const bool request_wedges,
                          const bool request_corners) :
    nx_(nx), ny_(ny), nz_(nz),
    x0_(x0), x1_(x1),
    y0_(y0), y1_(y1),
    z0_(z0), z1_(z1),
    nodes_per_face_(4), faces_per_cell_(6), nodes_per_cell_(8),
    faces_per_node_aug_(13), cells_per_node_aug_(9),
    Mesh(request_faces,request_edges,request_wedges,request_corners)
{
  Mesh::set_mesh_type(RECTANGULAR);
  if (gm != (JaliGeometry::GeometricModelPtr) NULL) 
    Mesh::set_geometric_model(gm);

  clear_internals_3d_();
  update_internals_3d_();
}


// Have to define this dummy routine because we are not able to
// eliminate the need in FrameworkTraits.cc which uses boost
// functionality extensively

Mesh_simple::Mesh_simple (double x0, double y0,
                          double x1, double y1,
                          int nx, int ny, 
                          const MPI_Comm& comm,
                          const JaliGeometry::GeometricModelPtr &gm,
                          const bool request_faces,
                          const bool request_edges,
                          const bool request_wedges,
                          const bool request_corners) 
{
  Exceptions::Jali_throw(Errors::Message("Simple mesh cannot generate 2D meshes"));
}
  

// x is the node coordinates in 1d.
// keeping faces_per_node_aug and cells_per_node_aug to reuse 3d SimpleMesh
// logic

Mesh_simple::Mesh_simple (std::vector<double> x,
                          const MPI_Comm& comm,
                          const JaliGeometry::GeometricModelPtr &gm,
                          const bool request_faces,
                          const bool request_edges,
                          const bool request_wedges,
                          const bool request_corners) :
    nx_(x.size()-1), ny_(-3), nz_(-3),
    coordinates_(x),
    nodes_per_face_(1), faces_per_cell_(2), nodes_per_cell_(2),
    faces_per_node_aug_(2), cells_per_node_aug_(3),
    Mesh(request_faces,request_edges,request_wedges,request_corners)
{
  set_space_dimension(1);
  set_cell_dimension(1);

  clear_internals_1d_();
  update_internals_1d_();
  //  Exceptions::Jali_throw(Errors::Message("Simple mesh cannot generate 1D meshes"));  
}
  

//--------------------------------------
// Constructor - Construct a new mesh from a subset of an existing mesh
//--------------------------------------

Mesh_simple::Mesh_simple (const Mesh *inmesh, 
                          const std::vector<std::string>& setnames, 
                          const Entity_kind setkind,
                          const bool flatten,
                          const bool extrude,
                          const bool request_faces,
                          const bool request_edges,
                          const bool request_wedges,
                          const bool request_corners)
{  
  Errors::Message mesg("Construction of new mesh from an existing mesh not yet implemented in the Simple mesh framework\n");
  Exceptions::Jali_throw(mesg);
}

Mesh_simple::Mesh_simple (const Mesh& inmesh, 
                          const std::vector<std::string>& setnames, 
                          const Entity_kind setkind,
                          const bool flatten,
                          const bool extrude,
                          const bool request_faces,
                          const bool request_edges,
                          const bool request_wedges,
                         const bool request_corners)
{  
  Errors::Message mesg("Construction of new mesh from an existing mesh not yet implemented in the Simple mesh framework\n");
  Exceptions::Jali_throw(mesg);
}

Mesh_simple::Mesh_simple (const Mesh& inmesh, 
                          const std::vector<int>& entity_id_list, 
                          const Entity_kind entity_kind,
                          const bool flatten,
                          const bool extrude,
                          const bool request_faces,
                          const bool request_edges,
                          const bool request_wedges,
                          const bool request_corners)
{  
  Errors::Message mesg("Construction of new mesh from an existing mesh not yet implemented in the Simple mesh framework\n");
  Exceptions::Jali_throw(mesg);
}


Mesh_simple::~Mesh_simple()
{
}

void Mesh_simple::clear_internals_3d_ ()
{

  coordinates_.resize(0);

  cell_to_face_.resize(0);
  cell_to_node_.resize(0);
  face_to_node_.resize(0);

  side_sets_.resize(0);

}

void Mesh_simple::clear_internals_1d_ ()
{

  cell_to_face_.resize(0);
  cell_to_node_.resize(0);
  face_to_node_.resize(0);

  side_sets_.resize(0);

}

void Mesh_simple::update_internals_1d_ () {
  num_cells_ = nx_;
  num_nodes_ = nx_+1;
  num_faces_ = num_nodes_;

  cell_to_face_.resize(faces_per_cell_*num_cells_);
  cell_to_face_dirs_.resize(faces_per_cell_*num_cells_);
  cell_to_node_.resize(nodes_per_cell_*num_cells_);
  face_to_node_.resize(nodes_per_face_*num_faces_);
  node_to_face_.resize(faces_per_node_aug_*num_nodes_); // 1 extra for num faces
  node_to_cell_.resize(cells_per_node_aug_*num_nodes_); // 1 extra for num cells
  face_to_cell_.resize(2*num_faces_);
  face_to_cell_.assign(2*num_faces_,-1);

  // loop over cells and initialize cell_to_node_
  for (int ix=0; ix<nx_; ix++)
  {
    int jstart=0;
    int istart = nodes_per_cell_ * cell_index_(ix);
    int ncell=0;

    cell_to_node_[istart]   = node_index_(ix);
    cell_to_node_[istart+1] = node_index_(ix+1);

    // +1 because of num cells in node_to_cell_
    jstart = cells_per_node_aug_ * node_index_(ix);
    ncell = node_to_cell_[jstart];
    node_to_cell_[jstart+1+ncell] = cell_index_(ix);
    (node_to_cell_[jstart])++;

    jstart = cells_per_node_aug_ * node_index_(ix+1);
    ncell = node_to_cell_[jstart];
    node_to_cell_[jstart+1+ncell] = cell_index_(ix);
    (node_to_cell_[jstart])++;
  }

  // loop over cells and initialize cell_to_face_
  for (int ix=0; ix<nx_; ix++)
  {
    int istart = faces_per_cell_ * cell_index_(ix);
    int jstart = 0;

    cell_to_face_[istart]  = yzface_index_(ix+1);
    cell_to_face_[istart+1]  = yzface_index_(ix);

    cell_to_face_dirs_[istart] = 1;
    cell_to_face_dirs_[istart+1] = -1;

    jstart = 2*yzface_index_(ix+1);
    face_to_cell_[jstart+1] = cell_index_(ix);

    jstart = 2*yzface_index_(ix);
    face_to_cell_[jstart] = cell_index_(ix);
  }

  // loop over faces and initialize face_to_node_
  // the yz faces
  for (int ix=0; ix<=nx_; ix++)
  {
    int istart = nodes_per_face_ * yzface_index_(ix);
    int jstart = 0;
    int nfaces = 0;

    face_to_node_[istart]   = node_index_(ix);

    // +1 because of num faces in node_to_face_
    jstart = faces_per_node_aug_*node_index_(ix);
    nfaces = node_to_face_[jstart];
    node_to_face_[jstart+1+nfaces] = xyface_index_(ix);
    (node_to_face_[jstart])++;
  }

  // populate entity ids arrays in the base class so that iterators work

  int lid;
  std::vector<int>::iterator it;

  Mesh::nodeids.resize(num_nodes_);
  lid = 0;
  it = Mesh::nodeids.begin();
  while (it != Mesh::nodeids.end()) {
    *it = lid++;
    ++it;
  }

  Mesh::faceids.resize(num_faces_);
  lid = 0;
  it = Mesh::faceids.begin();
  while (it != Mesh::faceids.end()) {
    *it = lid++;
    ++it;
  }

  Mesh::cellids.resize(num_cells_);
  lid = 0;
  it = Mesh::cellids.begin();
  while (it != Mesh::cellids.end()) {
    *it = lid++;
    ++it;
  }

}


void Mesh_simple::update_internals_3d_ ()
{

  num_cells_ = nx_ * ny_ * nz_;
  num_nodes_ = (nx_+1)*(ny_+1)*(nz_+1);
  num_faces_ = (nx_+1)*(ny_)*(nz_) + (nx_)*(ny_+1)*(nz_) + (nx_)*(ny_)*(nz_+1);

  coordinates_.resize(3*num_nodes_);

  double hx = (x1_ - x0_)/nx_;
  double hy = (y1_ - y0_)/ny_;
  double hz = (z1_ - z0_)/nz_;

  for (int iz=0; iz<=nz_; iz++)
    for (int iy=0; iy<=ny_; iy++)
      for (int ix=0; ix<=nx_; ix++)
      {
        int istart = 3*node_index_(ix,iy,iz);

        coordinates_[ istart ]     = x0_ + ix*hx;
        coordinates_[ istart + 1 ] = y0_ + iy*hy;
        coordinates_[ istart + 2 ] = z0_ + iz*hz;
      }

  cell_to_face_.resize(faces_per_cell_*num_cells_);
  cell_to_face_dirs_.resize(faces_per_cell_*num_cells_);
  cell_to_node_.resize(nodes_per_cell_*num_cells_);
  face_to_node_.resize(nodes_per_face_*num_faces_);
  node_to_face_.resize(faces_per_node_aug_*num_nodes_); // 1 extra for num faces
  node_to_cell_.resize(cells_per_node_aug_*num_nodes_); // 1 extra for num cells
  face_to_cell_.resize(2*num_faces_); 
  face_to_cell_.assign(2*num_faces_,-1); 
                                          

  // loop over cells and initialize cell_to_node_
  for (int iz=0; iz<nz_; iz++)
    for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<nx_; ix++)
      {
        int jstart=0;
        int istart = nodes_per_cell_ * cell_index_(ix,iy,iz);
        int ncell=0;

        cell_to_node_[istart]   = node_index_(ix,iy,iz);
        cell_to_node_[istart+1] = node_index_(ix+1,iy,iz);
        cell_to_node_[istart+2] = node_index_(ix+1,iy+1,iz);
        cell_to_node_[istart+3] = node_index_(ix,iy+1,iz);
        cell_to_node_[istart+4] = node_index_(ix,iy,iz+1);
        cell_to_node_[istart+5] = node_index_(ix+1,iy,iz+1);
        cell_to_node_[istart+6] = node_index_(ix+1,iy+1,iz+1);
        cell_to_node_[istart+7] = node_index_(ix,iy+1,iz+1);

        jstart = cells_per_node_aug_ * node_index_(ix,iy,iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = cells_per_node_aug_ * node_index_(ix+1,iy,iz); 
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = cells_per_node_aug_ * node_index_(ix+1,iy+1,iz); 
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = cells_per_node_aug_ * node_index_(ix,iy+1,iz); 
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

        jstart = cells_per_node_aug_ * node_index_(ix,iy,iz+1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        node_to_cell_[jstart]++;

        // 1 extra for num cells
        jstart = cells_per_node_aug_ * node_index_(ix+1,iy,iz+1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

         // 1 extra for num cells
        jstart = cells_per_node_aug_ * node_index_(ix+1,iy+1,iz+1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

         // 1 extra for num cells
        jstart = cells_per_node_aug_ * node_index_(ix,iy+1,iz+1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix,iy,iz);
        (node_to_cell_[jstart])++;

      }


  // loop over cells and initialize cell_to_face_
  for (int iz=0; iz<nz_; iz++)
    for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<nx_; ix++)
      {
        int istart = faces_per_cell_ * cell_index_(ix,iy,iz);
        int jstart = 0;

        cell_to_face_[istart]    = xzface_index_(ix,iy,iz);
        cell_to_face_[istart+1]  = yzface_index_(ix+1,iy,iz);
        cell_to_face_[istart+2]  = xzface_index_(ix,iy+1,iz);
        cell_to_face_[istart+3]  = yzface_index_(ix,iy,iz);
        cell_to_face_[istart+4]  = xyface_index_(ix,iy,iz);
        cell_to_face_[istart+5]  = xyface_index_(ix,iy,iz+1);

        cell_to_face_dirs_[istart]   = 1;
        cell_to_face_dirs_[istart+1] = 1;
        cell_to_face_dirs_[istart+2] = -1;
        cell_to_face_dirs_[istart+3] = -1;
        cell_to_face_dirs_[istart+4] = -1;
        cell_to_face_dirs_[istart+5] = 1;

        jstart = 2*xzface_index_(ix,iy,iz);
        face_to_cell_[jstart+1] = cell_index_(ix,iy,iz);

        jstart = 2*yzface_index_(ix+1,iy,iz);
        face_to_cell_[jstart+1] = cell_index_(ix,iy,iz);

        jstart = 2*xzface_index_(ix,iy+1,iz);
        face_to_cell_[jstart] = cell_index_(ix,iy,iz);

        jstart = 2*yzface_index_(ix,iy,iz);
        face_to_cell_[jstart] = cell_index_(ix,iy,iz);

        jstart = 2*xyface_index_(ix,iy,iz);
        face_to_cell_[jstart] = cell_index_(ix,iy,iz);

        jstart = 2*xyface_index_(ix,iy,iz+1);
        face_to_cell_[jstart+1] = cell_index_(ix,iy,iz);
      }


  // loop over faces and initialize face_to_node_
  // first we do the xy faces
  for (int iz=0; iz<=nz_; iz++)
    for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<nx_; ix++)  
      {
        int istart = nodes_per_face_ * xyface_index_(ix,iy,iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]   = node_index_(ix,iy,iz);
        face_to_node_[istart+1] = node_index_(ix+1,iy,iz);
        face_to_node_[istart+2] = node_index_(ix+1,iy+1,iz);
        face_to_node_[istart+3] = node_index_(ix,iy+1,iz);

        jstart = faces_per_node_aug_*node_index_(ix,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix+1,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix+1,iy+1,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix,iy+1,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

      }
  // then we do the xz faces
  for (int iz=0; iz<nz_; iz++)
    for (int iy=0; iy<=ny_; iy++)
      for (int ix=0; ix<nx_; ix++)  
      {
        int istart = nodes_per_face_ * xzface_index_(ix,iy,iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]   = node_index_(ix,iy,iz);
        face_to_node_[istart+1] = node_index_(ix+1,iy,iz);
        face_to_node_[istart+2] = node_index_(ix+1,iy,iz+1);
        face_to_node_[istart+3] = node_index_(ix,iy,iz+1);

        jstart = faces_per_node_aug_*node_index_(ix,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix+1,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix+1,iy,iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix,iy,iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

      }
  // finally we do the yz faces
  for (int iz=0; iz<nz_; iz++)
    for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<=nx_; ix++)  
      {
        int istart = nodes_per_face_ * yzface_index_(ix,iy,iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]   = node_index_(ix,iy,iz);
        face_to_node_[istart+1] = node_index_(ix,iy+1,iz);
        face_to_node_[istart+2] = node_index_(ix,iy+1,iz+1);
        face_to_node_[istart+3] = node_index_(ix,iy,iz+1);

        jstart = faces_per_node_aug_*node_index_(ix,iy,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix,iy+1,iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix,iy+1,iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix,iy,iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix,iy,iz);
        (node_to_face_[jstart])++;

      }

  // populate entity ids arrays in the base class so that iterators work

  int lid;
  std::vector<int>::iterator it;

  Mesh::nodeids.resize(num_nodes_);
  lid = 0;
  it = Mesh::nodeids.begin();
  while (it != Mesh::nodeids.end()) {
    *it = lid++;
    ++it;
  }

  Mesh::faceids.resize(num_faces_);
  lid = 0;
  it = Mesh::faceids.begin();
  while (it != Mesh::faceids.end()) {
    *it = lid++;
    ++it;
  }

  Mesh::cellids.resize(num_cells_);
  lid = 0;
  it = Mesh::cellids.begin();
  while (it != Mesh::cellids.end()) {
    *it = lid++;
    ++it;
  }

}


Parallel_type Mesh_simple::entity_get_ptype(const Entity_kind kind, 
					    const Entity_ID entid) const 
{
  return OWNED; // Its a serial code
}




// Get cell type
    
Jali::Cell_type Mesh_simple::cell_get_type(const Jali::Entity_ID cellid) const 
{
  return HEX;
}
        
    
Entity_ID Mesh_simple::GID(const Jali::Entity_ID lid, 
			      const Jali::Entity_kind kind) const
{
  return lid;  // Its a serial code
}



unsigned int Mesh_simple::num_entities (Jali::Entity_kind kind, 
					Jali::Parallel_type ptype) const
{
  switch (kind) {
    case Jali::FACE: 
      return (ptype != Jali::GHOST) ? num_faces_ : 0;
      break;
    case Jali::NODE:
      return (ptype != Jali::GHOST) ? num_nodes_ : 0;
      break;
    case Jali::CELL:
      return (ptype != Jali::GHOST) ? num_cells_ : 0;
      break;
    default:
      throw std::exception();
      break;
  }
}


void Mesh_simple::cell_get_faces_and_dirs_internal (const Jali::Entity_ID cellid,
                                           Jali::Entity_ID_List *faceids,
                                           std::vector<int> *cfacedirs,
                                           const bool ordered) const
{
  unsigned int offset = (unsigned int) faces_per_cell_*cellid;

  faceids->clear();
  if (cfacedirs) cfacedirs->clear();

  for (int i = 0; i < faces_per_cell_; i++) {
    faceids->push_back(*(cell_to_face_.begin()+offset));
    if (cfacedirs)
      cfacedirs->push_back(*(cell_to_face_dirs_.begin()+offset));
    offset++;
  }
}



void Mesh_simple::cell_get_nodes (Jali::Entity_ID cell, 
				  Jali::Entity_ID_List *nodeids) const
{
  unsigned int offset = (unsigned int) nodes_per_cell_*cell;

  nodeids->clear();

  for (int i = 0; i < nodes_per_cell_; i++) {
    nodeids->push_back(*(cell_to_node_.begin()+offset));
    offset++;
  }
}


void Mesh_simple::face_get_nodes (Jali::Entity_ID face, 
				  Jali::Entity_ID_List *nodeids) const
{
  unsigned int offset = (unsigned int) nodes_per_face_*face;

  nodeids->clear();

  for (int i = 0; i < nodes_per_face_; i++) {
    nodeids->push_back(*(face_to_node_.begin()+offset));
    offset++;
  }
}


// Cooordinate Getters
// -------------------

void Mesh_simple::node_get_coordinates (const Jali::Entity_ID local_node_id, 
					JaliGeometry::Point *ncoords) const
{
  unsigned int offset = (unsigned int) Mesh::space_dimension()*local_node_id;

  ncoords->set( Mesh::space_dimension(), &(coordinates_[offset]) );
}


void Mesh_simple::face_get_coordinates (Jali::Entity_ID local_face_id, 
					std::vector<JaliGeometry::Point> *fcoords) const
{
  Entity_ID_List node_indices(nodes_per_face_);

  face_get_nodes (local_face_id, &node_indices);

  fcoords->clear();

  JaliGeometry::Point node_coords(Mesh::space_dimension());
  for (std::vector<Entity_ID>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
  {
    node_get_coordinates ( *it, &node_coords);
    fcoords->push_back(node_coords);
  }

}

void Mesh_simple::cell_get_coordinates (Jali::Entity_ID local_cell_id, 
					std::vector<JaliGeometry::Point> *ccoords) const
{  
  std::vector<Entity_ID> node_indices(nodes_per_cell_);

  cell_get_nodes (local_cell_id, &node_indices);

  ccoords->clear();

  JaliGeometry::Point node_coords(Mesh::space_dimension());
  for (std::vector<Entity_ID>::iterator it = node_indices.begin(); 
       it != node_indices.end(); ++it)
  {      
    node_get_coordinates ( *it, &node_coords);
    ccoords->push_back(node_coords);
  }
}


void Mesh_simple::node_set_coordinates(const Jali::Entity_ID local_node_id, 
                                      const double *ncoord)
{
  int spdim = Mesh::space_dimension();
  unsigned int offset = (unsigned int) spdim*local_node_id;


  ASSERT(ncoord != NULL);
  
  std::vector<double>::iterator destination_begin = coordinates_.begin() + offset;
  for (int i = 0; i < spdim; i++) {
    *destination_begin = ncoord[i];
    destination_begin++;
  }
}

void Mesh_simple::node_set_coordinates(const Jali::Entity_ID local_node_id, 
                                       const JaliGeometry::Point ncoord)
{
  int spdim = Mesh::space_dimension();
  unsigned int offset = (unsigned int) spdim*local_node_id;

  std::vector<double>::iterator destination_begin = coordinates_.begin() + offset;
  for (int i = 0; i < spdim; i++) {
    *destination_begin = ncoord[i];
    destination_begin++;
  }
}



void Mesh_simple::node_get_cells (const Jali::Entity_ID nodeid, 
                                  const Jali::Parallel_type ptype,
                                  Jali::Entity_ID_List *cellids) const 
{
  unsigned int offset = (unsigned int) cells_per_node_aug_*nodeid;
  unsigned int ncells = node_to_cell_[offset];

  cellids->clear();

  for (int i = 0; i < ncells; i++) 
    cellids->push_back(node_to_cell_[offset+i+1]);
}
    


// Faces of type 'ptype' connected to a node
    
void Mesh_simple::node_get_faces (const Jali::Entity_ID nodeid, 
                                  const Jali::Parallel_type ptype,
                                  Jali::Entity_ID_List *faceids) const
{
  unsigned int offset = (unsigned int) faces_per_node_aug_*nodeid;
  unsigned int nfaces = node_to_face_[offset];

  faceids->clear();

  for (int i = 0; i < nfaces; i++) 
    faceids->push_back(node_to_face_[offset+i+1]);
}   


 
// Get faces of ptype of a particular cell that are connected to the
// given node

void Mesh_simple::node_get_cell_faces (const Jali::Entity_ID nodeid, 
                                       const Jali::Entity_ID cellid,
                                       const Jali::Parallel_type ptype,
                                       Jali::Entity_ID_List *faceids) const
{
  unsigned int offset = (unsigned int) faces_per_cell_*cellid;

  faceids->clear();

  for (int i = 0; i < faces_per_cell_; i++) {
    Entity_ID cellfaceid = cell_to_face_[offset+i];

    unsigned int offset2 = (unsigned int) nodes_per_face_*cellfaceid;

    Jali::Entity_ID_List fnodes;
    face_get_nodes(cellfaceid,&fnodes);
 
    for (int j = 0; j < nodes_per_face_; j++) {

      if (face_to_node_[offset2+j] == nodeid) {
	faceids->push_back(cellfaceid);
	break;
      }
    }
  }
}
    
// Cells connected to a face
    
void Mesh_simple::face_get_cells_internal (const Jali::Entity_ID faceid, 
                                         const Jali::Parallel_type ptype,
                                         Jali::Entity_ID_List *cellids) const
{
  unsigned int offset = (unsigned int) 2*faceid;

  cellids->clear();

  if (face_to_cell_[offset] != -1)
    cellids->push_back(face_to_cell_[offset]);
  if (face_to_cell_[offset+1] != -1)
    cellids->push_back(face_to_cell_[offset+1]);
}
    


// Same level adjacencies
//-----------------------

// Face connected neighboring cells of given cell of a particular ptype
// (e.g. a hex has 6 face neighbors)

// The order in which the cellids are returned cannot be
// guaranteed in general except when ptype = USED, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces

void Mesh_simple::cell_get_face_adj_cells(const Jali::Entity_ID cellid,
                                          const Jali::Parallel_type ptype,
                                          Jali::Entity_ID_List *fadj_cellids) const
{
  unsigned int offset = (unsigned int) faces_per_cell_*cellid;

  fadj_cellids->clear();

  for (int i = 0; i < faces_per_cell_; i++) {    
    Entity_ID faceid = cell_to_face_[offset];

    unsigned int foffset = (unsigned int) 2*faceid;

    int adjcell0 = face_to_cell_[foffset];
    if (adjcell0 != -1 && adjcell0 != cellid)
      fadj_cellids->push_back(adjcell0);    
    else {
      int adjcell1 = face_to_cell_[foffset+1];
      if (adjcell1 != -1 && adjcell1 != cellid)
	fadj_cellids->push_back(adjcell1);
    }

    offset++;
  }
}

// Node connected neighboring cells of given cell
// (a hex in a structured mesh has 26 node connected neighbors)
// The cells are returned in no particular order

void Mesh_simple::cell_get_node_adj_cells(const Jali::Entity_ID cellid,
                                          const Jali::Parallel_type ptype,
                                          Jali::Entity_ID_List *nadj_cellids) const
{
  unsigned int offset = (unsigned int) nodes_per_cell_*cellid;

  nadj_cellids->clear();
  
  for (int i = 0; i < nodes_per_cell_; i++) {
    Entity_ID nodeid = cell_to_node_[offset+i];

    unsigned int offset2 = (unsigned int) cells_per_node_aug_*nodeid;
    unsigned int ncell = node_to_cell_[offset2];

    for (int j = 0; j < nodes_per_cell_; j++) {
      Entity_ID nodecell = node_to_cell_[offset2+j];
      
      unsigned int found = 0;
      unsigned int size = nadj_cellids->size();
      for (int k = 0; k < size; k++) {
	if ((*nadj_cellids)[k] == nodecell) {
	  found = 1;
	  break;
	}
      }

      if (!found)
	nadj_cellids->push_back(nodecell);
    }
  }

}

    
unsigned int Mesh_simple::get_set_size (const char *setname,
					const Jali::Entity_kind kind,
					const Jali::Parallel_type ptype) const
{
  Entity_ID_List setents;
  get_set_entities(setname,kind,ptype,&setents);

  return setents.size();
}

unsigned int Mesh_simple::get_set_size (const std::string setname, 
					const Jali::Entity_kind kind,
					const Jali::Parallel_type ptype) const
{
  Entity_ID_List setents;
  get_set_entities(setname,kind,ptype,&setents);

  return setents.size();
}


void Mesh_simple::get_set_entities (const char *setname, 
				    const Jali::Entity_kind kind, 
				    const Jali::Parallel_type ptype, 
				    Jali::Entity_ID_List *setents) const
{
  std::string setname1(setname);
  get_set_entities(setname1,kind,ptype,setents);
}

void Mesh_simple::get_set_entities (const std::string setname, 
				    const Jali::Entity_kind kind, 
				    const Jali::Parallel_type ptype, 
				    Jali::Entity_ID_List *setents) const
{
  // we ignore ptype since this is a serial implementation

  JaliGeometry::GeometricModelPtr gm = Mesh::geometric_model();

  if (setents == NULL) {
    Errors::Message mesg("NULL setents list");
    Exceptions::Jali_throw(mesg);
  }
  setents->clear();
  
  switch (kind) {
    case Jali::FACE:
      {
      Entity_ID_List ss;

      // Does this set exist?

      int nss = side_sets_.size();  // number of side sets
      bool found = false;
      for (int i = 0; i < nss; i++) {
        if ((side_set_regions_[i])->name() == setname) {
          found = true;
          ss = side_sets_[i];
        }
      }

      if (!found) // create the side set from the region definition
        {
          JaliGeometry::RegionPtr rgn = gm->FindRegion(setname);
          
          if (rgn == NULL) 
            {
              std::cerr << "Geometric model has no region named " << setname << std::endl;
              std::cerr << "Cannot construct set by this name" << std::endl;
              throw std::exception();
            }
          
          if (rgn->dimension() != Mesh::cell_dimension()-1)
            {
              std::cerr << "Geometric model does not have a region named" << setname << "with the appropriate dimension" << std::endl;
              std::cerr << "Cannot construct set by this name" << std::endl;
              throw std::exception();
            }

          if (rgn->type() == JaliGeometry::BOX ||
              rgn->type() == JaliGeometry::PLANE) 
            {
              for (int ix=0; ix<nx_; ix++)
                for (int iz=0; iz<nz_; iz++)
                  {
                    int face;
                    std::vector<JaliGeometry::Point> fxyz;
                    bool inbox;

                    face = xzface_index_(ix,0,iz);
                    face_get_coordinates(face,&fxyz);

                    inbox = true;
                    if (rgn->type() == JaliGeometry::BOX)
                      {
                        JaliGeometry::Point cenpnt(Mesh::space_dimension());
                        for (int j = 0; j < nodes_per_face_; j++)
                          cenpnt += fxyz[j];
                        cenpnt /= static_cast<double>(nodes_per_face_);

                        if (!rgn->inside(cenpnt)) inbox = false;
                      }
                    else
                      {
                        for (int j = 0; j < nodes_per_face_; j++)
                          {
                            if (!rgn->inside(fxyz[j]))
                              {
                                inbox = false;
                                break;
                              }
                          }                    
                      }
                      
                    if (inbox)
                      ss.push_back(face);

                    face = xzface_index_(ix,ny_,iz);
                    face_get_coordinates(face,&fxyz);

                    inbox = true;
                    if (rgn->type() == JaliGeometry::BOX)
                      {
                        JaliGeometry::Point cenpnt(Mesh::space_dimension());
                        for (int j = 0; j < nodes_per_face_; j++)
                          cenpnt += fxyz[j];
                        cenpnt /= static_cast<double>(nodes_per_face_);

                        if (!rgn->inside(cenpnt)) inbox = false;
                      }
                    else
                      {
                        for (int j = 0; j < nodes_per_face_; j++)
                          {
                            if (!rgn->inside(fxyz[j]))
                              {
                                inbox = false;
                                break;
                              }
                          }                    
                      }
                      
                    if (inbox)
                      ss.push_back(face);
                  }


              for (int ix=0; ix<nx_; ix++)
                for (int iy=0; iy<ny_; iy++)
                  {
                    int face;
                    std::vector<JaliGeometry::Point> fxyz;
                    bool inbox;

                    face = xyface_index_(ix,iy,0);
                    face_get_coordinates(face,&fxyz);

                    inbox = true;
                    if (rgn->type() == JaliGeometry::BOX)
                      {
                        JaliGeometry::Point cenpnt(Mesh::space_dimension());
                        for (int j = 0; j < nodes_per_face_; j++)
                          cenpnt += fxyz[j];
                        cenpnt /= static_cast<double>(nodes_per_face_);

                        if (!rgn->inside(cenpnt)) inbox = false;
                      }
                    else
                      {
                        for (int j = 0; j < 4; j++)
                          {
                            if (!rgn->inside(fxyz[j]))
                              {
                                inbox = false;
                                break;
                              }
                          }                    
                      }
                      
                    if (inbox)
                      ss.push_back(face);

                    face = xyface_index_(ix,iy,nz_);
                    face_get_coordinates(face,&fxyz);

                    inbox = true;
                    if (rgn->type() == JaliGeometry::BOX)
                      {
                        JaliGeometry::Point cenpnt(Mesh::space_dimension());
                        for (int j = 0; j < nodes_per_face_; j++)
                          cenpnt += fxyz[j];
                        cenpnt /= static_cast<double>(nodes_per_face_);

                        if (!rgn->inside(cenpnt)) inbox = false;
                      }
                    else
                      {
                        for (int j = 0; j < nodes_per_face_; j++)
                          {
                            if (!rgn->inside(fxyz[j]))
                              {
                                inbox = false;
                                break;
                              }
                          }                    
                      }
                      
                    if (inbox)
                      ss.push_back(face);
                  }

              for (int iy=0; iy<ny_; iy++)
                for (int iz=0; iz<nz_; iz++)
                  {
                    int face;
                    std::vector<JaliGeometry::Point> fxyz;
                    bool inbox;

                    face = yzface_index_(0,iy,iz);
                    face_get_coordinates(face,&fxyz);

                    inbox = true;
                    if (rgn->type() == JaliGeometry::BOX)
                      {
                        JaliGeometry::Point cenpnt(Mesh::space_dimension());
                        for (int j = 0; j < nodes_per_face_; j++)
                          cenpnt += fxyz[j];
                        cenpnt /= static_cast<double>(nodes_per_face_);

                        if (!rgn->inside(cenpnt)) inbox = false;
                      }
                    else
                      {
                        for (int j = 0; j < nodes_per_face_; j++)
                          {
                            if (!rgn->inside(fxyz[j]))
                              {
                                inbox = false;
                                break;
                              }
                          }                    
                      }
                      
                    if (inbox)
                      ss.push_back(face);

                    face = yzface_index_(nx_,iy,iz);
                    face_get_coordinates(face,&fxyz);

                    inbox = true;
                    if (rgn->type() == JaliGeometry::BOX)
                      {
                        JaliGeometry::Point cenpnt(Mesh::space_dimension());
                        for (int j = 0; j < nodes_per_face_; j++)
                          cenpnt += fxyz[j];
                        cenpnt /= static_cast<double>(nodes_per_face_);

                        if (!rgn->inside(cenpnt)) inbox = false;
                      }
                    else
                      {
                        for (int j = 0; j < nodes_per_face_; j++)
                          {
                            if (!rgn->inside(fxyz[j]))
                              {
                                inbox = false;
                                break;
                              }
                          }                    
                      }
                      
                    if (inbox)
                      ss.push_back(face);
                  }

              if (Mesh::space_dimension() == 1) {
                  int face;
                  std::vector<JaliGeometry::Point> fx;
                  bool inbox;

                  face = yzface_index_(0);
                  face_get_coordinates(face,&fx);

                  inbox = true;
                  if (rgn->type() == JaliGeometry::BOX)
                  {
                    JaliGeometry::Point cenpnt(Mesh::space_dimension());
                    for (int j = 0; j < nodes_per_face_; j++)
                      cenpnt += fx[j];
                    cenpnt /= static_cast<double>(nodes_per_face_);

                    if (!rgn->inside(cenpnt)) inbox = false;
                  }
                  else
                  {
                    for (int j = 0; j < nodes_per_face_; j++)
                    {
                      if (!rgn->inside(fx[j]))
                      {
                        inbox = false;
                        break;
                      }
                    }                    
                  }
                      
                  if (inbox)
                    ss.push_back(face);

                  face = yzface_index_(nx_);
                  face_get_coordinates(face,&fx);

                  inbox = true;
                  if (rgn->type() == JaliGeometry::BOX)
                  {
                    JaliGeometry::Point cenpnt(Mesh::space_dimension());
                    for (int j = 0; j < nodes_per_face_; j++)
                      cenpnt += fx[j];
                    cenpnt /= static_cast<double>(nodes_per_face_);

                    if (!rgn->inside(cenpnt)) inbox = false;
                  }
                  else
                  {
                    for (int j = 0; j < nodes_per_face_; j++)
                    {
                      if (!rgn->inside(fx[j]))
                      {
                        inbox = false;
                        break;
                      }
                    }                    
                  }
                      
                  if (inbox)
                    ss.push_back(face);
              }

              side_sets_.push_back(ss);
              side_set_regions_.push_back(rgn);
            }
          else
            {
              std::cerr << "Region type not suitable/applicable for sidesets" << std::endl;
              throw std::exception();
            }
        }

      *setents = ss;
      break;
      }
    case Jali::CELL:
      {
      Entity_ID_List cs; // cell set
          
      int ncs = element_blocks_.size();
      bool found = false;
      for (int i = 0; i < ncs; i++) {
        if ((element_block_regions_[i])->name() == setname) {
          found = true;
          cs = element_blocks_[i];
        }
      }

      if (!found) // create the cell set from the region definition
        {
          JaliGeometry::RegionPtr rgn = gm->FindRegion(setname);
          
          if (rgn == NULL) 
            {
              std::cerr << "Geometric model has no region named " << setname << std::endl;
              std::cerr << "Cannot construct set by this name" << std::endl;
              throw std::exception();
            }
          
          if (rgn->dimension() != Mesh::cell_dimension())
            {
              std::cerr << "Geometric model does not have a region named" << setname << "with the appropriate dimension" << std::endl;
              std::cerr << "Cannot construct set by this name" << std::endl;
              throw std::exception();
            }

          if (rgn->type() == JaliGeometry::BOX || 
              rgn->type() == JaliGeometry::COLORFUNCTION)
            {
              for (int ix=0; ix<nx_; ix++)
                for (int iy=0; iy<ny_; iy++)
                  for (int iz=0; iz<nz_; iz++)
                    {
                      int cell = cell_index_(ix,iy,iz);
                      std::vector<JaliGeometry::Point> cxyz;
                      cell_get_coordinates(cell,&cxyz);

                      JaliGeometry::Point cenpnt(Mesh::space_dimension());

                      for (int j = 0; j < nodes_per_cell_; j++)
                        cenpnt += cxyz[j];
                      cenpnt /= static_cast<double>(nodes_per_cell_);

                      if (rgn->inside(cenpnt))
                        cs.push_back(cell);
                    }

              if (Mesh::space_dimension() == 1) {
                for (int ix=0; ix<nx_; ix++)
                {
                  int cell = cell_index_(ix);
                  std::vector<JaliGeometry::Point> cx;
                  cell_get_coordinates(cell,&cx);

                  JaliGeometry::Point cenpnt(Mesh::space_dimension());

                  for (int j = 0; j < nodes_per_cell_; j++)
                    cenpnt += cx[j];
                  cenpnt /= static_cast<double>(nodes_per_cell_);

                  if (rgn->inside(cenpnt))
                    cs.push_back(cell);
                }
                
              }

              element_blocks_.push_back(cs);
              element_block_regions_.push_back(rgn);
            }
          else
            {
              std::cerr << "Region type not suitable/applicable for cellsets" << std::endl;
              throw std::exception();
            }
        }

      *setents = cs;
      break;
      }
  case Jali::NODE:
    {
      Entity_ID_List ns;

      // Does this set exist?

      int nns = node_sets_.size();  // number of node sets
      bool found = false;
      for (int i = 0; i < nns; i++) {
        if ((node_set_regions_[i])->name() == setname) {
          found = true;
          ns = node_sets_[i];
        }
      }

      if (!found) // create the side set from the region definition
        {
          JaliGeometry::RegionPtr rgn = gm->FindRegion(setname);
          
          if (rgn == NULL) 
            {
              std::cerr << "Geometric model has no region named " << setname << std::endl;
              std::cerr << "Cannot construct set by this name" << std::endl;
              throw std::exception();
            }


          bool done=false;
          for (int ix=0; ix<nx_+1 && !done; ix++)
            for (int iy=0; iy<ny_+1 && !done; iy++)
              for (int iz=0; iz<nz_+1 && !done; iz++)
                {
                  int node = node_index_(ix,iy,iz);
                  JaliGeometry::Point nxyz;
                  node_get_coordinates(node,&nxyz);
                  
                  if (rgn->inside(nxyz)) {
                    ns.push_back(node);
                    
                    if (rgn->type() == JaliGeometry::POINT)
                      done = true;
                  }                   
                }

          if (Mesh::space_dimension() == 1) {
            bool done=false;
            for (int ix=0; ix<nx_+1 && !done; ix++)
            {
              int node = node_index_(ix);
              JaliGeometry::Point nx;
              node_get_coordinates(node,&nx);

              if (rgn->inside(nx)) {
                ns.push_back(node);

                if (rgn->type() == JaliGeometry::POINT)
                  done = true;
              }
            }
            
          }
              
          node_sets_.push_back(ns);
          node_set_regions_.push_back(rgn);
        }      
      
      *setents = ns;
      break;
    }
  default:
    // set type not recognized
    throw std::exception();
  }
}




} // close namespace Jali
