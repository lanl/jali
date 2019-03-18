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

#include "Mesh_simple.hh"

#include <algorithm>

#include "mpi.h"   // only for MPI_COMM_WORLD in Mesh constructor

#include "dbc.hh"
#include "errors.hh"

namespace Jali {

Mesh_simple::Mesh_simple(double x0, double y0, double z0,
                         double x1, double y1, double z1,
                         int nx, int ny, int nz,
                         const MPI_Comm& comm,
                         const JaliGeometry::GeometricModelPtr gm,
                         const bool request_faces,
                         const bool request_edges,
                         const bool request_sides,
                         const bool request_wedges,
                         const bool request_corners,
                         const int num_tiles_ini,
                         const int num_ghost_layers_tile,
                         const int num_ghost_layers_distmesh,
                         const bool boundary_ghosts_requested,
                         const Partitioner_type partitioner) :
    nx_(nx), ny_(ny), nz_(nz),
    x0_(x0), x1_(x1),
    y0_(y0), y1_(y1),
    z0_(z0), z1_(z1),
    nodes_per_face_(4), faces_per_cell_(6), nodes_per_cell_(8),
    faces_per_node_aug_(13), cells_per_node_aug_(9),
  Mesh(request_faces, request_edges, request_sides, request_wedges,
       request_corners, num_tiles_ini, num_ghost_layers_tile,
       num_ghost_layers_distmesh, boundary_ghosts_requested,
       partitioner, JaliGeometry::Geom_type::CARTESIAN, comm) {

  assert(!boundary_ghosts_requested);  // Cannot yet make boundary ghosts

  Mesh::set_mesh_type(Mesh_type::RECTANGULAR);
  if (gm != (JaliGeometry::GeometricModelPtr) NULL)
    Mesh::set_geometric_model(gm);

  clear_internals_3d_();
  update_internals_3d_();

  cache_extra_variables();

  if (Mesh::num_tiles_ini_)
    Mesh::build_tiles();
}


// Have to define this dummy routine because we are not able to
// eliminate the need in FrameworkTraits.cc which uses boost
// functionality extensively

Mesh_simple::Mesh_simple(double x0, double y0,
                         double x1, double y1,
                         int nx, int ny,
                         const MPI_Comm& comm,
                         const JaliGeometry::GeometricModelPtr &gm,
                         const bool request_faces,
                         const bool request_edges,
                         const bool request_sides,
                         const bool request_wedges,
                         const bool request_corners,
                         const int num_tiles_ini,
                         const int num_ghost_layers_tile,
                         const int num_ghost_layers_distmesh,
                         const bool boundary_ghosts_requested,
                         const Partitioner_type partitioner,
                         const JaliGeometry::Geom_type geom_type) {
  Exceptions::Jali_throw(Errors::Message("Simple mesh cannot generate 2D meshes"));
}

// x is the node coordinates in 1d.
// keeping faces_per_node_aug and cells_per_node_aug to reuse 3D SimpleMesh
// logic

Mesh_simple::Mesh_simple(const std::vector<double>& x,
                         const MPI_Comm& comm,
                         const JaliGeometry::GeometricModelPtr &gm,
                         const bool request_faces,
                         const bool request_edges,
                         const bool request_sides,
                         const bool request_wedges,
                         const bool request_corners,
                         const int num_tiles_ini,
                         const int num_ghost_layers_tile,
                         const int num_ghost_layers_distmesh,
                         const bool boundary_ghosts_requested,
                         const Partitioner_type partitioner,
                         const JaliGeometry::Geom_type geom_type) :
  nx_(x.size()-1), ny_(-3), nz_(-3),
  coordinates_(x),
  nodes_per_face_(1), faces_per_cell_(2), nodes_per_cell_(2),
  faces_per_node_aug_(2), cells_per_node_aug_(3),
  Mesh(request_faces, request_edges, request_sides, request_wedges,
       request_corners, num_tiles_ini, num_ghost_layers_tile,
       num_ghost_layers_distmesh, boundary_ghosts_requested,
       partitioner, geom_type, comm) {
  set_space_dimension(1);
  set_manifold_dimension(1);

  assert(!boundary_ghosts_requested);  // Cannot yet make boundary ghosts

  clear_internals_1d_();
  update_internals_1d_();

  cache_extra_variables();

  if (Mesh::num_tiles_ini_)
    Mesh::build_tiles();
}

//--------------------------------------
// Constructor - Construct a new mesh from a subset of an existing mesh
//--------------------------------------

Mesh_simple::Mesh_simple(const std::shared_ptr<Mesh> inmesh,
                         const std::vector<std::string>& setnames,
                         const Entity_kind setkind,
                         const bool flatten,
                         const bool extrude,
                         const bool request_faces,
                         const bool request_edges,
                         const bool request_wedges,
                         const bool request_sides,
                         const bool request_corners,
                         const int num_tiles_ini,
                         const int num_ghost_layers_tile,
                         const int num_ghost_layers_distmesh,
                         const bool boundary_ghosts_requested,
                         const Partitioner_type partitioner,
                         const JaliGeometry::Geom_type geom_type) {
  Errors::Message mesg("Construction of new mesh from an existing mesh not yet"
                       " implemented in the Simple mesh framework\n");
  Exceptions::Jali_throw(mesg);
}

Mesh_simple::Mesh_simple(const Mesh& inmesh,
                         const std::vector<std::string>& setnames,
                         const Entity_kind setkind,
                         const bool flatten,
                         const bool extrude,
                         const bool request_faces,
                         const bool request_edges,
                         const bool request_sides,
                         const bool request_wedges,
                         const bool request_corners,
                         const int num_tiles_ini,
                         const int num_ghost_layers_tile,
                         const int num_ghost_layers_distmesh,
                         const bool boundary_ghosts_requested,
                         const Partitioner_type partitioner,
                         const JaliGeometry::Geom_type geom_type) {
  Errors::Message mesg("Construction of new mesh from an existing mesh not yet"
                       " implemented in the Simple mesh framework\n");
  Exceptions::Jali_throw(mesg);
}

Mesh_simple::Mesh_simple(const Mesh& inmesh,
                         const std::vector<int>& entity_id_list,
                         const Entity_kind entity_kind,
                         const bool flatten,
                         const bool extrude,
                         const bool request_faces,
                         const bool request_edges,
                         const bool request_sides,
                         const bool request_wedges,
                         const bool request_corners,
                         const int num_tiles_ini,
                         const int num_ghost_layers_tile,
                         const int num_ghost_layers_distmesh,
                         const bool boundary_ghosts_requested,
                         const Partitioner_type partitioner,
                         const JaliGeometry::Geom_type geom_type) {
  Errors::Message mesg("Construction of new mesh from an existing mesh not yet"
                       " implemented in the Simple mesh framework\n");
  Exceptions::Jali_throw(mesg);
}


Mesh_simple::~Mesh_simple() { }

void Mesh_simple::clear_internals_3d_() {
  coordinates_.resize(0);

  cell_to_face_.resize(0);
  cell_to_node_.resize(0);
  face_to_node_.resize(0);

  side_sets_.resize(0);
}

void Mesh_simple::clear_internals_1d_() {
  cell_to_face_.resize(0);
  cell_to_node_.resize(0);
  face_to_node_.resize(0);

  side_sets_.resize(0);
}

void Mesh_simple::update_internals_1d_() {
  num_cells_ = nx_;
  num_nodes_ = nx_+1;
  num_faces_ = num_nodes_;

  cell_to_face_.resize(faces_per_cell_*num_cells_);
  cell_to_face_dirs_.resize(faces_per_cell_*num_cells_);
  cell_to_node_.resize(nodes_per_cell_*num_cells_);
  face_to_node_.resize(nodes_per_face_*num_faces_);
  // 1 extra for num faces
  node_to_face_.resize(faces_per_node_aug_*num_nodes_);
  // 1 extra for num cells
  node_to_cell_.resize(cells_per_node_aug_*num_nodes_);
  face_to_cell_.resize(2*num_faces_);
  face_to_cell_.assign(2*num_faces_, -1);

  // loop over cells and initialize cell_to_node_
  for (int ix = 0; ix < nx_; ++ix) {
    int jstart = 0;
    int istart = nodes_per_cell_ * cell_index_(ix);
    int ncell = 0;

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
  for (int ix = 0; ix < nx_; ++ix) {
    int istart = faces_per_cell_ * cell_index_(ix);
    int jstart = 0;

    cell_to_face_[istart]   = yzface_index_(ix);
    cell_to_face_[istart+1] = yzface_index_(ix+1);

    cell_to_face_dirs_[istart] = -1;
    cell_to_face_dirs_[istart+1] = 1;

    jstart = 2*yzface_index_(ix+1);
    face_to_cell_[jstart+1] = cell_index_(ix);

    jstart = 2*yzface_index_(ix);
    face_to_cell_[jstart] = cell_index_(ix);
  }

  // loop over faces and initialize face_to_node_
  // the yz faces
  for (int ix = 0; ix <= nx_; ++ix) {
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

  Mesh::nodeids_owned_.resize(num_nodes_);
  for (int i = 0; i < num_nodes_; ++i)
    nodeids_owned_[i] = i;
  Mesh::nodeids_ghost_.resize(0);
  Mesh::nodeids_all_ = Mesh::nodeids_owned_;

  if (Mesh::edges_requested) {
    int num_edges = num_cells_+1;  // edges are same as faces and nodes in 1D
    Mesh::edgeids_owned_.resize(num_edges);
    for (int i = 0; i < num_edges; ++i)
      edgeids_owned_[i] = i;
    Mesh::edgeids_ghost_.resize(0);
    Mesh::edgeids_all_ = Mesh::edgeids_owned_;
  } else {
    edgeids_owned_.resize(0);
    edgeids_ghost_.resize(0);
    edgeids_all_.resize(0);
  }

  if (Mesh::faces_requested) {
    Mesh::faceids_owned_.resize(num_faces_);
    for (int i = 0; i < num_faces_; ++i)
      faceids_owned_[i] = i;
    Mesh::faceids_ghost_.resize(0);
    Mesh::faceids_all_ = Mesh::faceids_owned_;
  } else {
    faceids_owned_.resize(0);
    faceids_ghost_.resize(0);
    faceids_all_.resize(0);
  }

  Mesh::cellids_owned_.resize(num_cells_);
  for (int i = 0; i < num_cells_; ++i)
    cellids_owned_[i] = i;
  Mesh::cellids_ghost_.resize(0);
  Mesh::cellids_all_ = Mesh::cellids_owned_;

}


void Mesh_simple::update_internals_3d_() {
  num_cells_ = nx_ * ny_ * nz_;
  num_nodes_ = (nx_+1)*(ny_+1)*(nz_+1);
  num_faces_ = (nx_+1)*(ny_)*(nz_) + (nx_)*(ny_+1)*(nz_) + (nx_)*(ny_)*(nz_+1);

  coordinates_.resize(3*num_nodes_);

  double hx = (x1_ - x0_)/nx_;
  double hy = (y1_ - y0_)/ny_;
  double hz = (z1_ - z0_)/nz_;

  for (int iz = 0; iz <= nz_; ++iz)
    for (int iy = 0; iy <= ny_; ++iy)
      for (int ix = 0; ix <= nx_; ++ix) {
        int istart = 3*node_index_(ix, iy, iz);

        coordinates_[ istart ]     = x0_ + ix*hx;
        coordinates_[ istart + 1 ] = y0_ + iy*hy;
        coordinates_[ istart + 2 ] = z0_ + iz*hz;
      }

  cell_to_face_.resize(faces_per_cell_*num_cells_);
  cell_to_face_dirs_.resize(faces_per_cell_*num_cells_);
  cell_to_node_.resize(nodes_per_cell_*num_cells_);
  face_to_node_.resize(nodes_per_face_*num_faces_);
  // 1 extra for num faces
  node_to_face_.resize(faces_per_node_aug_*num_nodes_);
  // 1 extra for num cells
  node_to_cell_.resize(cells_per_node_aug_*num_nodes_);
  face_to_cell_.resize(2*num_faces_);
  face_to_cell_.assign(2*num_faces_, -1);

  // loop over cells and initialize cell_to_node_
  for (int iz = 0; iz < nz_; ++iz)
    for (int iy = 0; iy < ny_; ++iy)
      for (int ix = 0; ix < nx_; ++ix) {
        int jstart = 0;
        int istart = nodes_per_cell_ * cell_index_(ix, iy, iz);
        int ncell = 0;

        cell_to_node_[istart]   = node_index_(ix, iy, iz);
        cell_to_node_[istart+1] = node_index_(ix+1, iy, iz);
        cell_to_node_[istart+2] = node_index_(ix+1, iy+1, iz);
        cell_to_node_[istart+3] = node_index_(ix, iy+1, iz);
        cell_to_node_[istart+4] = node_index_(ix, iy, iz+1);
        cell_to_node_[istart+5] = node_index_(ix+1, iy, iz+1);
        cell_to_node_[istart+6] = node_index_(ix+1, iy+1, iz+1);
        cell_to_node_[istart+7] = node_index_(ix, iy+1, iz+1);

        jstart = cells_per_node_aug_ * node_index_(ix, iy, iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart = cells_per_node_aug_ * node_index_(ix+1, iy, iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart = cells_per_node_aug_ * node_index_(ix+1, iy+1, iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart = cells_per_node_aug_ * node_index_(ix, iy+1, iz);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

        jstart = cells_per_node_aug_ * node_index_(ix, iy, iz+1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix, iy, iz);
        node_to_cell_[jstart]++;

        // 1 extra for num cells
        jstart = cells_per_node_aug_ * node_index_(ix+1, iy, iz+1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

         // 1 extra for num cells
        jstart = cells_per_node_aug_ * node_index_(ix+1, iy+1, iz+1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;

         // 1 extra for num cells
        jstart = cells_per_node_aug_ * node_index_(ix, iy+1, iz+1);
        ncell = node_to_cell_[jstart];
        node_to_cell_[jstart+1+ncell] = cell_index_(ix, iy, iz);
        (node_to_cell_[jstart])++;
      }


  // loop over cells and initialize cell_to_face_
  for (int iz = 0; iz < nz_; ++iz)
    for (int iy = 0; iy < ny_; ++iy)
      for (int ix = 0; ix < nx_; ++ix) {
        int istart = faces_per_cell_ * cell_index_(ix, iy, iz);
        int jstart = 0;

        cell_to_face_[istart]    = xzface_index_(ix, iy, iz);
        cell_to_face_[istart+1]  = yzface_index_(ix+1, iy, iz);
        cell_to_face_[istart+2]  = xzface_index_(ix, iy+1, iz);
        cell_to_face_[istart+3]  = yzface_index_(ix, iy, iz);
        cell_to_face_[istart+4]  = xyface_index_(ix, iy, iz);
        cell_to_face_[istart+5]  = xyface_index_(ix, iy, iz+1);

        cell_to_face_dirs_[istart]   = 1;
        cell_to_face_dirs_[istart+1] = 1;
        cell_to_face_dirs_[istart+2] = -1;
        cell_to_face_dirs_[istart+3] = -1;
        cell_to_face_dirs_[istart+4] = -1;
        cell_to_face_dirs_[istart+5] = 1;

        jstart = 2*xzface_index_(ix, iy, iz);
        face_to_cell_[jstart+1] = cell_index_(ix, iy, iz);

        jstart = 2*yzface_index_(ix+1, iy, iz);
        face_to_cell_[jstart+1] = cell_index_(ix, iy, iz);

        jstart = 2*xzface_index_(ix, iy+1, iz);
        face_to_cell_[jstart] = cell_index_(ix, iy, iz);

        jstart = 2*yzface_index_(ix, iy, iz);
        face_to_cell_[jstart] = cell_index_(ix, iy, iz);

        jstart = 2*xyface_index_(ix, iy, iz);
        face_to_cell_[jstart] = cell_index_(ix, iy, iz);

        jstart = 2*xyface_index_(ix, iy, iz+1);
        face_to_cell_[jstart+1] = cell_index_(ix, iy, iz);
      }


  // loop over faces and initialize face_to_node_
  // first we do the xy faces

  for (int iz = 0; iz <= nz_; ++iz)
    for (int iy = 0; iy < ny_; ++iy)
      for (int ix = 0; ix < nx_; ++ix) {
        int istart = nodes_per_face_ * xyface_index_(ix, iy, iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]   = node_index_(ix, iy, iz);
        face_to_node_[istart+1] = node_index_(ix+1, iy, iz);
        face_to_node_[istart+2] = node_index_(ix+1, iy+1, iz);
        face_to_node_[istart+3] = node_index_(ix, iy+1, iz);

        jstart = faces_per_node_aug_*node_index_(ix, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix+1, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix+1, iy+1, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix, iy+1, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;
      }
  // then we do the xz faces
  for (int iz = 0; iz < nz_; ++iz)
    for (int iy = 0; iy <= ny_; ++iy)
      for (int ix = 0; ix < nx_; ++ix) {
        int istart = nodes_per_face_ * xzface_index_(ix, iy, iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]   = node_index_(ix, iy, iz);
        face_to_node_[istart+1] = node_index_(ix+1, iy, iz);
        face_to_node_[istart+2] = node_index_(ix+1, iy, iz+1);
        face_to_node_[istart+3] = node_index_(ix, iy, iz+1);

        jstart = faces_per_node_aug_*node_index_(ix, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix+1, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix+1, iy, iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix, iy, iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;
      }
  // finally we do the yz faces
  for (int iz = 0; iz < nz_; ++iz)
    for (int iy = 0; iy < ny_; ++iy)
      for (int ix = 0; ix <= nx_; ++ix) {
        int istart = nodes_per_face_ * yzface_index_(ix, iy, iz);
        int jstart = 0;
        int nfaces = 0;

        face_to_node_[istart]   = node_index_(ix, iy, iz);
        face_to_node_[istart+1] = node_index_(ix, iy+1, iz);
        face_to_node_[istart+2] = node_index_(ix, iy+1, iz+1);
        face_to_node_[istart+3] = node_index_(ix, iy, iz+1);

        jstart = faces_per_node_aug_*node_index_(ix, iy, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix, iy+1, iz);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix, iy+1, iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;

        jstart = faces_per_node_aug_*node_index_(ix, iy, iz+1);
        nfaces = node_to_face_[jstart];
        node_to_face_[jstart+1+nfaces] = xyface_index_(ix, iy, iz);
        (node_to_face_[jstart])++;
      }

  // populate entity ids arrays in the base class so that iterators work

  int lid;
  std::vector<int>::iterator it;

  Mesh::nodeids_owned_.resize(num_nodes_);
  for (int i = 0; i < num_nodes_; ++i)
    nodeids_owned_[i] = i;
  Mesh::nodeids_ghost_.resize(0);
  Mesh::nodeids_all_ = Mesh::nodeids_owned_;

  // Simple mesh does not handle edges in 3D

  Mesh::edgeids_owned_.resize(0);
  Mesh::edgeids_ghost_.resize(0);
  Mesh::edgeids_all_.resize(0);

  if (Mesh::faces_requested) {
    Mesh::faceids_owned_.resize(num_faces_);
    for (int i = 0; i < num_faces_; ++i)
      faceids_owned_[i] = i;
    Mesh::faceids_ghost_.resize(0);
    Mesh::faceids_all_ = Mesh::faceids_owned_;
  } else {
    Mesh::faceids_owned_.resize(0);
    Mesh::faceids_ghost_.resize(0);
    Mesh::faceids_all_.resize(0);
  }

  Mesh::cellids_owned_.resize(num_cells_);
  for (int i = 0; i < num_cells_; ++i)
    cellids_owned_[i] = i;
  Mesh::cellids_ghost_.resize(0);
  Mesh::cellids_all_ = Mesh::cellids_owned_;

}


// Get cell type
Jali::Cell_type Mesh_simple::cell_get_type(const Jali::Entity_ID cellid) const {
  return Cell_type::HEX;
}

Entity_ID Mesh_simple::GID(const Jali::Entity_ID lid,
                           const Jali::Entity_kind kind) const {
  return lid;  // Its a serial code
}




void Mesh_simple::cell_get_faces_and_dirs_internal(const Jali::Entity_ID cellid,
                                                   Jali::Entity_ID_List
                                                   *faceids,
                                                   std::vector<dir_t> *cfacedirs,
                                                   const bool ordered) const {
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



void Mesh_simple::cell_get_nodes(Jali::Entity_ID cell,
                                 Jali::Entity_ID_List *nodeids) const {
  unsigned int offset = (unsigned int) nodes_per_cell_*cell;

  nodeids->clear();

  for (int i = 0; i < nodes_per_cell_; i++) {
    nodeids->push_back(*(cell_to_node_.begin()+offset));
    offset++;
  }
}


void Mesh_simple::face_get_nodes(Jali::Entity_ID face,
                                 Jali::Entity_ID_List *nodeids) const {
  unsigned int offset = (unsigned int) nodes_per_face_*face;

  nodeids->clear();

  for (int i = 0; i < nodes_per_face_; i++) {
    nodeids->push_back(*(face_to_node_.begin()+offset));
    offset++;
  }
}


// Cooordinate Getters
// -------------------

void Mesh_simple::node_get_coordinates(const Jali::Entity_ID local_node_id,
                                       JaliGeometry::Point *ncoords) const {
  unsigned int offset = (unsigned int) Mesh::space_dimension()*local_node_id;

  ncoords->set(Mesh::space_dimension(), &(coordinates_[offset]));
}
void Mesh_simple::node_get_coordinates(const Jali::Entity_ID local_node_id,
                                       std::array<double, 3> *ncoords) const {
  assert(Mesh::space_dimension() == 3);
  unsigned int offset = 3*local_node_id;
  (*ncoords)[0] = coordinates_[offset];
  (*ncoords)[1] = coordinates_[offset + 1];
  (*ncoords)[2] = coordinates_[offset + 2];
}
void Mesh_simple::node_get_coordinates(const Jali::Entity_ID local_node_id,
                                       double *ncoords) const {
  assert(Mesh::space_dimension() == 1);
  unsigned int offset = 1*local_node_id;
  *ncoords = coordinates_[offset];
}


void Mesh_simple::face_get_coordinates(Jali::Entity_ID local_face_id,
                                       std::vector<JaliGeometry::Point>
                                       *fcoords) const {
  Entity_ID_List node_indices(nodes_per_face_);

  face_get_nodes(local_face_id, &node_indices);

  fcoords->clear();

  JaliGeometry::Point node_coords(Mesh::space_dimension());
  for (std::vector<Entity_ID>::iterator it = node_indices.begin();
       it != node_indices.end(); ++it) {
    node_get_coordinates(*it, &node_coords);
    fcoords->push_back(node_coords);
  }
}

void Mesh_simple::cell_get_coordinates(Jali::Entity_ID local_cell_id,
                                       std::vector<JaliGeometry::Point>
                                       *ccoords) const {
  std::vector<Entity_ID> node_indices(nodes_per_cell_);

  cell_get_nodes(local_cell_id, &node_indices);

  ccoords->clear();

  JaliGeometry::Point node_coords(Mesh::space_dimension());
  for (std::vector<Entity_ID>::iterator it = node_indices.begin();
       it != node_indices.end(); ++it) {
    node_get_coordinates(*it, &node_coords);
    ccoords->push_back(node_coords);
  }
}


void Mesh_simple::node_set_coordinates(const Jali::Entity_ID local_node_id,
                                      const double *ncoord) {
  int spdim = Mesh::space_dimension();
  unsigned int offset = (unsigned int) spdim*local_node_id;


  ASSERT(ncoord != NULL);

  std::vector<double>::iterator destination_begin = coordinates_.begin()
      + offset;
  for (int i = 0; i < spdim; i++) {
    *destination_begin = ncoord[i];
    destination_begin++;
  }
}

void Mesh_simple::node_set_coordinates(const Jali::Entity_ID local_node_id,
                                       const JaliGeometry::Point ncoord) {
  int spdim = Mesh::space_dimension();
  unsigned int offset = (unsigned int) spdim*local_node_id;

  std::vector<double>::iterator destination_begin = coordinates_.begin()
      + offset;
  for (int i = 0; i < spdim; i++) {
    *destination_begin = ncoord[i];
    destination_begin++;
  }
}


void Mesh_simple::node_get_cells(const Jali::Entity_ID nodeid,
                                 const Jali::Entity_type ptype,
                                 Jali::Entity_ID_List *cellids) const {
  unsigned int offset = (unsigned int) cells_per_node_aug_*nodeid;
  unsigned int ncells = node_to_cell_[offset];

  cellids->clear();

  for (int i = 0; i < ncells; i++)
    cellids->push_back(node_to_cell_[offset+i+1]);
}


// Faces of type 'ptype' connected to a node
void Mesh_simple::node_get_faces(const Jali::Entity_ID nodeid,
                                 const Jali::Entity_type ptype,
                                 Jali::Entity_ID_List *faceids) const {
  unsigned int offset = (unsigned int) faces_per_node_aug_*nodeid;
  unsigned int nfaces = node_to_face_[offset];

  faceids->clear();

  for (int i = 0; i < nfaces; i++)
    faceids->push_back(node_to_face_[offset+i+1]);
}


// Get faces of ptype of a particular cell that are connected to the
// given node

void Mesh_simple::node_get_cell_faces(const Jali::Entity_ID nodeid,
                                      const Jali::Entity_ID cellid,
                                      const Jali::Entity_type ptype,
                                      Jali::Entity_ID_List *faceids) const {
  unsigned int offset = (unsigned int) faces_per_cell_*cellid;

  faceids->clear();

  for (int i = 0; i < faces_per_cell_; i++) {
    Entity_ID cellfaceid = cell_to_face_[offset+i];

    unsigned int offset2 = (unsigned int) nodes_per_face_*cellfaceid;

    Jali::Entity_ID_List fnodes;
    face_get_nodes(cellfaceid, &fnodes);

    for (int j = 0; j < nodes_per_face_; j++) {
      if (face_to_node_[offset2+j] == nodeid) {
        faceids->push_back(cellfaceid);
        break;
      }
    }
  }
}

// Cells connected to a face
void Mesh_simple::face_get_cells_internal(const Jali::Entity_ID faceid,
                                          const Jali::Entity_type ptype,
                                          Jali::Entity_ID_List *cellids) const {
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
                                          const Jali::Entity_type ptype,
                                          Jali::Entity_ID_List
                                          *fadj_cellids) const {
  unsigned int offset = (unsigned int) faces_per_cell_*cellid;

  fadj_cellids->clear();

  for (int i = 0; i < faces_per_cell_; i++) {
    Entity_ID faceid = cell_to_face_[offset];

    unsigned int foffset = (unsigned int) 2*faceid;

    int adjcell0 = face_to_cell_[foffset];
    if (adjcell0 != -1 && adjcell0 != cellid) {
      fadj_cellids->push_back(adjcell0);
    } else {
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
                                          const Jali::Entity_type ptype,
                                          Jali::Entity_ID_List
                                          *nadj_cellids) const {
  unsigned int offset = (unsigned int) nodes_per_cell_*cellid;

  nadj_cellids->clear();

  for (int i = 0; i < nodes_per_cell_; i++) {
    Entity_ID nodeid = cell_to_node_[offset+i];

    unsigned int offset2 = (unsigned int) cells_per_node_aug_*nodeid;
    unsigned int ncell = node_to_cell_[offset2];

    for (int j = 0; j < ncell; j++) {
      Entity_ID nodecell = node_to_cell_[offset2+j+1];
      if (nodecell == cellid) continue;

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

void
Mesh_simple::get_labeled_set_entities(const JaliGeometry::LabeledSetRegionPtr r,
                                      const Entity_kind kind,
                                      Entity_ID_List *owned_entities,
                                      Entity_ID_List *ghost_entities) const {
  // simple mesh does not read from any mesh file so it cannot have
  // pre-existing labeled sets

  owned_entities->clear();
  ghost_entities->clear();
}
                              





}  // close namespace Jali
