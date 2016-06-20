//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#include "Mesh.hh"

#ifdef Jali_HAVE_METIS
#include "metis.h"
#endif

#ifdef Jali_HAVE_ZOLTAN
#include "zoltan.h"
#endif

#include <math.h>
#include <cmath>
#include <vector>

#include "MeshTile.hh"

#include "Geometry.hh"
#include "dbc.hh"
#include "errors.hh"
#include "LabeledSetRegion.hh"


namespace Jali {

// Gather and cache type info for cells, faces, edges and nodes.
// The parallel type for other entities is derived

void Mesh::cache_type_info() const {
  cell_type.resize(num_cells());
  for (auto const& c : cells<Entity_type::PARALLEL_OWNED>())
    cell_type[c] = Entity_type::PARALLEL_OWNED;
  for (auto const& c : cells<Entity_type::PARALLEL_GHOST>())
    cell_type[c] = Entity_type::PARALLEL_GHOST;
  for (auto const& c : cells<Entity_type::BOUNDARY_GHOST>())
    cell_type[c] = Entity_type::BOUNDARY_GHOST;

  if (faces_requested) {
    face_type.resize(num_faces());
    for (auto const& f : faces<Entity_type::PARALLEL_OWNED>())
      face_type[f] = Entity_type::PARALLEL_OWNED;
    for (auto const& f : faces<Entity_type::PARALLEL_GHOST>())
      face_type[f] = Entity_type::PARALLEL_GHOST;
  }

  if (edges_requested) {
    edge_type.resize(num_edges());
    for (auto const& e : edges<Entity_type::PARALLEL_OWNED>())
      edge_type[e] = Entity_type::PARALLEL_OWNED;
    for (auto const& e : edges<Entity_type::PARALLEL_GHOST>())
      edge_type[e] = Entity_type::PARALLEL_GHOST;
  }

  node_type.resize(num_nodes());
  for (auto const& n : nodes<Entity_type::PARALLEL_OWNED>())
    node_type[n] = Entity_type::PARALLEL_OWNED;
  for (auto const& n : nodes<Entity_type::PARALLEL_GHOST>())
    node_type[n] = Entity_type::PARALLEL_GHOST;

  type_info_cached = true;
}

// Gather and cache cell to face connectivity info.
//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh class for further explanation

void Mesh::cache_cell2face_info() const {
  int ncells = num_cells<Entity_type::ALL>();
  cell_face_ids.resize(ncells);
  cell_face_dirs.resize(ncells);

  for (int c = 0; c < ncells; c++)
    cell_get_faces_and_dirs_internal(c, &(cell_face_ids[c]),
                                     &(cell_face_dirs[c]), false);

  cell2face_info_cached = true;
}


// Gather and cache face to cell connectivity info.
//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh class for further explanation

void Mesh::cache_face2cell_info() const {
  int nfaces = num_faces<Entity_type::ALL>();
  face_cell_ids.resize(nfaces);

  std::vector<Entity_ID> fcells;

  for (int f = 0; f < nfaces; f++) {
    face_get_cells_internal(f, Entity_type::ALL, &fcells);

    face_cell_ids[f].resize(2);

    for (int i = 0; i < fcells.size(); ++i)
      face_cell_ids[f][i] = fcells[i];
    for (int i = fcells.size(); i < 2; i++)
      face_cell_ids[f][i] = -1;
  }

  face2cell_info_cached = true;
}

// Gather and cache face to edge connectivity info.
//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh class for further explanation

void Mesh::cache_face2edge_info() const {
  int nfaces = num_faces<Entity_type::ALL>();
  face_edge_ids.resize(nfaces);
  face_edge_dirs.resize(nfaces);

  for (int f = 0; f < nfaces; ++f)
    face_get_edges_and_dirs_internal(f, &(face_edge_ids[f]),
                                     &(face_edge_dirs[f]), true);

  face2edge_info_cached = true;
}

// Gather and cache cell to edge connectivity info.
//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh class for further explanation

void Mesh::cache_cell2edge_info() const {
  int ncells = num_cells<Entity_type::ALL>();
  cell_edge_ids.resize(ncells);

  if (spacedim == 1) {
    for (auto const& c : cells())
      cell_get_nodes(c, &(cell_edge_ids[c]));   // edges are same as nodes
  } else if (spacedim == 2) {
    cell_2D_edge_dirs.resize(ncells);
    for (auto const& c : cells())
      cell_2D_get_edges_and_dirs_internal(c, &(cell_edge_ids[c]),
                                          &(cell_2D_edge_dirs[c]));
  } else if (spacedim == 3) {
    for (auto const& c : cells())
      cell_get_edges_internal(c, &(cell_edge_ids[c]));
  }

  cell2edge_info_cached = true;
}

// Gather and cache edge to node connectivity info

//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh class for further explanation

void Mesh::cache_edge2node_info() const {
  int nedges = num_edges<Entity_type::ALL>();
  edge_node_ids.resize(nedges);

  if (spacedim == 1) {
    // nodes and edges are identical in 1D - no need to go the derived
    // class for info
    for (auto const& e : edges())
      edge_node_ids[e][0] = edge_node_ids[e][1] = e;
  } else {
    for (auto const& e : edges())
      edge_get_nodes_internal(e, &(edge_node_ids[e][0]),
                              &(edge_node_ids[e][1]));
  }

  edge2node_info_cached = true;
}

// Gather and cache side information

void Mesh::cache_side_info() const {
  int ncells_owned = num_cells<Entity_type::PARALLEL_OWNED>();
  int ncells_ghost = num_cells<Entity_type::PARALLEL_GHOST>();
  int ncells_bndry_ghost = num_cells<Entity_type::BOUNDARY_GHOST>();
  int ncells = ncells_owned + ncells_ghost + ncells_bndry_ghost;

  cell_side_ids.resize(ncells);

  int nnodes_owned = num_nodes<Entity_type::PARALLEL_OWNED>();
  int nnodes_ghost = num_nodes<Entity_type::PARALLEL_GHOST>();
  int nnodes = nnodes_owned + nnodes_ghost;

  int num_sides_all = 0;
  int num_sides_owned = 0;
  int num_sides_ghost = 0;
  int num_sides_bndry_ghost = 0;

  if (celldim == 1) {  // in 1D there are always 2 sides per cell
    num_sides_all = 2*ncells;
    num_sides_owned = 2*ncells_owned;
    num_sides_ghost = 2*ncells_ghost;
    num_sides_bndry_ghost = 2*ncells_bndry_ghost;

    for (auto const & c : cells())
      cell_side_ids[c].reserve(2);
  } else {
    for (auto const & c : cells()) {
      std::vector<Entity_ID> cfaces;
      cell_get_faces(c, &cfaces);
      
      int numsides_in_cell = 0;
      for (auto const & f : cfaces) {
        std::vector<Entity_ID> fedges;
        std::vector<dir_t> fedirs;
        face_get_edges_and_dirs(f, &fedges, &fedirs);
        
        int nfedges = fedges.size();
        num_sides_all += nfedges;  // In 2D there is 1 edge per face
        numsides_in_cell += nfedges;
        
        if (cell_type[c] == Entity_type::PARALLEL_OWNED)
          num_sides_owned += nfedges;
        else if (cell_type[c] == Entity_type::PARALLEL_GHOST)
          num_sides_ghost += nfedges;
        else if (cell_type[c] == Entity_type::BOUNDARY_GHOST)
          num_sides_bndry_ghost += nfedges;
      }

      cell_side_ids[c].reserve(numsides_in_cell);
    }
  }

  sideids_owned_.resize(num_sides_owned);
  sideids_ghost_.resize(num_sides_ghost);
  sideids_boundary_ghost_.resize(num_sides_bndry_ghost);
  sideids_all_.resize(num_sides_all);
  side_edge_id.resize(num_sides_all, -1);
  side_edge_use.resize(num_sides_all, true);
  side_face_id.resize(num_sides_all, -1);
  side_cell_id.resize(num_sides_all, -1);
  side_opp_side_id.resize(num_sides_all, -1);
  //  side_node_ids.resize(num_sides_all, {-1, -1});  // intel 15.0.3 does not like this
  side_node_ids.resize(num_sides_all);

  if (celldim == 1) {
    int iall = 0, iown = 0, ighost = 0, ibndry = 0;
    for (auto const& c : cells()) {
      // always 2 sides per cell
      int sideid = 2*c;
      
      Entity_ID_List nodeids;
      cell_get_nodes(c, &nodeids);
      
      cell_side_ids[c].push_back(sideid);
      cell_side_ids[c].push_back(sideid+1);
      sideids_all_[iall++] = sideid;
      sideids_all_[iall++] = sideid+1;
      if (cell_type[c] == Entity_type::PARALLEL_OWNED) {
        sideids_owned_[iown++] = sideid;
        sideids_owned_[iown++] = sideid+1;
      } else if (cell_type[c] == Entity_type::PARALLEL_GHOST) {  // not yet ..
        sideids_ghost_[ighost++] = sideid;
        sideids_ghost_[ighost++] = sideid+1;
      } else if (cell_type[c] == Entity_type::BOUNDARY_GHOST) {  // not yet ..
        sideids_boundary_ghost_[ibndry++] = sideid;
        sideids_boundary_ghost_[ibndry++] = sideid+1;
      }

      // Sides are degenerate
      side_node_ids[sideid][0] =  nodeids[0];
      side_node_ids[sideid][1] =  nodeids[0];
      side_node_ids[sideid+1][0] = nodeids[1];
      side_node_ids[sideid+1][1] = nodeids[1];

      // edges are the same as NODES - CHANGED FROM PREVIOUS
      // IMPLEMENTATION - VERIFY THATS OK?
      side_edge_id[sideid  ] = nodeids[0];
      side_edge_id[sideid+1] = nodeids[1];
      // faces are the same as nodes
      side_face_id[sideid  ] = nodeids[0];
      side_face_id[sideid+1] = nodeids[1];
      side_cell_id[sideid  ] = c;
      side_cell_id[sideid+1] = c;
      // across face boundaries
      side_opp_side_id[sideid  ] = (sideid == 0) ? -1 : sideid-1;
      side_opp_side_id[sideid+1] = (sideid+1 == num_sides_all-1) ?
                                    -1 : sideid+2;
    }
  } else {
    int nedges = num_edges<Entity_type::ALL>();
    std::vector<std::vector<Entity_ID>> sides_of_edge(nedges);  // Temp. var.

    int sideid = 0;
    int iall = 0, iown = 0, ighost = 0, ibndry = 0;
    for (auto const & c : cells()) {
      std::vector<Entity_ID> cfaces;
      std::vector<dir_t> cfdirs;
      cell_get_faces_and_dirs(c, &cfaces, &cfdirs);
      
      Entity_ID_List::iterator itf = cfaces.begin();
      std::vector<dir_t>::iterator itfd = cfdirs.begin();
      while (itf != cfaces.end()) {
        Entity_ID f = *itf;
        int fdir = *itfd;  // -1/1
        
        Entity_ID_List fedges;
        std::vector<dir_t> fedirs;
        face_get_edges_and_dirs(f, &fedges, &fedirs);
        
        Entity_ID_List::iterator ite = fedges.begin();
        std::vector<dir_t>::iterator ited = fedirs.begin();
        while (ite != fedges.end()) {
          Entity_ID e = *ite;
          int edir = *ited;  // -1/1
          
          Entity_ID enodes[2];
          edge_get_nodes(e, &(enodes[0]), &(enodes[1]));

          if (celldim == 2) {  // 2D
            side_node_ids[sideid][0] = (fdir == 1) ? enodes[0] : enodes[1];
            side_node_ids[sideid][1] = (fdir == 1) ? enodes[1] : enodes[0];
            side_edge_use[sideid] = (fdir == 1) ? true : false;
          } else {  // 3D
            
            // In 3D, if the cell uses the face in the +ve dir (normal
            // pointing out) and the face is using the edge in the +ve
            // dir, then the side will use the edge in the -ve dir,
            // because we want the facet formed by side point 0, side
            // point 1 and face center for the side to have a normal
            // pointing into the cell
            
            int fdir01 = (fdir+1)/2;  // convert from -1/1 to 0/1
            int edir01 = (edir+1)/2;
            int dir = fdir01^edir01;  // gives 0 if both are same
            side_node_ids[sideid][0] = dir ? enodes[0] : enodes[1];
            side_node_ids[sideid][1] = dir ? enodes[1] : enodes[0];
            side_edge_use[sideid] = dir ? true : false;

          }
            
          side_edge_id[sideid] = e;
          side_face_id[sideid] = f;
          side_cell_id[sideid] = c;
          cell_side_ids[c].push_back(sideid);
          
          sideids_all_[iall++] = sideid;
          if (cell_type[c] == Entity_type::PARALLEL_OWNED)
            sideids_owned_[iown++] = sideid;
          else if (cell_type[c] == Entity_type::BOUNDARY_GHOST)
            sideids_boundary_ghost_[ibndry++] = sideid;
          else
            sideids_ghost_[ighost++] = sideid;
          
          // See if any of the other sides attached to the edge
          // shares the same edge and face but is in the
          // adjacent cell. This is called the opposite side
          
          for (auto const& s2 : sides_of_edge[e]) {
            if (side_edge_id[s2] == e && side_face_id[s2] == f &&
                side_cell_id[s2] != c) {
              side_opp_side_id[sideid] = s2;
              side_opp_side_id[s2] = sideid;
              break;
            }
          }
          sides_of_edge[e].push_back(sideid);

          sideid++;

          ++ite;
          ++ited;
        }  // while (ite != fedges.end())
        
        ++itf;
        ++itfd;
      }  // while (itf != cfaces.end())
    }  // for (c : cells())
  }  // if (celldim)

  side_info_cached = true;
}  // cache_side_info


// Gather and cache wedge information

void Mesh::cache_wedge_info() const {
  int ncells_owned = num_cells<Entity_type::PARALLEL_OWNED>();
  int ncells_ghost = num_cells<Entity_type::PARALLEL_GHOST>();
  int ncells_boundary_ghost = num_cells<Entity_type::BOUNDARY_GHOST>();
  int ncells = ncells_owned + ncells_ghost + ncells_boundary_ghost;

  int nnodes_owned = num_nodes<Entity_type::PARALLEL_OWNED>();
  int nnodes_ghost = num_nodes<Entity_type::PARALLEL_GHOST>();
  int nnodes = nnodes_owned + nnodes_ghost;

  int nsides_owned = num_sides<Entity_type::PARALLEL_OWNED>();
  int nsides_ghost = num_sides<Entity_type::PARALLEL_GHOST>();
  int nsides_boundary_ghost = num_sides<Entity_type::BOUNDARY_GHOST>();
  int nsides_all = nsides_owned + nsides_ghost + nsides_boundary_ghost;

  int num_wedges_all = 2*nsides_all;
  int num_wedges_owned = 2*nsides_owned;
  int num_wedges_ghost = 2*nsides_ghost;
  int num_wedges_boundary_ghost = 2*nsides_boundary_ghost;

  wedgeids_owned_.resize(num_wedges_owned);
  wedgeids_ghost_.resize(num_wedges_ghost);
  wedgeids_boundary_ghost_.resize(num_wedges_boundary_ghost);
  wedge_corner_id.resize(num_wedges_all, -1);  // filled when building corners

  
  int iown = 0, ighost = 0, ibndry = 0;
  for (auto const& s : sides()) {
    int wedgeid0 = 2*s;
    int wedgeid1 = 2*s + 1;
    int c = side_cell_id[s];
    if (cell_type[c] == Entity_type::PARALLEL_OWNED) {
      wedgeids_owned_[iown++] = wedgeid0;
      wedgeids_owned_[iown++] = wedgeid1;
    } else if (cell_type[c] == Entity_type::PARALLEL_GHOST) {
      wedgeids_ghost_[ighost++] = wedgeid0;
      wedgeids_ghost_[ighost++] = wedgeid1;
    } else if (cell_type[c] == Entity_type::BOUNDARY_GHOST) {
      wedgeids_boundary_ghost_[ibndry++] = wedgeid0;
      wedgeids_boundary_ghost_[ibndry++] = wedgeid1;
    }
  }

  wedgeids_all_.reserve(num_wedges_all);
  wedgeids_all_ = wedgeids_owned_;  // copy
  wedgeids_all_.insert(wedgeids_all_.end(), wedgeids_ghost_.begin(),
              wedgeids_ghost_.end());
  wedgeids_all_.insert(wedgeids_all_.end(), wedgeids_boundary_ghost_.begin(),
              wedgeids_boundary_ghost_.end());

  wedge_info_cached = true;
}  // cache_wedge_info


void Mesh::cache_corner_info() const {
  int ncells_owned = num_cells<Entity_type::PARALLEL_OWNED>();
  int ncells_ghost = num_cells<Entity_type::PARALLEL_GHOST>();
  int ncells_boundary_ghost = num_cells<Entity_type::BOUNDARY_GHOST>();
  int ncells = ncells_owned + ncells_ghost + ncells_boundary_ghost;

  int nnodes_owned = num_nodes<Entity_type::PARALLEL_OWNED>();
  int nnodes_ghost = num_nodes<Entity_type::PARALLEL_GHOST>();
  int nnodes = nnodes_owned + nnodes_ghost;

  cell_corner_ids.resize(ncells);
  node_corner_ids.resize(nnodes);

  int num_corners_all = 0;
  int num_corners_owned = 0;
  int num_corners_ghost = 0;
  int num_corners_boundary_ghost = 0;

  for (auto const& c : cells()) {
    std::vector<Entity_ID> cnodes;
    cell_get_nodes(c, &cnodes);
    cell_corner_ids[c].reserve(cnodes.size());

    num_corners_all += cnodes.size();  // as many corners as nodes in cell
    if (cell_type[c] == Entity_type::PARALLEL_OWNED)
      num_corners_owned += cnodes.size();
    else if (cell_type[c] == Entity_type::PARALLEL_GHOST)
      num_corners_ghost += cnodes.size();
    else if (cell_type[c] == Entity_type::BOUNDARY_GHOST)
      num_corners_boundary_ghost += cnodes.size();
  }

  cornerids_owned_.resize(num_corners_owned);
  cornerids_ghost_.resize(num_corners_ghost);
  cornerids_boundary_ghost_.resize(num_corners_boundary_ghost);
  corner_wedge_ids.resize(num_corners_all);

  int cornerid = 0;
  int iown = 0, ighost = 0, ibndry = 0;
  for (auto const& c : cells()) {
    std::vector<Entity_ID> cnodes;
    cell_get_nodes(c, &cnodes);

    std::vector<Entity_ID> cwedges;
    cell_get_wedges(c, &cwedges);

    for (auto const& n : cnodes) {
      cell_corner_ids[c].push_back(cornerid);
      node_corner_ids[n].push_back(cornerid);

      if (cell_type[c] == Entity_type::PARALLEL_OWNED)
        cornerids_owned_[iown++] = cornerid;
      else if (cell_type[c] == Entity_type::PARALLEL_GHOST)
        cornerids_ghost_[ighost++] = cornerid;
      else if (cell_type[c] == Entity_type::BOUNDARY_GHOST)
        cornerids_boundary_ghost_[ibndry++] = cornerid;

      for (auto const& w : cwedges) {
        Entity_ID n2 = wedge_get_node(w);
        if (n == n2) {
          corner_wedge_ids[cornerid].push_back(w);
          wedge_corner_id[w] = cornerid;
        }
      }  // for (w : cwedges)

      ++cornerid;
    }  // for (n : cnodes)
  }  // for (c : cells())

  cornerids_all_.reserve(num_corners_all);
  cornerids_all_ = cornerids_owned_;  // list copy
  cornerids_all_.insert(cornerids_all_.end(), cornerids_ghost_.begin(),
                        cornerids_ghost_.end());
  cornerids_all_.insert(cornerids_all_.end(), cornerids_boundary_ghost_.begin(),
                        cornerids_boundary_ghost_.end());

  corner_info_cached = true;
}  // cache_corner_info

void Mesh::update_geometric_quantities() {
  if (faces_requested) compute_face_geometric_quantities();
  if (edges_requested) compute_edge_geometric_quantities();
  compute_cell_geometric_quantities();
  if (sides_requested || wedges_requested) compute_side_geometric_quantities();
  if (corners_requested) compute_corner_geometric_quantities();
}


void Mesh::cache_extra_variables() {
  // Should be before side, wedge and corner info is processed
  cache_type_info();
    
  if (faces_requested) {
    cache_cell2face_info();
    cache_face2cell_info();
  }

  if (edges_requested) {
    cache_face2edge_info();
    cache_cell2edge_info();
    cache_edge2node_info();
  }

  if (sides_requested)
    cache_side_info();
  if (wedges_requested) {  // keep this order
    if (!side_info_cached) cache_side_info();
    cache_wedge_info();
  }
  if (corners_requested) {  // Keep this order
    if (!side_info_cached) cache_side_info();
    if (!wedge_info_cached) cache_wedge_info();
    cache_corner_info();
  }

  update_geometric_quantities();
}


// Partition the mesh on this compute node into submeshes or tiles

void Mesh::build_tiles() {

  std::vector<std::vector<int>> partitions;
  partitions.resize(num_tiles_ini_);

  std::cerr << "Calling partitioner " << partitioner_pref_ << "\n";
  get_partitioning(num_tiles_ini_, partitioner_pref_, &partitions);

  // Make the tiles - make mesh tile will also call a routine of this mesh
  // to add the tile to the mesh's list of tiles

  for (int i = 0; i < num_tiles_ini_; ++i)
    make_meshtile(*this, partitions[i], num_ghost_layers_tile_, faces_requested,
                  edges_requested, sides_requested, wedges_requested,
                  corners_requested);
}



// Initialize some arrays for storing master tile ID of entity. On
// this tile, the entity is OWNED. It can be a GHOST on any number of
// other entities

void Mesh::init_tiles() {
  tiles_initialized_ = true;
  node_master_tile_ID_.resize(num_nodes(), -1);
  if (edges_requested) edge_master_tile_ID_.resize(num_edges(), -1);
  if (faces_requested) face_master_tile_ID_.resize(num_faces(), -1);
  cell_master_tile_ID_.resize(num_cells(), -1);
}


// Add one tile to the mesh

void Mesh::add_tile(std::shared_ptr<MeshTile> const tile2add) {
  meshtiles.emplace_back(tile2add);
}

  
Entity_ID Mesh::entity_get_parent(const Entity_kind kind,
                                  const Entity_ID entid) const {
  Errors::Message mesg("Parent/daughter entities not"
                       " enabled in this framework.");
  Exceptions::Jali_throw(mesg);
}


unsigned int Mesh::cell_get_num_faces(const Entity_ID cellid) const {
#if JALI_CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //
  assert(cell2face_info_cached);

  return cell_face_ids[cellid].size();

#else

  // Non-cached version

  Entity_ID_List cfaceids;
  std::vector<dir_t> cfacedirs;

  cell_get_faces_and_dirs_internal(cellid, &cfaceids, &cfacedirs, false);

  return cfaceids.size();

#endif

}


void Mesh::cell_get_faces_and_dirs(const Entity_ID cellid,
                                   Entity_ID_List *faceids,
                                   std::vector<dir_t> *face_dirs,
                                   const bool ordered) const {
#if JALI_CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //
  assert(cell2face_info_cached);

  if (ordered) {
    cell_get_faces_and_dirs_internal(cellid, faceids, face_dirs, ordered);
  } else {
    Entity_ID_List &cfaceids = cell_face_ids[cellid];

    *faceids = cfaceids;  // copy operation

    if (face_dirs) {
      std::vector<dir_t> &cfacedirs = cell_face_dirs[cellid];
      *face_dirs = cfacedirs;  // copy operation
    }
  }

#else

  //
  // Non-cached version
  //

  cell_get_faces_and_dirs_internal(cellid, faceids, face_dirs, ordered);

#endif

}


// Cells connected to a face - cache the results the first time it
// is called and then return the cached results subsequently

void Mesh::face_get_cells(const Entity_ID faceid, const Entity_type ptype,
                          Entity_ID_List *cellids) const {
#if JALI_CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //

  assert(face2cell_info_cached);


  cellids->clear();

  switch (ptype) {
  case Entity_type::ALL:
    for (int i = 0; i < 2; i++) {
      Entity_ID c = face_cell_ids[faceid][i];
      if (c != -1) cellids->push_back(c);
    }
    break;
  case Entity_type::PARALLEL_OWNED:
    for (int i = 0; i < 2; i++) {
      Entity_ID c = face_cell_ids[faceid][i];
      if (c != -1 && cell_type[c] == Entity_type::PARALLEL_OWNED)
        cellids->push_back(c);
    }
    break;
  case Entity_type::PARALLEL_GHOST:
    for (int i = 0; i < 2; i++) {
      Entity_ID c = face_cell_ids[faceid][i];
      if (c != -1 && cell_type[c] == Entity_type::PARALLEL_GHOST)
        cellids->push_back(c);
    }
    break;
  }

#else

  //
  // Non-cached version
  //

  Entity_ID_List fcells;

  face_get_cells_internal(faceid, Entity_type::ALL, &fcells);

  cellids->clear();

  switch (ptype) {
  case Entity_type::ALL:
    for (int i = 0; i < fcells.size(); i++)
      cellids->push_back(fcells[i]);
    break;
  case Entity_type::PARALLEL_OWNED:
    for (int i = 0; i < fcells.size(); i++)
      if (cell_type[fcells[i]] == Entity_type::PARALLEL_OWNED)
        cellids->push_back(fcells[i]);
    break;
  case Entity_type::PARALLEL_GHOST:
    for (int i = 0; i < fcells.size(); i++)
      if (cell_type[fcells[i]] == Entity_type::PARALLEL_GHOST)
        cellids->push_back(fcells[i]);
    break;
  }

#endif
}


void Mesh::face_get_edges_and_dirs(const Entity_ID faceid,
                                   Entity_ID_List *edgeids,
                                   std::vector<dir_t> *edge_dirs,
                                   const bool ordered) const {
#if JALI_CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //

  assert(face2edge_info_cached);

  *edgeids = face_edge_ids[faceid];  // copy operation

  if (edge_dirs) {
    std::vector<dir_t> &fedgedirs = face_edge_dirs[faceid];
    *edge_dirs = fedgedirs;  // copy operation
  }


#else

  //
  // Non-cached version
  //

  face_get_edges_and_dirs_internal(faceid, edgeids, edge_dirs, ordered);

#endif
}


// Get the local ID of a face edge in a cell edge list

void Mesh::face_to_cell_edge_map(const Entity_ID faceid,
                                 const Entity_ID cellid,
                                 std::vector<int> *map) const {
#if JALI_CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //

  assert(face2edge_info_cached && cell2edge_info_cached);

  map->resize(face_edge_ids[faceid].size());
  for (int f = 0; f < face_edge_ids[faceid].size(); ++f) {
    Entity_ID fedge = face_edge_ids[faceid][f];

    for (int c = 0; c < cell_edge_ids[cellid].size(); ++c) {
      if (fedge == cell_edge_ids[cellid][c]) {
        (*map)[f] = c;
        break;
      }
    }
  }

#else

  Entity_ID_List fedgeids, cedgeids;
  std::vector<dir_t> fedgedirs;

  face_get_edges_and_dirs(faceid, &fedgeids, &fedgedirs, true);
  cell_get_edges(cellid, &cedgeids);

  map->resize(fedgeids.size(), -1);
  for (int f = 0; f < fedgeids.size(); ++f) {
    Entity_ID fedge = fedgeids[f];

    for (int c = 0; c < cedgeids.size(); ++c) {
      if (fedge == cedgeids[c]) {
        (*map)[f] = c;
        break;
      }
    }
  }

#endif
}


void Mesh::cell_get_edges(const Entity_ID cellid,
                          Entity_ID_List *edgeids) const {
#if JALI_CACHE_VARS != 0

  //
  // Cached version - turn off for profiling
  //

  assert(cell2edge_info_cached);

  Entity_ID_List &cedgeids = cell_edge_ids[cellid];

  *edgeids = cell_edge_ids[cellid];  // copy operation

#else

  //
  // Non-cached version
  //

  cell_get_edges_internal(cellid, edgeids);

#endif
}  // Mesh::cell_get_edges


void Mesh::cell_2D_get_edges_and_dirs(const Entity_ID cellid,
                                      Entity_ID_List *edgeids,
                                      std::vector<dir_t> *edgedirs) const {
#if JALI_CACHE_VARS != 0

  //
  // Cached version - turn off for profiling
  //

  assert(cell2edge_info_cached);

  *edgeids = cell_edge_ids[cellid];  // copy operation
  *edgedirs = cell_2D_edge_dirs[cellid];

#else

  //
  // Non-cached version
  //

  cell_2D_get_edges_and_dirs_internal(cellid, edgeids, edgedirs);

#endif
}  // Mesh::cell_get_edges_and_dirs


void Mesh::cell_get_sides(const Entity_ID cellid,
                           Entity_ID_List *sideids) const {
  assert(sides_requested);
  assert(side_info_cached);

  int nsides = cell_side_ids[cellid].size();
  sideids->resize(nsides);
  std::copy(cell_side_ids[cellid].begin(), cell_side_ids[cellid].end(),
            sideids->begin());
}


void Mesh::cell_get_wedges(const Entity_ID cellid,
                           Entity_ID_List *wedgeids) const {
  assert(wedges_requested);
  assert(side_info_cached);

  std::vector<Entity_ID>& csides = cell_side_ids[cellid];
  int nsides = csides.size();
  int nwedges = 2*nsides;
  wedgeids->resize(nwedges);
  for (int i = 0; i < nsides; ++i) {
    (*wedgeids)[2*i] = 2*csides[i];
    (*wedgeids)[2*i + 1] = 2*csides[i] + 1;
  }
}


void Mesh::cell_get_corners(const Entity_ID cellid,
                            Entity_ID_List *cornerids) const {
  assert(corners_requested);
  assert(corner_info_cached);

  int ncorners = cell_corner_ids[cellid].size();
  cornerids->resize(ncorners);
  std::copy(cell_corner_ids[cellid].begin(), cell_corner_ids[cellid].end(),
            cornerids->begin());
}


Entity_ID Mesh::cell_get_corner_at_node(const Entity_ID cellid,
                                        const Entity_ID nodeid) const {
  assert(corners_requested);
  assert(corner_info_cached);

  Entity_ID_List::iterator itc = cell_corner_ids[cellid].begin();
  while (itc != cell_corner_ids[cellid].end()) {
    int cornerid = *itc;
    if (corner_get_node(cornerid) == nodeid)
      return cornerid;
    else
      ++itc;
  }
  return -1;   // shouldn't come here unless node does not belong to cell
}


void Mesh::node_get_wedges(const Entity_ID nodeid, Entity_type ptype,
                           Entity_ID_List *wedgeids) const {
  assert(wedges_requested);
  assert(wedge_info_cached && corner_info_cached);

  wedgeids->clear();
  for (auto const& cn : node_corner_ids[nodeid]) {
    Entity_ID_List cnwedges = corner_wedge_ids[cn];
    for (auto const& w : cnwedges) {
      Entity_ID s = static_cast<Entity_ID>(w/2);
      Entity_ID c = side_cell_id[s];
      if (ptype == Entity_type::ALL || cell_type[c] == ptype)
        wedgeids->push_back(w);
    }
  }
}


void Mesh::node_get_corners(const Entity_ID nodeid, Entity_type ptype,
                            Entity_ID_List *cornerids) const {
  assert(corners_requested);
  assert(corner_info_cached);

  switch (ptype) {
    case Entity_type::ALL:
      cornerids->resize(node_corner_ids[nodeid].size());
      std::copy(node_corner_ids[nodeid].begin(), node_corner_ids[nodeid].end(),
                cornerids->begin());
      break;
    default:
      cornerids->clear();
      for (auto const& cn : node_corner_ids[nodeid]) {
        Entity_ID w0 = corner_wedge_ids[cn][0];
        Entity_ID s = static_cast<Entity_ID>(w0/2);
        Entity_ID c = side_cell_id[s];
        if (cell_type[c] == ptype)
          cornerids->push_back(cn);
      }
      break;
  }
}


int Mesh::compute_cell_geometric_quantities() const {
  int ncells = num_cells<Entity_type::ALL>();

  cell_volumes.resize(ncells);
  cell_centroids.resize(ncells);
  
  std::vector<double> zerovec(spacedim, 0.0);
  for (int c = 0; c < ncells; c++) {
    if (cell_type[c] == Entity_type::BOUNDARY_GHOST) {
      cell_volumes[c] = 0.0;
      cell_centroids[c].set(spacedim, &(zerovec[0]));
    } else {
      double volume;
      JaliGeometry::Point centroid(spacedim);
      
      compute_cell_geometry(c, &volume, &centroid);
      
      cell_volumes[c] = volume;
      cell_centroids[c] = centroid;
    }
  }

  cell_geometry_precomputed = true;
  return 1;
}  // Mesh::compute_cell_geometric_quantities



int Mesh::compute_face_geometric_quantities() const {
  int nfaces = num_faces<Entity_type::ALL>();

  face_areas.resize(nfaces);
  face_centroids.resize(nfaces);
  face_normal0.resize(nfaces);
  face_normal1.resize(nfaces);

  for (int i = 0; i < nfaces; i++) {
    double area;
    JaliGeometry::Point centroid(spacedim), normal0(spacedim),
        normal1(spacedim);

    // normal0 and normal1 are outward normals of the face with
    // respect to the cell0 and cell1 of the face. The natural normal
    // of the face points out of cell0 and into cell1. If one of these
    // cells do not exist, then the normal is the null vector.

    compute_face_geometry(i, &area, &centroid, &normal0, &normal1);

    face_areas[i] = area;
    face_centroids[i] = centroid;
    face_normal0[i] = normal0;
    face_normal1[i] = normal1;
  }

  face_geometry_precomputed = true;
  return 1;
}  // Mesh::compute_face_geometric_quantities



int Mesh::compute_edge_geometric_quantities() const {
  int nedges = num_edges<Entity_type::ALL>();

  edge_vectors.resize(nedges);
  edge_lengths.resize(nedges);

  for (int i = 0; i < nedges; i++) {
    double length;
    JaliGeometry::Point evector(spacedim), ecenter;

    compute_edge_geometry(i, &length, &evector, &ecenter);

    edge_lengths[i] = length;
    edge_vectors[i] = evector;
  }

  edge_geometry_precomputed = true;
  return 1;
}  // Mesh::compute_edge_geometric_quantities


int Mesh::compute_side_geometric_quantities() const {
  side_volumes.resize(num_sides());

  // Cannot use resize for the facet normals because we cannot tell
  // the initialization operator what the dimensionality of the points
  // will be, unless we can write a statement like this
  //
  //  wedge_facet_normals0.resize(nwedges, Point(3));
  //
  // It is also unclear that this is good to do because it will cause
  // a copy operator to be triggered for each element (?)

  side_outward_facet_normal.reserve(num_sides());
  side_mid_facet_normal.reserve(num_sides());

  JaliGeometry::Point outward_facet_normal(spacedim);
  JaliGeometry::Point mid_facet_normal(spacedim);
  for (auto const & s : sides()) {
    if (entity_get_type(Entity_kind::SIDE, s) == Entity_type::BOUNDARY_GHOST) {
      side_volumes[s] = 0.0;
      outward_facet_normal.set(0.0);
      mid_facet_normal.set(0.0);
    } else {
      compute_side_geometry(s, &(side_volumes[s]),
                            &(outward_facet_normal),
                            &(mid_facet_normal));
    }
    side_outward_facet_normal.push_back(outward_facet_normal);
    side_mid_facet_normal.push_back(mid_facet_normal);
  }

  side_geometry_precomputed = true;
  return 1;
}

int Mesh::compute_corner_geometric_quantities() const {
  corner_volumes.resize(num_corners());
  for (auto const & cn : corners())
    if (entity_get_type(Entity_kind::CORNER, cn) == Entity_type::BOUNDARY_GHOST)
      corner_volumes[cn] = 0.0;
    else
      compute_corner_geometry(cn, &(corner_volumes[cn]));
  corner_geometry_precomputed = true;
  return 1;
}

int Mesh::compute_cell_geometry(const Entity_ID cellid, double *volume,
                                JaliGeometry::Point *centroid) const {
  if (celldim == 3) {

    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine

    // General polyhedra always need to have an explicit face
    // representation - special elements like hexes can get away
    // without (but we have yet to put in the code for the standard
    // node ordering and computation for these special elements)

    Entity_ID_List faces;
    std::vector<unsigned int> nfnodes;
    std::vector<dir_t> fdirs;
    std::vector<JaliGeometry::Point> ccoords, cfcoords, fcoords;

    cell_get_faces_and_dirs(cellid, &faces, &fdirs);

    int nf = faces.size();
    nfnodes.resize(nf);

    for (int j = 0; j < nf; j++) {

      face_get_coordinates(faces[j], &fcoords);
      nfnodes[j] = fcoords.size();

      if (fdirs[j] == 1) {
        for (int k = 0; k < nfnodes[j]; k++)
          cfcoords.push_back(fcoords[k]);
      } else {
        for (int k = nfnodes[j]-1; k >=0; k--)
          cfcoords.push_back(fcoords[k]);
      }
    }

    cell_get_coordinates(cellid, &ccoords);

    JaliGeometry::polyhed_get_vol_centroid(ccoords, nf, nfnodes, cfcoords,
                                           volume, centroid);
    return 1;
  } else if (celldim == 2) {
    std::vector<JaliGeometry::Point> ccoords;

    cell_get_coordinates(cellid, &ccoords);

    JaliGeometry::Point normal(spacedim);

    JaliGeometry::polygon_get_area_centroid_normal(ccoords, volume, centroid,
                                                   &normal);

    return 1;
  } else if (celldim == 1) {
    std::vector<JaliGeometry::Point> ccoords;

    cell_get_coordinates(cellid, &ccoords);

    JaliGeometry::segment_get_vol_centroid(ccoords, geomtype,
                                           volume, centroid);
    return 1;
  }

  return 0;
}  // Mesh::compute_cell_geometry


int Mesh::compute_face_geometry(const Entity_ID faceid, double *area,
                                JaliGeometry::Point *centroid,
                                JaliGeometry::Point *normal0,
                                JaliGeometry::Point *normal1) const {
  JaliGeometry::Point_List fcoords;

  (*normal0).set(0.0L);
  (*normal1).set(0.0L);

  if (celldim == 3) {

    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine

    face_get_coordinates(faceid, &fcoords);

    JaliGeometry::Point normal(3);
    JaliGeometry::polygon_get_area_centroid_normal(fcoords, area, centroid,
                                                   &normal);

    Entity_ID_List cellids;
    face_get_cells(faceid, Entity_type::ALL, &cellids);

    for (int i = 0; i < cellids.size(); i++) {
      Entity_ID_List cellfaceids;
      std::vector<dir_t> cellfacedirs;
      dir_t dir = 1;

      cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);

      bool found = false;
      for (int j = 0; j < cellfaceids.size(); j++) {
        if (cellfaceids[j] == faceid) {
          found = true;
          dir = cellfacedirs[j];
          break;
        }
      }

      ASSERT(found);

      if (dir == 1)
        *normal0 = normal;
      else
        *normal1 = -normal;
    }

    return 1;
  } else if (celldim == 2) {

    if (spacedim == 2) {   // 2D mesh

      face_get_coordinates(faceid, &fcoords);

      JaliGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);

      *centroid = 0.5*(fcoords[0]+fcoords[1]);

      JaliGeometry::Point normal(evec[1], -evec[0]);

      Entity_ID_List cellids;
      face_get_cells(faceid, Entity_type::ALL, &cellids);

      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        std::vector<dir_t> cellfacedirs;
        dir_t dir = 1;

        cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);

        bool found = false;
        for (int j = 0; j < cellfaceids.size(); j++) {
          if (cellfaceids[j] == faceid) {
            found = true;
            dir = cellfacedirs[j];
            break;
          }
        }

        ASSERT(found);

        if (dir == 1)
          *normal0 = normal;
        else
          *normal1 = -normal;
      }

      return 1;
    } else {  // Surface mesh - cells are 2D, coordinates are 3D
      // edge normals are ambiguous for surface mesh
      // So we won't compute them

      face_get_coordinates(faceid, &fcoords);

      JaliGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);

      *centroid = 0.5*(fcoords[0]+fcoords[1]);

      Entity_ID_List cellids;
      face_get_cells(faceid, Entity_type::ALL, &cellids);

      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        std::vector<dir_t> cellfacedirs;
        dir_t dir = 1;

        cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);

        bool found = false;
        for (int j = 0; j < cellfaceids.size(); j++) {
          if (cellfaceids[j] == faceid) {
            found = true;
            dir = cellfacedirs[j];
            break;
          }
        }

        ASSERT(found);

        JaliGeometry::Point cvec = fcoords[0]-cell_centroids[cellids[i]];
        JaliGeometry::Point trinormal = cvec^evec;

        JaliGeometry::Point normal = evec^trinormal;

        double len = norm(normal);
        normal /= len;
        normal *= *area;

        if (dir == 1)
          *normal0 = normal;
        else
          *normal1 = normal;  // Note that we are not flipping the sign here
      }

      return 1;
    }

  } else if (celldim == 1) {
    face_get_coordinates(faceid, &fcoords);

    JaliGeometry::face1d_get_area(fcoords, geomtype, area);
    JaliGeometry::Point normal(spacedim);
    normal.set(*area);

    Entity_ID_List cellids;
    face_get_cells(faceid, Entity_type::ALL, &cellids);

    for (int i = 0; i < cellids.size(); i++) {
      Entity_ID_List cellfaceids;
      std::vector<dir_t> cellfacedirs;
      dir_t dir = 1;

      cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);

      bool found = false;
      for (int j = 0; j < cellfaceids.size(); j++) {
        if (cellfaceids[j] == faceid) {
          found = true;
          dir = cellfacedirs[j];
          break;
        }
      }

      ASSERT(found);

      if (dir == 1)
        *normal0 = normal;
      else
        *normal1 = -normal;
    }

    return 1;
  }

  return 0;
}  // Mesh::compute_face_geometry


int Mesh::compute_edge_geometry(const Entity_ID edgeid, double *edge_length,
                                JaliGeometry::Point *edge_vector,
                                JaliGeometry::Point *centroid) const {
  (*edge_vector).set(0.0L);
  *edge_length = 0.0;

  Entity_ID node0, node1;

  edge_get_nodes(edgeid, &node0, &node1);

  JaliGeometry::Point point0, point1;
  node_get_coordinates(node0, &point0);
  node_get_coordinates(node1, &point1);

  *edge_vector = point1 - point0;
  *edge_length = norm(*edge_vector);

  *centroid = 0.5*(point0+point1);

  return 0;
}  // Mesh::compute_edge_geometry


void Mesh::compute_side_geometry(Entity_ID const sideid,
                                 double *side_volume,
                                 JaliGeometry::Point *outward_facet_normal,
                                 JaliGeometry::Point *mid_facet_normal) const {
  if (celldim == 3) {
    std::vector<JaliGeometry::Point> scoords;

    // Get vertex coordinates of side
    //
    // These are always - node 0 coordinate, node 1 coordinate, face
    // center, cell center. The way we designated node 0 and
    // node 1, the triangle formed by node 0, node 1 and face center
    // will have a normal point inward.

    side_get_coordinates(sideid, &scoords);

    // vector from node 0 to node 1
    JaliGeometry::Point vec0 = scoords[1] - scoords[0];

    // vector from node 0 to face center
    JaliGeometry::Point vec1 = scoords[2] - scoords[0];

    // vector from node 0 to cell center
    JaliGeometry::Point vec2 = scoords[3] - scoords[0];

    JaliGeometry::Point cpvec = vec0^vec1;

    *side_volume = cpvec*vec2/6.0; 

    // Area weighted normal to the triangular facet formed by node 0
    // coordinate, node 1 coordinate and face center such that the
    // normal is pointing out of the wedge and cell

    *outward_facet_normal = -0.5*cpvec;   // note the -ve sign

    // Area weighted normal of triangular facet formed by edge center,
    // face center and zone center. This is common (with a sign change)
    // to the two wedges of the side

    JaliGeometry::Point ecen = edge_centroid(side_edge_id[sideid]);
    
    JaliGeometry::Point vec3 = scoords[2] - ecen;
    JaliGeometry::Point vec4 = scoords[3] - ecen;

    *mid_facet_normal = 0.5*(vec3^vec4);

  } else if (celldim == 2) {
    std::vector<JaliGeometry::Point> scoords;
    
    // Get vertex coordinates of side
    //
    // These are always - node coordinate, edge/face center, cell center

    side_get_coordinates(sideid, &scoords);

    // vector from node 0 to node 1
    JaliGeometry::Point vec0 = scoords[1]-scoords[0];

    // vector from node 0 to cell center
    JaliGeometry::Point vec1 = scoords[2]-scoords[0];

    // Area of side is 1/2 of the cross-product of vec0 and vec1
    
    *side_volume = norm(vec0^vec1)/2.0;
    
    // length weighted normal to the segment formed
    // by node coordinate and edge/face center

    *outward_facet_normal = JaliGeometry::Point(vec0[1], -vec0[0]);

    JaliGeometry::Point ecen = edge_centroid(side_edge_id[sideid]);

    JaliGeometry::Point vec3 = scoords[2] - ecen;

    *mid_facet_normal = JaliGeometry::Point(vec1[1], -vec1[0]);

  } else if (celldim == 1) {
    std::vector<JaliGeometry::Point> scoords;

    // Get vertex coordinates of side
    //
    // These are always - node coordinate and cell center

    side_get_coordinates(sideid, &scoords);

    // vector from node to cell center
    JaliGeometry::Point vec0 = scoords[1]-scoords[0];
    if (geomtype == JaliGeometry::Geom_type::SPHERICAL) {
      *side_volume = (4.0/3.0) * M_PI * fabs(pow(scoords[1][0], 3) -
                                             pow(scoords[0][0], 3));

      // should point out of zone and have area of face
      *outward_facet_normal = -vec0/JaliGeometry::norm(vec0) *
          4.0 * M_PI * pow(scoords[0][0], 2);
    } else if (geomtype == JaliGeometry::Geom_type::CARTESIAN) {
      *side_volume = JaliGeometry::norm(vec0);

      // should point out of zone and have area of face
      *outward_facet_normal = -vec0/JaliGeometry::norm(vec0);
    }
    
    // Flip sign depending on which node/edge of the cell the side is
    // associated with

    Entity_ID nodeid = side_node_ids[sideid][0];
    Entity_ID cellid = side_cell_id[sideid];

    Entity_ID_List cnodes;
    cell_get_nodes(cellid, &cnodes);
    if (nodeid == cnodes[1]) {  // have to reverse sign of volume and normal
      *side_volume = -(*side_volume);
      *outward_facet_normal = -(*outward_facet_normal);
    }

    // This is not really defined - nobody should try to use this
    mid_facet_normal->set(0.0);
  }
}  // Compute side geometry


void Mesh::compute_corner_geometry(const Entity_ID cornerid,
                                   double *volume) const {
  Entity_ID_List cwedges;
  corner_get_wedges(cornerid, &cwedges);

  *volume = 0;
  Entity_ID_List::iterator itw = cwedges.begin();
  while (itw != cwedges.end()) {
    Entity_ID w = *itw;
    *volume += wedge_volume(w);
    ++itw;
  }
}  // compute corner geometry

// Volume/Area of cell

double Mesh::cell_volume(const Entity_ID cellid, const bool recompute) const {
  if (recompute) {
    double volume;
    JaliGeometry::Point centroid(spacedim);
    compute_cell_geometry(cellid, &volume, &centroid);
    return volume;
  } else {
    return cell_volumes[cellid];
  }
}

// Area/length of face

double Mesh::face_area(const Entity_ID faceid, const bool recompute) const {
  ASSERT(faces_requested);
  ASSERT(face_geometry_precomputed);

  if (recompute) {
    double area;
    JaliGeometry::Point centroid(spacedim);
    JaliGeometry::Point normal0(spacedim), normal1(spacedim);
    compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
    return area;
  } else {
    return face_areas[faceid];
  }
}

// Length of an edge

double Mesh::edge_length(const Entity_ID edgeid, const bool recompute) const {
  ASSERT(edges_requested);
  ASSERT(edge_geometry_precomputed);

  if (recompute) {
    double length;
    JaliGeometry::Point vector(spacedim), centroid(spacedim);
    compute_edge_geometry(edgeid, &length, &vector, &centroid);
    return length;
  } else {
    return edge_lengths[edgeid];
  }
}

// Volume/Area of side

double Mesh::side_volume(const Entity_ID sideid, const bool recompute) const {
  ASSERT(sides_requested);
  ASSERT(side_geometry_precomputed);

  if (recompute) {
    double side_volume;
    JaliGeometry::Point outward_facet_normal, mid_facet_normal;
    compute_side_geometry(sideid, &side_volume, &outward_facet_normal,
                          &mid_facet_normal);
    return side_volume;
  } else {
    return side_volumes[sideid];
  }
}

// Volume/Area of wedge

double Mesh::wedge_volume(const Entity_ID wedgeid, const bool recompute) const {
  ASSERT(wedges_requested);
  ASSERT(side_geometry_precomputed);

  Entity_ID sideid = static_cast<Entity_ID>(wedgeid/2);

  if (recompute) {
    double side_volume;
    JaliGeometry::Point outward_facet_normal, mid_facet_normal;
    compute_side_geometry(sideid, &side_volume, &outward_facet_normal,
                          &mid_facet_normal);
    return side_volume/2.0;
  } else {
    return side_volumes[sideid]/2.0;
  }
}

// Corner volume

double Mesh::corner_volume(const Entity_ID cornerid,
                           const bool recompute) const {
  ASSERT(corners_requested);
  ASSERT(corner_geometry_precomputed);

  if (recompute) {
    double volume;
    compute_corner_geometry(cornerid, &volume);
    return volume;
  } else {
    return corner_volumes[cornerid];
  }
}

// Centroid of cell

JaliGeometry::Point Mesh::cell_centroid(const Entity_ID cellid,
                                        const bool recompute) const {
  ASSERT(cell_geometry_precomputed);

  if (recompute) {
    double volume;
    JaliGeometry::Point centroid(spacedim);
    compute_cell_geometry(cellid, &volume, &centroid);
    return centroid;
  } else {
    return cell_centroids[cellid];
  }
}

// Centroid of face

JaliGeometry::Point Mesh::face_centroid(const Entity_ID faceid,
                                        const bool recompute) const {
  ASSERT(faces_requested);
  ASSERT(face_geometry_precomputed);

  if (recompute) {
    double area;
    JaliGeometry::Point centroid(spacedim);
    JaliGeometry::Point normal0(spacedim), normal1(spacedim);
    compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
    return centroid;
  } else {
    return face_centroids[faceid];
  }
}

// Normal to face
// The vector is normalized and then weighted by the area of the face
//
// If recompute is TRUE, then the normal is recalculated using current
// face coordinates but not stored. (If the recomputed normal must be
// stored, then call recompute_geometric_quantities).
//
// If cellid is not specified, the normal is the natural normal of the face
// If cellid is specified, the normal is the outward normal with respect
// to the cell. In planar and solid meshes, the normal with respect to
// the cell on one side of the face is just the negative of the normal
// with respect to the cell on the other side. In general surfaces meshes,
// this will not be true at C1 discontinuities
//
// if cellid is specified, then orientation returns the direction of
// the natural normal of the face with respect to the cell (1 is
// pointing out of the cell and -1 pointing in)


JaliGeometry::Point Mesh::face_normal(const Entity_ID faceid,
                                      const bool recompute,
                                      const Entity_ID cellid,
                                      int *orientation) const {
  ASSERT(faces_requested);
  ASSERT(face_geometry_precomputed);

  JaliGeometry::Point normal0(spacedim);
  JaliGeometry::Point normal1(spacedim);

  if (recompute) {
    double area;
    JaliGeometry::Point centroid(spacedim);
    
    compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
  } else {
    normal0 = face_normal0[faceid];
    normal1 = face_normal1[faceid];
  }

  if (cellid == -1) {
    // Just the natural normal of the face
    // Since normal0 and normal1 are outward facing normals with respect
    // to their respective cells, we can return normal0 as is but have
    // to negate normal1.

    if (orientation)
      *orientation = 1;

    if (L22(normal0) != 0.0) {
      return normal0;
    } else {
      ASSERT(L22(normal1) != 0.0);
      return -normal1;
    }
  } else {
    Entity_ID_List faceids;
    std::vector<dir_t> face_dirs;

    cell_get_faces_and_dirs(cellid, &faceids, &face_dirs);

    int nf = faceids.size();
    bool found = false;
    dir_t dir = 1;
    for (int i = 0; i < nf; i++)
      if (faceids[i] == faceid) {
        dir = face_dirs[i];
        found = true;
        break;
      }

    ASSERT(found);

    if (orientation) *orientation = dir;
    if (dir == 1) {
      // ASSERT(L22(normal0) != 0.0);
      return normal0;              // Copy to output
    } else {
      // ASSERT(L22(normal1) != 0.0);
      return normal1;              // Copy to output
    }
  }

  return normal0;      // Copy to output
}


// Direction vector of edge

JaliGeometry::Point Mesh::edge_vector(const Entity_ID edgeid,
                                      const bool recompute,
                                      const Entity_ID pointid,
                                      int *orientation) const {
  ASSERT(edges_requested);
  ASSERT(edge_geometry_precomputed);

  JaliGeometry::Point evector(spacedim), ecenter(spacedim);
  JaliGeometry::Point& evector_ref = evector;  // to avoid extra copying

  if (recompute) {
    double length;
    compute_edge_geometry(edgeid, &length, &evector, &ecenter);
    // evector_ref already points to evector
  } else {
    evector_ref = edge_vectors[edgeid];
  }

  if (orientation) *orientation = 1;

  if (pointid == -1) {
    return evector_ref;
  } else {
    Entity_ID p0, p1;
    edge_get_nodes(edgeid, &p0, &p1);

    if (pointid == p0) {
      return evector_ref;
    } else {
      if (orientation) *orientation=-1;
      return -evector_ref;
    }
  }
}  // edge_vector


// Center/Centroid of edge

JaliGeometry::Point Mesh::edge_centroid(const Entity_ID edgeid) const {
  Entity_ID p0, p1;
  JaliGeometry::Point xyz0, xyz1;

  ASSERT(edges_requested);
  ASSERT(edge_geometry_precomputed);

  edge_get_nodes(edgeid, &p0, &p1);
  node_get_coordinates(p0, &xyz0);
  node_get_coordinates(p1, &xyz1);
  return (xyz0+xyz1)/2.0;
}


// Coordinates of a side
//
// If posvol_order = true, then the coordinates will be returned in an
// order that will result in a positive volume (in 3D this assumes
// that the computation for volume is done as (V01xV02).V03 where V0i
// is a vector from coordinate 0 to coordinate i of the tet). If
// posvol_order is false, the coordinates will be returned in a fixed
// order - in 1D, this is node point and cell center, in 2D, this is
// node point 0, node point 1, edge/face center, cell center and in
// 3D, this is node point 0, node point 1, edge center, face center,
// cell center
//
// By default the coordinates are returned in fixed order
// (posvol_order = false)

void Mesh::side_get_coordinates(const Entity_ID sideid,
                                std::vector<JaliGeometry::Point> *scoords,
                                bool posvol_order) const {
  ASSERT(sides_requested);

  if (celldim == 3) {
    scoords->resize(4);  // sides are tets in 3D cells
    Entity_ID n0 = side_get_node(sideid, 0);
    node_get_coordinates(n0, &((*scoords)[0]));

    Entity_ID n1 = side_get_node(sideid, 1);
    node_get_coordinates(n1, &((*scoords)[1]));

    Entity_ID f = side_get_face(sideid);
    (*scoords)[2] = face_centroid(f);

    Entity_ID c = side_get_cell(sideid);
    (*scoords)[3] = cell_centroid(c);

  } else if (celldim == 2) {

    scoords->resize(3);  // sides are tris in 2D cells
    Entity_ID n0 = side_get_node(sideid, 0);
    node_get_coordinates(n0, &((*scoords)[0]));

    Entity_ID n1 = side_get_node(sideid, 1);
    node_get_coordinates(n1, &((*scoords)[1]));

    Entity_ID c = side_get_cell(sideid);
    (*scoords)[2] = cell_centroid(c);

  } else if (celldim == 1) {

    scoords->resize(2);  // sides are segments in 1D cells
    Entity_ID n0 = side_get_node(sideid, 0);
    node_get_coordinates(n0, &((*scoords)[0]));

    Entity_ID c = side_get_cell(sideid);
    (*scoords)[1] = cell_centroid(c);
  }

  // In 2D and 3D, the natural order of coordinates of a side will
  // give a positive volume (unless the cell is too distorted). In 1D,
  // however, we have to check explicitly

  if (posvol_order && celldim == 1) {
    Entity_ID c = side_get_cell(sideid);
    Entity_ID_List cnodes;
    cell_get_nodes(c, &cnodes);
    if (side_get_node(sideid, 0) == cnodes[1])
      std::swap((*scoords)[0], (*scoords)[1]);
  }
}  // side_get_coordinates


// Coordinates of a wedge
//
// If posvol_order = true, then the coordinates will be returned
// in an order that will result in a positive volume (in 3D this assumes
// that the computation for volume is done as (V01xV02).V03 where V0i
// is a vector from coordinate 0 to coordinate i of the tet). If posvol_order
// is false, the coordinates will be returned in a fixed order - in 2D,
// this is node point, edge/face center, cell center and in 3D, this is
// node point, edge center, face center, cell center
//
// By default the coordinates are returned in fixed order
// (posvol_order = false)

void Mesh::wedge_get_coordinates(const Entity_ID wedgeid,
                                 std::vector<JaliGeometry::Point> *wcoords,
                                 bool posvol_order) const {
  ASSERT(wedges_requested);

  int np = celldim + 1;
  wcoords->resize(np);

  Entity_ID n = wedge_get_node(wedgeid);
  node_get_coordinates(n, &((*wcoords)[0]));

  if (celldim == 3) {
    Entity_ID e = wedge_get_edge(wedgeid);
    (*wcoords)[1] = edge_centroid(e);

    Entity_ID f = wedge_get_face(wedgeid);
    (*wcoords)[2] = face_centroid(f);
  }
  else if (celldim == 2) {
    Entity_ID f = wedge_get_face(wedgeid);
    (*wcoords)[1] = face_centroid(f);
  }

  Entity_ID c = wedge_get_cell(wedgeid);
  (*wcoords)[np-1] = cell_centroid(c);

  // If caller has requested that coordinates be ordered such that
  // they will give a +ve volume, we have to reverse the coordinate
  // order for odd numbered wedges. This follows by the convention
  // used for retrieving side coordinates and for the numbering of
  // wedges.

  if (posvol_order && wedgeid%2) {
    if (celldim == 1)
      std::swap((*wcoords)[0], (*wcoords)[1]);
    else
      std::swap((*wcoords)[1], (*wcoords)[2]);
  }
}  // wedge_get_coordinates


// side facet area-weighted normal
//
// side facet is part of the face between the cell of the side and the
// adjacent cell - facet normal points out from the cell into the
// adjacent cell face

JaliGeometry::Point Mesh::side_facet_normal(const int sideid,
                                            const bool recompute) const {
  assert(sides_requested);
  assert(side_geometry_precomputed);

  JaliGeometry::Point normal(spacedim);

  if (recompute) {
    double volume;
    JaliGeometry::Point outward_facet_normal(spacedim),
        mid_facet_normal(spacedim);
    
    compute_side_geometry(sideid, &volume, &outward_facet_normal,
                          &mid_facet_normal);
    return outward_facet_normal;
  } else {
    return side_outward_facet_normal[sideid];
  }
}  // side_facet_normal


// wedge facet area-weighted normal
//
// facet 0 is part of the face between the cell of the wedge and the
// adjacent cell - facet 0 normal points out from the cell into the
// adjacent cell face
//
// facet 1 is shared by the wedge and the other wedge on a common edge
// and face - facet 1 normal points in the same direction as the
// vector from from the node of the wedge to the edge center

JaliGeometry::Point Mesh::wedge_facet_normal(const int wedgeid,
                                             const unsigned int which_facet,
                                             const bool recompute) const {
  assert(wedges_requested);
  assert(side_geometry_precomputed);
  assert(which_facet == 0 || which_facet == 1);

  Entity_ID sideid = static_cast<Entity_ID>(wedgeid/2);

  if (recompute) {
    double volume;
    JaliGeometry::Point facet_normal(spacedim);
    JaliGeometry::Point outward_facet_normal(spacedim),
        mid_facet_normal(spacedim);
        
    compute_side_geometry(sideid, &volume, &outward_facet_normal,
                          &mid_facet_normal);
    if (which_facet == 0)
      return outward_facet_normal/2.0;  // wedge facet area = 1/2 side facet
    else
      return (wedgeid%2) ? -mid_facet_normal : mid_facet_normal;
  } else {
    if (which_facet == 0)
      return side_outward_facet_normal[sideid]/2.0;
    else
      return (wedgeid%2) ? -side_mid_facet_normal[sideid] :
          side_mid_facet_normal[sideid];
  }
}  // wedge_facet_normal


// Triangular facets describing a Corner in 3D

void
Mesh::corner_get_facetization(const Entity_ID cornerid,
                              std::vector<JaliGeometry::Point> *pointcoords,
                              std::vector<std::array<Entity_ID, 3>>
                              *facetpoints) const {
  assert(corners_requested);
  assert(corner_geometry_precomputed);

  Entity_ID_List cwedges;
  corner_get_wedges(cornerid, &cwedges);

  assert(celldim == 3);
  pointcoords->clear();
  pointcoords->reserve(4*(cwedges.size()));  // upper limit
  facetpoints->clear();
  facetpoints->reserve(2*(cwedges.size()));  // 2 facets per wedge will be on
  //                                         // boundary of the corner

  JaliGeometry::Point p(spacedim);

  std::vector< std::pair<Entity_ID, Entity_kind> > point_entity_list;

  int n = corner_get_node(cornerid);
  node_get_coordinates(n, &p);
  pointcoords->push_back(p);        // emplace_back when we switch to C++11
  point_entity_list.push_back(std::pair<Entity_ID, Entity_kind>(n, Entity_kind::NODE));

  int c = corner_get_cell(cornerid);
  JaliGeometry::Point ccen = cell_centroid(c);
  pointcoords->push_back(ccen);
  point_entity_list.push_back(std::pair<Entity_ID, Entity_kind>(c, Entity_kind::CELL));
  JaliGeometry::Point vec0 = ccen-p;

  Entity_ID_List::iterator itw = cwedges.begin();
  while (itw != cwedges.end()) {
    Entity_ID w = *itw;

    Entity_ID f = wedge_get_face(w);
    int idxe = 0;
    bool found = false;
    while (!found && idxe < point_entity_list.size()) {
      if (point_entity_list[idxe] ==
          std::pair<Entity_ID, Entity_kind>(f, Entity_kind::FACE))
        found = true;
      else
        idxe++;
    }
    if (!found) {
      point_entity_list.push_back(std::pair<Entity_ID, Entity_kind>(f, Entity_kind::FACE));
      pointcoords->push_back(face_centroid(f));
    }

    Entity_ID e = wedge_get_edge(w);
    int idxf = 0;
    found = false;
    while (!found && idxf < point_entity_list.size()) {
      if (point_entity_list[idxf] ==
          std::pair<Entity_ID, Entity_kind>(e, Entity_kind::EDGE))
        found = true;
      else
        idxf++;
    }
    if (!found) {
      point_entity_list.push_back(std::pair<Entity_ID, Entity_kind>(e, Entity_kind::EDGE));
      pointcoords->push_back(edge_centroid(e));
    }

    // Now record the facet 0 coords but only after checking that it is
    // pointing out of the wedge and corner
    // Note, idxf is the index for the edge; idxe is the index for the face

    JaliGeometry::Point vec1 = (*pointcoords)[idxf]-(*pointcoords)[0];
    JaliGeometry::Point vec2 = (*pointcoords)[idxe]-(*pointcoords)[0];
    JaliGeometry::Point crossvec = vec1^vec2;
    double volume = crossvec*vec0;

    if (volume > 0) {
      std::array<Entity_ID, 3> pointlist = {0, idxe, idxf};
      facetpoints->push_back(pointlist);
    } else {
      std::array<Entity_ID, 3> pointlist = {0, idxf, idxe};
      facetpoints->push_back(pointlist);
    }

    // record facet 1 points using the same test

    if (volume > 0) {
      std::array<Entity_ID, 3> pointlist = {1, idxf, idxe};
      facetpoints->push_back(pointlist);
    } else {
      std::array<Entity_ID, 3> pointlist = {1, idxe, idxf};
      facetpoints->push_back(pointlist);
    }

    ++itw;
  }  // while (itw != cwedges.end())
}  // corner get facetization for 3D


// "facets" (line segments) describing a corner in 2D

void
Mesh::corner_get_facetization(const Entity_ID cornerid,
                              std::vector<JaliGeometry::Point> *pointcoords,
                              std::vector< std::array<Entity_ID, 2> >
                              *facetpoints) const {
  assert(corners_requested);
  assert(corner_info_cached);

  Entity_ID_List cwedges;
  corner_get_wedges(cornerid, &cwedges);

  assert(celldim == 2);
  pointcoords->clear();
  pointcoords->reserve(4);  // upper limit
  facetpoints->clear();
  facetpoints->reserve(8);  // 2 facets per wedge (2 points per facet) will be
  //                        // on boundary of the corner

  JaliGeometry::Point p(spacedim);

  int n = corner_get_node(cornerid);
  node_get_coordinates(n, &p);
  pointcoords->push_back(p);        // emplace_back when we switch to C++11

  int c = corner_get_cell(cornerid);
  JaliGeometry::Point ccen = cell_centroid(c);
  JaliGeometry::Point vec0 = ccen-p;

  Entity_ID f0 = wedge_get_face(cwedges[0]);
  JaliGeometry::Point fcen0 = face_centroid(f0);
  Entity_ID f1 = wedge_get_face(cwedges[1]);
  JaliGeometry::Point fcen1 = face_centroid(f1);

  JaliGeometry::Point vec1 = fcen0-p;
  JaliGeometry::Point cpvec = vec1^vec0;
  if (cpvec[0] > 0) {
    pointcoords->push_back(fcen0);
    pointcoords->push_back(ccen);
    pointcoords->push_back(fcen1);
  } else {
    pointcoords->push_back(fcen1);
    pointcoords->push_back(ccen);
    pointcoords->push_back(fcen0);
  }

  facetpoints->push_back({{0, 1}});
  facetpoints->push_back({{1, 2}});
  facetpoints->push_back({{2, 3}});
  facetpoints->push_back({{3, 0}});

}  // corner get facetization for 2D


// "facets" (points - node and cell center) describing a corner in 1D :)

void
Mesh::corner_get_facetization(const Entity_ID cornerid,
                              std::vector<JaliGeometry::Point> *pointcoords,
                              std::vector< std::array<Entity_ID, 1> >
                              *facetpoints) const {
  assert(corners_requested);
  assert(corner_info_cached);

  // corner and wedge are the same in 1d
  Entity_ID_List cwedges;
  corner_get_wedges(cornerid, &cwedges);
  wedge_get_coordinates(cwedges[0], pointcoords);
  // ordering is because wedge_get_coordinates comes back in (node, cell) order
  // and node is the facet external to the zone and cell is the facet between
  // the wedges of the side, which is consistent with 2/3D.
  facetpoints->push_back({{0}});
  facetpoints->push_back({{1}});
}

// list of points describing a Corner in 2D in ccw direction and in no
// particular order in 3D

void
Mesh::corner_get_coordinates(const Entity_ID cornerid,
                             std::vector<JaliGeometry::Point>
                             *pointcoords) const {
  assert(corners_requested);
  assert(corner_info_cached);

  Entity_ID_List cwedges;
  corner_get_wedges(cornerid, &cwedges);

  pointcoords->clear();

  if (celldim == 1) {
    pointcoords->reserve(2);

    // 1D wedge coordinates are - node point and zone center; we always want the
    // positive volume order, hence the 'true' in the call
    // wedges are the same as corners

    wedge_get_coordinates(cwedges[0], pointcoords, true);
  } else if (celldim == 2) {
    pointcoords->reserve(4);

    // 2D wedge coordinates are - node point, edge center, zone center

    wedge_get_coordinates(cwedges[0], pointcoords);

    std::vector<JaliGeometry::Point> wpoints;
    wedge_get_coordinates(cwedges[1], &wpoints);

    pointcoords->push_back(wpoints[1]);

    // Make sure we got the coordinates in the right orientation

    JaliGeometry::Point vec0 = (*pointcoords)[1] - (*pointcoords)[0];
    JaliGeometry::Point vec1 = (*pointcoords)[3] - (*pointcoords)[0];
    JaliGeometry::Point cpvec = vec0^vec1;  // 2 times the vol of cwedges[0]

    // If the sign of cpvec[0] and wedge_volume(cwedges[0]) not the same, then
    // switch the order of the coordinates

    if (cpvec[0]*wedge_volume(cwedges[0]) < 0) {
      // reverse the order of the points
      std::swap((*pointcoords)[1], (*pointcoords)[3]);
    }
  } else if (celldim == 3) {
    pointcoords->reserve(4*(cwedges.size()));  // upper limit

    JaliGeometry::Point p(spacedim);

    std::vector< std::pair<Entity_ID, Entity_kind> > point_entity_list;

    int n = corner_get_node(cornerid);
    node_get_coordinates(n, &p);
    pointcoords->push_back(p);        // emplace_back when we switch to C++11
    point_entity_list.push_back(std::pair<Entity_ID, Entity_kind>(n, Entity_kind::NODE));

    int c = corner_get_cell(cornerid);
    JaliGeometry::Point ccen = cell_centroid(c);

    Entity_ID_List::iterator itw = cwedges.begin();
    while (itw != cwedges.end()) {
      Entity_ID w = *itw;

      Entity_ID f = wedge_get_face(w);
      int idxe = 0;
      bool found = false;
      while (!found && idxe < point_entity_list.size()) {
        if (point_entity_list[idxe] ==
            std::pair<Entity_ID, Entity_kind>(f, Entity_kind::FACE))
          found = true;
        else
          idxe++;
      }
      if (!found) {
        point_entity_list.push_back(std::pair<Entity_ID, Entity_kind>(f, Entity_kind::FACE));
        pointcoords->push_back(face_centroid(f));
      }

      Entity_ID e = wedge_get_edge(w);
      int idxf = 0;
      found = false;
      while (!found && idxf < point_entity_list.size()) {
        if (point_entity_list[idxf] ==
            std::pair<Entity_ID, Entity_kind>(e, Entity_kind::EDGE))
          found = true;
        else
          idxf++;
      }
      if (!found) {
        point_entity_list.push_back(std::pair<Entity_ID, Entity_kind>(e, Entity_kind::EDGE));
        pointcoords->push_back(edge_centroid(e));
      }

      ++itw;
    }

    pointcoords->push_back(ccen);
  }
}  // corner_get_points


Entity_type Mesh::entity_get_type(const Entity_kind kind,
                                     const Entity_ID entid) const {
  switch (kind) {
    case Entity_kind::CELL:
      return cell_type[entid];
      break;
    case Entity_kind::FACE:
      return faces_requested ? face_type[entid] :
          Entity_type::TYPE_UNKNOWN;
      break;
    case Entity_kind::EDGE:
      return edges_requested ? edge_type[entid] :
          Entity_type::TYPE_UNKNOWN;
      break;
    case Entity_kind::NODE:
      return node_type[entid];
      break;
    case Entity_kind::SIDE:
      if (sides_requested) {
        Entity_ID cellid = side_cell_id[entid];
        return cell_type[cellid];
      } else
        return Entity_type::TYPE_UNKNOWN;
      break;
    case Entity_kind::WEDGE:
      if (wedges_requested) {
        Entity_ID sideid = static_cast<int>(entid/2);
        Entity_ID cellid = side_cell_id[sideid];
        return cell_type[cellid];
      } else
        return Entity_type::TYPE_UNKNOWN;
      break;
    case Entity_kind::CORNER:
      if (corners_requested) {
        Entity_ID wedgeid = corner_wedge_ids[entid][0];
        Entity_ID sideid = static_cast<int>(wedgeid/2);
        Entity_ID cellid = side_cell_id[sideid];
        return cell_type[cellid];
      } else
        return Entity_type::TYPE_UNKNOWN;
      break;
    default:
      return Entity_type::TYPE_UNKNOWN;
  }
}


// Is there a set with this name and entity type

bool Mesh::valid_set_name(std::string name, Entity_kind kind) const {
  if (!geometric_model_) {
    Errors::Message mesg("Mesh sets not enabled because mesh was created"
                         " without reference to a geometric model");
    Exceptions::Jali_throw(mesg);
  }

  unsigned int gdim = geometric_model_->dimension();

  unsigned int ngr = geometric_model_->Num_Regions();
  for (int i = 0; i < ngr; i++) {
    JaliGeometry::RegionPtr rgn = geometric_model_->Region_i(i);

    unsigned int rdim = rgn->dimension();

    if (rgn->name() == name) {

      // For regions of type Color Function, the dimension
      // parameter is not guaranteed to be correct

      if (rgn->type() == JaliGeometry::Region_type::COLORFUNCTION) return true;

      // For regions of type Labeled set, extract some more info and verify

      if (rgn->type() == JaliGeometry::Region_type::LABELEDSET) {
        JaliGeometry::LabeledSetRegionPtr lsrgn =
            dynamic_cast<JaliGeometry::LabeledSetRegionPtr> (rgn);
        std::string entity_type = lsrgn->entity_str();

        if ((kind == Entity_kind::CELL && entity_type == "CELL") ||
            (kind == Entity_kind::FACE && entity_type == "FACE") ||
            (kind == Entity_kind::EDGE && entity_type == "EDGE") ||
            (kind == Entity_kind::NODE && entity_type == "NODE"))
          return true;
        else
          return false;
      }

      // If we are looking for a cell set the region has to be
      // of the same topological dimension as the cells or it
      // has to be a point region

      if (kind == Entity_kind::CELL && (rdim >= celldim || rdim == 0))
        return true;

      // If we are looking for a side set, the region has to be
      // one topological dimension less than the cells

      if (kind == Entity_kind::FACE && rdim >= celldim-1) return true;

      // If we are looking for a node set, the region can be of any
      // dimension upto the spatial dimension of the domain

      if (kind == Entity_kind::NODE) return true;
    }
  }
  return false;
}


bool Mesh::point_in_cell(const JaliGeometry::Point &p,
                         const Entity_ID cellid) const {
  std::vector<JaliGeometry::Point> ccoords;

  if (celldim == 3) {

    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine

    int nf;
    Entity_ID_List faces;
    std::vector<unsigned int> nfnodes;
    std::vector<dir_t> fdirs;
    std::vector<JaliGeometry::Point> cfcoords;

    cell_get_faces_and_dirs(cellid, &faces, &fdirs);

    nf = faces.size();
    nfnodes.resize(nf);

    for (int j = 0; j < nf; j++) {
      std::vector<JaliGeometry::Point> fcoords;

      face_get_coordinates(faces[j], &fcoords);
      nfnodes[j] = fcoords.size();

      if (fdirs[j] == 1) {
        for (int k = 0; k < nfnodes[j]; k++)
          cfcoords.push_back(fcoords[k]);
      } else {
        for (int k = nfnodes[j]-1; k >=0; k--)
          cfcoords.push_back(fcoords[k]);
      }
    }

    cell_get_coordinates(cellid, &ccoords);
    return JaliGeometry::point_in_polyhed(p, ccoords, nf, nfnodes, cfcoords);

  } else if (celldim == 2) {

    cell_get_coordinates(cellid, &ccoords);
    return JaliGeometry::point_in_polygon(p, ccoords);

  } else if (celldim == 1) {
    cell_get_coordinates(cellid, &ccoords);
    if (p[0]-ccoords[0][0] >= 0.0 &&
        ccoords[1][0] - p[0] >= 0.0) return true;
  }

  return false;
}


void
Mesh::get_partitioning(int const num_parts,
                       Partitioner_type const partitioner,
                       std::vector<std::vector<Entity_ID>> *partitions) {


  partitions->resize(num_parts);

  if (partitioner == Partitioner_type::METIS && space_dimension() != 1) {

#ifdef Jali_HAVE_METIS

    get_partitioning_with_metis(num_parts, partitions);

// #elif Jali_HAVE_ZOLTAN

//     std::cerr << "Mesh::get_partitioning() - " <<
//         "Preferred partitioner METIS not available. Using ZOLTAN.\n";

//     get_partitioning_with_zoltan(num_tiles_ini, partitions);

#else
    std::cerr << "Mesh::get_partitioning() - " <<
        "Preferred partitioner METIS not available. " <<
        "Partitioning by entity IDs.\n";

     get_partitioning_by_index_space(num_parts, partitions);

#endif
  } else if ((partitioner == Partitioner_type::ZOLTAN_GRAPH ||
              partitioner == Partitioner_type::ZOLTAN_RCB) &&
             space_dimension() != 1) {

// #ifdef Jali_HAVE_ZOLTAN

//     get_partitioning_with_zoltan(num_parts, partitions);

// #elif Jali_HAVE_METIS

//     std::cerr << "Mesh::get_partitioning() - " <<
//         "Preferred partitioner ZOLTAN not available. " <<
//         " Using METIS.\n";

//     get_partitioning_with_metis(num_parts, partitions);

// #else

    std::cerr << "Mesh::get_partitioning() - " <<
        "Requested partitioner ZOLTAN not available. " <<
        "Partitioning by entity IDs\n";

    get_partitioning_by_index_space(num_parts, partitions);

// #endif
  } else if (partitioner == Partitioner_type::BLOCK) {

    get_partitioning_by_blocks(num_parts, partitions);

  } else {
    if (space_dimension() != 1) {
      std::cerr << "Mesh::get_partitioning() - " <<
          "No valid partitioner or partioning scheme specified. " <<
          "Partitioning by entity IDs\n";
    }

    get_partitioning_by_index_space(num_parts, partitions);
  }
}


void Mesh::get_partitioning_by_index_space(int const num_parts,
                                 std::vector<std::vector<int>> *partitions) {
        
  std::cerr << "No partitioner defined - " <<
      "subdividing cell index space into equal parts\n";

  int ncells = num_cells<Entity_type::PARALLEL_OWNED>();
  int maxcells_per_tile =
      static_cast<int> (static_cast<double>(ncells)/num_parts + 0.5);
  int index = 0;

  int *ncells_per_tile = new int[num_parts];
  for (int i = 0; i < num_parts; ++i) {
    (*partitions)[i].resize(maxcells_per_tile);  // excludes MPI ghosts
    ncells_per_tile[i] = 0;

    // Owned cells in this partition
    for (int j = 0; j < maxcells_per_tile && index < ncells; ++j) {
      (*partitions)[i][j] = index++;
      (ncells_per_tile[i])++;
    }

    (*partitions)[i].resize(ncells_per_tile[i]);
  }
  delete [] ncells_per_tile;
}

#ifdef Jali_HAVE_METIS

void Mesh::get_partitioning_with_metis(int const num_parts,
                                  std::vector<std::vector<int>> *partitions) {

  std::cerr << "Partitioning mesh on each compute node into " << num_parts <<
      " parts using METIS\n";

  // Build the adjacency info

  int ncells_owned = num_cells<Entity_type::PARALLEL_OWNED>();
  int nfaces_owned = num_faces<Entity_type::PARALLEL_OWNED>();
  int ncells_all = num_cells<Entity_type::ALL>();
  int nfaces_all = num_faces<Entity_type::ALL>();
  idx_t *xadj = new idx_t[ncells_all+1];  // offset indicator into adjacency
  //                                      // array - overallocate
  idx_t *adjncy = new idx_t[2*nfaces_all];  // adjacency array (nbrs of a cell)
  //                                        // - overallocate
  xadj[0] = 0;  // starting offset
  int ipos = 0;
  int i = 0;
  Entity_ID_List nbrcells;
  for (auto const& c : cells<Entity_type::PARALLEL_OWNED>()) {
    cell_get_face_adj_cells(c, Entity_type::PARALLEL_OWNED, &nbrcells);
    for (auto const& nbr : nbrcells)
      adjncy[ipos++] = nbr;
    xadj[++i] = ipos;
  }
  if (ipos > 2*nfaces_all)
    std::cerr << "Mesh::get_partitioning - Wrote past end of adjncy array" <<
        "Expand adjncy array or expect crashes\n";

  // Partition the graph

  idx_t wtflag = 0;  // No weights
  idx_t *vwt = nullptr;
  idx_t *adjwt = nullptr;

  int numflag = 0;  // C style numbering of cells (nodes of the dual graph)
  idx_t ngraphvtx = ncells_owned;
  idx_t nparts = static_cast<idx_t>(num_parts);
  
  idx_t *idxpart = new idx_t[ngraphvtx];

  idx_t ncons = 1;  // Number of constraints
  idx_t *vsize = nullptr;
  real_t *tpwts = nullptr;
  real_t *ubvec = nullptr;
  idx_t nedgecut;
  
  idx_t metisopts[METIS_NOPTIONS];
  METIS_SetDefaultOptions(metisopts);  // Get default options from METIS
  metisopts[METIS_OPTION_NUMBERING] = 0;   // Indicate C style numbering

  if (num_parts <= 8)
    METIS_PartGraphRecursive(&ngraphvtx, &ncons, xadj, adjncy, vwt, vsize,
                             adjwt, &nparts, tpwts, ubvec,
                             metisopts, &nedgecut, idxpart);
  else
    METIS_PartGraphKway(&ngraphvtx, &ncons, xadj, adjncy, vwt, vsize,
                        adjwt, &nparts, tpwts, ubvec, metisopts,
                        &nedgecut, idxpart);

  for (int i = 0; i < ncells_owned; ++i) {
    int ipart = idxpart[i];
    (*partitions)[ipart].push_back(i);
  }

  delete [] xadj;
  delete [] adjncy;
  delete [] idxpart;
}

#endif

#ifdef Jali_HAVE_ZOLTAN

void Mesh::get_partitioning_with_zoltan(int const num_parts,
                                  std::vector<std::vector<int>> *partitions) {
}

#endif


// Get a partitioning of a mesh by recursive subdivision such that
// each partition is a block. Will fail for meshes that are not
// regular meshes with cells of regular spacing.

void Mesh::get_partitioning_by_blocks(int const num_parts,
                                      std::vector<std::vector<int>> *partitions) {
  // Figure out the domain limits

  std::vector<double> domainlimits(2*spacedim);
  for (int d = 0; d < spacedim; ++d) {
    domainlimits[2*d] = 1e20;
    domainlimits[2*d+1] = -1e20;
  }

  for (auto const& c : cells<Entity_type::PARALLEL_OWNED>()) {
    std::vector<JaliGeometry::Point> cellpnts;
    cell_get_coordinates(c, &cellpnts);
    for (auto const& p : cellpnts) {
      for (int d = 0; d < spacedim; ++d) {
        if (p[d] < domainlimits[2*d]) domainlimits[2*d] = p[d];
        if (p[d] > domainlimits[2*d+1]) domainlimits[2*d+1] = p[d];
      }
    }
  }
  
  // Try to figure out if this is a regular grid from which we
  // can decipher how many cells are in each direction

  std::vector<int> num_cells_in_dir(spacedim, 0);

  // First find the lower left most cell
 
  Entity_ID c0 = -1;
  bool found = false;
  for (auto const& c : cells<Entity_type::PARALLEL_OWNED>()) {
    std::vector<JaliGeometry::Point> cellpnts;
    cell_get_coordinates(c, &cellpnts);
    for (auto const& p : cellpnts) {
      double dist = 0.0;
      for (int d = 0; d < spacedim; d++)
        dist += (domainlimits[2*d] - p[d])*(domainlimits[2*d] - p[d]);
      if (dist < 1.0e-20) {
        found = true;
        c0 = c;
        break;
      }
    }
    if (found) break;
  }
  if (!found) {
    std::cerr << "This is not a regular mesh - cannot partition using " <<
        "Partitioner_type::BLOCK\n";
    return;
  }

  // Try to walk from cell to cell in each direction to the high
  // boundary in that direction

  for (int d = 0; d < spacedim; d++) {
    // direction in which we should be walking

    JaliGeometry::Point refdir(spacedim);
    refdir[d] = 1.0;
    
    Entity_ID c = c0;
    bool done = 0;
    num_cells_in_dir[d] = 1;
    while (!done) {
      Entity_ID_List cfaces;
      cell_get_faces(c, &cfaces);
      
      // Find a face whose normal w.r.t. the cell points in the
      // 'refdir' direction

      bool next_face_found = false;
      for (auto const& f : cfaces) {
        JaliGeometry::Point fnormal = face_normal(f, false, c);
        fnormal /= norm(fnormal);  // normalization
        double dp = fnormal*refdir;
        if (fabs(dp-1) < 1.0e-12) {
          next_face_found = true;
          
          // find the adjacent cell

          Entity_ID_List fcells;
          face_get_cells(f, Entity_type::PARALLEL_OWNED, &fcells);

          if (fcells.size() == 2) {
            c = (fcells[0] == c) ? fcells[1] : fcells[0];
            num_cells_in_dir[d]++;
            next_face_found = true;
            break;
          } else {
            // reached a boundary face - this had better be the high
            // boundary

            JaliGeometry::Point fcen = face_centroid(f);
            if (fabs(fcen[d]-domainlimits[2*d+1]) < 1.0e-12) {
              done = 1;
              break;
            } else {
              std::cerr << "Reached a boundary but it does not match the " <<
                  "high boundary of the domain in this direction! " <<
                  "Mesh may not be a regular mesh\n";
              return;
            }
          }
        }
      }

      if (!next_face_found && !done) {
        std::cerr << "Could not find the next face and cell along direction " <<
            d << ". Mesh may not be a regular mesh\n";
        return;
      }
    }

    if (!done)
      std::cerr << "Could not find a regular structure along direction " << d
                << "\n";
  }


  // Now get the partitioning blocks
  std::vector<std::array<double, 6>> blocklimits;
  std::vector<std::array<int, 3>> blocknumcells;

  block_partition_regular_mesh(spacedim, &(domainlimits[0]),
                               &(num_cells_in_dir[0]), num_parts,
                               &blocklimits, &blocknumcells);

  for (auto const& c : cells<Entity_type::PARALLEL_OWNED>()) {
    JaliGeometry::Point ccen = cell_centroid(c);
    for (int i = 0; i < num_parts; i++) {
      bool inside = true;
      for (int d = 0; d < spacedim; d++) {
        inside = (blocklimits[i][2*d] < ccen[d] &&
                  ccen[d] < blocklimits[i][2*d+1]);
        if (!inside) break;
      }
      if (inside) {
        ((*partitions)[i]).push_back(c);
        break;
      }
    }
  }

  // Sanity check - make sure each block got the number of cells it
  // was supposed to

  for (int i = 0; i < num_parts; i++) {
    int expected_count = 1;
    for (int d = 0; d < spacedim; d++)
      expected_count *= blocknumcells[i][d];
    if (((*partitions)[i]).size() != expected_count)
      std::cerr << "Partition " << i << " has fewer cells than expected\n";
  }
}

// @brief Get the partitioning of a regular mesh such that each
// partition is a rectangular block
//
// @param dim Dimension of problem - 1, 2 or 3
// @param domain 2*dim values for min/max of domain
//  (xmin, xmax, ymin, ymax, zmin, zmax)
// @param num_cells_in_dir  number of cells in each direction
// @param num_blocks_requested number of blocks requested
// @param blocklimits min/max limits for each block
// @param blocknumcells num cells in each direction for blocks
//
// Returns 1 if successful, 0 otherwise

int Mesh::block_partition_regular_mesh(int const dim,
                                       double const * const domain,
                                       int const * const num_cells_in_dir,
                                       int const num_blocks_requested,
                                       std::vector<std::array<double, 6>> *blocklimits,
                                       std::vector<std::array<int, 3>> *blocknumcells) {

  // Create local block arrays that are larger than the requested number of
  // blocks. The extra storage is for temporary blocks while we are subdividing

  std::vector<std::array<double, 6>> blimits(2*num_blocks_requested);
  std::vector<std::array<int, 3>> bnumcells(2*num_blocks_requested);
  int nblocks = 1;
  int nblocks_in_dir[3] = {1, 1, 1};
  int nblockcells_in_dir[3] = {0, 0, 0};
  for (int i = 0; i < dim; i++)
    nblockcells_in_dir[i] = num_cells_in_dir[i];

  if (num_blocks_requested > 1) {
    // First try bisection
    
    bool done = false;
    while (!done) {
      bool bisected = false;
      if (nblockcells_in_dir[0]%2 == 0) {  // even number of cells in x
        nblockcells_in_dir[0] /= 2;
        nblocks *= 2;
        nblocks_in_dir[0] *= 2;
        bisected = true;
        
        if (nblocks == num_blocks_requested) {
          done = 1;
          continue;
        }
      }
      if (dim > 1 && nblockcells_in_dir[1]%2 == 0) {  // even num of cells in y
        nblockcells_in_dir[1] /= 2;
        nblocks *= 2;
        nblocks_in_dir[1] *= 2;
        
        bisected = true;
        if (nblocks == num_blocks_requested) {
          done = 1;
          continue;
        }
      }
      if (dim > 2 && nblockcells_in_dir[2]%2 == 0) {  // even num of cells in z
        nblockcells_in_dir[2] /= 2;
        nblocks *= 2;
        nblocks_in_dir[2] *= 2;
        bisected = true;
        
        if (nblocks == num_blocks_requested) {
          done = 1;
          continue;
        }
      }
      if (!bisected) {
        // Reached a state where we cannot evenly subdivide the number
        // of elements in any one direction
        done = 1;
      }
    }
  }

  // Populate the block details

  double delta[3] = {0.0, 0.0, 0.0};
  for (int dir = 0; dir < dim; dir++)
    delta[dir] = (domain[2*dir+1]-domain[2*dir])/nblocks_in_dir[dir];

  int partnum = 0;
  for (int i = 0; i < nblocks_in_dir[0]; i++) {
    for (int j = 0; j < nblocks_in_dir[1]; j++) {
      for (int k = 0; k < nblocks_in_dir[2]; k++) {
        blimits[partnum][0] = domain[0] + i*delta[0];
        blimits[partnum][1] = (i < nblocks_in_dir[0]-1) ?
            domain[0] + (i+1)*delta[0] : domain[1];
        bnumcells[partnum][0] = nblockcells_in_dir[0];

        blimits[partnum][2] = domain[2] + j*delta[1];
        blimits[partnum][3] = (j < nblocks_in_dir[1]-1) ?
            domain[2] + (j+1)*delta[1] : domain[3];
        bnumcells[partnum][1] = nblockcells_in_dir[1];

        blimits[partnum][4] = domain[4] + k*delta[2];
        blimits[partnum][5] = (k < nblocks_in_dir[2]-1) ?
            domain[4] + (k+1)*delta[2] : domain[5];
        bnumcells[partnum][2] = nblockcells_in_dir[2];
        partnum++;
      }
    }
  }

  if (nblocks != num_blocks_requested) {

    // Start dividing the blocks as unevenly to get the number of
    // partitions we need
    
    bool done = false;
    while (!done) {
      for (int dir = 0; dir < dim; dir++) {  // split blocks in x, then y etc
        int nnewblocks = 0;
        for (int ib = 0; ib < nblocks; ib++) {
          if (bnumcells[ib][0] == 0 && bnumcells[ib][1] == 0 &&
              bnumcells[ib][2] == 0) continue;  // Blanked out block
          
          double diff = blimits[ib][2*dir+1] - blimits[ib][2*dir];
          double delta = diff/bnumcells[ib][dir];
          int ncells1 = bnumcells[ib][dir]/2;
          int ncells2 = bnumcells[ib][dir] - ncells1;
          if (ncells1 == 0 || ncells2 == 0) continue;
          
          // Make two new partitions at the end of the list
          
          // Copy the original block details over
          int ib1 = nblocks + nnewblocks;
          for (int dir1 = 0; dir1 < dim; dir1++) {
            bnumcells[ib1][dir1] = bnumcells[ib][dir1];
            blimits[ib1][2*dir1] = blimits[ib][2*dir1];
            blimits[ib1][2*dir1+1] = blimits[ib][2*dir1+1];
          }
          /* overwrite the data in the direction of the refinement */
          bnumcells[ib1][dir] = ncells1;
          blimits[ib1][2*dir] = blimits[ib][2*dir];
          blimits[ib1][2*dir+1] = blimits[ib][2*dir] + ncells1*delta;
          
          // copy the initial block details over
          int ib2 = nblocks + nnewblocks + 1;
          for (int dir1 = 0; dir1 < dim; dir1++) {
            bnumcells[ib2][dir1] = bnumcells[ib][dir1];
            blimits[ib2][2*dir1] = blimits[ib][2*dir1];
            blimits[ib2][2*dir1+1] = blimits[ib][2*dir1+1];
          }
          /* overwrite the data in the direction of the refinement */
          bnumcells[ib2][dir] = ncells2;
          blimits[ib2][2*dir] = blimits[ib][2*dir] + ncells1*delta;
          blimits[ib2][2*dir+1] = blimits[ib][2*dir+1];
          
          // Blank out the original block
          bnumcells[ib][0] = bnumcells[ib][1] =
              bnumcells[ib][2] = 0;
          nnewblocks += 2;
          
          // Check if we reached the requested number of blocks. Each
          // block that was split into two will cause the loss of one
          // block and gain of two new blocks, so the net gain is just
          // nnewblocks/2
          
          if (nblocks + nnewblocks/2 >= num_blocks_requested) {
            done = 1;
            break;
          }
        }
        
        // Squeeze out the blocks that were split and dummied out
        for (int ib = nblocks-1; ib >= 0; ib--) {
          if (bnumcells[ib][0] == 0 && bnumcells[ib][1] == 0 &&
              bnumcells[ib][2] == 0) {  // dummy block
            for (int ib1 = ib; ib1 < nblocks+nnewblocks-1; ib1++) {
              for (int dir1 = 0; dir1 < 3; dir1++) {
                bnumcells[ib1][dir1] = bnumcells[ib1+1][dir1];
                blimits[ib1][2*dir1] = blimits[ib1+1][2*dir1];
                blimits[ib1][2*dir1+1] = blimits[ib1+1][2*dir1+1];
              }
            }
            // Blank out the last block
            bnumcells[nblocks+nnewblocks-1][0] =
                bnumcells[nblocks+nnewblocks-1][1] =
                bnumcells[nblocks+nnewblocks-1][2] = 0;
            nblocks--;
          }
        }
        nblocks += nnewblocks;
        if (done)
          break;
      }
    }
  }  // if (nblocks != num_blocks_requested)

  blimits.resize(nblocks);
  bnumcells.resize(nblocks);

  blocklimits->resize(nblocks);
  blocknumcells->resize(nblocks);
  std::copy(blimits.begin(), blimits.end(), blocklimits->begin());
  std::copy(bnumcells.begin(), bnumcells.end(), blocknumcells->begin());

  return 1;
}

}  // close namespace Jali
