//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#include "Geometry.hh"
#include "dbc.hh"
#include "errors.hh"
#include "LabeledSetRegion.hh"

#include "Mesh.hh"
#include "MeshTile.hh"  // Always include MeshTile.hh after Mesh.hh


#define CACHE_VARS 1

namespace Jali {

// Forward declaration of make_meshtile function which the mesh class
// references - need to do this to avoid circular dependence between
// Mesh and MeshTile classes

std::shared_ptr<MeshTile> make_meshtile(Mesh const & parent_mesh,
                                        const std::vector<int> & cells,
                                        const bool request_faces,
                                        const bool request_edges,
                                        const bool request_wedges,
                                        const bool request_corners);

// Gather and cache cell to face connectivity info.
//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh calss for further explanation

void Mesh::cache_cell2face_info() const {

  int ncells = num_cells<Parallel_type::ALL>();
  cell_face_ids.resize(ncells);
  cell_face_dirs.resize(ncells);

  for (int c = 0; c < ncells; c++)
    cell_get_faces_and_dirs_internal(c, &(cell_face_ids[c]),
                                     &(cell_face_dirs[c]), false);

  cell2face_info_cached = true;
  faces_requested = true;
}


// Gather and cache face to cell connectivity info.
//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh calss for further explanation

void Mesh::cache_face2cell_info() const {

  int nfaces = num_faces<Parallel_type::ALL>();
  face_cell_ids.resize(nfaces);
  face_cell_ptype.resize(nfaces);

  std::vector<Entity_ID> fcells;

  for (int f = 0; f < nfaces; f++) {
    face_get_cells_internal(f, Parallel_type::ALL, &fcells);

    face_cell_ids[f].resize(2);
    face_cell_ptype[f].resize(2);

    for (int i = 0; i < fcells.size(); i++) {
      int c = fcells[i];
      face_cell_ids[f][i] = c;
      face_cell_ptype[f][i] = entity_get_ptype(Entity_kind::CELL, c);
    }
    for (int i = fcells.size(); i < 2; i++) {
      face_cell_ids[f][i] = -1;
      face_cell_ptype[f][i] = Parallel_type::PTYPE_UNKNOWN;
    }
  }

  face2cell_info_cached = true;
  faces_requested = true;
}

// Gather and cache face to edge connectivity info.
//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh calss for further explanation

void Mesh::cache_face2edge_info() const {

  int nfaces = num_faces<Parallel_type::ALL>();
  face_edge_ids.resize(nfaces);
  face_edge_dirs.resize(nfaces);

  for (int f = 0; f < nfaces; f++) {
    Entity_ID_List fedgeids;
    std::vector<int> fedgedirs;

    face_get_edges_and_dirs_internal(f, &(face_edge_ids[f]),
                                     &(face_edge_dirs[f]), true);
  }

  face2edge_info_cached = true;
  faces_requested = true;
  edges_requested = true;
}

// Gather and cache cell to edge connectivity info.
//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh calss for further explanation

void Mesh::cache_cell2edge_info() const {

  int ncells = num_cells<Parallel_type::ALL>();
  cell_edge_ids.resize(ncells);

  if (spacedim == 2) {
    cell_2D_edge_dirs.resize(ncells);
    for (int c = 0; c < ncells; c++)
      cell_2D_get_edges_and_dirs_internal(c, &(cell_edge_ids[c]),
                                          &(cell_2D_edge_dirs[c]));
  }
  else
    for (int c = 0; c < ncells; c++)
      cell_get_edges_internal(c, &(cell_edge_ids[c]));

  cell2edge_info_cached = true;
}


// Gather and cache wedge information

void Mesh::cache_wedge_info() const {

  if (wedge_info_cached) return;

  int ncells_owned = num_cells<Parallel_type::OWNED>();
  int ncells_ghost = num_cells<Parallel_type::GHOST>();
  int ncells = ncells_owned + ncells_ghost;

  cell_wedge_ids.resize(ncells);

  int nnodes_owned = num_nodes<Parallel_type::OWNED>();
  int nnodes_ghost = num_nodes<Parallel_type::GHOST>();
  int nnodes = nnodes_owned + nnodes_ghost;

  node_wedge_ids.resize(nnodes);

  int num_wedges_all = 0;
  int num_wedges_owned = 0;
  int num_wedges_ghost = 0;

  for (auto const & c : cells()) {
    Parallel_type ptype = entity_get_ptype(Entity_kind::CELL,c);

    std::vector<Entity_ID> cfaces;
    cell_get_faces(c, &cfaces);

    int numwedges_in_cell = 0;
    for (auto const & f : cfaces) {
      std::vector<Entity_ID> fedges;
      std::vector<int> fedirs;
      face_get_edges_and_dirs(f, &fedges, &fedirs);

      num_wedges_all += 2*fedges.size();  // In 2D there will be 1 edge per face
      numwedges_in_cell += 2*fedges.size();

      if (ptype == Parallel_type::OWNED) 
        num_wedges_owned += 2*fedges.size();
      else
        num_wedges_ghost += 2*fedges.size();
    }

    cell_wedge_ids[c].reserve(numwedges_in_cell);
  }

  wedgeids_owned_.resize(num_wedges_owned);
  wedgeids_ghost_.resize(num_wedges_ghost);
  wedgeids_all_.resize(num_wedges_all);
  wedge_node_id.resize(num_wedges_all, -1);
  wedge_edge_id.resize(num_wedges_all, -1);
  wedge_face_id.resize(num_wedges_all, -1);
  wedge_cell_id.resize(num_wedges_all, -1);
  wedge_corner_id.resize(num_wedges_all, -1);  // filled when building corners
  wedge_parallel_type.resize(num_wedges_all, Parallel_type::OWNED);
  wedge_adj_wedge_id.resize(num_wedges_all, -1);
  wedge_opp_wedge_id.resize(num_wedges_all, -1);
  wedge_posvol_flag.resize(num_wedges_all, true);

  int wedgeid = 0;
  int iall = 0, iown = 0, ighost = 0;
  for (auto const & c : cells()) {
    std::vector<Entity_ID> cfaces;
    std::vector<int> cfdirs;
    cell_get_faces_and_dirs(c, &cfaces, &cfdirs);

    Entity_ID_List::iterator itf = cfaces.begin();
    std::vector<int>::iterator itfd = cfdirs.begin();
    while (itf != cfaces.end()) {
      Entity_ID f = *itf;
      int fdir = *itfd;

      Entity_ID_List fedges;
      std::vector<int> fedirs;
      face_get_edges_and_dirs(f, &fedges, &fedirs);

      Entity_ID_List::iterator ite = fedges.begin();
      std::vector<int>::iterator ited = fedirs.begin();
      while (ite != fedges.end()) {
        Entity_ID e = *ite;
        int edir = *ited;

        int dir = fdir*edir;
        Entity_ID enode[2], ewedge[2];
        edge_get_nodes(e, &(enode[0]), &(enode[1]));
        for (int i = 0; i < 2; i++) {
          Entity_ID n = enode[i];
          wedge_node_id[wedgeid] = n;
          wedge_edge_id[wedgeid] = e;
          wedge_face_id[wedgeid] = f;
          wedge_cell_id[wedgeid] = c;
          wedge_parallel_type[wedgeid] = entity_get_ptype(Entity_kind::CELL, c);
          cell_wedge_ids[c].push_back(wedgeid);
          node_wedge_ids[n].push_back(wedgeid);

          wedgeids_all_[iall++] = wedgeid;
          if (wedge_parallel_type[wedgeid] == Parallel_type::OWNED)
            wedgeids_owned_[iown++] = wedgeid;
          else
            wedgeids_ghost_[ighost++] = wedgeid;


          // Whethter the fixed ordering of coordinates (n,e,f,c)
          // gives a +ve or -ve volume depends also on cell dimension
          // (2D or 3D). If all dirs (cell-to-face) and (face-to-edge)
          // are +ve, the triangle formed by n,e,c where n is node
          // point 0 of edge, will give a +ve area. On the other hand, the tet
          // formed by n,e,f,c will give a -ve area because triangle
          // n,e,f will point out of the cell.  Here n is node point
          // 0 of the edge, e is the edge center, f is the face center
          // and c is the cell center

          if (cell_dimension() == 2)
            wedge_posvol_flag[wedgeid] = !((i == 0)^(dir == 1));
          else if (cell_dimension() == 3)
            wedge_posvol_flag[wedgeid] = !((i == 0)^(dir == -1));

          ewedge[i] = wedgeid;

          // See if any of the other wedges attached to the node
          // shares the same edge, face and node but is in the
          // adjacent cell. This is called the opposite wedge

          Entity_ID_List::iterator itw = node_wedge_ids[n].begin();
          bool found = false;
          while (!found && itw != node_wedge_ids[n].end()) {
            Entity_ID w2 = *itw;
            if (w2 != wedgeid &&
                wedge_node_id[w2] == n && wedge_edge_id[w2] == e &&
                wedge_face_id[w2] == f && wedge_cell_id[w2] != c) {
              found = true;
              wedge_opp_wedge_id[wedgeid] = w2;
              wedge_opp_wedge_id[w2] = wedgeid;
            }
            ++itw;
          }
          wedgeid++;
        }

        // The two wedges associated with this edge, face and cell
        // are called adjacent wedges

        wedge_adj_wedge_id[ewedge[0]] = ewedge[1];
        wedge_adj_wedge_id[ewedge[1]] = ewedge[0];

        ++ite;
        ++ited;
      }  // while (ite != fedges.end())

      ++itf;
      ++itfd;
    }  // while (itf != cfaces.end())
  }  // for (int c = 0;....)

  wedge_info_cached = true;
}  // cache_wedge_info


void Mesh::cache_corner_info() const {

  if (corner_info_cached) return;

  if (!wedge_info_cached)
    cache_wedge_info();

  int ncells_owned = num_cells<Parallel_type::OWNED>();
  int ncells_ghost = num_cells<Parallel_type::GHOST>();
  int ncells = ncells_owned + ncells_ghost;

  int nnodes_owned = num_nodes<Parallel_type::OWNED>();
  int nnodes_ghost = num_nodes<Parallel_type::GHOST>();
  int nnodes = nnodes_owned + nnodes_ghost;
  
  cell_corner_ids.resize(ncells);
  node_corner_ids.resize(nnodes);

  int num_corners_all = 0;
  int num_corners_owned = 0;
  int num_corners_ghost = 0;

  for (int c = 0; c < ncells; c++) {
    std::vector<Entity_ID> cnodes;
    cell_get_nodes(c, &cnodes);
    cell_corner_ids[c].reserve(cnodes.size());

    num_corners_all += cnodes.size();  // as many corners as there are nodes in cell
    Parallel_type ptype = entity_get_ptype(Entity_kind::CELL, c);
    if (ptype == Parallel_type::OWNED)
      num_corners_owned += cnodes.size();
    else
      num_corners_ghost += cnodes.size();
  }

  cornerids_owned_.resize(num_corners_owned);
  cornerids_ghost_.resize(num_corners_ghost);
  cornerids_all_.resize(num_corners_all);
  corner_wedge_ids.resize(num_corners_all);
  corner_node_id.resize(num_corners_all, -1);
  corner_cell_id.resize(num_corners_all, -1);
  corner_parallel_type.resize(num_corners_all, Parallel_type::ALL);

  int cornerid = 0;
  int iall = 0, iown = 0, ighost = 0;
  for (int c = 0; c < ncells; c++) {
    std::vector<Entity_ID> cnodes;
    cell_get_nodes(c, &cnodes);

    Entity_ID_List::iterator itn = cnodes.begin();
    while (itn != cnodes.end()) {
      Entity_ID n = *itn;

      cell_corner_ids[c].push_back(cornerid);
      node_corner_ids[n].push_back(cornerid);
      corner_node_id[cornerid] = n;
      corner_cell_id[cornerid] = c;
      corner_parallel_type[cornerid] = entity_get_ptype(Entity_kind::CELL, c);

      cornerids_all_[iall++] = cornerid;
      if (corner_parallel_type[cornerid] == Parallel_type::OWNED)
        cornerids_owned_[iown++] = cornerid;
      else
        cornerids_ghost_[ighost++] = cornerid;

      std::vector<Entity_ID> nwedges;
      node_get_wedges(n, Parallel_type::ALL, &nwedges);
      Entity_ID_List::iterator itw = nwedges.begin();
      while (itw != nwedges.end()) {
        Entity_ID w = *itw;
        Entity_ID c2 = wedge_get_cell(w);
        if (c == c2) {
          corner_wedge_ids[cornerid].push_back(w);
          wedge_corner_id[w] = cornerid;
        }

        ++itw;
      }  // while (itw != nwedges.end())

      ++cornerid;
      ++itn;
    }  // while (itn != cnodes.end())
  }  // for (int c = ...)

  corner_info_cached = true;
}  // cache_corner_info


void Mesh::cache_extra_variables() {
  if (faces_requested) {
    cache_cell2face_info();
    cache_face2cell_info();
  }
  if (edges_requested) {
    cache_face2edge_info();
    cache_cell2edge_info();
  }

  if (wedges_requested) cache_wedge_info();
  if (corners_requested) cache_corner_info();

  if (faces_requested) compute_face_geometric_quantities();
  if (edges_requested) compute_edge_geometric_quantities();
  compute_cell_geometric_quantities();
  if (corners_requested) compute_corner_geometric_quantities();
  if (wedges_requested) compute_wedge_geometric_quantities();
}


// Partition the mesh on this compute node into submeshes or tiles
// Each tile contains only cell indices that are owned by it - no
// halo/ghost cell indices

void Mesh::build_tiles() {

  std::vector<std::vector<int>> partitions(num_tiles_);

  // #ifdef HAVE_METIS
  //  get_partitioning_with_metis(num_tiles, &partitions);
  // #elif HAVE_ZOLTAN
  //  get_partitioning_with_zoltan(num_tiles, &partitions);
  // #else
  std::cerr << "No partitioner defined - " <<
      "subdividing cell index space into equal parts\n";
  
  int ncells = num_cells<Parallel_type::OWNED>();
  int ncells_per_tile =
      static_cast<int> (static_cast<double>(ncells)/num_tiles_ + 0.5);
  int index = 0;
  for (int i = 0; i < num_tiles_; ++i) {
    partitions[i].resize(ncells_per_tile);  // excludes MPI ghosts

    // Owned cells in this partition
    for (int j = 0; j < ncells_per_tile && j < ncells; ++j)
      partitions[i][j] = index++;
  }

  // Put any leftover cells in the last partition

  for (int j = index; j < ncells; ++j)
    partitions[num_tiles_-1][j] = j;

  // #endif

  // Make the tiles and store shared pointers to them

  for (int i = 0; i < num_tiles_; ++i)
    meshtiles[i] = make_meshtile(*this, partitions[i],
                                 faces_requested, edges_requested,
                                 wedges_requested, corners_requested);
}



Entity_ID Mesh::entity_get_parent(const Entity_kind kind, const Entity_ID entid)
    const {
  Errors::Message mesg("Parent/daughter entities not enabled in this framework.");
  Exceptions::Jali_throw(mesg);
}


unsigned int Mesh::cell_get_num_faces(const Entity_ID cellid) const {

#if CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //
  if (!cell2face_info_cached) cache_cell2face_info();

  return cell_face_ids[cellid].size();

#else

  // Non-cached version

  Entity_ID_List cfaceids;
  std::vector<int> cfacedirs;

  cell_get_faces_and_dirs_internal(cellid, &cfaceids, &cfacedirs, false);

  return cfaceids.size();

#endif

}


void Mesh::cell_get_faces_and_dirs(const Entity_ID cellid,
                                   Entity_ID_List *faceids,
                                   std::vector<int> *face_dirs,
                                   const bool ordered) const {

#if CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //

  if (!cell2face_info_cached) cache_cell2face_info();

  if (ordered) {
    cell_get_faces_and_dirs_internal(cellid, faceids, face_dirs, ordered);
  } else {
    Entity_ID_List &cfaceids = cell_face_ids[cellid];

    *faceids = cfaceids;  // copy operation

    if (face_dirs) {
      std::vector<int> &cfacedirs = cell_face_dirs[cellid];
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

void Mesh::face_get_cells(const Entity_ID faceid, const Parallel_type ptype,
                          Entity_ID_List *cellids) const {

#if CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //

  if (!face2cell_info_cached) cache_face2cell_info();


  cellids->clear();

  switch (ptype) {
  case Parallel_type::ALL:
    for (int i = 0; i < 2; i++)
      if (face_cell_ptype[faceid][i] != Parallel_type::PTYPE_UNKNOWN)
        cellids->push_back(face_cell_ids[faceid][i]);
    break;
  case Parallel_type::OWNED:
    for (int i = 0; i < 2; i++)
      if (face_cell_ptype[faceid][i] == Parallel_type::OWNED)
        cellids->push_back(face_cell_ids[faceid][i]);
    break;
  case Parallel_type::GHOST:
    for (int i = 0; i < 2; i++)
      if (face_cell_ptype[faceid][i] == Parallel_type::GHOST)
        cellids->push_back(face_cell_ids[faceid][i]);
    break;
  }

#else

  //
  // Non-cached version
  //

  Entity_ID_List fcells;

  face_get_cells_internal(faceid, Parallel_type::ALL, &fcells);

  cellids->clear();

  switch (ptype) {
  case Parallel_type::ALL:
    for (int i = 0; i < fcells.size(); i++)
      if (entity_get_ptype(Entity_kind::CELL, fcells[i]) !=
          Parallel_type::PTYPE_UNKNOWN)
        cellids->push_back(fcells[i]);
    break;
  case Parallel_type::OWNED:
    for (int i = 0; i < fcells.size(); i++)
      if (entity_get_ptype(Entity_kind::CELL, fcells[i]) ==
          Parallel_type::OWNED)
        cellids->push_back(fcells[i]);
    break;
  case Parallel_type::GHOST:
    for (int i = 0; i < fcells.size(); i++)
      if (entity_get_ptype(Entity_kind::CELL, fcells[i]) ==
          Parallel_type::GHOST)
        cellids->push_back(fcells[i]);
    break;
  }

#endif

}


void Mesh::face_get_edges_and_dirs(const Entity_ID faceid,
                                   Entity_ID_List *edgeids,
                                   std::vector<int> *edge_dirs,
                                   const bool ordered) const {

#if CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //

  if (!face2edge_info_cached) cache_face2edge_info();

  *edgeids = face_edge_ids[faceid];  // copy operation

  if (edge_dirs) {
    std::vector<int> &fedgedirs = face_edge_dirs[faceid];
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

#if CACHE_VARS != 0

  //
  // Cached version - turn off for profiling or to save memory
  //

  if (!face2edge_info_cached) cache_face2edge_info();
  if (!cell2edge_info_cached) cache_cell2edge_info();

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
  std::vector<int> fedgedirs;

  face_get_edges_and_dirs(faceid,&fedgeids,&fedgedirs,true);
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


void Mesh::cell_get_edges (const Entity_ID cellid, 
                           Entity_ID_List *edgeids) const {


#if CACHE_VARS != 0

  //
  // Cached version - turn off for profiling
  //

  if (!cell2edge_info_cached) cache_cell2edge_info();

  Entity_ID_List &cedgeids = cell_edge_ids[cellid];

  *edgeids = cell_edge_ids[cellid]; // copy operation

#else

  //
  // Non-cached version
  //

  cell_get_edges_internal(cellid, edgeids);

#endif

} // Mesh::cell_get_edges


void Mesh::cell_2D_get_edges_and_dirs (const Entity_ID cellid,
                                       Entity_ID_List *edgeids,
                                       std::vector<int> *edgedirs) const {


#if CACHE_VARS != 0

  //
  // Cached version - turn off for profiling
  //

  if (!cell2edge_info_cached) cache_cell2edge_info();

  *edgeids = cell_edge_ids[cellid];  // copy operation
  *edgedirs = cell_2D_edge_dirs[cellid];

#else

  //
  // Non-cached version
  //

  cell_2D_get_edges_and_dirs_internal(cellid, edgeids, edgedirs);

#endif

}  // Mesh::cell_get_edges_and_dirs


void Mesh::cell_get_wedges(const Entity_ID cellid,
                           Entity_ID_List *wedgeids) const {
  assert(wedges_requested);
  if (!wedge_info_cached) cache_wedge_info();

  int nwedges = cell_wedge_ids[cellid].size();
  wedgeids->resize(nwedges);
  std::copy(cell_wedge_ids[cellid].begin(), cell_wedge_ids[cellid].end(),
            wedgeids->begin());
}


void Mesh::cell_get_corners(const Entity_ID cellid,
                            Entity_ID_List *cornerids) const {
  assert(corners_requested);
  if (!corner_info_cached) cache_corner_info();

  int ncorners = cell_corner_ids[cellid].size();
  cornerids->resize(ncorners);
  std::copy(cell_corner_ids[cellid].begin(), cell_corner_ids[cellid].end(),
            cornerids->begin());
}


Entity_ID Mesh::cell_get_corner_at_node(const Entity_ID cellid,
                                        const Entity_ID nodeid) const {
  assert(corners_requested);
  if (!corner_info_cached) cache_corner_info();

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

void Mesh::corner_get_wedges(const Entity_ID cornerid,
                             Entity_ID_List *wedgeids) const {
  assert(corners_requested);
  if (!corner_info_cached) cache_corner_info();

  int nwedges = corner_wedge_ids[cornerid].size();
  (*wedgeids).resize(nwedges);
  std::copy(corner_wedge_ids[cornerid].begin(),
            corner_wedge_ids[cornerid].end(), wedgeids->begin());
}


void Mesh::node_get_wedges(const Entity_ID nodeid, Parallel_type ptype,
                           Entity_ID_List *wedgeids) const {

  assert(wedges_requested);
  if (!wedge_info_cached) cache_wedge_info();

  switch (ptype) {
    case Parallel_type::ALL:
      wedgeids->resize(node_wedge_ids[nodeid].size());
      std::copy(node_wedge_ids[nodeid].begin(), node_wedge_ids[nodeid].end(),
                wedgeids->begin());
      break;
    default:
      wedgeids->clear();
      Entity_ID_List::iterator it = node_wedge_ids[nodeid].begin();
      while (it != node_wedge_ids[nodeid].end()) {
        Entity_ID w = *it;
        if (wedge_parallel_type[w] == ptype)
          wedgeids->push_back(w);
        ++it;
      }
      break;
  }

}


void Mesh::node_get_corners(const Entity_ID nodeid, Parallel_type ptype,
                            Entity_ID_List *cornerids) const {

  assert(corners_requested);
  if (!corner_info_cached) cache_corner_info();

  switch (ptype) {
    case Parallel_type::ALL:
      cornerids->resize(node_corner_ids[nodeid].size());
      std::copy(node_corner_ids[nodeid].begin(), node_corner_ids[nodeid].end(),
                cornerids->begin());
      break;
    default:
      cornerids->clear();
      Entity_ID_List::iterator it = node_corner_ids[nodeid].begin();
      while (it != node_corner_ids[nodeid].end()) {
        Entity_ID cn = *it;
        if (corner_parallel_type[cn] == ptype)
          cornerids->push_back(cn);
        ++it;
      }
      break;
  }

}


int Mesh::compute_cell_geometric_quantities() const {

  int ncells = num_cells<Parallel_type::ALL>();

  cell_volumes.resize(ncells);
  cell_centroids.resize(ncells);
  for (int i = 0; i < ncells; i++) {
    double volume;
    JaliGeometry::Point centroid(spacedim);

    compute_cell_geometry(i, &volume, &centroid);

    cell_volumes[i] = volume;
    cell_centroids[i] = centroid;
  }

  cell_geometry_precomputed = true;

  return 1;

}  // Mesh::compute_cell_geometric_quantities



int Mesh::compute_face_geometric_quantities() const {

  if (space_dimension() == 3 && cell_dimension() == 2) {
    // need cell centroids to compute normals

    if (!cell_geometry_precomputed)
      compute_cell_geometric_quantities();
  }

  int nfaces = num_faces<Parallel_type::ALL>();

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

  int nedges = num_edges<Parallel_type::ALL>();

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


int Mesh::compute_wedge_geometric_quantities() const {

  wedge_volumes.resize(num_wedges());

  // Cannot use resize for the facet normals because we cannot tell
  // the initialization operator what the dimensionality of the points
  // will be, unless we can write a statement like this
  //
  //  wedge_facet_normals0.resize(nwedges, Point(3));
  //
  // It is also unclear that this is good to do because it will cause
  // a copy operator to be triggered for each element (?)

  wedge_facet_normals0.reserve(num_wedges());
  wedge_facet_normals1.reserve(num_wedges());

  for (auto const & w : wedges()) {
    JaliGeometry::Point facet_normal0(spacedim), facet_normal1(spacedim);

    compute_wedge_geometry(w, &(wedge_volumes[w]),
                           &(facet_normal0), &(facet_normal1));
    wedge_facet_normals0.push_back(facet_normal0);
    wedge_facet_normals1.push_back(facet_normal1);
  }

  wedge_geometry_precomputed = true;

  return 1;
}

int Mesh::compute_corner_geometric_quantities() const {
  corner_volumes.resize(num_corners());
  for (auto const & c : corners())
    compute_corner_geometry(c, &(corner_volumes[c]));

  corner_geometry_precomputed = true;
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
    std::vector<int> fdirs;
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

    JaliGeometry::polyhed_get_vol_centroid(ccoords, nf, nfnodes,
            cfcoords, volume, centroid);
    return 1;

  } else if (celldim == 2) {

    std::vector<JaliGeometry::Point> ccoords;

    cell_get_coordinates(cellid, &ccoords);

    JaliGeometry::Point normal(spacedim);

    JaliGeometry::polygon_get_area_centroid_normal(ccoords, volume, centroid,
                                                     &normal);

    return 1;
  }

  return 0;
}  // Mesh::compute_cell_geometry


int Mesh::compute_face_geometry(const Entity_ID faceid, double *area,
        JaliGeometry::Point *centroid, JaliGeometry::Point *normal0,
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
    face_get_cells(faceid, Parallel_type::ALL, &cellids);

    for (int i = 0; i < cellids.size(); i++) {
      Entity_ID_List cellfaceids;
      std::vector<int> cellfacedirs;
      int dir = 1;

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
      face_get_cells(faceid, Parallel_type::ALL, &cellids);

      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        std::vector<int> cellfacedirs;
        int dir = 1;

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
      face_get_cells(faceid, Parallel_type::ALL, &cellids);

      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        std::vector<int> cellfacedirs;
        int dir = 1;

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


void Mesh::compute_wedge_geometry(Entity_ID const wedgeid,
                                 double *wedge_volume,
                                 JaliGeometry::Point *facet_normal0,
                                 JaliGeometry::Point *facet_normal1) const {

  // First record the flag which indicates if the default ordering
  // gives us positive volumes for the wedges if we assume n, e, c for
  // 2D and n, e, f, c for 3D (n - node point, e - edge center, f - face
  // center, c - cell center. This flag was recorded during the
  // construction of wedges using purely topological info and so is robust

  bool posvol = wedge_posvol_flag[wedgeid];

  if (celldim == 3) {

    std::vector<JaliGeometry::Point> wcoords;

    // Get vertex coordinates of wedge
    //
    // These are always - node coordinate, edge center, face center,
    // cell center

    wedge_get_coordinates(wedgeid, &wcoords);

    // vector from edge center to node
    JaliGeometry::Point vec0 = wcoords[0]-wcoords[1];

    // vector from edge center to face center
    JaliGeometry::Point vec1 = wcoords[2]-wcoords[1];

    // vector from edge center to cell center
    JaliGeometry::Point vec2 = wcoords[3]-wcoords[1];


    // Area weighted normal to the triangular facet formed by node
    // coordinate, edge center and face center such that the normal is
    // pointing out of the wedge and cell

    *facet_normal0 = posvol ? 0.5*(vec0^vec1) : -0.5*(vec0^vec1);

    // Volume of wedge adjusted for +ve volume

    *wedge_volume = posvol ? (vec1^vec0)*vec2/6.0 : -(vec1^vec0)*vec2/6.0;

    // Compute normal of triangular facet formed by edge center, face
    // center and zone center

    *facet_normal1 = posvol ? 0.5*(vec1^vec2) : -0.5*(vec1^vec2);

  } else if (celldim == 2) {

    std::vector<JaliGeometry::Point> wcoords;

    // Get vertex coordinates of wedge
    //
    // These are always - node coordinate, edge/face center, cell center

    wedge_get_coordinates(wedgeid, &wcoords);

    // vector from edge/face center to node
    JaliGeometry::Point vec0 = wcoords[0]-wcoords[1];

    // vector from cell center to edge/face center
    JaliGeometry::Point vec1 = wcoords[2]-wcoords[1];

    // length weighted normal to the segment formed
    // by node coordinate and edge/face center

    JaliGeometry::Point normal(JaliGeometry::Point(-vec0[1], vec0[0]));

    *facet_normal0 = posvol ? normal : -normal;

    // Area of wedge is 1/2 of the cross-product of vec1 and vec0

    JaliGeometry::Point cpvec = (vec1^vec0)/2.0;
    *wedge_volume = posvol ? cpvec[0] : -cpvec[0];

    // length weighted normal to the segment formed
    // by cell center and edge/face center

    normal.set(vec1[1], -vec1[0]);

    // Adjust sign of the normal to ensure that its pointing out
    // from the node.

    *facet_normal1 = posvol ? normal : -normal;
  } else if (celldim == 1) {
    // unclear what the definitions of these normals are
  }

}  // Compute wedge geometry


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

  if (!cell_geometry_precomputed) {
    compute_cell_geometric_quantities();
    return cell_volumes[cellid];
  } else {
    if (recompute) {
      double volume;
      JaliGeometry::Point centroid(spacedim);
      compute_cell_geometry(cellid, &volume, &centroid);
      return volume;
    }
    else
      return cell_volumes[cellid];
  }
}

// Area/length of face

double Mesh::face_area(const Entity_ID faceid, const bool recompute) const {

  ASSERT(faces_requested);

  if (!face_geometry_precomputed) {
    compute_face_geometric_quantities();
    return face_areas[faceid];
  } else {
    if (recompute) {
      double area;
      JaliGeometry::Point centroid(spacedim);
      JaliGeometry::Point normal0(spacedim), normal1(spacedim);
      compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
      return area;
    }
    else
      return face_areas[faceid];
  }
}

// Length of an edge

double Mesh::edge_length(const Entity_ID edgeid, const bool recompute) const {

  ASSERT(edges_requested);

  if (!edge_geometry_precomputed) {
    compute_edge_geometric_quantities();
    return edge_lengths[edgeid];
  } else {
    if (recompute) {
      double length;
      JaliGeometry::Point vector(spacedim), centroid(spacedim);
      compute_edge_geometry(edgeid, &length, &vector, &centroid);
      return length;
    }
    else
      return edge_lengths[edgeid];
  }
}

// Volume/Area of wedge

double Mesh::wedge_volume(const Entity_ID wedgeid, const bool recompute) const {

  if (!wedge_geometry_precomputed) {
    compute_wedge_geometric_quantities();
    return wedge_volumes[wedgeid];
  } else {
    if (recompute) {
      double volume;
      JaliGeometry::Point facet_normal0, facet_normal1;
      compute_wedge_geometry(wedgeid, &volume, &facet_normal0, &facet_normal1);
      return volume;
    }
    else
      return wedge_volumes[wedgeid];
  }
}

// Corner volume

double Mesh::corner_volume(const Entity_ID cornerid, const bool recompute) const {
  double corner_volume = 0.0;

  if (!corner_geometry_precomputed) {
    compute_corner_geometric_quantities();
    return corner_volumes[cornerid];
  } else {
    if (recompute) {
      double volume;
      compute_corner_geometry(cornerid, &volume);
      return volume;
    }
    else
      return corner_volumes[cornerid];
  }
}

// Centroid of cell

JaliGeometry::Point Mesh::cell_centroid(const Entity_ID cellid,
                                        const bool recompute) const {

  if (!cell_geometry_precomputed) {
    compute_cell_geometric_quantities();
    return cell_centroids[cellid];
  } else {
    if (recompute) {
      double volume;
      JaliGeometry::Point centroid(spacedim);
      compute_cell_geometry(cellid, &volume, &centroid);
      return centroid;
    }
    else
      return cell_centroids[cellid];
  }

}

// Centroid of face

JaliGeometry::Point Mesh::face_centroid(const Entity_ID faceid,
                                        const bool recompute) const {

  ASSERT(faces_requested);

  if (!face_geometry_precomputed) {
    compute_face_geometric_quantities();
    return face_centroids[faceid];
  } else {
    if (recompute) {
      double area;
      JaliGeometry::Point centroid(spacedim);
      JaliGeometry::Point normal0(spacedim), normal1(spacedim);
      compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
      return centroid;
    }
    else
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


JaliGeometry::Point Mesh::face_normal (const Entity_ID faceid,
                                         const bool recompute,
                                         const Entity_ID cellid,
                                         int *orientation) const {

  ASSERT(faces_requested);

  JaliGeometry::Point normal0(spacedim);
  JaliGeometry::Point normal1(spacedim);

  if (!face_geometry_precomputed) {
    compute_face_geometric_quantities();

    normal0 = face_normal0[faceid];
    normal1 = face_normal1[faceid];
  } else {
    if (recompute) {
      double area;
      JaliGeometry::Point centroid(spacedim);

      compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
    } else {
      normal0 = face_normal0[faceid];
      normal1 = face_normal1[faceid];
    }
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
    std::vector<int> face_dirs;

    cell_get_faces_and_dirs(cellid, &faceids, &face_dirs);

    int nf = faceids.size();
    bool found = false;
    int dir = 1;
    for (int i = 0; i < nf; i++)
      if (faceids[i] == faceid) {
        dir = face_dirs[i];
        found = true;
        break;
      }

    ASSERT(found);

    if (orientation) *orientation = dir;
    if (dir == 1) {
      ASSERT(L22(normal0) != 0.0);
      return normal0;              // Copy to output
    } else {
      ASSERT(L22(normal1) != 0.0);
      return normal1;              // Copy to output
    }
  }

  return normal0;      // Copy to output
}


// Direction vector of edge

JaliGeometry::Point Mesh::edge_vector (const Entity_ID edgeid,
                                         const bool recompute,
                                         const Entity_ID pointid,
                                         int *orientation) const {

  ASSERT(edges_requested);

  JaliGeometry::Point evector(spacedim), ecenter(spacedim);
  JaliGeometry::Point& evector_ref = evector;  // to avoid extra copying

  if (!edge_geometry_precomputed)
    compute_edge_geometric_quantities();

  if (recompute) {
    double length;
    compute_edge_geometry(edgeid, &length, &evector, &ecenter);
    // evector_ref already points to evector
  }
  else
    evector_ref = edge_vectors[edgeid];

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

  edge_get_nodes(edgeid, &p0, &p1);
  node_get_coordinates(p0, &xyz0);
  node_get_coordinates(p1, &xyz1);
  return (xyz0+xyz1)/2.0;
}


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

  if (celldim == 3) {
    wcoords->resize(4);  // wedges are tets in 3D cells
    Entity_ID n = wedge_get_node(wedgeid);
    node_get_coordinates(n, &((*wcoords)[0]));

    Entity_ID e = wedge_get_edge(wedgeid);
    (*wcoords)[1] = edge_centroid(e);

    Entity_ID f = wedge_get_face(wedgeid);
    (*wcoords)[2] = face_centroid(f);

    Entity_ID c = wedge_get_cell(wedgeid);
    (*wcoords)[3] = cell_centroid(c);
  } else if (celldim == 2) {
    wcoords->resize(3);  // wedges are tris in 2D cells
    Entity_ID n = wedge_get_node(wedgeid);
    node_get_coordinates(n, &((*wcoords)[0]));

    Entity_ID f = wedge_get_face(wedgeid);
    (*wcoords)[1] = face_centroid(f);

    Entity_ID c = wedge_get_cell(wedgeid);
    (*wcoords)[2] = cell_centroid(c);
  }

  // If caller has requested that coordinates be ordered such that
  // they will give a +ve volume AND the wedge has been flag as not
  // giving a +ve volume with its natural coordinate order (node
  // point, edge center, face center, cell center), then swap the edge
  // and face centers (in 3D) or face and cell centers (in 2D)

  if (posvol_order && !wedge_posvol_flag[wedgeid])
    std::swap((*wcoords)[1], (*wcoords)[2]);

} // wedge_get_coordinates


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
                                             const int which_facet,
                                             const bool recompute) const {

  assert(wedges_requested);

  JaliGeometry::Point normal(spacedim);

  if (!wedge_geometry_precomputed) {
    compute_wedge_geometric_quantities();

    return which_facet ?
        wedge_facet_normals1[wedgeid] : wedge_facet_normals0[wedgeid];
  } else {
    if (recompute) {
      double volume;
      JaliGeometry::Point facet_normal0(spacedim), facet_normal1(spacedim);

      compute_wedge_geometry(wedgeid, &volume, &facet_normal0, &facet_normal1);
      return which_facet ? facet_normal1 : facet_normal0;
    }
    else
      return which_facet ?
          wedge_facet_normals1[wedgeid] : wedge_facet_normals0[wedgeid];
  }

} // wedge_facet_normal


// Triangular facets describing a Corner in 3D

void
Mesh::corner_get_facetization(const Entity_ID cornerid,
                              std::vector<JaliGeometry::Point> *pointcoords,
                              std::vector<std::array<Entity_ID, 3>> *facetpoints) const {

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
      if (point_entity_list[idxe] == std::pair<Entity_ID, Entity_kind>(f, Entity_kind::FACE))
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
      if (point_entity_list[idxf] == std::pair<Entity_ID, Entity_kind>(e, Entity_kind::EDGE))
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
Mesh::corner_get_facetization (const Entity_ID cornerid,
                               std::vector<JaliGeometry::Point> *pointcoords,
                               std::vector<std::array<Entity_ID, 2>> *facetpoints) const {

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
  facetpoints->push_back({{2, 4}});

}  // corner get facetization for 3D



// list of points describing a Corner in 2D in ccw direction and in no
// particular order in 3D

void
Mesh::corner_get_coordinates(const Entity_ID cornerid, 
                             std::vector<JaliGeometry::Point> *pointcoords)
    const {


  Entity_ID_List cwedges;
  corner_get_wedges(cornerid, &cwedges);

  pointcoords->clear();

  if (celldim == 2) {
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
    JaliGeometry::Point vec0 = ccen-(*pointcoords)[0];

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
  }

}  // corner_get_points



// Is there a set with this name and entity type

bool Mesh::valid_set_name(std::string name, Entity_kind kind) const
{

  if (!geometric_model_) {
    Errors::Message mesg("Mesh sets not enabled because mesh was created without reference to a geometric model");
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

      if (rgn->type() == JaliGeometry::COLORFUNCTION) return true;

      // For regions of type Labeled set, extract some more info and verify

      if (rgn->type() == JaliGeometry::LABELEDSET) {
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

      if (kind == Entity_kind::CELL && (rdim >= celldim || rdim == 0)) return true;

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


bool Mesh::point_in_cell(const JaliGeometry::Point &p, const Entity_ID cellid)
    const {

  std::vector<JaliGeometry::Point> ccoords;

  if (celldim == 3) {

    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine

    int nf;
    Entity_ID_List faces;
    std::vector<unsigned int> nfnodes;
    std::vector<int> fdirs;
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

  }

  return false;
}


std::string Mesh::cell_type_to_name(const Cell_type type) {

  switch (type) {
    case TRI:
      return "triangle";
    case QUAD:
      return "quad";
    case POLYGON:
      return "polygon";
    case TET:
      return "tetrahedron";
    case PYRAMID:
      return "pyramid";
    case PRISM:
      return "prism";
    case HEX:
      return "hexahedron";
    case POLYHED:
      return "polyhedron";
    default:
      return "unknown";
  }
}

}  // close namespace Jali
