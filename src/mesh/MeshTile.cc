//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#include <vector>
#include <algorithm>
#include <memory>

#include "mpi.h"

#include "MeshDefs.hh"
#include "Mesh.hh"
#include "MeshTile.hh"

namespace Jali {

/*! 
  @brief Constructor for MeshTile
  
  Mesh tiles are small groupings of mesh cells on a compute node and
  are typically obtained by partitioning the full mesh on a compute
  node. They are lightweight structures that are merely lists of
  cells that are in a partition. Mesh tiles include ghost nodes at
  the compute node level (i.e. MPI level).
  
  Meshtiles are a way for the work on a mesh to be broken up into
  manageable, parallely executable chunks.  
*/


// Constructor

MeshTile::MeshTile(Mesh const & parent_mesh,
                   std::vector<Entity_ID> const meshcells,
                   bool const request_faces, bool const request_edges,
                   bool const request_wedges, bool const request_corners) :
    mesh_(parent_mesh) {

  int ncells = meshcells.size(), nowned, nghost;
  cellids_all_.resize(ncells);
  cellids_owned_.resize(ncells);
  cellids_ghost_.resize(ncells);

  ncells = 0; nowned = 0; nghost = 0;
  for (auto const & c : meshcells) {
    Parallel_type ptype = mesh_.entity_get_ptype(CELL, c);
    if (ptype == OWNED) {
      cellids_owned_[nowned++] = c;
    }
    else if (ptype == GHOST) {
      cellids_ghost_[nghost++] = c;
    }
    else 
      std::cerr << "Invalid parallel type for mesh entity " << ptype <<
          std::endl;
    cellids_all_[ncells++] = c;
  }
  cellids_owned_.resize(nowned);
  cellids_ghost_.resize(nghost);

  // Make a list of nodeids in the tile taking care to place the
  // off-processor ghost nodes at the end

  for (auto const & c : meshcells) {
    std::vector<int> cnodes;
    mesh_.cell_get_nodes(c, &cnodes);
    for (auto n : cnodes) {
      if (mesh_.entity_get_ptype(NODE, n) == OWNED) {
        if (std::find(nodeids_owned_.begin(), nodeids_owned_.end(), n) ==
            nodeids_owned_.end())
          nodeids_owned_.emplace_back(n);
      }
      else {
        if (std::find(nodeids_ghost_.begin(), nodeids_ghost_.end(), n) ==
            nodeids_ghost_.end())
          nodeids_ghost_.emplace_back(n);
      }
    }
  }
  nodeids_all_ = nodeids_owned_;
  nodeids_all_.insert(nodeids_all_.end(),
                      nodeids_ghost_.begin(), nodeids_ghost_.end());

  // Make a list of faces, edges, wedges and corners similarly if requested

  if (request_faces) {    
    for (auto c : meshcells) {
      std::vector<int> cfaces;
      mesh_.cell_get_faces(c, &cfaces);
      for (auto f : cfaces) {
        if (mesh_.entity_get_ptype(FACE, f) == OWNED) {
          if (std::find(faceids_owned_.begin(), faceids_owned_.end(), f) ==
              faceids_owned_.end())
            faceids_owned_.emplace_back(f);
        }
        else {
          if (std::find(faceids_ghost_.begin(), faceids_ghost_.end(), f) ==
              faceids_ghost_.end())
            faceids_ghost_.emplace_back(f);
        }
      }
    }
    faceids_all_ = faceids_owned_;
    faceids_all_.insert(faceids_all_.end(),
                        faceids_ghost_.begin(), faceids_ghost_.end());
  }

  if (request_edges) {
    for (auto c : meshcells) {
      std::vector<int> cedges;
      mesh_.cell_get_edges(c, &cedges);
      for (auto e : cedges) {
        if (mesh_.entity_get_ptype(EDGE, e) == OWNED) {
          if (std::find(edgeids_owned_.begin(), edgeids_owned_.end(), e) ==
              edgeids_owned_.end())
            edgeids_owned_.emplace_back(e);
        }
        else {
          if (std::find(edgeids_ghost_.begin(), edgeids_ghost_.end(), e) ==
              edgeids_ghost_.end())
            edgeids_ghost_.emplace_back(e);
        }
      }
    }
    edgeids_all_ = edgeids_owned_;
    edgeids_all_.insert(edgeids_all_.end(),
                        edgeids_ghost_.begin(), edgeids_ghost_.end());
  }

  if (request_wedges) {
    for (auto c : meshcells) {
      std::vector<int> cwedges;
      mesh_.cell_get_wedges(c, &cwedges);
      for (auto w : cwedges) {
        if (mesh_.entity_get_ptype(WEDGE, w) == OWNED) {
          if (std::find(wedgeids_owned_.begin(), wedgeids_owned_.end(), w) ==
              wedgeids_owned_.end())
            wedgeids_owned_.emplace_back(w);
        }
        else {
          if (std::find(wedgeids_ghost_.begin(), wedgeids_ghost_.end(), w) ==
              wedgeids_ghost_.end())
            wedgeids_ghost_.emplace_back(w);
        }
      }
    }
    wedgeids_all_ = wedgeids_owned_;
    wedgeids_all_.insert(wedgeids_all_.end(),
                         wedgeids_all_.begin(), wedgeids_all_.end());
  }

  if (request_corners) {
    for (auto c : meshcells) {
      std::vector<int> ccorners;
      mesh_.cell_get_corners(c, &ccorners);
      for (auto cn : ccorners) {
        if (mesh_.entity_get_ptype(CORNER, cn) == OWNED) {
          if (std::find(cornerids_owned_.begin(), cornerids_owned_.end(), cn) ==
              cornerids_owned_.end())
            cornerids_owned_.emplace_back(cn);
        }
        else {
          if (std::find(cornerids_ghost_.begin(), cornerids_ghost_.end(), cn) ==
              cornerids_ghost_.end())
            cornerids_ghost_.emplace_back(cn);
        }
      }
    }
    cornerids_all_ = cornerids_owned_;
    cornerids_all_.insert(cornerids_all_.end(),
                          cornerids_ghost_.begin(), cornerids_ghost_.end());
  }
}  // MeshTile::MeshTile

// Standalone function to make a tile and return a pointer to it so
// that Mesh.hh can use a forward declaration of MeshTile and this
// function to create new tiles

std::shared_ptr<MeshTile> make_meshtile(Mesh const & parent_mesh,
                                        const std::vector<Entity_ID> & cells,
                                        const bool request_faces,
                                        const bool request_edges,
                                        const bool request_wedges,
                                        const bool request_corners) {
  return std::make_shared<MeshTile>(parent_mesh, cells,
                                    request_faces, request_edges,
                                    request_wedges, request_corners);
}



} // end namespace Jali

