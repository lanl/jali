//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#include <mpi.h>

#include <vector>
#include <algorithm>
#include <memory>

#include "MeshDefs.hh"
#include "Mesh.hh"
#include "MeshTile.hh"

namespace Jali {

/*! 
  @brief Constructor for MeshTile
  
  Mesh tiles are small groupings of mesh cells on a compute node and
  are typically obtained by partitioning the full mesh on a compute
  node. They are lightweight structures that are merely lists of
  cells that are in a partition.
  
  Meshtiles are a way for the work on a mesh to be broken up into
  manageable, parallely executable chunks.  

  Meshtiles can have a layer of halo or ghost cells for ease of
  computations.  Note that ghost cells of a meshtile may or may not be
  an MPI ghost (Parallel_type::GHOST)
*/


// Constructor - should we make this private and only allow Mesh to
// call it as a friend? MeshTile can send a reference to itself to the
// parent_mesh so that it can be added to the list of tiles

MeshTile::MeshTile(Mesh& parent_mesh,
                   std::vector<Entity_ID> const& meshcells_owned,
                   int const num_halo_layers,
                   bool const request_faces, bool const request_edges,
                   bool const request_wedges, bool const request_corners) :
    mesh_(parent_mesh),
    mytileid_(parent_mesh.tiles().size()) {

  cellids_owned_ = meshcells_owned;
  cellids_all_ = meshcells_owned;

  // Build up halos if requested

  if (num_halo_layers > 0) {
    for (int i = 0; i < num_halo_layers; ++i) {
      Entity_ID_List next_halo_layer;

      for (auto const& c : cellids_all_) {
        Entity_ID_List nbrs;
        mesh_.cell_get_node_adj_cells(c, Parallel_type::ALL, &nbrs);
        
        for (auto const& cnbr : nbrs) {
          // Check if the neighbor is already in the all cells list or
          // in this halo

          if (std::find(next_halo_layer.begin(), next_halo_layer.end(), cnbr)
              != next_halo_layer.end()) continue;
          if (std::find(cellids_all_.begin(), cellids_all_.end(), cnbr) !=
              cellids_all_.end()) continue;

          // Add neighbor to the next halo layer
          next_halo_layer.push_back(cnbr);
        }
      }

      cellids_all_.insert(cellids_all_.end(),
                          next_halo_layer.begin(), next_halo_layer.end());
      cellids_ghost_.insert(cellids_ghost_.end(),
                            next_halo_layer.begin(), next_halo_layer.end());
    }
  }


  for (auto const& c : cellids_owned_)
    mesh_.set_master_tile_ID_of_cell(c, mytileid_);

  // Make a list of nodeids in the tile

  for (auto const& c : cellids_owned_) {
    std::vector<int> cnodes;
    mesh_.cell_get_nodes(c, &cnodes);
    for (auto const& n : cnodes) {
      int tileid = mesh_.master_tile_ID_of_node(n);
      if (tileid == -1) {
        // Node not yet in any tile - so this tile owns it
        nodeids_owned_.emplace_back(n);
        mesh_.set_master_tile_ID_of_node(n, mytileid_);
      } else {
        // If node is owned by another tile put it in the ghost list
        if (tileid != mytileid_ &&
            std::find(nodeids_ghost_.begin(), nodeids_ghost_.end(), n) ==
            nodeids_ghost_.end())
          nodeids_ghost_.emplace_back(n);
      }
    }
  }

  for (auto const& c : cellids_ghost_) {
    std::vector<int> cnodes;
    mesh_.cell_get_nodes(c, &cnodes);
    for (auto const& n : cnodes) {
      int tileid = mesh_.master_tile_ID_of_node(n);  // may be -1 (unassigned)
      if (tileid != mytileid_ &&
          std::find(nodeids_ghost_.begin(), nodeids_ghost_.end(), n) ==
          nodeids_ghost_.end())
        nodeids_ghost_.emplace_back(n);
    }
  }

  nodeids_all_ = nodeids_owned_;
  nodeids_all_.insert(nodeids_all_.end(),
                      nodeids_ghost_.begin(), nodeids_ghost_.end());

  // Make a list of faces similarly if requested

  if (request_faces) {

    for (auto const& c : cellids_owned_) {
      std::vector<int> cfaces;
      mesh_.cell_get_faces(c, &cfaces);
      for (auto const& f : cfaces) {
        int tileid = mesh_.master_tile_ID_of_face(f);
        if (tileid == -1) {
          // face not yet in any tile - so this tile owns it
          faceids_owned_.emplace_back(f);
          mesh_.set_master_tile_ID_of_face(f, mytileid_);
        } else {
          // If face is owned by another tile put it in the ghost list
          if (tileid != mytileid_ &&
              std::find(faceids_ghost_.begin(), faceids_ghost_.end(), f) ==
              faceids_ghost_.end())
            faceids_ghost_.emplace_back(f);
        }
      }
    }

    for (auto const& c : cellids_ghost_) {
      std::vector<int> cfaces;
      mesh_.cell_get_faces(c, &cfaces);
      for (auto const& f : cfaces) {
        int tileid = mesh_.master_tile_ID_of_face(f);  // may be -1 (unassigned)
        if (tileid != mytileid_ &&
            std::find(faceids_ghost_.begin(), faceids_ghost_.end(), f) ==
            faceids_ghost_.end())
          faceids_ghost_.emplace_back(f);
      }
    }

    faceids_all_ = faceids_owned_;
    faceids_all_.insert(faceids_all_.end(),
                        faceids_ghost_.begin(), faceids_ghost_.end());
  }

  // Make a list of faces similarly if requested

  if (request_edges) {
    for (auto const& c : cellids_owned_) {
      std::vector<int> cedges;
      mesh_.cell_get_edges(c, &cedges);
      for (auto const& e : cedges) {
        int etileid = mesh_.master_tile_ID_of_edge(e);
        if (etileid == -1) {
          // Edge not yet in any tile - so this tile owns it
          edgeids_owned_.emplace_back(e);
          mesh_.set_master_tile_ID_of_edge(e, mytileid_);
        } else {
          // If edge is owned by another tile put it in the ghost list
          if (etileid != mytileid_ &&
              std::find(edgeids_ghost_.begin(), edgeids_ghost_.end(), e) ==
              edgeids_ghost_.end())
            edgeids_ghost_.emplace_back(e);
        }
      }
    }

    for (auto const& c : cellids_ghost_) {
      std::vector<int> cedges;
      mesh_.cell_get_edges(c, &cedges);
      for (auto const& e : cedges) {
        int tileid = mesh_.master_tile_ID_of_edge(e);  // may be -1 (unassigned)
        if (tileid != mytileid_ &&
            std::find(edgeids_ghost_.begin(), edgeids_ghost_.end(), e) ==
            edgeids_ghost_.end())
          edgeids_ghost_.emplace_back(e);
      }
    }

    edgeids_all_ = edgeids_owned_;
    edgeids_all_.insert(edgeids_all_.end(),
                        edgeids_ghost_.begin(), edgeids_ghost_.end());
  }

  if (request_wedges) {
    for (auto const& c : cellids_owned_) {
      std::vector<int> cwedges;
      mesh_.cell_get_wedges(c, &cwedges);
      for (auto const& w : cwedges)
        wedgeids_owned_.emplace_back(w);
    }

    for (auto const& c : cellids_ghost_) {
      std::vector<int> cwedges;
      mesh_.cell_get_wedges(c, &cwedges);
      for (auto const& w : cwedges)
        wedgeids_ghost_.emplace_back(w);
    }

    wedgeids_all_ = wedgeids_owned_;
    wedgeids_all_.insert(wedgeids_all_.end(),
                         wedgeids_ghost_.begin(), wedgeids_ghost_.end());
  }

  if (request_corners) {
    for (auto const& c : cellids_owned_) {
      std::vector<int> ccorners;
      mesh_.cell_get_corners(c, &ccorners);
      for (auto const& cn : ccorners)
        cornerids_owned_.emplace_back(cn);
    }

    for (auto const& c : cellids_ghost_) {
      std::vector<int> ccorners;
      mesh_.cell_get_corners(c, &ccorners);
      for (auto const& cn : ccorners)
        cornerids_ghost_.emplace_back(cn);
    }

    cornerids_all_ = cornerids_owned_;
    cornerids_all_.insert(cornerids_all_.end(),
                          cornerids_ghost_.begin(), cornerids_ghost_.end());
  }
}  // MeshTile::MeshTile

// Standalone function to make a tile and return a pointer to it so
// that Mesh.hh can use a forward declaration of MeshTile and this
// function to create new tiles

std::shared_ptr<MeshTile> make_meshtile(Mesh& parent_mesh,
                                        std::vector<Entity_ID> const& cells,
                                        int const num_halo_layers,
                                        bool const request_faces,
                                        bool const request_edges,
                                        bool const request_wedges,
                                        bool const request_corners) {
  if (parent_mesh.num_tiles() == 0)
    parent_mesh.init_tiles();

  // This is a less than optimal use of the make_shared function since
  // it involves two memory allocations but I am not able to do it in
  // the optimal way since the MeshTile constructor is private (to
  // ensure mesh tiles are constructed correctly through this function
  // and proper mesh<-->meshtile links are established) - even when I
  // make the std::make_shared_ptr class a friend of MeshTile, it
  // complains that the constructor is private

  auto tile =  std::make_shared<MeshTile>(parent_mesh, cells,
                                          num_halo_layers,
                                          request_faces,
                                          request_edges,
                                          request_wedges,
                                          request_corners);
  
  parent_mesh.add_tile(tile);
  return tile;
}

}  // end namespace Jali

