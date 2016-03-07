//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#ifndef _JALI_MESHTILE_H_
#define _JALI_MESHTILE_H_

#include <vector>
#include <algorithm>

#include "mpi.h"

#include "MeshDefs.hh"
#include "Mesh.hh"

namespace Jali {

/*! 
  @class MeshTile "MeshTile.hh"
  @brief Encapsulates notion of a meshtile (a group of cells)
  
  Mesh tiles are small groupings of mesh cells on a compute node and
  are typically obtained by partitioning the full mesh on a compute
  node. They are lightweight structures that are merely lists of cells
  that are in a partition. Mesh tiles include ghost nodes at the
  compute node level (i.e. MPI level) BUT TYPICALLY ARE NOT EXPECTED
  TO INCLUDE GHOST CELLS (although they can).
  
  Meshtiles are a way for the work on a mesh to be broken up into
  manageable, parallely executable chunks.
  
  **** IMPORTANT NOTE ABOUT CONSTANTNESS OF THIS CLASS ****
  Instantiating a const version of this class only guarantees that
  the underlying mesh topology and geometry does not change (the
  public interfaces conforms strictly to this definition). However,
  for purposes of memory savings we use lazy initialization and
  caching of face data, edge data, geometry quantities, columns etc.,
  which means that these data may still change. We also cannot
  initialize the cached quantities in the constructor since they
  depend on initialization of data structures in the derived class -
  however, the base class gets constructed before the derived class
  gets constructed so it is not possible without more obscure
  acrobatics. This is why some of the caching data declarations are
  declared with the keyword 'mutable' and routines that modify the
  mutable data are declared with a constant qualifier.
  
*/

class MeshTile {

 public:

  /// @brief Constructor

  MeshTile(Mesh const & parent_mesh, std::vector<Entity_ID> const meshcells,
           bool const request_faces = true, bool const request_edges = false,
           bool const request_wedges = false,
           bool const request_corners = false);

  /// @brief Copy Constructor - deleted

  MeshTile(MeshTile const &meshtile_in) = delete;

  /// @brief Assignment operator - deleted

  MeshTile & operator=(MeshTile const &meshtile_in) = delete;

  /// @brief Destructor

  ~MeshTile() {}

  /// @brief The mesh that this tile belongs to

  Mesh const & mesh() {
    return mesh_;
  }

  //
  // General mesh tile information
  // -------------------------
  //

  /*! 
    @brief Number of nodes of parallel type
    @tparam ptype Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> unsigned int num_nodes() const;

  /*! 
    @brief Number of edges of parallel type
    @tparam ptype Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> unsigned int num_edges() const;

  /*! 
    @brief Number of faces of parallel type
    @tparam ptype Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> unsigned int num_faces() const;

  /*! 
    @brief Number of wedges of parallel type
    @tparam ptype Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> unsigned int num_wedges() const;

  /*! 
    @brief Number of corners of parallel type
    @tparam ptype Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> unsigned int num_corners() const;

  /*! 
    @brief Number of nodes of parallel type
    @tparam ptype Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> unsigned int num_cells() const;

  /*! 
    @brief Node list
    @tparam ptype   Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL>
  std::vector<Entity_ID> const & nodes() const;

  /*! 
    @brief Edge list
    @tparam ptype   Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> std::vector<Entity_ID>
  const & edges() const;

  /*! 
    @brief Face list
    @tparam ptype   Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> std::vector<Entity_ID>
  const & faces() const;

  /*! 
    @brief Wedge list
    @tparam ptype   Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> std::vector<Entity_ID>
  const & wedges() const;

  /*! 
    @brief Corner list
    @tparam ptype   Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> std::vector<Entity_ID>
  const & corners() const;

  /*! 
    @brief Cell list
    @tparam ptype   Parallel type (OWNED, GHOST, ALL)
  */
  template<Parallel_type ptype = ALL> std::vector<Entity_ID>
  const & cells() const;

 private:

  // Data

  Mesh const & mesh_;

  std::vector<Entity_ID> nodeids_owned_, nodeids_ghost_, nodeids_all_;
  std::vector<Entity_ID> edgeids_owned_, edgeids_ghost_, edgeids_all_;
  std::vector<Entity_ID> faceids_owned_, faceids_ghost_, faceids_all_;
  std::vector<Entity_ID> wedgeids_owned_, wedgeids_ghost_, wedgeids_all_;
  std::vector<Entity_ID> cornerids_owned_, cornerids_ghost_, cornerids_all_;
  std::vector<Entity_ID> cellids_owned_, cellids_ghost_, cellids_all_;
  std::vector<Entity_ID> dummy_list_;

  // Make the State class a friend so that it can access protected
  // methods for retrieving and storing mesh fields

  friend class State;

}; // End class MeshTile



// Default implementation of num_nodes prints an error and returns
// 0 while the meaningful implementation is in specialized routines
template<Parallel_type ptype> inline
unsigned int MeshTile::num_nodes() const {
  std::cerr << "MeshTile::num_nodes() - " <<
      "Meaningless to query for list of nodes of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_nodes<OWNED>() const {
  return nodeids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_nodes<GHOST>() const {
  return nodeids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_nodes<ALL>() const {
  return num_nodes<OWNED>() + num_nodes<GHOST>();
}

template<Parallel_type ptype> inline
unsigned int MeshTile::num_edges() const {
  std::cerr << "MeshTile::num_edges() - " <<
      "Meaningless to query for list of edges of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_edges<OWNED>() const {
  return edgeids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_edges<GHOST>() const {
  return edgeids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_edges<ALL>() const {
  return num_edges<OWNED>() + num_edges<GHOST>();
}

template<Parallel_type ptype> inline
unsigned int MeshTile::num_faces() const {
  std::cerr << "MeshTile::num_faces() - " <<
      "Meaningless to query for list of faces of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_faces<OWNED>() const {
  return faceids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_faces<GHOST>() const {
  return faceids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_faces<ALL>() const {
  return num_faces<OWNED>() + num_faces<GHOST>();
}

template<Parallel_type ptype> inline
unsigned int MeshTile::num_wedges() const {
  std::cerr << "MeshTile::num_wedges() - " <<
      "Meaningless to query for list of wedges of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_wedges<OWNED>() const {
return wedgeids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_wedges<GHOST>() const {
  return wedgeids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_wedges<ALL>() const {
  return num_wedges<OWNED>() + num_wedges<GHOST>();
}

template<Parallel_type ptype> inline
unsigned int MeshTile::num_corners() const {
  std::cerr << "MeshTile::num_corners() - " <<
      "Meaningless to query for list of corners of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_corners<OWNED>() const {
  return cornerids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_corners<GHOST>() const {
  return cornerids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_corners<ALL>() const {
  return num_corners<OWNED>() + num_corners<GHOST>();
}

template<Parallel_type ptype> inline
unsigned int MeshTile::num_cells() const {
  std::cerr << "MeshTile::num_cells() - " <<
      "Meaningless to query for list of cells of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_cells<OWNED>() const {
  return cellids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_cells<GHOST>() const {
  return cellids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_cells<ALL>() const {
  return num_cells<OWNED>() + num_cells<GHOST>();
}


// entity lists (default implementation prints error message -
// meaningful values returned through template specialization)

template<Parallel_type ptype> inline
const std::vector<Entity_ID> & MeshTile::nodes() const {
  std::cerr << "MeshTile::nodes() - " <<
      "Meaningless to query for list of nodes of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::nodes<OWNED>() const {
  return nodeids_owned_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::nodes<GHOST>() const {
  return nodeids_ghost_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::nodes<ALL>() const {
  return nodeids_all_;
}

template<Parallel_type ptype> inline
const std::vector<Entity_ID> & MeshTile::edges() const {
  std::cerr << "MeshTile::edges() - " <<
      "Meaningless to query for list of edges of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::edges<OWNED>() const {
  return edgeids_owned_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::edges<GHOST>() const {
  return edgeids_ghost_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::edges<ALL>() const {
  return edgeids_all_;
}

template<Parallel_type ptype> inline
const std::vector<Entity_ID> & MeshTile::faces() const {
  std::cerr << "MeshTile::faces() - " <<
      "Meaningless to query for list of faces of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::faces<OWNED>() const {
  return faceids_owned_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::faces<GHOST>() const {
  return faceids_ghost_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::faces<ALL>() const {
  return faceids_all_;
}


template<Parallel_type ptype> inline
const std::vector<Entity_ID> & MeshTile::wedges() const {
  std::cerr << "MeshTile::wedges() - " <<
      "Meaningless to query for list of wedges of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::wedges<OWNED>() const {
  return wedgeids_owned_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::wedges<GHOST>() const {
  return wedgeids_ghost_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::wedges<ALL>() const {
  return wedgeids_all_;
}


template<Parallel_type ptype> inline
const std::vector<Entity_ID> & MeshTile::corners() const {
  std::cerr << "MeshTile::corners() - " <<
      "Meaningless to query for list of corners of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::corners<OWNED>() const {
  return cornerids_owned_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::corners<GHOST>() const {
  return cornerids_ghost_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::corners<ALL>() const {
  return cornerids_all_;
}


template<Parallel_type ptype> inline
const std::vector<Entity_ID> & MeshTile::cells() const {
  std::cerr << "MeshTile::cells() - " <<
      "Meaningless to query for list of cells of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::cells<OWNED>() const {
  return cellids_owned_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::cells<GHOST>() const {
  return cellids_ghost_;
}
template<> inline
const std::vector<Entity_ID> & MeshTile::cells<ALL>() const {
  return cellids_all_;
}

// @brief MeshTile factory
//
// Standalone function to make a tile and return a pointer to it so
// that Mesh.hh can use a forward declaration of MeshTile and this
// function to create new tiles

std::shared_ptr<MeshTile> make_meshtile(Mesh const & parent_mesh,
                                        const std::vector<Entity_ID> & cells,
                                        const bool request_faces,
                                        const bool request_edges,
                                        const bool request_wedges,
                                        const bool request_corners);


} // end namespace Jali





#endif /* _JALI_MESHTILE_H_ */
