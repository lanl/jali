//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#ifndef _JALI_MESHTILE_H_
#define _JALI_MESHTILE_H_

#include <iostream>

#include <vector>
#include <algorithm>
#include <memory>

#include "mpi.h"

#include "MeshDefs.hh"

namespace Jali {

/*! 
  @class MeshTile "MeshTile.hh"
  @brief Encapsulates notion of a meshtile (a group of cells)
  
  Mesh tiles are small groupings of mesh cells on a compute node and
  are typically obtained by partitioning the full mesh on a compute
  node. They are lightweight structures that are merely lists of cells
  that are in a partition. 
  
  Meshtiles are a way for the work on a mesh to be broken up into
  manageable, parallely executable chunks.
  
  Meshtiles can have a layer of halo or ghost cells for ease of
  computations.  Note that ghost cells of a meshtile may or may not be
  an MPI ghost (Entity_type::PARALLEL_GHOST).  MPI ghost entities will not
  have a master tile ID (since they do not belong to any tile on this
  processor)

  **** IMPORTANT NOTE ABOUT CONSTANTNESS OF THIS CLASS ****
  Instantiating a const version of this class only guarantees that
  the underlying mesh topology and geometry does not change (the
  public interface conforms strictly to this definition). However,
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

// forward declaration of the mesh class

class Mesh;

class MeshTile {

 public:

  /// @brief Constructor
  // 
  // NOTE - should we make this private and only allow Mesh to
  // call it as a friend? MeshTile can send a reference to itself to the
  // parent_mesh so that it can be added to the list of tiles. I ran into
  // C++ trouble when trying to do this so I will need C++ guru help
  

  MeshTile(Mesh& parent_mesh,
           std::vector<Entity_ID> const& meshcells_owned,
           int const num_halo_layers = 0,
           bool const request_faces = true,
           bool const request_edges = false,
           bool const request_sides = false,
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

  /// @brief IDentifier for the tile

  int ID() const {
    return mytileid_;
  }

  //
  // General mesh tile information
  // -------------------------
  //

  /*! 
    @brief Number of entities of a particular kind and parallel type
    @param kind Entity_kind of the entities (CELL, NODE, WEDGE etc)
    @param parallel_type Entity_type of entities (PARALLEL_OWNED,
    PARALLEL_GHOST, ALL)
  */

  unsigned int num_entities(Entity_kind kind, Entity_type parallel_type) const;

  /*! 
    @brief Number of nodes of parallel type
    @tparam ptype Parallel type (Entity_type::PARALLEL_OWNED,
                                 Entity_type::PARALLEL_GHOST,
                                 Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL>
  unsigned int num_nodes() const;

  /*! 
    @brief Number of edges of parallel type
    @tparam ptype Parallel type (Entity_type::PARALLEL_OWNED,
                                 Entity_type::PARALLEL_GHOST,
                                 Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL>
  unsigned int num_edges() const;

  /*! 
    @brief Number of faces of parallel type
    @tparam ptype Parallel type (Entity_type::PARALLEL_OWNED,
                                 Entity_type::PARALLEL_GHOST,
                                 Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL>
  unsigned int num_faces() const;

  /*! 
    @brief Number of sides of parallel type
    @tparam ptype Parallel type (Entity_type::PARALLEL_OWNED,
                                 Entity_type::PARALLEL_GHOST,
                                 Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL>
  unsigned int num_sides() const;

  /*! 
    @brief Number of wedges of parallel type
    @tparam ptype Parallel type (Entity_type::PARALLEL_OWNED,
                                 Entity_type::PARALLEL_GHOST,
                                 Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL>
  unsigned int num_wedges() const;

  /*! 
    @brief Number of corners of parallel type
    @tparam ptype Parallel type (Entity_type::PARALLEL_OWNED,
                                 Entity_type::PARALLEL_GHOST,
                                 Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL>
  unsigned int num_corners() const;

  /*! 
    @brief Number of nodes of parallel type
    @tparam ptype Parallel type (Entity_type::PARALLEL_OWNED,
                                 Entity_type::PARALLEL_GHOST,
                                 Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL>
  unsigned int num_cells() const;

  /*! 
    @brief Node list
    @tparam ptype   Parallel type (Entity_type::PARALLEL_OWNED,
                                   Entity_type::PARALLEL_GHOST,
                                   Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL>
  std::vector<Entity_ID> const & nodes() const;

  /*! 
    @brief Edge list
    @tparam ptype   Parallel type (Entity_type::PARALLEL_OWNED,
                                   Entity_type::PARALLEL_GHOST,
                                   Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL> std::vector<Entity_ID>
  const & edges() const;

  /*! 
    @brief Face list
    @tparam ptype   Parallel type (Entity_type::PARALLEL_OWNED,
                                   Entity_type::PARALLEL_GHOST,
                                   Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL> std::vector<Entity_ID>
  const & faces() const;

  /*! 
    @brief Side list
    @tparam ptype   Parallel type (Entity_type::PARALLEL_OWNED,
                                   Entity_type::PARALLEL_GHOST,
                                   Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL> std::vector<Entity_ID>
  const & sides() const;

  /*! 
    @brief Wedge list
    @tparam ptype   Parallel type (Entity_type::PARALLEL_OWNED,
                                   Entity_type::PARALLEL_GHOST,
                                   Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL> std::vector<Entity_ID>
  const & wedges() const;

  /*! 
    @brief Corner list
    @tparam ptype   Parallel type (Entity_type::PARALLEL_OWNED,
                                   Entity_type::PARALLEL_GHOST,
                                   Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL> std::vector<Entity_ID>
  const & corners() const;

  /*! 
    @brief Cell list
    @tparam ptype   Parallel type (Entity_type::PARALLEL_OWNED,
                                   Entity_type::PARALLEL_GHOST,
                                   Entity_type::ALL)
  */
  template<Entity_type ptype = Entity_type::ALL> std::vector<Entity_ID>
  const & cells() const;


  //! Get list of tile entities of type 'kind' and 'ptype' in set ('setname')

  void get_set_entities(const Set_Name setname,
                        const Entity_kind kind,
                        const Entity_type ptype,
                        Entity_ID_List *entids) const;

 private:

  void get_nodes_of_set(const Set_Name setname, const Entity_type ptype,
                        Entity_ID_List *entids) const;
  void get_edges_of_set(const Set_Name setname, const Entity_type ptype,
                        Entity_ID_List *entids) const;
  void get_faces_of_set(const Set_Name setname, const Entity_type ptype,
                        Entity_ID_List *entids) const;
  void get_sides_of_set(const Set_Name setname, const Entity_type ptype,
                        Entity_ID_List *entids) const;
  void get_wedges_of_set(const Set_Name setname, const Entity_type ptype,
                        Entity_ID_List *entids) const;
  void get_corners_of_set(const Set_Name setname, const Entity_type ptype,
                        Entity_ID_List *entids) const;
  void get_cells_of_set(const Set_Name setname, const Entity_type ptype,
                        Entity_ID_List *entids) const;

  // Data

  Mesh& mesh_;

  unsigned int const mytileid_;

  Entity_ID_List nodeids_owned_, nodeids_ghost_, nodeids_all_;
  Entity_ID_List edgeids_owned_, edgeids_ghost_, edgeids_all_;
  Entity_ID_List faceids_owned_, faceids_ghost_, faceids_all_;
  Entity_ID_List sideids_owned_, sideids_ghost_, sideids_all_;
  Entity_ID_List wedgeids_owned_, wedgeids_ghost_, wedgeids_all_;
  Entity_ID_List cornerids_owned_, cornerids_ghost_, cornerids_all_;
  Entity_ID_List cellids_owned_, cellids_ghost_, cellids_all_;
  Entity_ID_List dummy_list_;

  

  // Make the State class a friend so that it can access protected
  // methods for retrieving and storing mesh fields

  friend class State;


};  // End class MeshTile



// Default implementation of num_nodes prints an error and returns
// 0 while the meaningful implementation is in specialized routines
template<Entity_type ptype> inline
unsigned int MeshTile::num_nodes() const {
  std::cerr << "MeshTile::num_nodes() - " <<
      "Meaningless to query for list of nodes of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_nodes<Entity_type::PARALLEL_OWNED>() const {
  return nodeids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_nodes<Entity_type::PARALLEL_GHOST>() const {
  return nodeids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_nodes<Entity_type::ALL>() const {
  return (num_nodes<Entity_type::PARALLEL_OWNED>() +
          num_nodes<Entity_type::PARALLEL_GHOST>());
}

template<Entity_type ptype> inline
unsigned int MeshTile::num_edges() const {
  std::cerr << "MeshTile::num_edges() - " <<
      "Meaningless to query for list of edges of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_edges<Entity_type::PARALLEL_OWNED>() const {
  return edgeids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_edges<Entity_type::PARALLEL_GHOST>() const {
  return edgeids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_edges<Entity_type::ALL>() const {
  return (num_edges<Entity_type::PARALLEL_OWNED>() +
          num_edges<Entity_type::PARALLEL_GHOST>());
}

template<Entity_type ptype> inline
unsigned int MeshTile::num_faces() const {
  std::cerr << "MeshTile::num_faces() - " <<
      "Meaningless to query for list of faces of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_faces<Entity_type::PARALLEL_OWNED>() const {
  return faceids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_faces<Entity_type::PARALLEL_GHOST>() const {
  return faceids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_faces<Entity_type::ALL>() const {
  return (num_faces<Entity_type::PARALLEL_OWNED>() +
          num_faces<Entity_type::PARALLEL_GHOST>());
}

template<Entity_type ptype> inline
unsigned int MeshTile::num_sides() const {
  std::cerr << "MeshTile::num_sides() - " <<
      "Meaningless to query for list of sides of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_sides<Entity_type::PARALLEL_OWNED>() const {
return sideids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_sides<Entity_type::PARALLEL_GHOST>() const {
  return sideids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_sides<Entity_type::ALL>() const {
  return (num_sides<Entity_type::PARALLEL_OWNED>() +
          num_sides<Entity_type::PARALLEL_GHOST>());
}

template<Entity_type ptype> inline
unsigned int MeshTile::num_wedges() const {
  std::cerr << "MeshTile::num_wedges() - " <<
      "Meaningless to query for list of wedges of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_wedges<Entity_type::PARALLEL_OWNED>() const {
return wedgeids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_wedges<Entity_type::PARALLEL_GHOST>() const {
  return wedgeids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_wedges<Entity_type::ALL>() const {
  return (num_wedges<Entity_type::PARALLEL_OWNED>() +
          num_wedges<Entity_type::PARALLEL_GHOST>());
}

template<Entity_type ptype> inline
unsigned int MeshTile::num_corners() const {
  std::cerr << "MeshTile::num_corners() - " <<
      "Meaningless to query for list of corners of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_corners<Entity_type::PARALLEL_OWNED>() const {
  return cornerids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_corners<Entity_type::PARALLEL_GHOST>() const {
  return cornerids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_corners<Entity_type::ALL>() const {
  return (num_corners<Entity_type::PARALLEL_OWNED>() +
          num_corners<Entity_type::PARALLEL_GHOST>());
}

template<Entity_type ptype> inline
unsigned int MeshTile::num_cells() const {
  std::cerr << "MeshTile::num_cells() - " <<
      "Meaningless to query for list of cells of kind parallel type " <<
      ptype << "\n";
  return 0;
}
template<> inline
unsigned int MeshTile::num_cells<Entity_type::PARALLEL_OWNED>() const {
  return cellids_owned_.size();
}
template<> inline
unsigned int MeshTile::num_cells<Entity_type::PARALLEL_GHOST>() const {
  return cellids_ghost_.size();
}
template<> inline
unsigned int MeshTile::num_cells<Entity_type::ALL>() const {
  return (num_cells<Entity_type::PARALLEL_OWNED>() +
          num_cells<Entity_type::PARALLEL_GHOST>());
}



inline
unsigned int MeshTile::num_entities(const Entity_kind kind,
                                   const Entity_type ptype) const {
  switch (kind) {
    case Entity_kind::NODE:
      switch (ptype) {
        case Entity_type::PARALLEL_OWNED:
          return num_nodes<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_nodes<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_nodes<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::EDGE:
      switch (ptype) {
        case Entity_type::PARALLEL_OWNED:
          return num_edges<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_edges<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_edges<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::FACE:
      switch (ptype) {
        case Entity_type::PARALLEL_OWNED:
          return num_faces<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_faces<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_faces<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::SIDE:
      switch (ptype) {
        case Entity_type::PARALLEL_OWNED:
          return num_sides<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_sides<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_sides<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::WEDGE:
      switch (ptype) {
        case Entity_type::PARALLEL_OWNED:
          return num_wedges<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_wedges<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_wedges<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::CORNER:
      switch (ptype) {
        case Entity_type::PARALLEL_OWNED:
          return num_corners<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_corners<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_corners<Entity_type::ALL>();
        default: return 0;
      }
    case Entity_kind::CELL:
      switch (ptype) {
        case Entity_type::PARALLEL_OWNED:
          return num_cells<Entity_type::PARALLEL_OWNED>();
        case Entity_type::PARALLEL_GHOST:
          return num_cells<Entity_type::PARALLEL_GHOST>();
        case Entity_type::ALL:
          return num_cells<Entity_type::ALL>();
        default: return 0;
      }
    default:
      return 0;
  }
}




// entity lists (default implementation prints error message -
// meaningful values returned through template specialization)

template<Entity_type ptype> inline
const std::vector<Entity_ID> & MeshTile::nodes() const {
  std::cerr << "MeshTile::nodes() - " <<
      "Meaningless to query for list of nodes of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::nodes<Entity_type::PARALLEL_OWNED>() const {
  return nodeids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::nodes<Entity_type::PARALLEL_GHOST>() const {
  return nodeids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::nodes<Entity_type::ALL>() const {
  return nodeids_all_;
}

template<Entity_type ptype> inline
const std::vector<Entity_ID> & MeshTile::edges() const {
  std::cerr << "MeshTile::edges() - " <<
      "Meaningless to query for list of edges of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::edges<Entity_type::PARALLEL_OWNED>() const {
  return edgeids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::edges<Entity_type::PARALLEL_GHOST>() const {
  return edgeids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::edges<Entity_type::ALL>() const {
  return edgeids_all_;
}

template<Entity_type ptype> inline
const std::vector<Entity_ID> & MeshTile::faces() const {
  std::cerr << "MeshTile::faces() - " <<
      "Meaningless to query for list of faces of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::faces<Entity_type::PARALLEL_OWNED>() const {
  return faceids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::faces<Entity_type::PARALLEL_GHOST>() const {
  return faceids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::faces<Entity_type::ALL>() const {
  return faceids_all_;
}


template<Entity_type ptype> inline
const std::vector<Entity_ID> & MeshTile::sides() const {
  std::cerr << "MeshTile::sides() - " <<
      "Meaningless to query for list of sides of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::sides<Entity_type::PARALLEL_OWNED>() const {
  return sideids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::sides<Entity_type::PARALLEL_GHOST>() const {
  return sideids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::sides<Entity_type::ALL>() const {
  return sideids_all_;
}


template<Entity_type ptype> inline
const std::vector<Entity_ID> & MeshTile::wedges() const {
  std::cerr << "MeshTile::wedges() - " <<
      "Meaningless to query for list of wedges of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::wedges<Entity_type::PARALLEL_OWNED>() const {
  return wedgeids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::wedges<Entity_type::PARALLEL_GHOST>() const {
  return wedgeids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::wedges<Entity_type::ALL>() const {
  return wedgeids_all_;
}


template<Entity_type ptype> inline
const std::vector<Entity_ID>&
MeshTile::corners() const {
  std::cerr << "MeshTile::corners() - " <<
      "Meaningless to query for list of corners of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::corners<Entity_type::PARALLEL_OWNED>() const {
  return cornerids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::corners<Entity_type::PARALLEL_GHOST>() const {
  return cornerids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::corners<Entity_type::ALL>() const {
  return cornerids_all_;
}


template<Entity_type ptype> inline
const std::vector<Entity_ID> & MeshTile::cells() const {
  std::cerr << "MeshTile::cells() - " <<
      "Meaningless to query for list of cells of parallel type " <<
      ptype << "\n";
  return dummy_list_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::cells<Entity_type::PARALLEL_OWNED>() const {
  return cellids_owned_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::cells<Entity_type::PARALLEL_GHOST>() const {
  return cellids_ghost_;
}
template<> inline
const std::vector<Entity_ID>&
MeshTile::cells<Entity_type::ALL>() const {
  return cellids_all_;
}

// @brief MeshTile factory
//
// Standalone function to make a tile and return a pointer to it so
// that Mesh.hh can use a forward declaration of MeshTile and this
// function to create new tiles

std::shared_ptr<MeshTile> make_meshtile(Mesh& parent_mesh,
                                        std::vector<Entity_ID> const& cells,
                                        int const num_halo_layers,
                                        bool const request_faces,
                                        bool const request_edges,
                                        bool const request_sides,
                                        bool const request_wedges,
                                        bool const request_corners);


}  // end namespace Jali





#endif /* _JALI_MESHTILE_H_ */
