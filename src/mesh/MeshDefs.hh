//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//
/*!
 * @file   MeshDefs.hh
 * @brief  Various definitions needed by Mesh
 * @author Originally by Bill Perkins (PNNL), Enhanced by Rao Garimella (LANL)
 *
 */

#ifndef _MeshDefs_hh_
#define _MeshDefs_hh_

#include <vector>
#include <string>
#include <cstdint>

namespace Jali {

// Necessary typedefs and enumerations

typedef int Entity_ID;
typedef int Set_ID;
typedef std::string Set_Name;
typedef std::vector<Entity_ID> Entity_ID_List;
typedef std::vector<Set_ID> Set_ID_List;

// Mesh Type

enum class Mesh_type {
  RECTANGULAR,   // Equivalent of structured but can't use i,j,k notation
  GENERAL        // general unstructured
};

// Cells (aka zones/elements) are the highest dimension entities in a mesh
// Nodes (aka vertices) are lowest dimension entities in a mesh
// Faces in a 3D mesh are 2D entities, in a 2D mesh are 1D entities
// BOUNDARY_FACE is a special type of entity that is need so that process
// kernels can define composite vectors (see src/data_structures) on
// exterior boundary faces of the mesh only
//
// Wedges are special subcell entities that are a simplicial
// decomposition of cell. In 3D, a wedge is tetrahedron formed by one
// point of the edge, the midpoint of the edge, the "center" of the
// face and the "center" of the cell volume. In 2D, a wedge is a
// triangle formed by an end-point of the edge, the mid-point of the
// edge and the center of the cell. In 1D, (IS THIS CORRECT?), wedges
// are lines, that are formed by the endpoint of the cell and the
// midpoint of the cell. There are two wedges associated with an edge
// of cell face in 3D.
//
// Corners are also subcell entities that are associated uniquely with
// a node of a cell. Each corner is the union of all the wedges incident
// upon that node in the cell
//
// Facets are the boundary entity between two wedges in adjacent
// cells. In 3D, a facet is a triangular subface of the cell face
// shared by two wedges in adjacent cells. In 2D, a facet is half of
// an edge that is shared by two wedges in adjacent cells
//


enum class Entity_kind : std::int8_t {
  ALL_KIND = -3,
  ANY_KIND = -2,
  UNKNOWN_KIND = -1,
  NODE = 0,
  EDGE = 1,
  FACE = 2,
  CELL = 3,
  SIDE = 4,
  WEDGE = 5,
  CORNER = 6,
  FACET = 7,
  BOUNDARY_FACE = 8
};

const int NUM_ENTITY_KINDS = 9;  // Don't want to count the 3 catch-all types

// Return a string describing the entity kind that can be printed out
inline
std::string Entity_kind_string(Entity_kind kind) {
  static const std::string entity_kind_string[9] =
      {"Entity_kind::NODE", "Entity_kind::EDGE", "Entity_kind::FACE",
       "Entity_kind::CELL", "Entity_kind::SIDE", "Entity_kind::WEDGE",
       "Entity_kind::CORNER", "Entity_kind::FACET",
       "Entity_kind::BOUNDARY_FACE"};

  int ikind = static_cast<int>(kind);
  return (ikind >= 0 && ikind < NUM_ENTITY_KINDS) ?
      entity_kind_string[ikind] : "";
}

// Output operator for Entity_kind
inline
std::ostream& operator<<(std::ostream& os, const Entity_kind& kind) {
  os << " " << Entity_kind_string(kind) << " ";
  return os;
}

// Check if Entity_kind is valid

inline
bool entity_valid_kind(const Entity_kind kind) {
  return (kind >= Entity_kind::NODE && kind <= Entity_kind::CELL);
}




// Parallel status of entity

enum class Parallel_type :  std::uint8_t {
  PTYPE_UNKNOWN = 0,  // Initializer
  OWNED = 1,         // Owned by this processor
  GHOST = 2,         // Owned by another processor
  ALL  = 3           // Parallel_type::OWNED + Parall_type::Parallel_type::GHOST
};

const int NUM_PARALLEL_TYPES = 4;

// Return a string describing the entity kind that can be printed out
inline
std::string Parallel_type_string(Parallel_type const ptype) {
  static const std::string parallel_type_str[4] =
      {"Parallel_type::PTYPE_UNKNOWN", "Parallel_type::OWNED",
       "Parallel_type::GHOST", "Parallel_type::ALL"};

  int iptype = static_cast<int>(ptype);
  return (iptype >= 0 && iptype < NUM_PARALLEL_TYPES) ?
      parallel_type_str[iptype] : "";
}


// Output operator for Parallel_type
inline
std::ostream& operator<<(std::ostream& os, const Parallel_type& ptype) {
  os << " " << Parallel_type_string(ptype) << " ";
  return os;
}

// Check if Parallel_type is valid

inline
bool entity_valid_ptype(const Parallel_type ptype) {
  return (ptype >= Parallel_type::OWNED && ptype <= Parallel_type::ALL);
}




// Standard element types and catchall (POLYGON/POLYHED)

enum class Cell_type : std::uint8_t {
  CELLTYPE_UNKNOWN = 0,
  TRI = 1,
  QUAD,
  POLYGON,
  TET,
  PRISM,
  PYRAMID,
  HEX,
  POLYHED                // Polyhedron
};

const int NUM_CELL_TYPES = 9;

// Return an string description for each cell type
inline
std::string Cell_type_string(const Cell_type ctype) {
  static std::string cell_type_str[9] =
      {"Cell_type::CELLTYPE_UNKNOWN", "Cell_type::TRI", "Cell_type::QUAD",
       "Cell_type::POLYGON", "Cell_type::TET", "Cell_type::PRISM",
       "Cell_type::PYRAMID", "Cell_type::HEX", "Cell_type::POLYHED"};

  int ictype = static_cast<int>(ctype);
  return (ictype >= 0 && ictype < NUM_CELL_TYPES) ?
      cell_type_str[ictype] : "";
}

// Output operator for Cell_type
inline
std::ostream& operator<<(std::ostream& os, const Cell_type& ctype) {
  os << " " << Cell_type_string(ctype) << " ";
  return os;
}


// Check if Cell_type is valid
inline
bool cell_valid_type(const Cell_type type) {
  return (type >= Cell_type::TRI && type <= Cell_type::POLYHED);
}


// Types of partitioners (partitioning scheme bundled into the name)

enum class Partitioner_type : std::uint8_t {
    INDEX,
    METIS,
    ZOLTAN_GRAPH,
    ZOLTAN_RCB
};
constexpr int NUM_PARTITIONER_TYPES = 4;

// Return an string description for each partitioner type
inline
std::string Partitioner_type_string(const Partitioner_type partitioner_type) {
  static std::string partitioner_type_str[NUM_PARTITIONER_TYPES] =
      {"Partitioner_type::INDEX", "Partitioner_type::METIS",
       "Partitioner_type::ZOLTAN_GRAPH", "Partitioner_type::ZOLTAN_RCB"};

  int iptype = static_cast<int>(partitioner_type);
  return (iptype >= 0 && iptype < NUM_PARTITIONER_TYPES) ?
      partitioner_type_str[iptype] : "";
}

// Output operator for Partitioner_type
inline
std::ostream& operator<<(std::ostream& os,
                         const Partitioner_type& partitioner_type) {
  os << " " << Partitioner_type_string(partitioner_type) << " ";
  return os;
}

// Types of partitioning algorithms - Add as needed in the format METIS_RCB etc.

enum class Partitioning_scheme {DEFAULT};

}  // close namespace Jali



#endif
