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

/*!
 * @file   MeshDefs.hh
 * @brief  Various definitions needed by Mesh
 * @author Originally by Bill Perkins (PNNL), Enhanced by Rao Garimella (LANL)
 *
 */

#ifndef _MeshDefs_hh_
#define _MeshDefs_hh_

#include <iostream>
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
typedef std::int8_t dir_t;

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




// Entity type - includes UNKNOWN_TYPE, PARALLEL_OWNED, PARALLEL_GHOST,
// ALL, BOUNDARY_GHOST, DELETED

enum class Entity_type : std::int8_t {
    TYPE_UNKNOWN = -1,
    DELETED = 0,
    PARALLEL_OWNED = 1,    // Owned by this processor
    PARALLEL_GHOST = 2,    // Owned by another processor
    BOUNDARY_GHOST = 3,    // Ghost/Virtual entity on boundary
    ALL  = 4      // PARALLEL_OWNED + PARALLEL_GHOST + BOUNDARY_GHOST
};

const int NUM_ENTITY_TYPES = 6;

// Check if Entity_type is valid

inline
bool valid_entity_type(const Entity_type type) {
  int itype = static_cast<int>(type);
  return (itype >= 0 && itype < NUM_ENTITY_TYPES);
}

// Return a string describing the entity kind that can be printed out
inline
std::string Entity_type_string(Entity_type const type) {
  static const std::string type_str[6] =
      {"Entity_type::TYPE_UNKNOWN", "Entity_type::DELETED",
       "Entity_type::PARALLEL_OWNED", "Entity_type::PARALLEL_GHOST",
       "Entity_type::BOUNDARY_GHOST", "Entity_type::ALL"};

  int itype = static_cast<int>(type);
  return valid_entity_type(type) ? type_str[itype] : "";
}


// Output operator for Entity_type
inline
std::ostream& operator<<(std::ostream& os, const Entity_type& ptype) {
  os << " " << Entity_type_string(ptype) << " ";
  return os;
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
    BLOCK,
    METIS,
    ZOLTAN_GRAPH,
    ZOLTAN_RCB
};
constexpr int NUM_PARTITIONER_TYPES = 5;
constexpr Partitioner_type PARTITIONER_DEFAULT = Partitioner_type::METIS;

// Return an string description for each partitioner type
inline
std::string Partitioner_type_string(const Partitioner_type partitioner_type) {
  static std::string partitioner_type_str[NUM_PARTITIONER_TYPES] =
      {"Partitioner_type::INDEX", "Partitioner_type::BLOCK",
       "Partitioner_type::METIS",
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

}  // close namespace Jali



#endif
