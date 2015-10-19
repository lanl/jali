// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: MeshDefs.hh
// -------------------------------------------------------------
/**
 * @file   MeshDefs.hh
 * @author William A. Perkins
 * @date Mon May  2 13:03:23 2011
 * 
 * @brief  Various definitions needed by Mesh
 * 
 * 
 */
// -------------------------------------------------------------
// Created May  2, 2011 by William A. Perkins
// Last Change: Mon May  2 13:03:23 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _MeshDefs_hh_
#define _MeshDefs_hh_

#include <vector>
#include <string>

namespace Jali {

// Necessary typedefs and enumerations
      
typedef int Entity_ID; 
typedef int Set_ID;
typedef std::string Set_Name;
typedef std::vector<Entity_ID> Entity_ID_List;
typedef std::vector<Set_ID> Set_ID_List;

// Mesh Type

enum Mesh_type
{
  RECTANGULAR,   // Equivalent of structured but can't use i,j,k notation
  GENERAL        // general unstructured
};

// Mesh Geometry Type

enum Geom_type
{
  CARTESIAN,
  CYLINDRICAL,
  SPHERICAL
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

    
enum Entity_kind 
{
  ALL_KIND = -3,
  ANY_KIND = -2,
  UNKNOWN_KIND = -1,
  NODE = 0,
  EDGE,
  FACE,
  CELL,
  WEDGE,
  CORNER,
  FACET,
  BOUNDARY_FACE,
  NUM_ENTITY_KINDS
};

// Check if Entity_kind is valid
inline 
bool entity_valid_kind (const Entity_kind kind) {
  return (kind >= NODE && kind <= CELL);
}

// Parallel status of entity 
    
enum Parallel_type 
{
  PTYPE_UNKNOWN = 0, // Initializer
  OWNED = 1,         // Owned by this processor
  GHOST = 2,         // Owned by another processor
  ALL  = 3           // OWNED + GHOST
};

// Check if Parallel_type is valid

inline 
bool entity_valid_ptype (const Parallel_type ptype) {
  return (ptype >= OWNED && ptype <= ALL);
}
    
// Standard element types and catchall (POLYGON/POLYHED)

enum Cell_type 
{
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
    
// Check if Cell_type is valid
inline 
bool cell_valid_type (const Cell_type type) {
  return (type >= TRI && type <= POLYHED); 
}

} // close namespace Jali



#endif
