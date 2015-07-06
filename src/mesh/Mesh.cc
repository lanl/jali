
#include "Geometry.hh"
#include "dbc.hh"
#include "errors.hh"
#include "LabeledSetRegion.hh"

#include "Mesh.hh"

#define CACHE_VARS 1

using namespace std;

namespace Jali
{


// Gather and cache cell to face connectivity info.
//
// Method is declared constant because it is not modifying the mesh
// itself; rather it is modifying mutable data structures - see
// declaration of Mesh calss for further explanation

void Mesh::cache_cell2face_info() const {

  int ncells = num_entities(CELL,ALL);
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

  int nfaces = num_entities(FACE,ALL);
  face_cell_ids.resize(nfaces);
  face_cell_ptype.resize(nfaces);
  
  std::vector<Entity_ID> fcells;
  
  for (int f = 0; f < nfaces; f++) {      
    face_get_cells_internal(f, ALL, &fcells);
    
    face_cell_ids[f].resize(2);
    face_cell_ptype[f].resize(2);
    
    for (int i = 0; i < fcells.size(); i++) {
      int c = fcells[i];
      face_cell_ids[f][i] = c;
      face_cell_ptype[f][i] = entity_get_ptype(CELL,c);
    }
    for (int i = fcells.size(); i < 2; i++) {
      face_cell_ids[f][i] = -1;
      face_cell_ptype[f][i] = PTYPE_UNKNOWN;
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

  int nfaces = num_entities(FACE,ALL);
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

  int ncells = num_entities(CELL,ALL);
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


Entity_ID Mesh::entity_get_parent(const Entity_kind kind, const Entity_ID entid) const {
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

  if (ordered)
    cell_get_faces_and_dirs_internal(cellid, faceids, face_dirs, ordered);
  else {
    Entity_ID_List &cfaceids = cell_face_ids[cellid];

    *faceids = cfaceids; // copy operation

    if (face_dirs) {
      std::vector<int> &cfacedirs = cell_face_dirs[cellid];
      *face_dirs = cfacedirs; // copy operation
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

void Mesh::face_get_cells (const Entity_ID faceid, const Parallel_type ptype,
			   Entity_ID_List *cellids) const {

#if CACHE_VARS != 0

  // 
  // Cached version - turn off for profiling or to save memory
  //
  
  if (!face2cell_info_cached) cache_face2cell_info();


  cellids->clear();
  
  switch (ptype) {
  case ALL:
    for (int i = 0; i < 2; i++)
      if (face_cell_ptype[faceid][i] != PTYPE_UNKNOWN)
	cellids->push_back(face_cell_ids[faceid][i]);
    break;
  case OWNED:
    for (int i = 0; i < 2; i++)
      if (face_cell_ptype[faceid][i] == OWNED) 
	cellids->push_back(face_cell_ids[faceid][i]);
    break;
  case GHOST:
    for (int i = 0; i < 2; i++)
      if (face_cell_ptype[faceid][i] == GHOST)
	cellids->push_back(face_cell_ids[faceid][i]);
    break;
  }
  
#else

  // 
  // Non-cached version
  //
  
  Entity_ID_List fcells;
  
  face_get_cells_internal(faceid, ALL, &fcells);
  
  cellids->clear();
  
  switch (ptype) {
  case ALL:
    for (int i = 0; i < fcells.size(); i++)
      if (entity_get_ptype(CELL,fcells[i]) != PTYPE_UNKNOWN)
	cellids->push_back(fcells[i]);
    break;
  case OWNED:
    for (int i = 0; i < fcells.size(); i++)
      if (entity_get_ptype(CELL,fcells[i]) == OWNED) 
	cellids->push_back(fcells[i]);
    break;
  case GHOST:
    for (int i = 0; i < fcells.size(); i++)
      if (entity_get_ptype(CELL,fcells[i]) == GHOST)
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

  *edgeids = face_edge_ids[faceid]; // copy operation

  if (edge_dirs) {
    std::vector<int> &fedgedirs = face_edge_dirs[faceid];
    *edge_dirs = fedgedirs; // copy operation
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

  face_get_edges_and_dirs(faceid, &fedgeids, &fedgedirs, true);
  cell_get_edges(cellid, &cedgeids);

  map->resize(fedgeids.size(),-1);
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

  *edgeids = cell_edge_ids[cellid]; // copy operation
  *edgedirs = cell_2D_edge_dirs[cellid];

#else

  // 
  // Non-cached version
  //

  cell_2D_get_edges_and_dirs_internal(cellid, edgeids, edgedirs);

#endif

} // Mesh::cell_get_edges_and_dirs





int Mesh::compute_cell_geometric_quantities() const {

  int ncells = num_entities(CELL,ALL);

  cell_volumes.resize(ncells);
  cell_centroids.resize(ncells);
  for (int i = 0; i < ncells; i++) {
    double volume;
    JaliGeometry::Point centroid(spacedim);

    compute_cell_geometry(i,&volume,&centroid);

    cell_volumes[i] = volume;
    cell_centroids[i] = centroid;
  }

  cell_geometry_precomputed = true;

  return 1;

} // Mesh::compute_cell_geometric_quantities



int Mesh::compute_face_geometric_quantities() const {

  if (space_dimension() == 3 && cell_dimension() == 2) {
    // need cell centroids to compute normals 

    if (!cell_geometry_precomputed)
      compute_cell_geometric_quantities();
  }

  int nfaces = num_entities(FACE,ALL);

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

    compute_face_geometry(i,&area,&centroid,&normal0,&normal1);

    face_areas[i] = area;
    face_centroids[i] = centroid;
    face_normal0[i] = normal0;
    face_normal1[i] = normal1;
  }

  face_geometry_precomputed = true;

  return 1;

} // Mesh::compute_face_geometric_quantities



int Mesh::compute_edge_geometric_quantities() const {

  int nedges = num_entities(EDGE,ALL);

  edge_vectors.resize(nedges);
  edge_lengths.resize(nedges);

  for (int i = 0; i < nedges; i++) {
    double length;
    JaliGeometry::Point evector(spacedim);

    compute_edge_geometry(i,&length,&evector);

    edge_lengths[i] = length;
    edge_vectors[i] = evector;
  }

  edge_geometry_precomputed = true;

  return 1;

} // Mesh::compute_edge_geometric_quantities



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

    cell_get_faces_and_dirs(cellid,&faces,&fdirs);

    int nf = faces.size();
    nfnodes.resize(nf);

    for (int j = 0; j < nf; j++) {

      face_get_coordinates(faces[j],&fcoords);
      nfnodes[j] = fcoords.size();

      if (fdirs[j] == 1) {
        for (int k = 0; k < nfnodes[j]; k++)
          cfcoords.push_back(fcoords[k]);
      }
      else {
        for (int k = nfnodes[j]-1; k >=0; k--)
          cfcoords.push_back(fcoords[k]);
      }
    }

    cell_get_coordinates(cellid,&ccoords);

    JaliGeometry::polyhed_get_vol_centroid(ccoords,nf,nfnodes,
            cfcoords,volume,
            centroid);
    return 1;
  }
  else if (celldim == 2) {

    std::vector<JaliGeometry::Point> ccoords;

    cell_get_coordinates(cellid,&ccoords);

    JaliGeometry::Point normal(spacedim);

    JaliGeometry::polygon_get_area_centroid_normal(ccoords,volume,centroid,
						     &normal);

    return 1;
  }

  return 0;
} // Mesh::compute_cell_geometry


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

    face_get_coordinates(faceid,&fcoords);

    JaliGeometry::Point normal(3);
    JaliGeometry::polygon_get_area_centroid_normal(fcoords,area,centroid,&normal);

    Entity_ID_List cellids;
    face_get_cells(faceid, ALL, &cellids);

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
  }
  else if (celldim == 2) {

    if (spacedim == 2) {   // 2D mesh

      face_get_coordinates(faceid,&fcoords);

      JaliGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);

      *centroid = 0.5*(fcoords[0]+fcoords[1]);

      JaliGeometry::Point normal(evec[1],-evec[0]);

      Entity_ID_List cellids;
      face_get_cells(faceid, ALL, &cellids);

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
    }
    else {  // Surface mesh - cells are 2D, coordinates are 3D

      // edge normals are ambiguous for surface mesh
      // So we won't compute them

      face_get_coordinates(faceid,&fcoords);

      JaliGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);

      *centroid = 0.5*(fcoords[0]+fcoords[1]);

      Entity_ID_List cellids;
      face_get_cells(faceid, ALL, &cellids);

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
          *normal1 = normal; // Note that we are not flipping the sign here
      }

      return 1;
    }

  }

  return 0;

} // Mesh::compute_face_geometry


int Mesh::compute_edge_geometry(const Entity_ID edgeid, double *edge_length,
				JaliGeometry::Point *edge_vector) const {

  (*edge_vector).set(0.0L);
  *edge_length = 0.0;

  Entity_ID node0, node1;

  edge_get_nodes(edgeid,&node0,&node1);

  JaliGeometry::Point point0, point1;
  node_get_coordinates(node0,&point0);
  node_get_coordinates(node1,&point1);

  *edge_vector = point1 - point0;
  *edge_length = norm(*edge_vector);

  return 0;

} // Mesh::compute_face_geometry



// Volume/Area of cell

double Mesh::cell_volume (const Entity_ID cellid, const bool recompute) const {

  if (!cell_geometry_precomputed) {
    compute_cell_geometric_quantities();
    return cell_volumes[cellid];
  }
  else {
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
  }
  else {
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
  }
  else {
    if (recompute) {
      double length;
      JaliGeometry::Point vector(spacedim);
      compute_edge_geometry(edgeid, &length, &vector);
      return length;
    }
    else
      return edge_lengths[edgeid];
  }
}

// Centroid of cell

JaliGeometry::Point Mesh::cell_centroid (const Entity_ID cellid, 
					   const bool recompute) const {

  if (!cell_geometry_precomputed) {
    compute_cell_geometric_quantities();
    return cell_centroids[cellid];
  }
  else {
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

JaliGeometry::Point Mesh::face_centroid (const Entity_ID faceid, const bool recompute) const {

  ASSERT(faces_requested);

  if (!face_geometry_precomputed) {
    compute_face_geometric_quantities();
    return face_centroids[faceid];
  }
  else {
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
  }
  else {
    if (recompute) {
      double area;
      JaliGeometry::Point centroid(spacedim);
      
      compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
    }
    else {
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
 
    if (L22(normal0) != 0.0)
      return normal0;
    else {
      ASSERT(L22(normal1) != 0.0);
      return -normal1;
    }
  }
  else {
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
    }
    else {
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

  JaliGeometry::Point evector(spacedim);
  JaliGeometry::Point& evector_ref = evector; // to avoid extra copying

  if (!edge_geometry_precomputed)
    compute_edge_geometric_quantities();

  if (recompute) {
    double length;
    compute_edge_geometry(edgeid, &length, &evector);
    // evector_ref already points to evector
  }
  else
    evector_ref = edge_vectors[edgeid];

  if (orientation) *orientation = 1;

  if (pointid == -1) 
    return evector_ref;
  else {
    Entity_ID p0, p1;
    edge_get_nodes(edgeid, &p0, &p1);

    if (pointid == p0)
      return evector_ref;
    else {
      if (orientation) *orientation=-1;
      return -evector_ref;
    }
  }

} // edge_vector


// Get set ID given the name of the set - return 0 if no match is found

// Set_ID Mesh::set_id_from_name(const std::string setname) const
// {
//   if (!geometric_model_) return 0;

//   unsigned int ngr = geometric_model_->Num_Regions();
//   for (int i = 0; i < ngr; i++) {
//     JaliGeometry::RegionPtr rgn = geometric_model_->Region_i(i);

//     if (rgn->name() == setname)
//       return rgn->id();
//   }

//   return 0;
// }

// Get set name given the ID of the set - return 0 if no match is found

// std::string Mesh::set_name_from_id(const int setid) const
// {
//   std::string nullname("");

//   if (!geometric_model_) return nullname;

//   unsigned int ngr = geometric_model_->Num_Regions();
//   for (int i = 0; i < ngr; i++) {
//     JaliGeometry::RegionPtr rgn = geometric_model_->Region_i(i);

//     if (rgn->id() == setid)
//       return rgn->name();
//   }

//   return 0;
// }



// Is there a set with this id and entity type

// bool Mesh::valid_set_id(Set_ID id, Entity_kind kind) const
// {

//   if (!geometric_model_) return false;

//   unsigned int gdim = geometric_model_->dimension();

//   unsigned int ngr = geometric_model_->Num_Regions();
//   for (int i = 0; i < ngr; i++) {
//     JaliGeometry::RegionPtr rgn = geometric_model_->Region_i(i);

//     unsigned int rdim = rgn->dimension();

//     if (rgn->id() == id) {

//       // For regions of type Labeled Set and Color Function, the
//       // dimension parameter is not guaranteed to be correct

//       if (rgn->type() == JaliGeometry::LABELEDSET ||
//           rgn->type() == JaliGeometry::COLORFUNCTION) return true;

//       // If we are looking for a cell set the region has to be
//       // of the same topological dimension as the cells

//       if (kind == CELL && rdim == celldim) return true;

//       // If we are looking for a side set, the region has to be
//       // one topological dimension less than the cells

//       if (kind == FACE && rdim == celldim-1) return true;

//       // If we are looking for a node set, the region can be of any
//       // dimension upto the spatial dimension of the domain

//       if (kind == NODE) return true;

//     }
//   }

//   return false;
// }



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
        JaliGeometry::LabeledSetRegionPtr lsrgn = dynamic_cast<JaliGeometry::LabeledSetRegionPtr> (rgn);
        std::string entity_type = lsrgn->entity_str();
        
        if ((kind == CELL && entity_type == "CELL") ||
            (kind == FACE && entity_type == "FACE") ||
            (kind == EDGE && entity_type == "EDGE") ||
            (kind == NODE && entity_type == "NODE"))
          return true;
        else
          return false;
      }

      // If we are looking for a cell set the region has to be
      // of the same topological dimension as the cells or it
      // has to be a point region

      if (kind == CELL && (rdim >= celldim || rdim == 0)) return true;

      // If we are looking for a side set, the region has to be
      // one topological dimension less than the cells

      if (kind == FACE && rdim >= celldim-1) return true;

      // If we are looking for a node set, the region can be of any
      // dimension upto the spatial dimension of the domain

      if (kind == NODE) return true;

    }
  }

  return false;
}


bool Mesh::point_in_cell(const JaliGeometry::Point &p, const Entity_ID cellid) const
{
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

    cell_get_faces_and_dirs(cellid,&faces,&fdirs);

    nf = faces.size();
    nfnodes.resize(nf);

    for (int j = 0; j < nf; j++) {
      std::vector<JaliGeometry::Point> fcoords;

      face_get_coordinates(faces[j],&fcoords);
      nfnodes[j] = fcoords.size();

      if (fdirs[j] == 1) {
        for (int k = 0; k < nfnodes[j]; k++)
          cfcoords.push_back(fcoords[k]);
      }
      else {
        for (int k = nfnodes[j]-1; k >=0; k--)
          cfcoords.push_back(fcoords[k]);
      }
    }

    cell_get_coordinates(cellid,&ccoords);

    return JaliGeometry::point_in_polyhed(p,ccoords,nf,nfnodes,cfcoords);

  }
  else if (celldim == 2) {

    cell_get_coordinates(cellid,&ccoords);

    return JaliGeometry::point_in_polygon(p,ccoords);

  }

  return false;
}


std::string Mesh::cell_type_to_name (const Cell_type type)
{

  switch (type)
  {
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

} // close namespace Jali
