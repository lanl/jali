//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

#ifndef _JALI_MESHSET_H_
#define _JALI_MESHSET_H_

#include "Mesh.hh"

namespace Jali {
  
class MeshSet {
  
 public:
  // Constructor - empty meshset
  
  MeshSet(Mesh const * const mesh, std::string const name, 
          Entity_kind entity_kind, bool with_reverse_map=true) : 
      mesh_(mesh), name_(name), entity_kind_(entity_kind),
      has_reverse_map_(with_reverse_map)
  {}

  // Constructor with an initial set of entities

  MeshSet(Mesh const * const mesh, std::string const name, 
          Entity_kind entity_kind, unsigned int const nentities, 
          unsigned int const * const entities, bool with_reverse_map=true) :
      mesh_(mesh), name_(name), entity_kind_(entity_kind), 
      has_reverse_map_(with_reverse_map)
  {
    meshset_to_mesh_map_.resize(nentities,-1);
    std::copy(entities, entities+nentities, meshset_to_mesh_map_.begin());

    if (with_reverse_map) {
      int nmeshents = mesh->num_entities(entity_kind,Jali::ALL);
      mesh_to_meshset_map_.resize(nmeshents,-1);
      for (int i = 0; i < nentities; ++i) {
        int mesh_entity = meshset_to_mesh_map_[i];
        mesh_to_meshset_map_[mesh_entity] = i;
      }
    }
  }

  // Copy constructor (deleted - don't want two meshsets with the same name)

  MeshSet(MeshSet const & meshset_in) = delete;

  // Assignment operator (deleted - don't want two meshsets with the same name)
  
  MeshSet & operator=(MeshSet const & meshset_in) = delete;
    
  // Destructor

  ~MeshSet() {};

  // What mesh does this meshset belong to

  Mesh const * const mesh() const {return mesh_;}

  // What is the name of this meshset?

  std::string  name() const {return name_;}

  // What kind of entities does this meshset contain?

  Entity_kind entity_kind() const {return entity_kind_;}

  // Number of entities in this meshset?

  unsigned int size() const {return meshset_to_mesh_map_.size();}

  // Add mesh entities to the meshset

  void add_entities(int const nadd, int const * const entities2add) {

    int cursize = meshset_to_mesh_map_.size();
    meshset_to_mesh_map_.resize(cursize+nadd,-1);

    if (has_reverse_map_) {
      for (int i = 0; i < nadd; ++i) {
        if (mesh_to_meshset_map_[entities2add[i]] == -1) {

          meshset_to_mesh_map_[cursize] = entities2add[i];
          mesh_to_meshset_map_[entities2add[i]] = cursize;
          cursize++;
          
        }
      }
    }
    else {
      for (int i = 0; i < nadd; ++i) {
        if (std::find(meshset_to_mesh_map_.begin(),meshset_to_mesh_map_.end(),
                      entities2add[i]) == meshset_to_mesh_map_.end()) {
          
          meshset_to_mesh_map_[cursize] = entities2add[i];
          cursize++;
          
        }
      }
    }
  }

  // Remove mesh entities from the meshset

  void remove_entities(int const nrem, int const * const  entities2rem) {
    int cursize = meshset_to_mesh_map_.size();

    for (int i = 0; i < nrem; ++i) {
      int local_index = mesh_to_meshset_map_[entities2rem[i]];

      if (local_index == -1) { // entity to remove is not part of this subset
        continue;
      }
      else if (local_index == cursize-1) {
        // This is the last element - nothing to do because when we
        // resize to the smaller size, this will be deleted

        cursize--;
      }
      else {
        // Deleting an element in the middle. Using an 'erase' op for
        // this is expensive. Since the order of elements does not
        // matter, we will intead take the last element of the set and
        // put it in this element's place instead of actually deleting
        // it. Simultaneously update the entry in the reverse map from
        // the mesh to meshset
        
        int last_mesh_entity = meshset_to_mesh_map_[cursize-1];
        meshset_to_mesh_map_[local_index] = last_mesh_entity;
        if (has_reverse_map_)
          mesh_to_meshset_map_[last_mesh_entity] = local_index;
        cursize--;
      }
    }

    // Since this is a shortening of the vector, no reallocation is done
    // and little cost is incurred

    meshset_to_mesh_map_.resize(cursize);
  }

  // The local ID in this meshset of a mesh entity

  int meshset_entity(int const mesh_entity_id) const {    
#ifdef DEBUG
    assert(mesh_entity_id < mesh_->num_entities(entitykind_,Jali::ALL));
    if (!has_reverse_map_)
      std::cerr << "Jali::MeshSet " << name_ << 
          " constructed without requesting a reverse (mesh->meshset) map" << 
          std::endl;
    assert(has_reverse_map_);
#endif
    return has_reverse_map_ ? mesh_to_meshset_map_[mesh_entity_id] : -1;
  }

  // The mesh entity of an entity in the meshset

  int mesh_entity(int const meshset_entity_id) const {
#ifdef DEBUG
    assert(meshset_entity_id < meshset_to_mesh_map_.size());
#endif
    return meshset_to_mesh_map_[meshset_entity_id];
  }
  
  
  
 private:
  
  Mesh const * const mesh_;
  std::string const name_;
  Entity_kind const entity_kind_;
  bool has_reverse_map_;
  std::vector<int> meshset_to_mesh_map_;
  std::vector<int> mesh_to_meshset_map_;
  
}; // class MeshSet

} // namespace Jali  

// _JALI_MESHSET_H_
#endif 
