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

#include "MeshSet.hh"

#include <mpi.h>

#include <vector>
#include <algorithm>
#include <memory>
#include <cassert>


#include "MeshDefs.hh"
#include "Mesh.hh"

namespace Jali {

/*! 
  @brief Constructor for MeshSet
  
  Mesh sets are lists of mesh entities of a particular kind (typically
  CELL, FACE or NODE)

*/


// Constructor - should we make this private and only allow Mesh to
// call it as a friend? MeshSet can send a reference to itself to the
// parent_mesh so that it can be added to the list of sets

MeshSet::MeshSet(std::string const& name,
                 Mesh& parent_mesh,
                 Entity_kind const& kind,
                 std::vector<Entity_ID> const& owned_entities,
                 std::vector<Entity_ID> const& ghost_entities,
                 bool build_reverse_map) :
    mesh_(parent_mesh),
    name_(name),
    kind_(kind),
    entityids_owned_(owned_entities),
    entityids_ghost_(ghost_entities),
    have_reverse_map_(build_reverse_map) {

  entityids_all_ = entityids_owned_;
  entityids_all_.insert(entityids_all_.end(), entityids_ghost_.begin(),
                        entityids_ghost_.end());

  if (build_reverse_map) {
    have_reverse_map_ = true;
    mesh2subset_.resize(mesh_.num_entities(kind, Entity_type::ALL), -1);
    int nall = entityids_all_.size();
    for (int i = 0; i < nall; ++i)
      mesh2subset_[entityids_all_[i]] = i;
  }
}  // MeshSet::MeshSet


// Add entity to meshset (no check for duplicates)

void MeshSet::add_entity(Entity_ID const& mesh_entity) {
  Entity_type etype = mesh_.entity_get_type(kind_, mesh_entity);

  if (etype == Entity_type::PARALLEL_OWNED) {
    int nowned_old = entityids_owned_.size();

    entityids_owned_.push_back(mesh_entity);
    entityids_all_.insert(entityids_all_.begin()+nowned_old, mesh_entity);

    if (have_reverse_map_) mesh2subset_[mesh_entity] = nowned_old;

  } else if (etype == Entity_type::PARALLEL_GHOST) {

    entityids_ghost_.push_back(mesh_entity);

    entityids_all_.push_back(mesh_entity);
    if (have_reverse_map_) mesh2subset_[mesh_entity] = entityids_all_.size()-1;

  } else
    return;  // Doesn't make sense to add any other type like BOUNDARY_GHOST

}


// Remove entity from meshset (PREFERABLY USE rem_entities)

void MeshSet::rem_entity(Entity_ID const& mesh_entity) {
  Entity_type etype = mesh_.entity_get_type(kind_, mesh_entity);
    int size;
    if (etype == Entity_type::PARALLEL_OWNED) {
      size = entityids_owned_.size();
      auto const& it = std::find(entityids_owned_.begin(),
                                 entityids_owned_.end(), mesh_entity);
      if (it != entityids_owned_.end()) {
        *it = entityids_owned_[size-1];  // replace mesh_entity with last entry
        entityids_owned_.resize(size-1);
      }
    } else if (etype == Entity_type::PARALLEL_GHOST) {
      size = entityids_ghost_.size();
      auto const& it = std::find(entityids_ghost_.begin(),
                                 entityids_ghost_.end(), mesh_entity);
      if (it != entityids_ghost_.end()) {
        *it = entityids_ghost_[size-1];  // replace mesh_entity with last entry
        entityids_ghost_.resize(size-1);
      }
    }
    else
      return;  // No other type like BOUNDARY_GHOST can be part of set

    size = entityids_all_.size();
    auto const& it = std::find(entityids_all_.begin(),
                               entityids_all_.end(), mesh_entity);
    if (it != entityids_all_.end()) {
      *it = entityids_all_[size-1];  // replace mesh_entity with last entry
      entityids_all_.resize(size-1);
    }

    if (have_reverse_map_) mesh2subset_[mesh_entity] = -1;
  }


// Add a group of entities to meshset (no check for duplicates)

void MeshSet::add_entities(std::vector<Entity_ID> const& in_entities) {
  int nowned_old = entityids_owned_.size();
  int nghost_old = entityids_ghost_.size();
  int nall = in_entities.size();
  int nowned = 0;
  int nghost = 0;
  for (auto const& mesh_entity : in_entities) {
    Entity_type etype = mesh_.entity_get_type(kind_, mesh_entity);
    if (etype == Entity_type::PARALLEL_OWNED)
      nowned++;
    else if (etype == Entity_type::PARALLEL_GHOST)
      nghost++;
  }
  if (nall == nowned)
    entityids_owned_.insert(entityids_owned_.end(), in_entities.begin(),
                            in_entities.end());
  else if (nall == nghost)
    entityids_owned_.insert(entityids_ghost_.end(), in_entities.begin(),
                            in_entities.end());
  else
    for (auto const& mesh_entity : in_entities)
      add_entity(mesh_entity);
  
  // entityids_all should always have owned entities first and ghost
  // entities last - so we can't just put in entities at the end -
  // we have to make a new list with the updated owned and ghost
  // entities and concatenate them (if there are any ghosts)
  
  if (entityids_ghost_.size()) {
    std::vector<Entity_ID> tmp_list(entityids_owned_);
    tmp_list.insert(tmp_list.end(), entityids_ghost_.begin(),
                    entityids_ghost_.end());
    entityids_all_.swap(tmp_list);
  } else
    entityids_all_.insert(entityids_all_.begin(), in_entities.begin(),
                          in_entities.end());

  if (have_reverse_map_) {  // have to update mesh to subset map
    // new owned entities should have gone in after old owned entities
    int start = nowned_old;
    for (int i = start; i < start + nowned; i++)
      mesh2subset_[entityids_all_[i]] = i;
    
    // new ghost entities should have gone in after old+new owned
    // entities and old ghost entities
    start = nowned_old + nowned + nghost_old;
    for (int i = start; i < start + nghost; i++)
      mesh2subset_[entityids_all_[i]] = i;
  }
}


// Remove a group of entities from a subset

void MeshSet::rem_entities(std::vector<Entity_ID> const& in_entities) {
  // Replace each to-be-removed entry with an entry from the end of
  // the list.  Then shrink the lists.
  
  int size = entityids_owned_.size();
  int ndel = 0;
  for (auto const& mesh_entity : in_entities) {
    Entity_type etype = mesh_.entity_get_type(kind_, mesh_entity);
    if (etype == Entity_type::PARALLEL_OWNED) {
      auto const& it = std::find(entityids_owned_.begin(),
                                 entityids_owned_.end(), mesh_entity);
      if (it != entityids_owned_.end()) {
        *it = entityids_owned_[size-ndel-1];
        ndel++;
      }
    }
  }
  entityids_owned_.resize(size-ndel);
  
  size = entityids_ghost_.size();
  ndel = 0;
  for (auto const& mesh_entity : in_entities) {
    auto const& it = std::find(entityids_ghost_.begin(),
                               entityids_ghost_.end(), mesh_entity);
    if (it != entityids_ghost_.end()) {
      *it = entityids_ghost_[size-ndel-1];
      ndel++;
    }
  }
  entityids_ghost_.resize(size-ndel);
  
  size = entityids_all_.size();
  ndel = 0;
  for (auto const& mesh_entity : in_entities) {
    auto const& it = std::find(entityids_all_.begin(),
                               entityids_all_.end(), mesh_entity);
    if (it != entityids_all_.end()) {
      *it = entityids_all_[size-ndel-1];
      ndel++;
    }
  }
  entityids_all_.resize(size-ndel);

  if (have_reverse_map_) {
    int nall = entityids_all_.size();
    for (int i = 0; i < nall; i++)
      mesh2subset_[entityids_all_[i]] = i;
  }
}

// Standalone function to make a set and return a pointer to it so
// that Mesh.hh can use a forward declaration of MeshSet and this
// function to create new sets

std::shared_ptr<MeshSet> make_meshset(std::string const& name,
                                      Mesh& parent_mesh,
                                      Entity_kind const& kind,
                                      Entity_ID_List const& entityids_owned,
                                      Entity_ID_List const& entityids_ghost,
                                      bool build_reverse_map) {

  // This is a less than optimal use of the make_shared function since
  // it involves two memory allocations but I am not able to do it in
  // the optimal way since the MeshSet constructor is private (to
  // ensure mesh sets are constructed correctly through this function
  // and proper mesh<-->meshset links are established) - even when I
  // make the std::make_shared_ptr class a friend of MeshSet, it
  // complains that the constructor is private

  auto set =  std::make_shared<MeshSet>(name, parent_mesh, kind,
                                        entityids_owned, entityids_ghost,
                                        build_reverse_map);
                                        
  parent_mesh.add_set(set);
  return set;
}


// Union of two or more mesh sets

std::shared_ptr<MeshSet>
merge(std::vector<std::shared_ptr<MeshSet>> const& inpsets, bool temporary) {
  assert(inpsets.size());
  std::shared_ptr<MeshSet> set0 = inpsets[0];

  if (inpsets.size() == 1)
    return set0;

  for (auto const& set : inpsets) {
    if (set == set0) continue;
    assert(&(set0->mesh_) == &(set->mesh_));
    assert(set0->kind_ == set->kind_);
  }

  
  // Add elements that are in any of the sets to result
  
  Entity_ID_List owned_list = set0->entityids_owned_;
  int maxownsize = 0;
  for (auto const& set : inpsets)
    maxownsize += set->entityids_owned_.size();
  owned_list.reserve(maxownsize);

  for (auto const& set : inpsets) {
    if (set == set0) continue;
    for (auto const& ent : set->entityids_owned_)
      if (std::find(owned_list.begin(), owned_list.end(), ent) ==
          owned_list.end())
        owned_list.push_back(ent);
  }
  
  
  Entity_ID_List ghost_list = set0->entityids_ghost_;
  int maxghostsize = 0;
  for (auto const& set : inpsets)
    maxghostsize += set->entityids_ghost_.size();
  ghost_list.reserve(maxghostsize);

  for (auto const& set : inpsets) {
    if (set == set0) continue;
    for (auto const& ent : set->entityids_ghost_)
      if (std::find(ghost_list.begin(), ghost_list.end(), ent) ==
          ghost_list.end())
        ghost_list.push_back(ent);
  }
  
  
  std::string newname = "(" + set0->name_ + ")";
  for (auto const& set : inpsets) {
    if (set == set0) continue;
    newname += "_PLUS_(" + set->name_ + ")";
  }
  
  // If either of these sets has the reverse map, then the result has it too
  bool build_reverse_map = set0->mesh2subset_.size() ? true : false;
  
  // If the set is temporary, we don't need to call make_meshset and
  // add it to the mesh
  if (temporary)
    return std::make_shared<MeshSet>(newname, set0->mesh_, set0->kind_,
                                     owned_list, ghost_list,
                                     build_reverse_map);
  else
    return make_meshset(newname, set0->mesh_, set0->kind_,
                        owned_list, ghost_list,
                        build_reverse_map);
}

// Subtraction of mesh sets from the first mesh set

std::shared_ptr<MeshSet>
subtract(std::shared_ptr<MeshSet> const& set0,
         std::vector<std::shared_ptr<MeshSet>> const& subtractsets,
         bool temporary) {

  for (auto const& set : subtractsets) {
    assert(&(set0->mesh_) == &(set->mesh_));
    assert(set0->kind_ == set->kind_);
  }

  // Make a temporary union of all the sets to be subtracted
  std::shared_ptr<MeshSet> setunion = merge(subtractsets, true);

  // Add elements that are in set0 and but not in the rest of the sets
  // to the result
  
  Entity_ID_List owned_list;
  owned_list.reserve(set0->entityids_owned_.size());
  for (auto const& ent : set0->entityids_owned_) {
    if (std::find(setunion->entityids_owned_.begin(),
                  setunion->entityids_owned_.end(),
                  ent) == setunion->entityids_owned_.end())
      owned_list.push_back(ent);
  }
  
  Entity_ID_List ghost_list;
  ghost_list.reserve(set0->entityids_ghost_.size());
  for (auto const& ent : set0->entityids_ghost_) {
    if (std::find(setunion->entityids_ghost_.begin(),
                  setunion->entityids_ghost_.end(),
                  ent) == setunion->entityids_ghost_.end())
      ghost_list.push_back(ent);
  }
  
  std::string newname = "(" + set0->name_ + ")_MINUS_(" + setunion->name_ + ")";

  // If either of these sets has the reverse map, then the result has it too
  bool build_reverse_map = set0->mesh2subset_.size() ? true : false;
  
  // If the set is temporary, we don't need to call make_meshset and
  // add it to the mesh
  if (temporary)
    return std::make_shared<MeshSet>(newname, set0->mesh_, set0->kind_,
                                     owned_list, ghost_list, build_reverse_map);
  else
    return make_meshset(newname, set0->mesh_, set0->kind_,
                        owned_list, ghost_list, build_reverse_map);
}


// Intersection of multiple mesh sets

std::shared_ptr<MeshSet>
intersect(std::vector<std::shared_ptr<MeshSet>> const& inpsets,
          bool temporary) {
  std::shared_ptr<MeshSet> set0 = inpsets[0];
  for (auto const& set : inpsets) {
    if (set == set0) continue;
    assert(&(set0->mesh_) == &(set->mesh_));
    assert(set0->kind_ == set->kind_);
  }

  // Add elements that are in every set to the result
  
  Entity_ID_List owned_list;
  owned_list.reserve(set0->entityids_owned_.size());
  for (auto const& ent : set0->entityids_owned_) {
    bool present_in_all = true;
    for (auto const& set : inpsets) {
      if (set == set0) continue;
      if (std::find(set->entityids_owned_.begin(),
                    set->entityids_owned_.end(),
                    ent) == set->entityids_owned_.end()) {
        present_in_all = false;
        break;
      }
      if (!present_in_all) break;
    }
    if (present_in_all)
      owned_list.push_back(ent);
  }
  
  Entity_ID_List ghost_list;
  ghost_list.reserve(set0->entityids_ghost_.size());
  for (auto const& ent : set0->entityids_ghost_) {
    bool present_in_all = true;
    for (auto const& set : inpsets) {
      if (set == set0) continue;
      if (std::find(set->entityids_ghost_.begin(),
                    set->entityids_ghost_.end(),
                    ent) == set->entityids_ghost_.end()) {
        present_in_all = false;
        break;
      }
      if (!present_in_all) break;
    }
    if (present_in_all)
      ghost_list.push_back(ent);
  }
  
  
  std::string newname = "(" + set0->name_ + ")";
  for (auto const& set : inpsets) {
    if (set == set0) continue;
    newname += "_INTERSECT_(" + set->name_ + ")";
  }
  
  bool build_reverse_map = set0->mesh2subset_.size() ? true : false;
  
  // If the set is temporary, we don't need to call make_meshset and
  // add it to the mesh
  if (temporary)
    return std::make_shared<MeshSet>(newname, set0->mesh_, set0->kind_,
                                     owned_list, ghost_list,
                                     build_reverse_map);
  else
    return make_meshset(newname, set0->mesh_, set0->kind_,
                        owned_list, ghost_list,
                        build_reverse_map);
}

// Complement of sets (all mesh entities not in the given sets)

std::shared_ptr<MeshSet>
complement(std::vector<std::shared_ptr<MeshSet>> inpsets, bool temporary) {
  std::shared_ptr<MeshSet> set0 = inpsets[0];
  int nent_owned = set0->mesh_.num_entities(set0->kind_,
                                            Entity_type::PARALLEL_OWNED);
  int nent_ghost = set0->mesh_.num_entities(set0->kind_,
                                            Entity_type::PARALLEL_GHOST);

  // Create a temporary union of the input sets

  std::shared_ptr<MeshSet> setunion = merge(inpsets, true);

  Entity_ID_List owned_list;
  owned_list.reserve(nent_owned - setunion->entityids_owned_.size());
  for (int ent = 0; ent < nent_owned; ent++)
    if (std::find(setunion->entityids_owned_.begin(),
                  setunion->entityids_owned_.end(),
                  ent) == setunion->entityids_owned_.end())
      owned_list.push_back(ent);
  
  Entity_ID_List ghost_list;
  ghost_list.reserve(nent_ghost - setunion->entityids_ghost_.size());
  for (int ent = 0; ent < nent_ghost; ent++)
    if (std::find(setunion->entityids_ghost_.begin(),
                  setunion->entityids_ghost_.end(),
                  ent) == setunion->entityids_ghost_.end())
      ghost_list.push_back(ent);
  
  std::string newname = "NOT_(" + setunion->name_ + ")";
  
  // If this set has the reverse map, then the result has it too
  bool build_reverse_map = set0->mesh2subset_.size() ? true : false;
  
  // If the set is temporary, we don't need to call make_meshset and
  // add it to the mesh
  if (temporary)
    return std::make_shared<MeshSet>(newname, set0->mesh_, set0->kind_,
                                     owned_list, ghost_list,
                                     build_reverse_map);
  else
    return make_meshset(newname, set0->mesh_, set0->kind_,
                        owned_list, ghost_list,
                        build_reverse_map);
}



}  // end namespace Jali

