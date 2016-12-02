//
// Copyright Los Alamos National Security, LLC 2009-2015
// All rights reserved. See Copyright notice in main directory
//

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
    entityids_ghost_(ghost_entities) {

  entityids_all_ = entityids_owned_;
  entityids_all_.insert(entityids_all_.end(), entityids_ghost_.begin(),
                        entityids_ghost_.end());

  if (build_reverse_map) {
    mesh2subset_.resize(mesh_.num_entities(kind, Entity_type::ALL), -1);
    int nall = entityids_all_.size();
    for (int i = 0; i < nall; ++i)
      mesh2subset_[entityids_all_[i]] = i;
  }
}  // MeshSet::MeshSet

// Standalone function to make a set and return a pointer to it so
// that Mesh.hh can use a forward declaration of MeshSet and this
// function to create new sets

std::shared_ptr<MeshSet> make_meshset(std::string const& name,
                                      Mesh& parent_mesh,
                                      Entity_kind const& kind,
                                      Entity_ID_List const& entityids_owned,
                                      Entity_ID_List const& entityids_ghost,
                                      bool build_reverse_map) {
  if (parent_mesh.num_sets() == 0)
    parent_mesh.init_sets();

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

