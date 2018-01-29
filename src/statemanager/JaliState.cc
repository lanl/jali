/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "JaliState.h"
#include "JaliStateVector.h"

namespace Jali {

// Add new material to the state

void State::add_material(std::string const& matname,
                         std::vector<int> const& matcells) {
  for (auto const& mset : material_cellsets_)
    if (matname == mset->name()) {
      std::cerr << "Material name \"" << matname << "\" already used\n";
      return;
    }
    
  // Find a cell set with this name in the mesh. If it doesn't
  // exist, create it
  std::shared_ptr<MeshSet> matset =
      mymesh_->find_meshset(matname, Entity_kind::CELL);
  if (!matset) {
    Entity_ID_List owned_entities, ghost_entities;
    bool with_reverse_map = true;
    matset = make_meshset(matname, *mymesh_, Jali::Entity_kind::CELL,
                          owned_entities, ghost_entities, with_reverse_map);
  }

  matset->add_entities(matcells);
    
  material_cellsets_.push_back(matset);
    
  // GO TO EACH MMStateVector AND ADD ENTRIES FOR THIS MATERIAL
  for (auto &sv : state_vectors_) {
    if (sv->get_type() == StateVector_type::MULTIVAL) {
      auto mmv = std::dynamic_pointer_cast<MMStateVector<double, Mesh>>(sv);
      mmv->add_material(matcells.size());
    }
  }
}

// Remove material from the state -- EXPENSIVE (CAN JUST MAKE SET EMPTY)
void State::rem_material(int m) {
  if (m >= material_cellsets_.size())
    return;

  std::shared_ptr<MeshSet> matset = material_cellsets_[m];
  if (!matset) return;

  material_cellsets_.erase(material_cellsets_.begin()+m);

  // GO TO EACH MMStateVector AND REMOVE ENTRIES FOR THIS MATERIAL
  for (auto &sv : state_vectors_) {
    if (sv->get_type() == StateVector_type::MULTIVAL) {
      auto mmv = std::dynamic_pointer_cast<MMStateVector<double, Mesh>>(sv);
      mmv->rem_material(m);
    }
  }
}

/// Add cells to a material

void State::add_cells_to_material(int m, std::vector<int> const& cells) {
  assert(m < material_cellsets_.size());
  material_cellsets_[m]->add_entities(cells);
    
  for (auto & sv : state_vectors_) {
    if (sv->get_type() == StateVector_type::MULTIVAL) {
      auto mmv = std::dynamic_pointer_cast<MMStateVector<double, Mesh>>(sv);
      int size = mmv->size(m);
      mmv->resize(m, size+cells.size(), 0.0);
    }
  }
}
  
/// Remove cells from a material
  
void State::rem_cells_from_material(int m, std::vector<int> const& cells) {

  // NOT IMPLEMENTED - THERE ARE MAJOR CONSIDERATIONS OF EFFICIENCY
  // IN HOW WE DO THIS - IF WE REMOVE THESE ENTRIES FROM THE LIST
  // THEN WE HAVE TO REMOVE ENTRIES FROM ALL STATE VECTORS FROM
  // CORRESPONDING TO THIS MATERIAL AS WELL. BESIDES MESHSETS DO
  // FUNKY THINGS WITH REMOVING CELLS AND THAT HAS TO BE TAKEN
  // INTO ACCOUNT (See MeshSet::rem_entities)

  throw std::runtime_error("rem_cells_from_material NOT IMPLEMENTED");

}



//! \brief Add a state vectors from the mesh
//! Initialize a state vectors in the statemanager from mesh field data

void State::init_from_mesh() {

  int num;
  std::vector<std::string> varnames, vartypes;

  for (int ikind = 0; ikind < NUM_ENTITY_KINDS; ikind++) {
    Entity_kind kind = (Entity_kind) ikind;
    if (kind != Entity_kind::NODE && kind != Entity_kind::FACE &&
        kind != Entity_kind::CELL) continue;
    
    mymesh_->get_field_info(kind, &num, &varnames, &vartypes);
    if (!num) continue;

    int spacedim = mymesh_->space_dimension();
    
    int nent = mymesh_->num_entities(kind, Entity_type::ALL);
    
    for (int i = 0; i < num; i++) {
      if (vartypes[i] == "INT") {
        int *data = new int[nent];
        mymesh_->get_field(varnames[i], kind, data);
        Jali::StateVector<int, Mesh> & sv = add(varnames[i], mymesh_,
                                                kind,
                                                Entity_type::ALL, data);
      } else if (vartypes[i] == "DOUBLE") {
        double *data = new double[nent];
        mymesh_->get_field(varnames[i], kind, data);
        Jali::StateVector<double, Mesh> & sv = add(varnames[i], mymesh_,
                                                         kind,
                                                         Entity_type::ALL,
                                                         data);
      } else if (vartypes[i] == "VECTOR") {
        if (spacedim == 2) {
          std::array<double, 2> *data = new std::array<double, 2>[nent];
          mymesh_->get_field(varnames[i], kind, data);
          Jali::StateVector<std::array<double, 2>, Mesh> & sv =
              add(varnames[i], mymesh_, kind, Entity_type::ALL, data);
        } else if (spacedim == 3) {
          std::array<double, 3> *data = new std::array<double, 3>[nent];
          mymesh_->get_field(varnames[i], kind, data);
          Jali::StateVector<std::array<double, 3>, Mesh> & sv =
              add(varnames[i], mymesh_, kind, Entity_type::ALL, data);
        }
      } else if (vartypes[i] == "TENSOR") {  // assumes symmetric tensors
        if (spacedim == 2) {  // lower half & diagonal of 2x2 tensor
          std::array<double, 3> *data = new std::array<double, 3>[nent];
          mymesh_->get_field(varnames[i], kind, data);
          Jali::StateVector<std::array<double, 3>, Mesh> & sv =
              add(varnames[i], mymesh_, kind, Entity_type::ALL, data);
        } else if (spacedim == 3) {  // lower half & diagonal of 3x3 tensor
          std::array<double, 6> *data = new std::array<double, 6>[nent];
          mymesh_->get_field(varnames[i], kind, data);
          Jali::StateVector<std::array<double, 6>, Mesh> & sv =
              add(varnames[i], mymesh_, kind, Entity_type::ALL, data);
        }
      }  // TENSOR
    }  // for each field on entity kind
  }  // for each entity kind

}  // init_from_mesh


//! \brief Export field data to mesh
//! Export data from state vectors to mesh fields - Since the statevector is
//! templated, we have to go through case by case to see if the type matches
//! any of the types that the mesh can receive

void State::export_to_mesh() {

  State::const_iterator it = cbegin();

  while (it != cend()) {
    const std::shared_ptr<BaseStateVector> vec = *it;
    std::string name = vec->name();
    Entity_kind entity_kind = vec->entity_kind();
    bool status = false;

    if (vec->get_data_type() == typeid(double))
      status = mymesh_->store_field(name, entity_kind, (double *)vec->get_raw_data());
    else if (vec->get_data_type() == typeid(int))
      status = mymesh_->store_field(name, entity_kind, (int *)vec->get_raw_data());
    else if (vec->get_data_type() == typeid(std::array<double, 2>))
      status = mymesh_->store_field(name, entity_kind,
                                    (std::array<double, 2> *) vec->get_raw_data());
    else if (vec->get_data_type() == typeid(std::array<double, 3>))
      status = mymesh_->store_field(name, entity_kind,
                                    (std::array<double, 3> *) vec->get_raw_data());
    else if (vec->get_data_type() == typeid(std::array<double, 6>))
      status = mymesh_->store_field(name, entity_kind,
                                    (std::array<double, 6> *) vec->get_raw_data());
    

    if (!status)
      std::cerr << "Could not export vector " << name << " to mesh file\n";
      
    ++it;
  }
}


//! Print all state vectors

std::ostream & operator<<(std::ostream & os, State const & s) {
  State::const_iterator it = s.cbegin();
  while (it != s.cend()) {
    const std::shared_ptr<BaseStateVector> vec = *it;
    vec->print(os);
    ++it;
  }
}

}  // namespace Jali
