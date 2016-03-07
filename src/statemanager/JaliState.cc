/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "JaliState.h"

namespace Jali {

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
    
    int nent = mymesh_->num_entities(kind, Parallel_type::ALL);
    
    for (int i = 0; i < num; i++) {
      if (vartypes[i] == "INT") {
        int *data = new int[nent];
        mymesh_->get_field(varnames[i], kind, data);
        Jali::StateVector<int> & sv = add(varnames[i], kind, data);
      } else if (vartypes[i] == "DOUBLE") {
        double *data = new double[nent];
        mymesh_->get_field(varnames[i], kind, data);
        Jali::StateVector<double> & sv = add(varnames[i], kind, data);          
      } else if (vartypes[i] == "VECTOR") {
        if (spacedim == 2) {
          std::array<double, 2> *data = new std::array<double, 2>[nent];
          mymesh_->get_field(varnames[i], kind, data);
          Jali::StateVector<std::array<double, 2>> & sv =
              add(varnames[i], kind, data);
        } else if (spacedim == 3) {
          std::array<double, 3> *data = new std::array<double, 3>[nent];
          mymesh_->get_field(varnames[i], kind, data);
          Jali::StateVector<std::array<double, 3>> & sv =
              add(varnames[i], kind, data);
        }
      } else if (vartypes[i] == "TENSOR") {  // assumes symmetric tensors
        if (spacedim == 2) {  // lower half & diagonal of 2x2 tensor
          std::array<double, 3> *data = new std::array<double, 3>[nent];
          mymesh_->get_field(varnames[i], kind, data);
          Jali::StateVector<std::array<double, 3>> & sv =
              add(varnames[i], kind, data);
        } else if (spacedim == 3) { // lower half & diagonal of 3x3 tensor
          std::array<double, 6> *data = new std::array<double, 6>[nent];
          mymesh_->get_field(varnames[i], kind, data);
          Jali::StateVector<std::array<double, 6>> & sv =
              add(varnames[i], kind, data);
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
    Entity_kind on_what = vec->on_what();
    bool status = false;

    if (vec->get_type() == typeid(double))
      status = mymesh_->store_field(name, on_what, (double *)vec->get_data());
    else if (vec->get_type() == typeid(int))
      status = mymesh_->store_field(name, on_what, (int *)vec->get_data());
    else if (vec->get_type() == typeid(std::array<double, 2>))
      status = mymesh_->store_field(name, on_what,
                                    (std::array<double, 2> *) vec->get_data());
    else if (vec->get_type() == typeid(std::array<double, 3>))
      status = mymesh_->store_field(name, on_what,
                                    (std::array<double, 3> *)vec->get_data());
    else if (vec->get_type() == typeid(std::array<double, 6>))
      status = mymesh_->store_field(name, on_what,
                                    (std::array<double,6> *)vec->get_data());
    

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
