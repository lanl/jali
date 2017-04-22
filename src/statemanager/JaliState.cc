/*
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.

Copyright 2017. Los Alamos National Security, LLC. This software was
produced under U.S. Government contract DE-AC52-06NA25396 for Los
Alamos National Laboratory (LANL), which is operated by Los Alamos
National Security, LLC for the U.S. Department of Energy. The
U.S. Government has rights to use, reproduce, and distribute this
software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so
as not to confuse it with the version available from LANL.
 
Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

1.  Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3.  Neither the name of Los Alamos National Security, LLC, Los Alamos
National Laboratory, LANL, the U.S. Government, nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


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

    if (vec->get_type() == typeid(double))
      status = mymesh_->store_field(name, entity_kind, (double *)vec->get_raw_data());
    else if (vec->get_type() == typeid(int))
      status = mymesh_->store_field(name, entity_kind, (int *)vec->get_raw_data());
    else if (vec->get_type() == typeid(std::array<double, 2>))
      status = mymesh_->store_field(name, entity_kind,
                                    (std::array<double, 2> *) vec->get_raw_data());
    else if (vec->get_type() == typeid(std::array<double, 3>))
      status = mymesh_->store_field(name, entity_kind,
                                    (std::array<double, 3> *) vec->get_raw_data());
    else if (vec->get_type() == typeid(std::array<double, 6>))
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
