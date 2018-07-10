/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <algorithm>
#include <typeinfo>

#include "JaliStateVector.h"

#include "Mesh.hh"    // jali mesh header
#include "MeshSet.hh"
#include "JaliState.h"

namespace Jali {

int state_get_num_materials(std::weak_ptr<State> state) {
  if (!state.expired()) {
    std::shared_ptr<State> state_shared_ptr = state.lock();
    return state_shared_ptr->num_materials();
  }
  return 0;
}

std::shared_ptr<MeshSet> state_get_material_set(std::weak_ptr<State> state,
                                                 int matindex) {
  if (!state.expired()) {
    std::shared_ptr<State> state_shared_ptr = state.lock();
    return state_shared_ptr->material_set(matindex);
  }
  return nullptr;
}

}  // namespace Jali
