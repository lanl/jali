/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "JaliState.h"

namespace Jali {

//! Print all state vectors

std::ostream & operator<<(std::ostream & os, State const & s) {
  State::const_iterator it = s.cbegin();
  while (it != s.cend()) {
    const std::shared_ptr<BaseStateVector> vec = *it;
    vec->print(os);
    ++it;
  }
}

} // namespace Portage
