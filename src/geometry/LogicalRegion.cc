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
/**
 * @file   LogicalRegion.cc
 * @author Rao Garimella
 * @date
 *
 * @brief  Implementation of Labeled Set Region class which derives its
 *         definition from named set of mesh entities in a mesh file
 *
 *
 */

#include "LogicalRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace JaliGeometry {

// -------------------------------------------------------------
//  class LogicalRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// LogicalRegion:: constructors / destructor
// -------------------------------------------------------------
LogicalRegion::LogicalRegion(const std::string name,
                             const unsigned int id,
                             const std::string operation_str,
                             const std::vector<std::string> region_names,
                             const LifeCycle_type lifecycle)
    : Region(name,id,3,lifecycle), operation_(Bool_type::NOBOOLEAN),
    region_names_(region_names)
{
  // Region dimension is set arbitrarily as 3 since the set of
  // entities in the mesh will determine the dimension

  if (operation_str == "Complement")
    operation_ = Bool_type::COMPLEMENT;
  else if (operation_str == "Union")
    operation_ = Bool_type::UNION;
  else if (operation_str == "Intersect")
    operation_ = Bool_type::INTERSECT;
  else if (operation_str == "Subtract")
    operation_ = Bool_type::SUBTRACT;
  else {
    Errors::Message mesg("Unknown logical operation type requested on regions");
    Exceptions::Jali_throw(mesg);
  }

}

LogicalRegion::LogicalRegion(const std::string name,
                             const unsigned int id,
                             const Bool_type bool_op_type,
                             const std::vector<std::string> region_names,
                             const LifeCycle_type lifecycle)
    : Region(name,id,3,lifecycle), operation_(bool_op_type),
    region_names_(region_names)
{
}

LogicalRegion::LogicalRegion(const LogicalRegion& old)
  : Region(old)
{
  // empty
}

LogicalRegion::~LogicalRegion(void)
{
  // empty
}


// -------------------------------------------------------------
// LogicalRegion::inside
// -------------------------------------------------------------
bool
LogicalRegion::inside(const Point& p) const
{
  Errors::Message mesg("In/out check not implemented for logical regions \n because the check may not be implemented for one of its component regions");
  Exceptions::Jali_throw(mesg);
  return false;
}

} // namespace JaliGeometry
