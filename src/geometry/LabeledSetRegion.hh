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
 * @file   LabeledSetRegion.hh
 * @author Rao Garimella
 * @date
 *
 * @brief  Declaration of Labeled Set Region class which derives
 *         its definition from a named set of mesh entities in
 *         a mesh file
 *
 *
 */

#ifndef _LabeledSetRegion_hh_
#define _LabeledSetRegion_hh_

#include "Region.hh"

  namespace JaliGeometry {

// -------------------------------------------------------------
//  class LabeledSetRegion
// -------------------------------------------------------------
/// A region defined by a set of mesh entities in a mesh file
///
/// Strictly speaking, we should tie this region class to a particular
/// mesh or mesh file but that cause a circular dependency of meshes
/// on regions and of labeled set regions on meshes. We will rely on the
/// fact that when a mesh is created specifying a geometric model, it
/// will create mesh entity sets based on the labeled sets in that
/// geometric model.
///
/// If we need to change this behavior, then we can make a forward
/// declaration of Jali::Mesh, make the Mesh class a friend, add
/// a mesh variable to this class and have a protected method to set
/// the mesh

class LabeledSetRegion : public Region {
public:

  /// Default constructor

  LabeledSetRegion(const std::string name,
                   const unsigned int id,
                   const std::string entity_str,
                   const std::string file,
                   const std::string format,
                   const std::string label,
                   const LifeCycle_type lifecycle = LifeCycle_type::PERMANENT);

  /// Protected copy constructor to avoid unwanted copies.
  LabeledSetRegion(const LabeledSetRegion& old);

  /// Destructor
  ~LabeledSetRegion(void);

  // Type of the region
  inline Region_type type() const { return Region_type::LABELEDSET; }

  // Label in the file
  inline std::string label() const { return label_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

  inline std::string entity_str() const { return entity_str_; }

protected:
  const std::string entity_str_; // what kind of entities make up this set
  const std::string file_; // which file are we supposed to read it from
  const std::string format_; // format of the file
  const std::string label_; // Label used to identify set in the file (may be different from name)
};

/// A smart pointer to LabeledSetRegion instances
// typedef Teuchos::RCP<LabeledSetRegion> LabeledSetRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef LabeledSetRegion *LabeledSetRegionPtr;

} // namespace JaliGeometry


#endif
