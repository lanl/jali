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
 * @file   GeometricModel.cc
 * @author Rao V. Garimella
 * @date Mon Aug  1 10:05:25 2011
 *
 * @brief
 *
 *
 */

#include "GeometricModel.hh"
#include "errors.hh"

namespace JaliGeometry {

// -------------------------------------------------------------
//  class GeometricModel
// -------------------------------------------------------------

// -------------------------------------------------------------
// GeometricModel:: constructors / destructor
// -------------------------------------------------------------

// Constructor

  GeometricModel::GeometricModel(const unsigned int dim) :
    topo_dimension_(dim)
{
  if (dim != 2 && dim != 3) {
    std::cerr << "Only 2D and 3D domains are supported" << std::endl;
    throw std::exception();
  }
  Regions.clear();
}

// Copy constructor

GeometricModel::GeometricModel(const GeometricModel& old) :
  topo_dimension_(old.dimension())
{
  int i, nr;

  nr = old.Num_Regions();

  for (i = 0; i < nr; i++) {
    RegionPtr r = old.Region_i(i);
    Regions.push_back(r);
  }
}


// GeometricModel::GeometricModel(const unsigned int dim,
//                                Teuchos::ParameterList gm_params,
//                                const Epetra_MpiComm *comm,
//                                const VerboseObject *verbobj) :
//   topo_dimension_(dim), verbosity_obj_(verbobj)
// {

//   if (dim != 2 && dim != 3) {
//     Errors::Message mesg("Only 2D and 3D domains are supported");
//     Exceptions::Jali_throw(mesg);
//   }

//   const int region_id_offset = 59049; // arbitrary number to avoid clashing
//                                       // with IDs of LabeledSet regions
//   int ngregions = 0; // Number of regions

//   // Go through the parameter list and populate the geometric model with regions

//   for (Teuchos::ParameterList::ConstIterator i = gm_params.begin(); i != gm_params.end(); i++)
//     {
//       if (gm_params.isSublist(gm_params.name(i)))
//         {

//           // Region name - get that from parameter list

//           std::string region_name = gm_params.name(i);

//           // Region ID - our internal numerical identifier

//           unsigned int region_id = ++ngregions+region_id_offset;


//           // Extract sublist specifying region

//           const Teuchos::ParameterList &reg_spec = gm_params.sublist(gm_params.name(i));

//           // While the XML file does not prevent it, there should only be one
//           // specification of what the region looks like

//           unsigned int k = 0;
//           for (Teuchos::ParameterList::ConstIterator j = reg_spec.begin(); j != reg_spec.end(); j++, k++)
//             {

//               if (k > 1)
//                 {
//                   std::stringstream sstream;
//                   sstream << "ERROR: Region " << region_name <<
//                     " described in multiple ways";
//                   Errors::Message mesg(sstream.str());
//                   Exceptions::Jali_throw(mesg);
//                 }



//               // Shape of the region

//               std::string shape = reg_spec.name(j);


//               // Create the region

//               Jali::JaliGeometry::RegionPtr regptr =
//                 RegionFactory(region_name, region_id, reg_spec, dim, comm,
//                               verbosity_obj());


//               // Add it to the geometric model

//               Regions.push_back(regptr);

//             }

//         }
//       else
//         {
//           Errors::Message mesg("Error: Improper region specification");
//           Exceptions::Jali_throw(mesg);
//         }
//     }
// }


// Destructor

GeometricModel::~GeometricModel(void)
{

  // If a geometric model is deleted, we will not delete all the
  // the regions added to it because someone else may be holding
  // on to a pointer to the regions. For now, the top level routine
  // deleting the geometric model has to delete the regions first
  // to prevent a memory leak

  // Once we can get RegionFactory to work with Reference Counted
  // Pointers we can remove this comment

  Regions.clear();
}


// Constructor with Region List

GeometricModel::GeometricModel(const unsigned int dim,
                               const std::vector<RegionPtr>& in_Regions) :
  topo_dimension_(dim)
{
  Regions.clear();
  Regions = in_Regions;
}


// Add a Region

void GeometricModel::Add_Region(const RegionPtr& regptr)
{
  if (topo_dimension_ < regptr->dimension()) {
    Errors::Message mesg("Topological dimension of geometric model less than that of the region");
    Exceptions::Jali_throw(mesg);
  }

  Regions.push_back(regptr);
}


// Number of Regions

int GeometricModel::Num_Regions(void) const
{
  return Regions.size();
}


// Get the i'th Region

RegionPtr GeometricModel::Region_i(const int i) const
{
  return Regions[i];
}


// Get a region by its ID

RegionPtr GeometricModel::FindRegion(const int id) const
{

  // FIXME: Can't get this to compile
  // std::vector<RegionPtr>::iterator r;
  // for (r = Regions.begin(); r != Regions.end(); r++)
  //   {
  //     if (r->id() == id)
  //       return *r;
  //   }

  for (int i = 0; i < Regions.size(); i++)
    {
      RegionPtr r = Regions[i];
      if (r->id() == id)
        return r;
    }
  return NULL;
}


// Get a region by its name

RegionPtr GeometricModel::FindRegion(const std::string name) const
{

  // FIXME: Can't get this to compile
  //  std::vector<RegionPtr>::iterator r;
  //  for (r = Regions.begin(); r != Regions.end(); r++)
  //    {
  //      if (r->name() == name)
  //        return *r;
  //    }

  for (int i = 0; i < Regions.size(); i++)
    {
      RegionPtr r = Regions[i];
      if (r->name() == name)
        return r;
    }
  return NULL;
}


// Check if regions cover the domain extents.
// This will work perfectly for domains with rectangular regions
// but not so for other types of regions

bool GeometricModel::Rough_Check_Tiling(void) const
{
  return true;
}


} // namespace JaliGeometry
