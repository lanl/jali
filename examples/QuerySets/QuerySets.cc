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


#include <iostream>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Geometry.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

// Once we resurrect RegionFactory, you won't have to include the
// headers for the individual types of regions

#include "LabeledSetRegion.hh"
#include "BoxRegion.hh"
#include "PlaneRegion.hh"
#include "PointRegion.hh"
#include "LogicalRegion.hh"

using namespace Jali;
using namespace JaliGeometry;

// Fire up Jali, read in a mesh, query entity sets based on various
// geometric regions. The mesh has 3x3x3 hexahedral cells and the
// domain has three element blocks with 9 elements each - a bottom
// layer, middle layer and top layer. The Exodus II file being read in
// already has some of these cell sets tagged. It also has some side
// sets (face sets) tagged.

int main(int argc, char *argv[]) {

  std::string exofilename = "hex_3x3x3_sets.exo";

  // Some sets that we want to query

  std::string cellsetnames[3] = {"Bottom LabeledSet", "Bottom+Middle Box",
                                 "Middle Layer"};
  int ncellsets = 3;
  Entity_ID_List expected_cells[3] = {{0, 1, 2, 3, 4, 5, 6, 7, 8},
                                      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                       12, 13, 14, 15, 16, 17},
                                      {9, 10, 11, 12, 13, 14, 15, 16, 17}};

  std::string facesetnames[3] = {"Face 101", "ZLO FACE Plane", "YLO FACE Box"};
  int nfacesets = 3;
  Entity_ID_List expected_faces[3] = {{4, 19, 9, 32, 23, 14, 36, 27, 40},
                                      {4, 9, 14, 19, 23, 27, 32, 36, 40},
                                      {0, 6, 11, 42, 47, 51, 75, 80, 84}};
  

  std::string nodesetnames[1] = {"INTERIOR XY PLANE"};
  int nnodesets = 1;
  Entity_ID_List expected_nodes[1] = {{16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                                       26, 27, 28, 29, 30, 31}};
  

  // Jali depends on MPI

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;



  //------------------------------------------------------------------------
  // CREATE GEOMETRIC REGIONS AND ADD TO GEOMETRIC MODEL
  // -----------------------------------------------------------------------

  // Make the various geometric regions and add them to the geometric model

  std::vector<RegionPtr> gregions;
  
  // NOTE: Exodus can tag sets of elements using Element Blocks or
  // Element Sets.  The difference is that elements may belong to only
  // one element block but can belong to any number of element
  // sets. Jali will look for element blocks with the given label
  // first and if not found, it will look for element sets. For that
  // reason, a mesh should not have element blocks and element sets
  // with the same ID

  // Labeled set - an Exodus element block with ID 10000. The second
  // argument is a set ID and can be any unique number.

  gregions.emplace_back(new LabeledSetRegion(cellsetnames[0], 1, "CELL",
                                             exofilename, "Exodus II", "10000"));

  // Box Region

  Point boxlo(0.0, 0.0, 0.0), boxhi(1.0, 1.0, 0.66);
  gregions.emplace_back(new BoxRegion(cellsetnames[1], 2, boxlo, boxhi));

  // Logical region from two other regions - Extract cells of the
  // bottom layer by subtracting the middle layer from the
  // bottom+middle layer

  std::vector<std::string> setnames;
  setnames.push_back(cellsetnames[1]);  // Bottom+Middle box
  setnames.push_back(cellsetnames[0]);  // Bottom Labeled Set
  gregions.emplace_back(new LogicalRegion(cellsetnames[2], 99, "Subtract",
                                          setnames));
  

  // Labeled set - an Exodus side set (face set) with ID 101.

  gregions.emplace_back(new LabeledSetRegion(facesetnames[0], 5, "FACE",
                                             exofilename, "Exodus II", "101"));

  // Plane region - designed to capture faces on the bottom surface of the mesh

  Point planepoint(0.0, 0.0, 0.0), planenormal(0.0, 0.0, 1.0);
  gregions.emplace_back(new PlaneRegion(facesetnames[1], 6, planepoint,
                                        planenormal));

  // Box region - degenerate box to capture faces on the left surface
  // of the mesh

  boxlo.set(0.0, 0.0, 0.0);
  boxhi.set(1.0, 0.0, 1.0);
  gregions.emplace_back(new BoxRegion(facesetnames[2], 7, boxlo, boxhi));


  // Interior XY plane - will be used to query nodes on a plane

  planepoint.set(0.0, 0.0, 1.0/3.0);
  planenormal.set(0.0, 0.0, 1.0);
  gregions.emplace_back(new PlaneRegion(nodesetnames[0], 8, planepoint,
                                        planenormal));


  // Create a geometric model of spatial dimension 3

  GeometricModelPtr gm(new GeometricModel(3, gregions));


  //---------------------------------------------------------------------------
  // DONE CREATING REGIONS AND ADDING TO GEOMETRIC MODEL
  //---------------------------------------------------------------------------



  //---------------------------------------------------------------------------
  // CREATE MESH OBJECT BY READING FROM FILE
  //---------------------------------------------------------------------------
  
  // Create a mesh factory object - this object has methods for
  // specifying the preference of mesh frameworks and unified
  // interfaces for instantiating a mesh object of a particular
  // framework type

  MeshFactory mesh_factory(comm);

  // Specify that MSTK is the preferred mesh framework. Currently Jali is
  // compiled only with MSTK support

  std::shared_ptr<Mesh> mymesh;  // Pointer to a mesh object
  bool parallel_file = false;
  if (framework_available(MSTK) &&
      framework_reads(MSTK, parallel_file, ExodusII)) {

    mesh_factory.framework(MSTK);
  
    // Read in an exodus file. Request faces (in addition to the
    // default cells and nodes)

    mesh_factory.included_entities(Entity_kind::FACE);

    // MAKE SURE TO SPECIFY THE GEOMETRIC MODEL
    
    mesh_factory.geometric_model(gm);

    mymesh = mesh_factory(exofilename);
  }

  //---------------------------------------------------------------------------
  // CREATE MESH OBJECT BY READING FROM FILE
  //---------------------------------------------------------------------------
  


  // Print out the spatial dimension of the mesh

  std::cerr << "Spatial dimension of mesh: " << mymesh->space_dimension() <<
      std::endl;
  
  // Print out the number of cells in the mesh

  std::cerr << "Number of mesh cells: " <<
      mymesh->num_cells<Entity_type::ALL>() << std::endl;

  // Print out the number of nodes in the mesh

  std::cerr << "Number of mesh nodes: " <<
      mymesh->num_nodes<Entity_type::ALL>() << std::endl;


  

  //-------------------------------------------------------------------------
  // Query the various sets
  //-------------------------------------------------------------------------

  for (int i = 0; i < ncellsets; ++i) {
    int ncells = mymesh->get_set_size(cellsetnames[i], Entity_kind::CELL,
                                      Entity_type::PARALLEL_OWNED);
    if (ncells != expected_cells[i].size()) {
      std::cerr << "Wrong number of cells in cell set " << cellsetnames[i] <<
          "(Expected " << expected_cells[i].size() << " Got " << ncells <<
          ")" << std::endl;
      continue;
    }

    Entity_ID_List cellids;
    mymesh->get_set_entities(cellsetnames[i], Entity_kind::CELL,
                             Entity_type::PARALLEL_OWNED, &cellids);
    int j = 0;
    for (auto const & c : cellids) {
      if (expected_cells[i][j] != c) {
        std::cerr << "Mismatch in expected and retrieved cell for cell set " <<
            cellsetnames[i] << std::endl;
        break;
      }
      j++;
    }

    std::cerr << "Retrieved cell set " << cellsetnames[i] << " containing " <<
        ncells << " cells  -- ";
    for (auto const & c : cellids)
      std::cerr << c << " ";
    std::cerr << std::endl << std::endl;
  }

  for (int i = 0; i < nfacesets; ++i) {
    int nfaces = mymesh->get_set_size(facesetnames[i], Entity_kind::FACE,
                                      Entity_type::PARALLEL_OWNED);
    if (nfaces != expected_faces[i].size()) {
      std::cerr << "Wrong number of faces in face set " << facesetnames[i] <<
          "(Expected " << expected_faces[i].size() << " Got " << nfaces <<
          ")" << std::endl;
      continue;
    }

    Entity_ID_List faceids;
    mymesh->get_set_entities(facesetnames[i], Entity_kind::FACE,
                             Entity_type::PARALLEL_OWNED, &faceids);
    int j = 0;
    for (auto const & f : faceids) {
      if (expected_faces[i][j] != f) {
        std::cerr << "Mismatch in expected and retrieved face for face set " <<
            facesetnames[i] << std::endl;
        break;
      }
      j++;
    }

    std::cerr << "Retrieved face set " << facesetnames[i] << " containing " <<
        nfaces << " faces  -- ";
    for (auto const & f : faceids)
      std::cerr << f << " ";
    std::cerr << std::endl << std::endl;
  }

  for (int i = 0; i < nnodesets; ++i) {
    int nnodes = mymesh->get_set_size(nodesetnames[i], Entity_kind::NODE,
                                      Entity_type::PARALLEL_OWNED);
    if (nnodes != expected_nodes[i].size()) {
      std::cerr << "Wrong number of nodes in node set " << nodesetnames[i] <<
          "(Expected " << expected_nodes[i].size() << " Got " << nnodes <<
          ")" << std::endl;
      continue;
    }

    Entity_ID_List nodeids;
    mymesh->get_set_entities(nodesetnames[i], Entity_kind::NODE,
                             Entity_type::PARALLEL_OWNED, &nodeids);
    int j = 0;
    for (auto const & n : nodeids) {
      if (expected_nodes[i][j] != n) {
        std::cerr << "Mismatch in expected and retrieved node for node set " <<
            nodesetnames[i] << std::endl;
        break;
      }
      j++;
    }

    std::cerr << "Retrieved node set " << nodesetnames[i] << " containing " <<
        nnodes << " nodes  -- ";
    for (auto const & n : nodeids)
      std::cerr << n << " ";
    std::cerr << std::endl << std::endl;
  }


  // Clean up and exit

  MPI_Finalize();

}

