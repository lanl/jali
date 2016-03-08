#include <iostream>

#include "mpi.h"

#include "dbc.hh"
#include "errors.hh"
#include "Geometry.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"

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
  Entity_ID_List expected_cells[3] = {{0,1,2,3,4,5,6,7,8},
                                      {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17},
                                      {9,10,11,12,13,14,15,16,17}};

  std::string facesetnames[3] = {"Face 101", "ZLO FACE Plane", "YLO FACE Box"};
  int nfacesets = 3;
  Entity_ID_List expected_faces[3] = {{4,19,9,32,23,14,36,27,40},
                                      {4,9,14,19,23,27,32,36,40},
                                      {0,6,11,42,47,51,75,80,84}};
  

  std::string nodesetnames[1] = {"INTERIOR XY PLANE"};
  int nnodesets = 1;
  Entity_ID_List expected_nodes[1] = {{16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}};
  

  // Jali depends on MPI 

  MPI_Init(&argc,&argv);
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

  gregions.emplace_back(new LabeledSetRegion(cellsetnames[0],1,"CELL",
                                             exofilename,"Exodus II","10000"));

  // Box Region

  Point boxlo(0.0,0.0,0.0), boxhi(1.0,1.0,0.66);
  gregions.emplace_back(new BoxRegion(cellsetnames[1],2,boxlo,boxhi));

  // Logical region from two other regions - Extract cells of the
  // bottom layer by subtracting the middle layer from the
  // bottom+middle layer

  std::vector<std::string> setnames;
  setnames.push_back(cellsetnames[1]);  // Bottom+Middle box
  setnames.push_back(cellsetnames[0]);  // Bottom Labeled Set
  gregions.emplace_back(new LogicalRegion(cellsetnames[2],99,"Subtract",
                                          setnames));
  

  // Labeled set - an Exodus side set (face set) with ID 101.

  gregions.emplace_back(new LabeledSetRegion(facesetnames[0],5,"FACE",
                                             exofilename,"Exodus II","101"));

  // Plane region - designed to capture faces on the bottom surface of the mesh

  Point planepoint(0.0,0.0,0.0), planenormal(0.0,0.0,1.0);
  gregions.emplace_back(new PlaneRegion(facesetnames[1],6,planepoint,
                                        planenormal));

  // Box region - degenerate box to capture faces on the left surface of the mesh

  boxlo.set(0.0,0.0,0.0);
  boxhi.set(1.0,0.0,1.0);
  gregions.emplace_back(new BoxRegion(facesetnames[2],7,boxlo,boxhi));


  // Interior XY plane - will be used to query nodes on a plane

  planepoint.set(0.0,0.0,1.0/3.0);
  planenormal.set(0.0,0.0,1.0);
  gregions.emplace_back(new PlaneRegion(nodesetnames[0],8,planepoint,planenormal));


  // Create a geometric model of spatial dimension 3

  GeometricModelPtr gm(new GeometricModel(3,gregions));


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

  FrameworkPreference pref;
  pref.push_back(MSTK);

  std::unique_ptr<Mesh> mymesh;  // Pointer to a mesh object
  if (framework_available(MSTK)) {  // check if framework is available
    mesh_factory.preference(pref);  
  
    // Read in an exodus file. Request faces, edges, wedges and
    // corners (true, true, true, true). MAKE SURE TO SPECIFY THE
    // GEOMETRIC MODEL

    mymesh = mesh_factory(exofilename,gm,true,true,true,true);
  }

  //---------------------------------------------------------------------------
  // CREATE MESH OBJECT BY READING FROM FILE
  //---------------------------------------------------------------------------
  


  // Print out the spatial dimension of the mesh

  std::cerr << "Spatial dimension of mesh: " << mymesh->space_dimension() << std::endl;
  
  // Print out the number of cells in the mesh

  std::cerr << "Number of mesh cells: " <<
      mymesh->num_cells<Parallel_type::ALL>() << std::endl;

  // Print out the number of nodes in the mesh

  std::cerr << "Number of mesh nodes: " <<
      mymesh->num_nodes<Parallel_type::ALL>() << std::endl;


  

  //-------------------------------------------------------------------------
  // Query the various sets
  //-------------------------------------------------------------------------

  for (int i = 0; i < ncellsets; ++i) {
    int ncells = mymesh->get_set_size(cellsetnames[i], Entity_kind::CELL,
                                      Parallel_type::OWNED);
    if (ncells != expected_cells[i].size()) {
      std::cerr << "Wrong number of cells in cell set " << cellsetnames[i] <<
          "(Expected " << expected_cells[i].size() << " Got " << ncells << 
          ")" << std::endl;
      continue;
    }

    Entity_ID_List cellids;
    mymesh->get_set_entities(cellsetnames[i], Entity_kind::CELL,
                             Parallel_type::OWNED, &cellids);
    for (int j = 0; j < ncells; ++j) {
      if (expected_cells[i][j] != cellids[j]) {
        std::cerr << "Mismatch in expected and retrieved cell for cell set " <<
            cellsetnames[i] << std::endl;
        break;
      }
    }

    std::cerr << "Retrieved cell set " << cellsetnames[i] << " containing " <<
        ncells << " cells  -- ";
    for (int j = 0; j < ncells; ++j)
      std::cerr << cellids[j] << " ";
    std::cerr << std::endl << std::endl;
  }

  for (int i = 0; i < nfacesets; ++i) {
    int nfaces = mymesh->get_set_size(facesetnames[i], Entity_kind::FACE,
                                      Parallel_type::OWNED);
    if (nfaces != expected_faces[i].size()) {
      std::cerr << "Wrong number of faces in face set " << facesetnames[i] <<
          "(Expected " << expected_faces[i].size() << " Got " << nfaces << 
          ")" << std::endl;
      continue;
    }

    Entity_ID_List faceids;
    mymesh->get_set_entities(facesetnames[i], Entity_kind::FACE,
                             Parallel_type::OWNED, &faceids);
    for (int j = 0; j < nfaces; ++j) {
      if (expected_faces[i][j] != faceids[j]) {
        std::cerr << "Mismatch in expected and retrieved face for face set " <<
            facesetnames[i] << std::endl;
        break;
      }
    }

    std::cerr << "Retrieved face set " << facesetnames[i] << " containing " <<
        nfaces << " faces  -- ";
    for (int j = 0; j < nfaces; ++j)
      std::cerr << faceids[j] << " ";
    std::cerr << std::endl << std::endl;
  }

  for (int i = 0; i < nnodesets; ++i) {
    int nnodes = mymesh->get_set_size(nodesetnames[i], Entity_kind::NODE,
                                      Parallel_type::OWNED);
    if (nnodes != expected_nodes[i].size()) {
      std::cerr << "Wrong number of nodes in node set " << nodesetnames[i] <<
          "(Expected " << expected_nodes[i].size() << " Got " << nnodes << 
          ")" << std::endl;
      continue;
    }

    Entity_ID_List nodeids;
    mymesh->get_set_entities(nodesetnames[i], Entity_kind::NODE,
                             Parallel_type::OWNED, &nodeids);
    for (int j = 0; j < nnodes; ++j) {
      if (expected_nodes[i][j] != nodeids[j]) {
        std::cerr << "Mismatch in expected and retrieved node for node set " <<
            nodesetnames[i] << std::endl;
        break;
      }
    }

    std::cerr << "Retrieved node set " << nodesetnames[i] << " containing " <<
        nnodes << " nodes  -- ";
    for (int j = 0; j < nnodes; ++j)
      std::cerr << nodeids[j] << " ";
    std::cerr << std::endl << std::endl;
  }


  // Clean up and exit

  MPI_Finalize();

}

