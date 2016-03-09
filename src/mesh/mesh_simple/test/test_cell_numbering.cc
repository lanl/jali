
#include "UnitTest++.h"

#include "../Mesh_simple.hh"

// #include "State.hpp"

TEST(NUMBERING) {


  double expnodecoords[8][3] = {{0.0,0.0,0.0},
                                {1.0,0.0,0.0},
                                {0.0,1.0,0.0},
                                {1.0,1.0,0.0},
                                {0.0,0.0,1.0},
                                {1.0,0.0,1.0},
                                {0.0,1.0,1.0},
                                {1.0,1.0,1.0}};
  int expcellnodes[8] = {0,1,3,2,4,5,7,6};
  int expfacenodes[6][4] = {{0,1,3,2},
                            {4,5,7,6},
                            {0,1,5,4},
                            {2,3,7,6},
                            {0,2,6,4},
                            {1,3,7,5}};

  int expcellfaces[6] = {2,5,3,4,0,1};
  int expfacedirs[6] = {1,1,-1,-1,-1,1};
    



  // Create a single-cell mesh;
  Jali::Mesh *mesh(new Jali::Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, MPI_COMM_WORLD));

  //  State S(1,mesh);

  //  std::string gmvfile = "out.gmv";
  //  S.write_gmv(gmvfile);

  // Write node coordinates

  JaliGeometry::Point x;
  for (Jali::Entity_ID j = 0; j < 8; ++j) {
    mesh->node_get_coordinates(j, &x);
    CHECK_ARRAY_EQUAL(expnodecoords[j],x,3);
  }

  // Write face-node connectivity
  Jali::Entity_ID_List fnode;
  for (Jali::Entity_ID j = 0; j < 6; ++j) {
    mesh->face_get_nodes(j, &fnode);
    CHECK_EQUAL(4,fnode.size());
    CHECK_ARRAY_EQUAL(expfacenodes[j],fnode,fnode.size());
  }

  // Write cell-node connectivity
  Jali::Entity_ID_List cnode;
  for (Jali::Entity_ID j = 0; j < 1; ++j) {
    mesh->cell_get_nodes(j, &cnode);
    CHECK_EQUAL(8,cnode.size());
    CHECK_ARRAY_EQUAL(expcellnodes,cnode,cnode.size());
  }

  // Write cell face-node connectivity
  //  Jali::Entity_ID cface[6];
  //  int fdir[6];
  Jali::Entity_ID_List cface;
  std::vector<int> fdir;
  mesh->cell_get_faces_and_dirs(0,&cface,&fdir);
  CHECK_ARRAY_EQUAL(expcellfaces,cface,6);
  CHECK_ARRAY_EQUAL(expfacedirs,fdir,6);

}
