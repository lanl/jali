#include <vector>

#include "UnitTest++.h"

#include "../Mesh_simple.hh"

// #include "State.hpp"

TEST(NUMBERING) {
  double expnodecoords[8][3] = {{0.0, 0.0, 0.0},
                                {1.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0},
                                {1.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0},
                                {1.0, 0.0, 1.0},
                                {0.0, 1.0, 1.0},
                                {1.0, 1.0, 1.0}};
  int expcellnodes[8] = {0, 1, 3, 2, 4, 5, 7, 6};
  int expfacenodes[6][4] = {{0, 1, 3, 2},
                            {4, 5, 7, 6},
                            {0, 1, 5, 4},
                            {2, 3, 7, 6},
                            {0, 2, 6, 4},
                            {1, 3, 7, 5}};

  int expcellfaces[6] = {2, 5, 3, 4, 0, 1};
  int expfacedirs[6] = {1, 1, -1, -1, -1, 1};



  // Create a single-cell mesh;
  Jali::Mesh *mesh(new Jali::Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                                         1, 1, 1,
                                         MPI_COMM_WORLD));

  //  State S(1,mesh);

  //  std::string gmvfile = "out.gmv";
  //  S.write_gmv(gmvfile);

  // Write node coordinates

  JaliGeometry::Point x;
  for (Jali::Entity_ID j = 0; j < 8; ++j) {
    mesh->node_get_coordinates(j, &x);
    CHECK_ARRAY_EQUAL(expnodecoords[j], x, 3);
  }

  // Write face-node connectivity
  Jali::Entity_ID_List fnode;
  for (Jali::Entity_ID j = 0; j < 6; ++j) {
    mesh->face_get_nodes(j, &fnode);
    CHECK_EQUAL(4, fnode.size());
    CHECK_ARRAY_EQUAL(expfacenodes[j], fnode, fnode.size());
  }

  // Write cell-node connectivity
  Jali::Entity_ID_List cnode;
  for (Jali::Entity_ID j = 0; j < 1; ++j) {
    mesh->cell_get_nodes(j, &cnode);
    CHECK_EQUAL(8, cnode.size());
    CHECK_ARRAY_EQUAL(expcellnodes, cnode, cnode.size());
  }

  // Write cell face-node connectivity
  //  Jali::Entity_ID cface[6];
  //  int fdir[6];
  Jali::Entity_ID_List cface;
  std::vector<int> fdir;
  mesh->cell_get_faces_and_dirs(0, &cface, &fdir);
  CHECK_ARRAY_EQUAL(expcellfaces, cface, 6);
  CHECK_ARRAY_EQUAL(expfacedirs, fdir, 6);
}


TEST(NUMBERING_1D) {
  // Make a 1d cartesian mesh from 0.0 to 1.0, with 2 nodes and 1 cell
  double expnodecoords[2][1] = {{0.0},
                                {1.0}};
  // indices of nodes for cell
  int expcellnodes[2] = {0, 1};
  // indices of nodes for faces
  int expfacenodes[2][1] = {{0}, {1}};

  // indices of faces for cell
  int expcellfaces[2] = {0, 1};
  // directions for faces of cell
  int expfacedirs[2] = {-1, 1};



  // Create a single-cell mesh;
  std::vector<double> node_coords = {0.0, 1.0};
  Jali::Mesh *mesh(new Jali::Mesh_simple(node_coords,
                                         MPI_COMM_WORLD));

  // Write node coordinates

  JaliGeometry::Point x;
  for (Jali::Entity_ID j = 0; j < 2; ++j) {
    mesh->node_get_coordinates(j, &x);
    CHECK_ARRAY_EQUAL(expnodecoords[j], x, 1);
  }

  // Write face-node connectivity
  Jali::Entity_ID_List fnode;
  for (Jali::Entity_ID j = 0; j < 2; ++j) {
    mesh->face_get_nodes(j, &fnode);
    CHECK_EQUAL(1, fnode.size());
    CHECK_ARRAY_EQUAL(expfacenodes[j], fnode, fnode.size());
  }

  // Write cell-node connectivity
  Jali::Entity_ID_List cnode;
  for (Jali::Entity_ID j = 0; j < 1; ++j) {
    mesh->cell_get_nodes(j, &cnode);
    CHECK_EQUAL(2, cnode.size());
    CHECK_ARRAY_EQUAL(expcellnodes, cnode, cnode.size());
  }

  // Write cell face-node connectivity
  Jali::Entity_ID_List cface;
  std::vector<int> fdir;
  mesh->cell_get_faces_and_dirs(0, &cface, &fdir);
  CHECK_ARRAY_EQUAL(expcellfaces, cface, 2);
  CHECK_ARRAY_EQUAL(expfacedirs, fdir, 2);
}


