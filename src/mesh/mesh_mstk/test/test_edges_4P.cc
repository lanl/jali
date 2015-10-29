#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

// Test edge functions in 2D

TEST(MSTK_EDGES_2D)
{

  int rank, size;

  int initialized;
  MPI_Initialized(&initialized);
  
  if (!initialized)
    MPI_Init(NULL,NULL);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  CHECK_EQUAL(4,size);

  // Generate a 4x4 quad mesh distributed over four processors
  
  bool request_faces = true, request_edges = true;
  Jali::Mesh *mesh(new Jali::Mesh_MSTK(0.0,0.0,2.0,1.0,4,4,
                                                   MPI_COMM_WORLD,
                                                   NULL,request_faces,
                                                   request_edges));

  // Check that we get the expected number of edges

  int ne_owned = mesh->num_entities(Jali::EDGE,
				    Jali::OWNED);
  int ne_all = mesh->num_entities(Jali::EDGE,
				  Jali::ALL);

  // This assumes a symmetric partitioning - not always the case with
  // ZOLTAN graph partitioning

  CHECK_EQUAL(24,ne_all);

  // In 2D, faces and edges are the same - so face global IDs and edge
  // global IDs for a cell must match

  int nc_owned = mesh->num_entities(Jali::CELL,
				    Jali::OWNED);

  for (int c = 0; c < nc_owned; ++c) {
    Jali::Entity_ID_List cedges, cfaces, fedges;
    std::vector<int> cfdirs, fedirs, cedirs;    

    mesh->cell_get_edges(c,&cedges);
    mesh->cell_get_faces_and_dirs(c,&cfaces,&cfdirs);

    for (int e = 0; e < cedges.size(); ++e) {
      CHECK_EQUAL(mesh->GID(cedges[e],Jali::EDGE), 
		  mesh->GID(cfaces[e],Jali::FACE));
    }


    for (int f = 0; f < cfaces.size(); ++f) {
      mesh->face_get_edges_and_dirs(cfaces[f],&fedges,&fedirs);

      CHECK_EQUAL(1,fedges.size()); // face is same as edge in 2D
      CHECK_EQUAL(1,fedirs[0]); // direction is always 1
      
      // check the face-edges to cell-edges map

      std::vector<int> map;

      mesh->face_to_cell_edge_map(cfaces[f],c,&map);

      for (int e = 0; e < fedges.size(); ++e)
	CHECK_EQUAL(fedges[e],cedges[map[e]]);
    }
  }

  // owing to how we constructed the mesh, the length of horizontal edges 
  // should be 0.5 and vertical edges 0.25

  for (int e = 0; e < ne_owned; ++e) {
    JaliGeometry::Point evec(2);
    double elen;

    evec = mesh->edge_vector(e);
    elen = mesh->edge_length(e);
    if (evec[1] == 0.0) {
      CHECK_EQUAL(0.5,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
    else if (evec[0] == 0.0) {
      CHECK_EQUAL(0.25,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
  }

}



// Test edge functions in 3D

TEST(MSTK_EDGES_3D)
{

  int rank, size;

  int initialized;
  MPI_Initialized(&initialized);
  
  if (!initialized)
    MPI_Init(NULL,NULL);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  CHECK_EQUAL(4,size);

  // Generate a 4x4x4 quad mesh distributed over four processors
  
  bool request_faces = true, request_edges = true;
  Jali::Mesh *mesh(new Jali::Mesh_MSTK(0.0,0.0,0.0,2.0,1.0,4.0,4,4,4,MPI_COMM_WORLD,NULL,request_faces,request_edges));

  // How many owned and used edges are there?

  int ne_owned = mesh->num_entities(Jali::EDGE,Jali::OWNED);
  int ne_all = mesh->num_entities(Jali::EDGE,Jali::ALL);

  // Check that we got a non-zero number

  CHECK(ne_owned != 0);
  CHECK(ne_all != 0);  


  // Go through the cells and retrieve their edges to make sure it
  // works correctly. Also, get the faces of the cells and the edges
  // of these faces and do additional checks

  int nc_owned = mesh->num_entities(Jali::CELL,Jali::OWNED);

  for (int c = 0; c < nc_owned; ++c) {
    Jali::Entity_ID_List cedges, cfaces, fedges;
    std::vector<int> cfdirs, fedirs;    

    mesh->cell_get_edges(c,&cedges);
    mesh->cell_get_faces_and_dirs(c,&cfaces,&cfdirs);

    for (int f = 0; f < cfaces.size(); ++f) {
      mesh->face_get_edges_and_dirs(cfaces[f],&fedges,&fedirs);

      // check the face-edges to cell-edges map

      std::vector<int> map;
      mesh->face_to_cell_edge_map(cfaces[f],c,&map);

      for (int e = 0; e < fedges.size(); ++e)
	CHECK_EQUAL(fedges[e],cedges[map[e]]);
    }
  }

  // owing to how we constructed the mesh, the length of x-direction
  // should be 0.5, y-direction edges should
  // 0.25 and z-direction edges should be 1.0

  for (int e = 0; e < ne_owned; ++e) {
    JaliGeometry::Point evec(2);
    double elen;

    evec = mesh->edge_vector(e);
    elen = mesh->edge_length(e);
    if (evec[0] != 0.0  && evec[1] == 0.0) {
      CHECK_EQUAL(0.5,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
    else if (evec[0] == 0.0 && evec[1] != 0.0) {
      CHECK_EQUAL(0.25,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
    else  {
      CHECK_EQUAL(1.0,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
  }

}

