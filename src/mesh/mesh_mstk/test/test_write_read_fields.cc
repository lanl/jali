#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

// Special class just for testing protected functions

namespace Jali {
class Mesh_MSTK_Test_Protected : public Mesh_MSTK {
 public:

  // Simplified constructor from filename

  Mesh_MSTK_Test_Protected(const char *filename, const MPI_Comm& incomm,
                           int spatial_dimension) :
      Mesh_MSTK(filename, incomm, spatial_dimension) {
  }

  // Simplified constructor that generates a 3D mesh

  Mesh_MSTK_Test_Protected(const double x0, const double y0, const double z0,
                           const double x1, const double y1, const double z1,
                           const int nx, const int ny, const int nz,
                           const MPI_Comm& comm) :
      Mesh_MSTK(x0,y0,z0,x1,y1,z1,nx,ny,nz,comm) {
  }


  // Destructor

  ~Mesh_MSTK_Test_Protected() {};

  // Write the mesh out to an exodus file

  void write_to_exodus_file(const std::string exodusfilename) const {
    Mesh_MSTK::write_to_exodus_file(exodusfilename);
  }

  // pass through functions for accessing protected functions

  void get_field_info(Entity_kind on_what, int *num,
                       std::vector<std::string> *varnames,
                       std::vector<std::string> *vartypes) const {
    Mesh_MSTK::get_field_info(on_what, num, varnames, vartypes);
  }

  template<class T>
  bool get_field(std::string field_name, Entity_kind on_what, T *data) const {
    return Mesh_MSTK::get_field(field_name, on_what, data);
  }

  bool store_field(std::string field_name, Entity_kind on_what, int *data) {
    return Mesh_MSTK::store_field(field_name, on_what, data);
  }

  bool store_field(std::string field_name, Entity_kind on_what, double *data) {
    return Mesh_MSTK::store_field(field_name, on_what, data);
  }

  template<long unsigned int N>
  bool store_field(std::string field_name, Entity_kind on_what,
                   std::array<double,N> *data) {
    return Mesh_MSTK::store_field(field_name, on_what, data);
  }

}; // end class Mesh_MSTK_Test_Protected
} // end namespace Jali


// Serial only for now

TEST(MSTK_WRITE_READ_FIELDS)
{

  std::string filename="hex_3x3x3_att.exo";

  // Generate a mesh consisting of 3x3x3 elements

  Jali::Mesh_MSTK_Test_Protected *outmesh =
      new Jali::Mesh_MSTK_Test_Protected(0.0,0.0,0.0,1.0,1.0,1.0,3,3,3,
                                         MPI_COMM_WORLD);

  int nv = 64;
  int nc = 27;

  // Create a field of doubles on cells

  double *cellval1_out = new double[nc];
  for (int i = 0; i < nc; i++)
    cellval1_out[i] = 2.0 + i*1.5; // make data general

  bool status;
  status = outmesh->store_field("cellval1",Jali::CELL,cellval1_out);
  CHECK(status);

  // Create a field of doubles on cells

  double *cellval2_out = new double[nc];
  for (int i = 0; i < nc; i++)
    cellval2_out[i] = 2 + i*2; // make data general

  status = outmesh->store_field("cellval2",Jali::CELL,cellval2_out);
  CHECK(status);


  // Create a field of std::array<double,3> on nodes

  std::array<double,3> *nodevec_out = new std::array<double,3>[nv];
  for (int i = 0; i < nv; i++)
    for (int j = 0; j < 3; j++)
      nodevec_out[i][j] = 0.4*i+0.1*j;

  status = outmesh->store_field("nodevec",Jali::NODE,nodevec_out);
  CHECK(status);

  // Create a field of doubles on nodes

  double *nodeval_out = new double[nv];
  for (int i = 0; i < nv; i++)
    nodeval_out[i] = 25.0*i;

  status = outmesh->store_field("nodeval",Jali::NODE,nodeval_out);
  CHECK(status);


  // Write the mesh out

  outmesh->write_to_exodus_file(filename);





  // Now read the same mesh back in

  Jali::Mesh_MSTK_Test_Protected *inmesh =
      new Jali::Mesh_MSTK_Test_Protected(filename.c_str(),MPI_COMM_WORLD,3);
  CHECK(inmesh);

  int space_dim = inmesh->space_dimension();
  CHECK_EQUAL(nv,inmesh->num_entities(Jali::NODE,Jali::ALL));
  CHECK_EQUAL(nc,inmesh->num_entities(Jali::CELL,Jali::ALL));

  // Query the number of mesh fields on each entity type - 2 on cells,
  // 1 on nodes

  int ncellfields = 0;
  std::vector< std::string> cellvarnames, cellvartypes;
  inmesh->get_field_info(Jali::CELL, &ncellfields, &cellvarnames, &cellvartypes);
  int nfound = 0;
  for (int i = 0; i < ncellfields; i++) {
    if (cellvarnames[i] == "cellval1") {
      CHECK_EQUAL("DOUBLE",cellvartypes[i]);
      nfound++;
    }
    else if (cellvarnames[i] == "cellval2") {
      CHECK_EQUAL("DOUBLE",cellvartypes[i]);
      nfound++;
    }
  }
  CHECK_EQUAL(2,nfound);

  int nnodefields = 0;
  std::vector<std::string> nodevarnames,nodevartypes;
  inmesh->get_field_info(Jali::NODE, &nnodefields, &nodevarnames, &nodevartypes);

  nfound = 0;
  for (int i = 0; i < nnodefields; i++) {
    if (nodevarnames[i] == "nodeval") {
      CHECK_EQUAL("DOUBLE",nodevartypes[i]);
      nfound++;
    }
    else if (nodevarnames[i] == "nodevec") {
      CHECK_EQUAL("VECTOR",nodevartypes[i]);
      nfound++;
    }
  }
  CHECK_EQUAL(2,nfound);

  double *cellval1_in = new double[nc];
  status = inmesh->get_field("cellval1",Jali::CELL,cellval1_in);
  CHECK(status);
  CHECK_ARRAY_EQUAL(cellval1_out,cellval1_in,nc);

  double *cellval2_in = new double[nc];
  status = inmesh->get_field("cellval2",Jali::CELL,cellval2_in);
  CHECK(status);
  CHECK_ARRAY_EQUAL(cellval2_out,cellval2_in,nc);

  std::array<double,3> *nodevec_in = new std::array<double,3>[nv];
  status = inmesh->get_field("nodevec",Jali::NODE,nodevec_in);
  for (int i = 0; i < nv; i++)
    CHECK_ARRAY_EQUAL(nodevec_out[i],nodevec_in[i],3);

  double *nodeval_in = new double [nv];
  status = inmesh->get_field("nodeval",Jali::NODE,nodeval_in);
  CHECK_ARRAY_EQUAL(nodeval_out,nodeval_in,nv);

  delete outmesh;
  delete inmesh;
  delete [] cellval1_out;
  delete [] cellval2_out;
  delete [] nodevec_out;
  delete [] nodeval_out;
  delete [] cellval1_in;
  delete [] cellval2_in;
  delete [] nodevec_in;
  delete [] nodeval_in;

}

