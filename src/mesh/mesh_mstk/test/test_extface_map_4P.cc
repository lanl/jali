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
#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"
#include "MeshAudit.hh"

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"


TEST(MSTK_EXTFACE_MAP_4P)
{

  int i, j, k, err, nc, nf, nv;
  std::vector<Jali::Entity_ID> faces(6), nodes(8);
  std::vector<int> facedirs(6);
  std::vector<Jali::JaliGeometry::Point> ccoords(8), fcoords(4);

  Teuchos::RCP<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
			

  int rank, size;

  int initialized;
  MPI_Initialized(&initialized);

  if (!initialized)
    MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  CHECK_EQUAL(4, size);

  Teuchos::RCP<Jali::Mesh> mesh(new Jali::Mesh_MSTK("test/hex_3x3x3_sets.exo",
                                                    comm.get()));

  Epetra_Map face_map(mesh->face_map(false));
  Epetra_Map extface_map(mesh->exterior_face_map());

  Epetra_Import all_to_extface_importer = mesh->exterior_face_importer();

  for (int f = extface_map.MinLID(); f <= extface_map.MaxLID(); f++)
    {
      int gid = extface_map.GID(f);
      int f2 = face_map.LID(gid); // f2 is local face id in face_map

      CHECK_EQUAL(face_map.GID(f2), gid);

      Jali::Entity_ID_List fcells;
      mesh->face_get_cells(f2, Jali::Entity_type::PARALLEL_OWNED, &fcells);
      CHECK_EQUAL(1, fcells.size());
    }

  Epetra_Vector allvec(face_map);
  Epetra_Vector bdryvec(extface_map);

  // Insert the GlobalID of each face offsetted by 3 into the allvec

  for (int f = face_map.MinLID(); f < face_map.MaxLID(); f++)
      allvec[f] = face_map.GID(f)+3;

  bdryvec.Import(allvec, all_to_extface_importer, Insert);

  // Check if the importer got the right values from allvec into bdryvec
  // by checking if the values in the bdryvec minus the offset correspond
  // to the correct global IDs.

  for (int f = extface_map.MinLID(); f < extface_map.MaxLID(); f++)
    CHECK_EQUAL(extface_map.GID(f), bdryvec[f]-3);

}

