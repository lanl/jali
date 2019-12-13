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

#ifndef JALI_STATE_VECTOR_H_
#define JALI_STATE_VECTOR_H_

#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <algorithm>
#include <typeinfo>

#include "Mesh.hh"    // jali mesh header

namespace Jali {

enum class StateVector_type {UNIVAL, MULTIVAL};
enum class Data_layout {CELL_CENTRIC, MATERIAL_CENTRIC};

// Forward declaration of State class and some functions to resolve
// circular dependency (cannot include JaliState.h or use methods of
// the State class). The functions are defined in JaliStateVector.cc

class State;
int state_get_num_materials(std::weak_ptr<State> state);
std::shared_ptr<MeshSet> state_get_material_set(std::weak_ptr<State> state,
                                                 int matindex);

/*!
  @class StateVectorBase jali_state_vector.h
  @brief StateVectorBase provides a base class for state vectors on meshes, mesh tiles or mesh subsets

  A StateVector class inherits from this class and is templated on
  specific types of data. The DomainType class has to be of type mesh
  or a type that is part of a mesh and therefore, can answer the
  question "What is your mesh" and "How many entities of a particular kind
  (CELL, NODE, etc) do you have" e.g. MeshTile or MeshSet?
*/

class StateVectorBase {
 public:

  /*!
    @brief Constructor with a name
    @param name   Name of the StateVector
    @param state  The state manager holding this vector (can be nullptr)
    @param kind   Kind of entities it is defined on (CELL, NODE, etc.)
    @param type   Type of entities it is defined on (PARALLEL_OWNED, etc.)
  */
  
  explicit StateVectorBase(std::string name,
                           std::shared_ptr<State> state,
                           Entity_kind kind,
                           Entity_type type) :
      mystate_(state), myname_(name), entity_kind_(kind), entity_type_(type) {}


  //! Destructor

  virtual ~StateVectorBase() {}

  //! Virtual methods

  virtual std::ostream & print(std::ostream & os) const {
    os << "Print not implemented for data type of StateVector\n";
    return os;
  }
  virtual size_t size() const = 0;

  virtual const std::type_info& data_type() = 0;
  virtual StateVector_type type() = 0;

  //! Query Metadata

  std::string name() const { return myname_; }

  /// what kind of entity does it live on (CELL, WEDGE, NODE)?

  Entity_kind entity_kind() const { return entity_kind_; }

  /// What type of entity does it live on (PARALLEL_OWNED, PARALLEL_GHOST or
  /// BOUNDARY_GHOST or ALL)?

  Entity_type entity_type() const { return entity_type_; }

 protected:
  std::string myname_;
  Entity_kind entity_kind_;
  Entity_type entity_type_;

  // Circular dependency between state vectors and state manager.
  // JaliState, the state manager, will contain a bunch of shared
  // pointers to BaseVectors. But some types of derived state vectors
  // like MultiStateVectors need to know which JaliState object they
  // belong to so that they inquire about the number of materials, the
  // size and membership of material sets, etc. So keep a weak_ptr to
  // JaliState

  std::weak_ptr<State> mystate_;
};


///////////////////////////////////////////////////////////////////////////////

/*!
  @class UniStateVectorBase jali_state_vector.h

  @brief UniStateVectorBase   Base vector for univalued state vector
  defined on entities in a mesh, mesh tile or mesh subset. This vector
  stores only metadata and not the actual data so that some general
  functions can be called on it without knowing the type of data stored

  @tparam DomainType  Mesh, Mesh Tile or Mesh Subset
*/

template <class DomainType = Mesh>
class UniStateVectorBase : public StateVectorBase {
 public:

  //! Default constructor
  UniStateVectorBase() : StateVectorBase("UninitializedVector", nullptr,
                                         Entity_kind::UNKNOWN_KIND,
                                         Entity_type::TYPE_UNKNOWN),
                         mydomain_(nullptr)  {}

  /*!
    @brief Meaningful constructor with data and a string identifier
    @param name            String identifier of vector
    @param state           State manager holding the vector (can be nullptr)
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
  */
  
  UniStateVectorBase(std::string name,
                     std::shared_ptr<DomainType> domain,
                     std::shared_ptr<State> state,
                     Entity_kind kind,
                     Entity_type type) :
      StateVectorBase(name, state, kind, type), mydomain_(domain) {}


  /*! 
    @brief Copy constructor (copies metadata)

    Since mystate_ is a weak_ptr, we have to lock it to get a shared_ptr
    to send to the UniStateVectorBase constructor
  */

  UniStateVectorBase(UniStateVectorBase const & in_vector) :
      StateVectorBase(in_vector.myname_, in_vector.mystate_.lock(),
                      in_vector.entity_kind_, in_vector.entity_type_),
      mydomain_(in_vector.mydomain_)
  {}

  /*!
    @brief Assignment operator (copies metadata)
  */

  UniStateVectorBase & operator=(UniStateVectorBase const & in_vector) {
    UniStateVectorBase::myname_ = in_vector.myname_;
    UniStateVectorBase::mystate_ = in_vector.mystate_;
    UniStateVectorBase::entity_kind_ = in_vector.entity_kind_;
    UniStateVectorBase::entity_type_ = in_vector.entity_type_;
    mydomain_ = in_vector.mydomain_;

    return *this;
  }

  /// Destructor
  
  ~UniStateVectorBase() {}

  /// Domain on which UniStateVectorBase is defined (Mesh, MeshTile,
  /// MeshSubset)
  
  std::shared_ptr<DomainType> domain() const {return mydomain_;}

  /// Underlying mesh regardless of what type of domain UniStateVector is
  /// defined on. We have to return a reference to the mesh rather
  /// than a shared pointer because if the DomainType is a MeshTile,
  /// it only has a mesh reference not a pointer to the mesh

  Mesh & mesh() const { return get_mesh_of_domain(mydomain_); }

  /// Get the type of state vector

  StateVector_type type() {return StateVector_type::UNIVAL;}

  /// Number of entries in the vector
  virtual size_t size() const = 0;

  /// Change number of entries in the vector
  virtual void resize(size_t new_size) = 0;

  /// Clear the vector -> number of entries will become 0
  virtual void clear() = 0;

 protected:
  std::shared_ptr<DomainType> mydomain_;

  const Mesh & get_mesh_of_domain(std::shared_ptr<MeshTile> meshtile) const {
     return meshtile->mesh();
  }
  Mesh & get_mesh_of_domain(std::shared_ptr<Mesh> mesh) const {
    return *mesh;
  }
};  // UniStateVectorBase


///////////////////////////////////////////////////////////////////////////////



/*!
  @class UniStateVector jali_state_vector.h
  @brief UniStateVector stores univalued state data for entities in a mesh, mesh tile or mesh subset (one value per entity)

  Provides some limited functionality of a std::vector while adding
  some additional meta-data like the mesh associated with this data.

  @tparam T           Data type (int, double, some_custom_type)
  @tparam DomainType  Mesh, Mesh Tile or Mesh Subset
*/

template <class T, class DomainType = Mesh>
class UniStateVector : public UniStateVectorBase<DomainType> {
 public:

  //! Default constructor - not to be used
  UniStateVector() : UniStateVectorBase<DomainType>(), mydata_(nullptr)
  {}

  
  /*!
    @brief Constructor with array data
    @param name            Name of vector
    @param state           State manager holding the vector (can be nullptr)
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
    @param data            Pointer to array data to be used to initialize vector (optional)

    Since data is optional, this can be used to create an uninitialized vector
  */
  
  UniStateVector(std::string name,
              std::shared_ptr<DomainType> domain,
              std::shared_ptr<State> state,
              Entity_kind kind,
              Entity_type type,
              T const * const data = nullptr) :
      UniStateVectorBase<DomainType>(name, domain, state, kind, type) {

    int num = domain->num_entities(kind, type);
    if (data == nullptr)
      mydata_ = std::make_shared<std::vector<T>>(num);
    else
      mydata_ = std::make_shared<std::vector<T>>(data, data+num);
  }


  /*!
    @brief Meaningful constructor with uniform initializer
    @param name            Name of vector
    @param state           State manager holding the vector (can be nullptr)
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
    @param initval         Value to which all elements should be initialized to 
  */
  
  UniStateVector(std::string name,
                 std::shared_ptr<DomainType> domain,
                 std::shared_ptr<State> state,
                 Entity_kind kind,
                 Entity_type type,
                 T initval) :
      UniStateVectorBase<DomainType>(name, domain, state, kind, type) {

    int num = domain->num_entities(kind, type);
    mydata_ = std::make_shared<std::vector<T>>(num, initval);
  }


  /*! 
    @brief Copy constructor - DEEP COPY OF DATA

    Copy constructor creates a new vector and copies the meta data of
    the UniStateVector over. Additionally, it copies all of the vector
    data from the source vector to the new vector.  Modification of one
    vector's data has no effect on the other.

    Since mystate_ is a weak_ptr, we have to lock it to get a shared_ptr
    to send to the UniStateVectorBase constructor
  */

  UniStateVector(UniStateVector const & in_vector) :
      UniStateVectorBase<DomainType>(in_vector.myname_,
                                     in_vector.mydomain_,
                                     in_vector.mystate_.lock(),
                                     in_vector.entity_kind_,
                                     in_vector.entity_type_) {

    mydata_ = std::make_shared<std::vector<T>>((in_vector.mydata_)->begin(),
                                               (in_vector.mydata_)->end());
  }

  /*!
    @brief Assignment operator
  
    Assignment operator does a shallow copy of the metadata and a
    shared_ptr to the data. So modification of one state vector's
    data will result in modification of the other's data as well
  */

  UniStateVector & operator=(UniStateVector const & in_vector) {
    StateVectorBase::myname_ = in_vector.myname_;
    StateVectorBase::mystate_ = in_vector.mystate_;
    StateVectorBase::entity_kind_ = in_vector.entity_kind_;
    StateVectorBase::entity_type_ = in_vector.entity_type_;
    UniStateVectorBase<DomainType>::mydomain_ = in_vector.mydomain_;

    mydata_ = in_vector.mydata_;  // shared_ptr counter will increment

    return *this;
  }

  /// Destructor
  
  ~UniStateVector() {}

  /// Get the raw data

  T *get_raw_data() { return &((*mydata_)[0]); }

  /// Get the raw data

  T const *get_raw_data() const { return &((*mydata_)[0]); }

  /// Get a shared pointer to the data

  std::shared_ptr<T> get_data() { return mydata_; }

  /// Type of data

  const std::type_info& data_type() {
    const std::type_info& ti = typeid(T);
    return ti;
  }

  //! Subset of std::vector functionality. We can add others as needed

  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

  iterator begin() { return mydata_->begin(); }
  iterator end() { return mydata_->end(); }
  const_iterator cbegin() const { return mydata_->cbegin(); }
  const_iterator cend() const { return mydata_->cend(); }

  typedef T& reference;
  typedef T const& const_reference;
  reference operator[](int i) { return (*mydata_)[i]; }
  const_reference operator[](int i) const { return (*mydata_)[i]; }

  size_t size() const { return mydata_->size(); }
  void resize(size_t newsize) { mydata_->resize(newsize); }
  void resize(size_t newsize, T val) { mydata_->resize(newsize, val); }

  void clear() {mydata_->clear();}

  //! Output the data

  std::ostream& print(std::ostream& os) const {
    os << "\n";
    os << "Vector \"" << StateVectorBase::myname_ << "\" on entity kind " <<
        StateVectorBase::entity_kind_ << " :\n";
    os << size() << " elements\n";

    for (const_iterator it = cbegin(); it != cend(); it++)
      os << (*it) << "\n";
    os << std::endl;  // flush the output

    return os;
  }

 private:
  std::shared_ptr<std::vector<T>> mydata_;
};  // UniStateVector

//! Send UniStateVector to output stream

template <class T, class DomainType>
std::ostream & operator<<(std::ostream & os,
                          UniStateVector<T, DomainType> const & sv) {
  return sv.print(os);
}

//! Send a std::array to output stream

template <class T, std::size_t N>
std::ostream & operator<<(std::ostream & os, const std::array<T, N>& arr) {
  std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(os, " "));
  return os;
}



///////////////////////////////////////////////////////////////////////////////


/*!
  @class MultiStateVectorBase jali_state_vector.h

  @brief MultiStateVectorBase Base class storing meta-data for a
  vector that is designed to store multi-material or multi-phase state
  data on entities in a mesh or a mesh tile

  @tparam DomainType  Mesh or Mesh Tile 
*/

template <class DomainType = Mesh>
class MultiStateVectorBase : public StateVectorBase {
  // Clearly we cannot have a MultiStateVector comprised of StateVectors
  // defined on MeshSets to itself be defined on MeshSets!
  // Valid types for DomainType are Mesh and MeshTile

  // Cannot compile
  //  static_assert(DomainType != MeshSet,
  //    "MultiStateVector cannot be defined on MeshSets");


 public:

  //! Default constructor
  MultiStateVectorBase() : StateVectorBase("UninitializedVector", nullptr,
                                            Entity_kind::UNKNOWN_KIND,
                                            Entity_type::TYPE_UNKNOWN) {}

  /*!
    @brief Meaningful constructor
    @param name            name of vector
    @param domain          Domain on which vector is defined
    @param state           State manager holding the vector (CANNOT BE nullptr)
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
  */
  
  MultiStateVectorBase(std::string name,
                       std::shared_ptr<DomainType> domain,
                       std::shared_ptr<State> state,
                       Entity_kind kind,
                       Entity_type type) :
      StateVectorBase(name, state, kind, type), mydomain_(domain) {}



  /*! 
    @brief Copy constructor

    Copy constructor creates a new vector and copies the meta data of
    the StateVector over.

    Since mystate_ is a weak_ptr, we have to lock it to get a shared_ptr
    to send to the UniStateVectorBase constructor
  */

  MultiStateVectorBase(MultiStateVectorBase const & in_vector) :
      StateVectorBase(in_vector.myname_, in_vector.mystate_.lock(),
                      in_vector.entity_kind_, in_vector.entity_type_),
      mydomain_(in_vector.mydomain_) {}

  /*!
    @brief Assignment operator
  
    Assignment operator does a shallow copy of the metadata
  */

  MultiStateVectorBase & operator=(MultiStateVectorBase const & in_vector) {
    StateVectorBase::myname_ = in_vector.myname_;
    StateVectorBase::mystate_ = in_vector.mystate_;
    StateVectorBase::entity_kind_ = in_vector.entity_kind_;
    StateVectorBase::entity_type_ = in_vector.entity_type_;
    mydomain_ = in_vector.mydomain_;

    return *this;
  }


  /// Destructor
  
  ~MultiStateVectorBase() {}

  /// Domain on which MultiStateVector is defined on (Mesh, MeshTile)
  
  std::shared_ptr<DomainType> domain() const { return mydomain_; }

  /// Underlying mesh regardless of what type of domain StateVector is
  /// defined on. We have to return a reference to the mesh rather
  /// than a shared pointer because if the DomainType is a MeshTile,
  /// it only has a mesh reference not a pointer to the mesh

  Mesh & mesh() const { return get_mesh_of_domain(mydomain_); }

  /// Get the type of state vector

  StateVector_type type() { return StateVector_type::MULTIVAL; }

  /// Size (Number of materials) of a multi-material vector
  virtual size_t size() const = 0;

  /// Size of a particular material array
  virtual size_t size(int m) const = 0;

  /// Resize a particular material array
  virtual void resize(int m, size_t newsize) = 0;

  /// Clear out all the data
  virtual void clear() = 0;

  /// Clear out data for a material
  virtual void clear(int m) = 0;

  /// Add a material and its entries to the vector
  virtual void add_material(int ncells) = 0;

  // Remove a material and its entries from the vector
  virtual void rem_material(int m) = 0;

  //! Output the data (but only if it is arithmetic type)
  // DISABLED UNTIL WE CAN ENABLE IT ONLY FOR THOSE TYPES THAT CAN BE STREAMED

  // typename std::enable_if<std::is_arithmetic<T>::value, std::ostream&>::type
  // print(std::ostream& os) const {
  //   os << "\n";
  //   os << "MultiMaterial vector \"" << myname_ << "\" on entity kind " << 
  //       entity_kind_ << " :\n";

  //   int nmats = mydata_->size();
  //   for (int m = 0; m < nmats; m++) {
  //     os << "Material " << m << ": " << size(m) << " elements\n";
  //     for (const_iterator it = cbegin(m); it != cend(m); it++)
  //       os << (*it) << "\n";
  //     os << std::endl;  // flush the output
  //   }
  //   return os;
  // }

 protected:
  std::shared_ptr<DomainType> mydomain_;

  const Mesh & get_mesh_of_domain(std::shared_ptr<MeshTile> meshtile) const {
     return meshtile->mesh();
  }
  Mesh & get_mesh_of_domain(std::shared_ptr<Mesh> mesh) const {
    return *mesh;
  }
};  // MultiStateVectorBase



///////////////////////////////////////////////////////////////////////////////



/*!
  @class MultiStateVector jali_state_vector.h
  @brief MultiStateVector stores multi-material or multi-phase state data for
  entities in a mesh or a mesh tile

  Templated class for multi-material/multi-phase state vectors, i.e.,
  a single entity like a cell might contain multiple values of the
  particular state. Allows data to be accessed using a (entityID,
  materialID) operator. MultiStateVectors can be associated with a mesh
  or a mesh tile but as far as we can see, it does not make sense to
  associate it with a meshset.
  t
  @tparam DomainType  Mesh or Mesh Tile 
*/

template <class T, class DomainType = Mesh>
class MultiStateVector : public MultiStateVectorBase<DomainType> {
 public:

  //! Default constructor
  MultiStateVector() :
      MultiStateVectorBase<DomainType>("UninitializedVector",
                                       nullptr,
                                       nullptr,
                                       Entity_kind::UNKNOWN_KIND,
                                       Entity_type::TYPE_UNKNOWN) {}

  /*!
    @brief Meaningful constructor with data
    @param name            Name identifier of vector
    @param domain          Domain on which vector is defined
    @param state           State manager holding the vector (CANNOT BE nullptr)
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
    @param data            Pointer to array data to be used to initialize vector (optional)

    If provided, data is a 2D rectangular array which is ncells by
    nmaterials in dimension. However, the state manager will be
    queried to see which materials are in which cells and only those
    values will be read.
  */
  
  MultiStateVector(std::string name,
                   std::shared_ptr<DomainType> domain,
                   std::shared_ptr<State> state,
                   Entity_kind kind,
                   Entity_type type,
                   Data_layout layout = Jali::Data_layout::MATERIAL_CENTRIC,
                   T const * const * data = nullptr) :
      MultiStateVectorBase<DomainType>(name, domain, state, kind, type) {
    assert(state != nullptr);
    allocate();
    if (data) assign(layout, data);
  }


  /*!
    @brief Meaningful constructor with uniform initializer
    @param name            String identifier of vector
    @param domain          Domain on which vector is defined 
    @param state           State manager holding the vector (CANNOT BE nullptr)
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
    @param matflag         Pointer to 2D rectangular boolean array indicating if a material is present in a cell or not (optional)
    @param data            Pointer to 2D rectangular array to be used to initialize vector (optional)
    
    If provided, data is a 2D rectangular array which is ncells by
    nmaterials in dimension. However, the state manager will be
    queried to see which materials are in which cells and only those
    values will be read.
  */
  
  MultiStateVector(std::string name,
                   std::shared_ptr<DomainType> domain,
                   std::shared_ptr<State> state,
                   Entity_kind kind,
                   Entity_type type,
                   T initval) :
      MultiStateVectorBase<DomainType>(name, domain, state, kind, type) {
    assert(state != nullptr);
    allocate();
    assign(initval);
  }
  
  /*! 
    @brief Copy constructor - DEEP COPY OF DATA

    Copy constructor creates a new vector and copies the meta data of
    the StateVector over. Additionally, it copies all of the vector
    data from the source vector to the new vector.  Modification of one
    vector's data has no effect on the other.

    Since mystate_ is a weak_ptr, we have to lock it to get a shared_ptr
    to send to the UniStateVectorBase constructor
  */

  MultiStateVector(MultiStateVector const & in_vector) :
      MultiStateVectorBase<DomainType>(in_vector.myname_,
                                       in_vector.mydomain_,
                                       in_vector.mystate_.lock(),
                                       in_vector.entity_kind_,
                                       in_vector.entity_type_) {
    mydata_ =
    std::make_shared<std::vector<std::vector<T>>>((in_vector.mydata_)->begin(),
                                                  (in_vector.mydata_)->end());
  }

  /*!
    @brief Assignment operator
  
    Assignment operator does a shallow copy of the metadata and a
    shared_ptr to the data. So modification of one state vector's
    data will result in modification of the other's data as well
  */

  MultiStateVector & operator=(MultiStateVector const & in_vector) {
    StateVectorBase::myname_ = in_vector.myname_;
    StateVectorBase::mystate_ = in_vector.mystate_;
    StateVectorBase::entity_kind_ = in_vector.entity_kind_;
    StateVectorBase::entity_type_ = in_vector.entity_type_;
    MultiStateVectorBase<DomainType>::mydomain_ = in_vector.mydomain_;

    mydata_ = in_vector.mydata_;  // shared_ptr counter will increment

    return *this;
  }


  /*!
    @brief Allocate space for storing multi-material data in material-centric form
  */

  void allocate() {
    int nummats = state_get_num_materials(StateVectorBase::mystate_);
    mydata_ = std::make_shared<std::vector<std::vector<T>>>(nummats);
    for (int m = 0; m < nummats; m++) {
      // get entities in the material set 'm'
      std::shared_ptr<MeshSet> mset =
          state_get_material_set(StateVectorBase::mystate_, m);
      int numents = mset->entities().size();
      (*mydata_)[m].resize(numents);
    }
  }


  /*!
    @brief Assign 2D array data to a multi-material state vector

    @param layout Cell centric (first index is cells) or material centric
    @param data  Data in rectangular format (number of cells in domain by number of materials in the problem)

    The presence/absence of a material in a cell is independently
    determined based on material sets in the mesh
  */

  void assign(Data_layout layout, T const * const * const data) {
    int nummats = state_get_num_materials(StateVectorBase::mystate_);
    mydata_->resize(nummats);
    
    for (int m = 0; m < nummats; m++) {
      // get entities in the material set 'm'
      std::shared_ptr<MeshSet> mset =
          state_get_material_set(StateVectorBase::mystate_, m);
      
      std::vector<int> const& entities = mset->entities();
      int numents = entities.size();
      (*mydata_)[m].resize(numents);
      if (data) {
        for (int i = 0; i < numents; i++) {
          int c = entities[i];   // cell of material set
          // rectangular to compact storage
          (*mydata_)[m][i] =
              (layout == Data_layout::CELL_CENTRIC) ? data[c][m] : data[m][c];
        }
      }
    }
  }

  /*!
    @brief Assign a number to each entry in the multi-material state vector

    @param initval  Value to assign to entries

    The first index of the 'data' array is the cell index and and the
    second is the material index. The presence/absence of a material
    in a cell is independently determined based on material sets in
    the mesh
  */

  void assign(T initval) {
    int nummats = state_get_num_materials(StateVectorBase::mystate_);
    mydata_->resize(nummats);
    
    for (int m = 0; m < nummats; m++) {
      // get entities in the material set 'm'
      std::shared_ptr<MeshSet> mset =
          state_get_material_set(StateVectorBase::mystate_, m);
      
      std::vector<int> const& entities = mset->entities();
      int numents = entities.size();
      (*mydata_)[m].resize(numents);
      for (int i = 0; i < numents; i++) {
        int c = entities[i];   // cell of material set
        (*mydata_)[m][i] = initval;   // rectangular to compact storage
      }
    }
  }

  /// Destructor
  
  ~MultiStateVector() {}

  /// Get the raw data (NOT USEFUL)

  T *get_raw_data() { return &((*mydata_)[0]); }

  /// Get the raw data for a material

  T *get_raw_data(int m) { return &((*mydata_)[m][0]); }

  /// Get the raw data for a material

  T const *get_raw_data(int m) const { return &((*mydata_)[m][0]); }

  /// Get a shared ptr to the data

  std::shared_ptr<std::vector<std::vector<T>>> get_data() { return mydata_; }

  /// Get a reference to the data for one material

  std::vector<T>& get_matdata(int m) { return (*mydata_)[m]; }

  /// Get a reference to the data for one material

  std::vector<T> const& get_matdata(int m) const { return (*mydata_)[m]; }

  /// Type of data

  const std::type_info& data_type() {
    const std::type_info& ti = typeid(T);
    return ti;
  }

  //! Subset of std::vector functionality. We can add others as needed

  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

  iterator begin(int m) { return (*mydata_)[m].begin(); }
  iterator end(int m) { return (*mydata_)[m].end(); }
  const_iterator cbegin(int m) const { return (*mydata_)[m].cbegin(); }
  const_iterator cend(int m) const { return (*mydata_)[m].cend(); }

  
  /// @brief Value of field for a material 'm' in a cell 'c'
  /// @param i Cell ID if layout = CELL_CENTRIC, Material index
  ///          if layout = MATERIAL_CENTRIC
  /// @param j Material index (like 0th material, 2nd material etc.,
  ///          not material 51965 or Steel) if layout is CELL_CENTRIC,
  ///          cell index if layout = MATERIAL_CENTRIC
  /// @param layout Way in which we want to retrieve the data (CELL_CENTRIC
  //                means the first index is the cell index and the second is
  //                material index, MATERIAL_CENTRIC is the reverse)
  ///
  ///
  /// If an entry for this material in this cell does not exist, then
  /// the operator will throw an exception
  ///
  /// This could be quite inefficient since we are looking up entries
  /// in the meshset on the one hand and then switching and looking up
  /// entries in the data vector on the other hand. There is also the
  /// cost of the 'if' check and branch misprediction cost. We can
  /// avoid the 'if' check if we offset everything by 1 and access
  /// mydata_[m][cloc+1], with mydata_[m][0] being T(0.0)
  ///
  /// This could be more efficient if 'c' were a local index in the
  /// material - then we wouldn't have to check for the local index in
  /// the material

  T operator()(int i, int j,
               Data_layout layout = Data_layout::MATERIAL_CENTRIC) const {
    int m = (layout == Data_layout::MATERIAL_CENTRIC) ? i : j;
    int c = (layout == Data_layout::MATERIAL_CENTRIC) ? j : i;
    assert(m < state_get_num_materials(StateVectorBase::mystate_));

    std::shared_ptr<MeshSet> mset =
        state_get_material_set(StateVectorBase::mystate_, m);
    int cloc = mset->index_in_set(c);
    if (cloc != -1)
      return (*mydata_)[m][cloc];
    else
      return T(0);
  }
    
  typedef T& reference;

  /// @brief Reference of field for a material 'm' in a cell 'c'.
  /// @param i Cell ID if layout = CELL_CENTRIC, Material index
  ///          if layout = MATERIAL_CENTRIC
  /// @param j Material index (like 0th material, 2nd material etc.,
  ///          not material 51965 or Steel) if layout is CELL_CENTRIC,
  ///          cell index if layout = MATERIAL_CENTRIC
  /// @param layout Way in which we want to retreieve the data (CELL_CENTRIC
  //                means the first index is the cell index and the second is
  //                material index, MATERIAL_CENTRIC is the reverse)
   
  reference operator()(int i, int j,
                       Data_layout layout = Data_layout::MATERIAL_CENTRIC) {
    int m = (layout == Data_layout::MATERIAL_CENTRIC) ? i : j;
    int c = (layout == Data_layout::MATERIAL_CENTRIC) ? j : i;
    assert(m < state_get_num_materials(StateVectorBase::mystate_));

    assert(m < state_get_num_materials(StateVectorBase::mystate_));
    std::shared_ptr<MeshSet> mset =
        state_get_material_set(StateVectorBase::mystate_, m);
    int cloc = mset->index_in_set(c);
    if (cloc == -1)
      throw std::runtime_error("Cell does not contain material. Add it in the statemanager");
    
    return (*mydata_)[m][cloc];
  }

  /// Size (Number of materials) of a multi-material vector
  size_t size() const { return mydata_->size(); }

  /// Size of a particular material array
  size_t size(int m) const { return (*mydata_)[m].size(); }

  /// Resize a particular material array
  void resize(int m, size_t newsize) { (*mydata_)[m].resize(newsize); }

  /// Resize a particular material array and initialize new elements to val
  void resize(int m, size_t newsize, T val) {
    (*mydata_)[m].resize(newsize, val);
  }

  /// Clear out all the data
  void clear() { mydata_->clear(); }

  /// Clear out data for a material
  void clear(int m) { (*mydata_)[m].clear(); }

  /// Add a material and its entries to the vector
  void add_material(int ncells) {
    size_t nmats = mydata_->size();
    mydata_->resize(nmats+1);
    (*mydata_)[nmats].resize(ncells);
  }

  // Remove a material and its entries from the vector
  void rem_material(int m) {
    mydata_->erase(mydata_->begin()+m);
  }

  //! Output the data (but only if it is arithmetic type)
  // DISABLED UNTIL WE CAN ENABLE IT ONLY FOR THOSE TYPES THAT CAN BE STREAMED

  // typename std::enable_if<std::is_arithmetic<T>::value, std::ostream&>::type
  // print(std::ostream& os) const {
  //   os << "\n";
  //   os << "MultiMaterial vector \"" << myname_ << "\" on entity kind " << 
  //       entity_kind_ << " :\n";

  //   int nmats = mydata_->size();
  //   for (int m = 0; m < nmats; m++) {
  //     os << "Material " << m << ": " << size(m) << " elements\n";
  //     for (const_iterator it = cbegin(m); it != cend(m); it++)
  //       os << (*it) << "\n";
  //     os << std::endl;  // flush the output
  //   }
  //   return os;
  // }

 private:
  // Should we flatten this for read access and rework it when the
  // state manager adds new cells to material sets (same as adding
  // materials to cells)? May be needed for accelerators

  std::shared_ptr<std::vector<std::vector<T>>> mydata_;
};  // MultiStateVector


//! Send MultiStateVector to output stream

// template <class T, class DomainType>
// std::ostream & operator<<(std::ostream & os,
//                          MultiStateVector<T, DomainType> const & sv) {
//  return sv.print(os);
//}

}  // namespace Jali


#endif  // JALI_STATE_VECTOR_H_
