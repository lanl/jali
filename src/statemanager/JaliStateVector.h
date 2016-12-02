/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

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

/*!
  @class BaseStateVector jali_state_vector.h
  @brief BaseStateVector provides a base class for state vectors on meshes, mesh tiles or mesh subsets
  @param name   Name of the state vector

  A StateVector class inherits from this class and is templated on
  specific types of data. The DomainType class has to be of type mesh
  or a type that is part of a mesh and therefore, can answer the
  question "What is your mesh" and "How many entities of a particular kind
  (CELL, NODE, etc) do you have"?
*/

class BaseStateVector {
 public:

  //! Constructor
  
  explicit BaseStateVector(std::string const name,
                           Entity_kind const kind,
                           Entity_type const type) :
      myname_(name), entity_kind_(kind), entity_type_(type) {}

  //! Constructor using int/enum as identifier instead of string

  BaseStateVector(int const identifier,
                  Entity_kind const kind,
                  Entity_type const type) :
      myname_(int_to_string(identifier)), entity_kind_(kind),
      entity_type_(type) {}

  //! Destructor

  virtual ~BaseStateVector() {}

  //! Virtual methods

  virtual std::ostream & print(std::ostream & os) const = 0;
  virtual void* get_raw_data() = 0;
  virtual int size() const = 0;
  virtual const std::type_info& get_type() = 0;

  //! Convert enum to string for identifying state vectors. Uses ~
  //! (which should be forbidden in user-defined names) to avoid
  //! potential collisions with string names if both are used.
  static std::string int_to_string(int const identifier) {
    return ("~" + std::to_string(static_cast<long long>(identifier)));
  }

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
};





/*!
  @class StateVector jali_state_vector.h
  @brief StateVector stores state data for entities in a mesh, mesh tile or mesh subset

  Templated class for state vectors with specific types.
  Provides some limited functionality of a std::vector while adding
  some additional meta-data like the mesh associated with this data.
  @tparam DomainType  Mesh, Mesh Tile or Mesh Subset (coming soon)
*/

template <class T, class DomainType = Mesh>
class StateVector : public BaseStateVector {
 public:

  //! Default constructor
  StateVector() : BaseStateVector("UninitializedVector",
                                  Entity_kind::UNKNOWN_KIND,
                                  Entity_type::TYPE_UNKNOWN) {}

  /*!
    @brief Meaningful constructor with data and a string identifier
    @param name            String identifier of vector
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
    @param data            Pointer to array data to be used to initialize vector (optional)
  */
  
  StateVector(std::string const name, std::shared_ptr<DomainType> domain,
              Entity_kind const kind, Entity_type const type,
              T const * const data = nullptr) :
      BaseStateVector(name, kind, type), mydomain_(domain) {

    int num = mydomain_->num_entities(kind, type);
    if (data == nullptr)
      mydata_ = std::make_shared<std::vector<T>>(num);
    else
      mydata_ = std::make_shared<std::vector<T>>(data, data+num);
  }

  /*!
    @brief Meaningful constructor with data and a integer identifier
    @param name            String identifier of vector
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
    @param data            Pointer to array data to be used to initialize vector (optional)
  */
  
  StateVector(int const identifier, std::shared_ptr<DomainType> domain,
              Entity_kind const kind, Entity_type const type,
              T const * const data = nullptr) :
      BaseStateVector(identifier, kind, type), mydomain_(domain) {

    int num = mydomain_->num_entities(kind, type);
    if (data == nullptr)
      mydata_ = std::make_shared<std::vector<T>>(num);
    else
      mydata_ = std::make_shared<std::vector<T>>(data, data+num);
  }

  /*!
    @brief Meaningful constructor with uniform initializer and a string identifier
    @param name            String identifier of vector
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
    @param initval         Value to which all elements should be initialized to 
  */
  
  StateVector(std::string const name, std::shared_ptr<DomainType> domain,
              Entity_kind const kind, Entity_type const type,
              T const& initval) :
      BaseStateVector(name, kind, type), mydomain_(domain) {

    int num = mydomain_->num_entities(kind, type);
    mydata_ = std::make_shared<std::vector<T>>(num, initval);
  }

  /*!
    @brief Meaningful constructor with data and a integer identifier
    @param name            String identifier of vector
    @param kind            What kind of entity in the Domain does data live on
    @param type            What type of entity data lives on (PARALLEL_OWNED, PARALLEL_GHOST, etc)
    @param initval         Value to which all elements should be initialized to 
  */
  
  StateVector(int const identifier, std::shared_ptr<DomainType> domain,
              Entity_kind const kind, Entity_type const type,
              T const& initval) :
      BaseStateVector(identifier, kind, type), mydomain_(domain) {
    
    int num = mydomain_->num_entities(kind, type);
    mydata_ = std::make_shared<std::vector<T>>(num, initval);
  }

  /*! 
    @brief Copy constructor - DEEP COPY OF DATA

    Copy constructor creates a new vector and copies the meta data of
    the StateVector over. Additionally, it copies all of the vector
    data from the source vector to the new vector.  Modification of one
    vector's data has no effect on the other.
  */

  StateVector(StateVector const & in_vector) :
      BaseStateVector(in_vector.myname_, in_vector.entity_kind_,
                      in_vector.entity_type_),
      mydomain_(in_vector.mydomain_) {
    
    mydata_ = std::make_shared<std::vector<T>>((in_vector.mydata_)->begin(),
                                               (in_vector.mydata_)->end());
  }

  /*!
    @brief Assignment operator
  
    Assignment operator does a shallow copy of the metadata and a
    shared_ptr to the data. So modification of one state vector's
    data will result in modification of the other's data as well
  */

  StateVector & operator=(StateVector const & in_vector) {
    BaseStateVector::myname_ = in_vector.myname_;
    BaseStateVector::entity_kind_ = in_vector.entity_kind_;
    BaseStateVector::entity_type_ = in_vector.entity_type_;
    mydomain_ = in_vector.mydomain_;
    mydata_ = in_vector.mydata_;  // shared_ptr counter will increment
  }

  /// Destructor
  
  ~StateVector() {}

  /// Domain on which StateVector is defined on (Mesh, MeshTile, MeshSubset)
  
  std::shared_ptr<DomainType> domain() const { return mydomain_; }

  /// Underlying mesh regardless of what type of domain StateVector is
  /// defined on. We have to return a reference to the mesh rather
  /// than a shared pointer because if the DomainType is a MeshTile,
  /// it only has a mesh reference not a pointer to the mesh

  Mesh & mesh() const { return get_mesh_of_domain(mydomain_); }

  /// Get the raw data

  void* get_raw_data() { return (void*)(&((*mydata_)[0])); }

  /// Get a shared pointer to the data

  std::shared_ptr<T> get_data() { return mydata_; }

  /// Type of data

  const std::type_info& get_type() {
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

  int size() const {return mydata_->size();}
  void resize(size_t n, T val) { mydata_->resize(n, val); }

  void clear() { mydata_->clear(); }

  //! Output the data

  std::ostream& print(std::ostream& os) const {
    os << "\n";
    os << "Vector \"" << myname_ << "\" on entity kind " << entity_kind_ <<
        " :\n";
    os << size() << " elements\n";

    for (const_iterator it = cbegin(); it != cend(); it++)
      os << (*it) << "\n";
    os << std::endl;  // flush the output

    return os;
  }

 protected:
  std::shared_ptr<DomainType> mydomain_;
  std::shared_ptr<std::vector<T>> mydata_;

 private:
  Mesh & get_mesh_of_domain(std::shared_ptr<MeshTile> meshtile) const {
     return meshtile->mesh();
  }
  Mesh & get_mesh_of_domain(std::shared_ptr<Mesh> mesh) const {
    return *mesh;
  }
};


//! Send StateVector to output stream

template <class T, class DomainType>
std::ostream & operator<<(std::ostream & os,
                          StateVector<T, DomainType> const & sv) {
  return sv.print(os);
}


//! Send a std::array to output stream

template <class T, std::size_t N>
std::ostream & operator<<(std::ostream & os, const std::array<T, N>& arr) {
  std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(os, " "));
  return os;
}


}  // namespace Jali


#endif  // JALI_STATE_VECTOR_H_
