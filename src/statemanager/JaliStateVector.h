/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_STATE_VECTOR_H_
#define JALI_STATE_VECTOR_H_

/*!
  @class StateVector jali_state_vector.h
  @brief StateVector stores state data for entities in a mesh, mesh tile or mesh subset
*/

#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <algorithm>
#include <typeinfo>
#include <typeindex>

#include "Mesh.hh"    // jali mesh header

namespace Jali {

// Helper functions to get the mesh of a domain (Mesh, MeshTile or MeshSubset)

/*!
  @class BaseStateVector jali_state_vector.h
  @brief BaseStateVector provides a base class for state vectors on meshes, mesh tiles or mesh subsets
  @param name   Name of the state vector

  A StateVector class inherits from this class and is templated on
  specific types of data. The DomainType class has to be of type mesh
  or a type that is part of a mesh and therefore, can answer the
  question "What is your mesh" and "How many entities of type
  'on_what' do you have"?
*/

class BaseStateVector {
 public:

  //! Constructor
  
  explicit BaseStateVector(std::string const name,
                           Entity_kind const on_what,
                           Parallel_type const parallel_type) :
      myname_(name), on_what_(on_what), parallel_type_(parallel_type) {}

  //! Destructor

  virtual ~BaseStateVector() {}

  //! Virtual methods

  virtual std::ostream & print(std::ostream & os) const = 0;
  virtual void* get_raw_data() = 0;
  virtual int size() const = 0;
  virtual const std::type_info& get_type() = 0;

  /// Name of BaseStateVector

  std::string name() const { return myname_; }

  /// What type of entity does it live on (CELL, WEDGE, NODE)?

  Entity_kind on_what() const { return on_what_; }

  /// What type of parallel entity does it live on (OWNED, GHOST or ALL)?

  Parallel_type parallel_type() const { return parallel_type_; }

 protected:
  std::string myname_;
  Entity_kind on_what_;
  Parallel_type parallel_type_;
};


//! Templated class for state vectors with specific types.
//! Provides some limited functionality of a std::vector while adding
//! some additional meta-data like the mesh associated with this data.
//!  @param on_what  Data corresponds to what type of entity in the Domain
//!  @tparam DomainType  Pointer to Mesh, Mesh Tile or Mesh Subset (coming soon)

template <class T, class DomainType = Mesh>
class StateVector : public BaseStateVector {
 public:

  //! No Default constructor
  StateVector() : BaseStateVector("UninitializedVector",
                                  Entity_kind::UNKNOWN_KIND,
                                  Parallel_type::PTYPE_UNKNOWN) {}

  //! Meaningful constructor with data
  
  StateVector(std::string const name, std::shared_ptr<DomainType> domain,
              Entity_kind const on_what, Parallel_type const parallel_type,
              T const * const data) : 
      BaseStateVector(name, on_what, parallel_type), mydomain_(domain) {

    int num = mydomain_->num_entities(on_what_, parallel_type_);
    mydata_ = std::make_shared<std::vector<T>>(data, data+num);

    //    mymesh_ = get_mesh_of_domain(domain);
  }

  /*! 
    @brief Copy constructor - DEEP COPY OF DATA

    Copy constructor creates a new vector and copies the meta data of
    the StateVector over. Additionally, it copies all of the vector
    data from the source vector to the new vector.  Modification of one
    vector's data has no effect on the other.
  */

  StateVector(StateVector const & in_vector) :
      BaseStateVector(in_vector.myname_, in_vector.on_what_,
                      in_vector.parallel_type_),
      mydomain_(in_vector.mydomain_) {
    
    mydata_ = std::make_shared<std::vector<T>>((in_vector.mydata_)->begin(),
                                               (in_vector.mydata_)->end());
    //    mymesh_ = get_mesh_of_domain(in_vector.domain_);
  }

  /*!
    @brief Assignment operator
  
    Assignment operator does a shallow copy of the metadata and a
    shared_ptr to the data. So modification of one state vector's
    data will result in modification of the other's data as well
  */

  StateVector & operator=(StateVector const & in_vector) {
    BaseStateVector::myname_ = in_vector.myname_;
    BaseStateVector::on_what_ = in_vector.on_what_;
    BaseStateVector::parallel_type_ = in_vector.parallel_type_;
    mydomain_ = in_vector.mydomain_;
    mydata_ = in_vector.mydata_;  // shared_ptr counter will increment
  }

  /// Destructor
  
  ~StateVector() {}

  /// Domain on which StateVector is defined on (Mesh, MeshTile, MeshSubset)
  
  std::shared_ptr<DomainType> domain() const { return mydomain_; }

  /// Underlying mesh regardless of what type of domain StateVector is
  /// defined on

  std::shared_ptr<Mesh> mesh() const { return mymesh_; }

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
    os << "Vector \"" << myname_ << "\" on entity kind " << on_what_ <<
        " :\n";
    os << size() << " elements\n";

    for (const_iterator it = cbegin(); it != cend(); it++)
      os << (*it) << "\n";
    os << std::endl;  // flush the output

    return os;
  }

 protected:
  std::shared_ptr<DomainType> mydomain_;
  std::shared_ptr<Mesh> mymesh_;
  std::shared_ptr<std::vector<T>> mydata_;

 private:
  // template <DomainType>
  // std::shared_ptr<Mesh> get_mesh_of_domain(std::shared_ptr<DomainType> domain) {
  //   return domain->mesh();
  // }
};

// template<>
// std::shared_ptr<Mesh>
// StateVector::get_mesh_of_domain<Mesh>(std::shared_ptr<Mesh> mesh) {
//   return mesh;
// }


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
