/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_STATE_VECTOR_H_
#define JALI_STATE_VECTOR_H_

//!
//!  \class StateVector jali_state_vector.h
//!  \brief StateVector provides a mechanism to store state data for mesh entities
//!

#include <iostream>
#include <iterator>
#include <vector>
#include <typeinfo>
#include <string>
#include <algorithm>

#include "Mesh.hh"    // jali mesh header

namespace Jali {

//! Base class for state vectors.  Children inherit from this class
//! and hold specific types of data.

class BaseStateVector {
 public:


  //! Constructor

  BaseStateVector(std::string const name, Entity_kind const on_what,
                  const std::shared_ptr<Mesh> mesh) :
      myname_(name), on_what_(on_what), mymesh_(mesh) {}

  //! Constructor using int/enum as identifier instead of string

  BaseStateVector(int const identifier, Entity_kind const on_what,
                  const std::shared_ptr<Mesh> mesh) :
      myname_(int_to_string(identifier)), on_what_(on_what), mymesh_(mesh) {}

  //! Destructor

  virtual ~BaseStateVector() {}

  //! Virtual methods

  virtual std::ostream & print(std::ostream & os) const = 0;
  virtual void* get_data() = 0;
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
  Jali::Entity_kind on_what() const { return on_what_; }
  std::shared_ptr<Mesh> mesh() const { return mymesh_; }

 protected:

  Jali::Entity_kind on_what_;
  std::string myname_;
  std::shared_ptr<Mesh> mymesh_;
};


//! Templated class for state vectors with specific types.
//! Provides some limited functionality of a std::vector while adding
//! some additional meta-data like the mesh associated with this data.

template <class T>
class StateVector : public BaseStateVector {
 public:

  //! Default constructor

  StateVector() : BaseStateVector("NoName", UNKNOWN_KIND,
                                  std::shared_ptr<Mesh>()) {}

  //! Meaningful constructor with data

  StateVector(std::string const name, Entity_kind const on_what,
              const std::shared_ptr<Mesh> mesh, T* data) :
      BaseStateVector(name, on_what, mesh) {

    int num = mesh->num_entities(on_what, ALL);
    mydata_ = std::shared_ptr<std::vector<T>>(new std::vector<T>);
    mydata_->resize(num);
    std::copy(data, data+num, mydata_->begin());

  }

  //! Constructor with int/enum as identifier instead of string
  
  StateVector(int const identifier, Entity_kind const on_what,
              const std::shared_ptr<Mesh> mesh, T* data) :
      BaseStateVector(identifier, on_what, mesh) {

    int num = mesh->num_entities(on_what, ALL);
    mydata_ = std::shared_ptr<std::vector<T>>(new std::vector<T>);
    mydata_->resize(num);
    std::copy(data, data+num, mydata_->begin());

  }

  //! Copy constructor - DEEP COPY OF DATA

  //! Copy constructor creates a new vector and copies the meta data
  //! of the StateVector over. Additionally, it copies all of the
  //! vector data from the source vector to the new vector.
  //! Modification of one vector's data has no effect on the other.

  StateVector(StateVector const & in_vector) :
      BaseStateVector(in_vector.myname_, in_vector.on_what_, in_vector.mymesh_),
      mydata_(new std::vector<T>(in_vector.size())) {

    // deep copy of the data
    std::copy((in_vector.mydata_)->begin(), (in_vector.mydata_)->end(),
              mydata_->begin());
  }


  //! \brief Assignment operator

  //! Assignment operator does a shallow copy of the metadata and a
  //! shared_ptr to the data. So modification of one state vector's
  //! data will result in modification of the other's data as well

  StateVector & operator=(StateVector const & in_vector) {

    myname_ = in_vector.myname_;
    on_what_ = in_vector.on_what_;
    mymesh_ = in_vector.mymesh_;
    mydata_ = in_vector.mydata_;
  }

  //! Destructor

  ~StateVector() {}

  //! Get the raw data

  void* get_data() { return (void*)(&((*mydata_)[0])); }


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
  void resize(size_t n, T val) { mydata_->resize(n, val);}

  void clear() { mydata_->clear(); }

  //! Output the data

  std::ostream & print(std::ostream & os) const {

    os << "\n";
    os << "Vector \"" << myname_ << "\" on entity kind " << on_what_ << " :\n";
    os << size() << " elements\n";

    for (const_iterator it = cbegin(); it != cend(); it++)
      os << (*it) << "\n";
    os << std::endl;  // flush the output

    return os;
  }

 protected:

  std::shared_ptr<std::vector<T>> mydata_;
};


//! Send StateVector to output stream

template <class T>
std::ostream & operator<<(std::ostream & os, Jali::StateVector<T> const & sv) {
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
