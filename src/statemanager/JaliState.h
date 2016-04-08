/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_STATE_H_
#define JALI_STATE_H_

/*!
  @class State jali_state.h
  @brief State is a class that stores all of the state data associated
  with a mesh
*/

#include <iostream>
#include <vector>
#include <string>
#include <boost/iterator/permutation_iterator.hpp>

#include "Mesh.hh"    // Jali mesh header

#include "JaliStateVector.h"  // Jali-based state vector


namespace Jali {

class State {
 public:

  /// Constructor

  explicit State(const std::shared_ptr<Jali::Mesh> mesh) : mymesh_(mesh) {}

  /// Copy constructor (disabled)

  State(const State &) = delete;

  /// Assignment operator (disabled)

  State & operator=(const State &) = delete;

  /// Destructor

  ~State() {}

  /// Mesh

  std::shared_ptr<Jali::Mesh> mesh() {return mymesh_;}

  //! Typedefs for iterators for going through all the state vectors

  typedef
  std::vector<std::shared_ptr<BaseStateVector>>::iterator
  iterator;
  
  typedef
  std::vector<std::shared_ptr<BaseStateVector>>::const_iterator
  const_iterator;

  /// Iterators for going through all the state vectors

  iterator begin() {return state_vectors_.begin();}
  iterator end() {return state_vectors_.end();}
  const_iterator cbegin() const {return state_vectors_.begin();}
  const_iterator cend() const {return state_vectors_.end();}


  /// Typedefs for iterators for going through all the state vector name

  typedef std::vector<std::string>::iterator string_iterator;

  /// Iterators for vector names

  string_iterator names_begin() {return names_.begin();}
  string_iterator names_end()   {return names_.end();}


  /// Typedef for permutation iterators to allow iteration through only
  /// the state vectors on a specified entity

  typedef boost::permutation_iterator<
    std::vector<std::shared_ptr<BaseStateVector>>::iterator,
    std::vector<int>::iterator
    >
  permutation_type;


  /// Permutation iterators for iterating over state vectors on a
  /// specific entity type

  permutation_type entity_begin(Jali::Entity_kind entitykind) {
    const int ikind = static_cast<int>(entitykind);
    return boost::make_permutation_iterator(state_vectors_.begin(),
                                            entity_indexes_[ikind].begin());
  }
  permutation_type entity_end(Jali::Entity_kind entitykind) {
    const int ikind = static_cast<int>(entitykind);
    return boost::make_permutation_iterator(state_vectors_.begin(),
                                            entity_indexes_[ikind].end());
  }


  /// Typedef for permutation iterators to allow iteration through only
  /// the state vector _names_ on a specified entity

  typedef boost::permutation_iterator< std::vector<std::string>::iterator,
                                       std::vector<int>::iterator >
  string_permutation;


  /// Iterators for vector names of specific entity types

  string_permutation names_entity_begin(Jali::Entity_kind entitykind) {
    const int ikind = static_cast<int>(entitykind);
    return boost::make_permutation_iterator(names_.begin(),
                                            entity_indexes_[ikind].begin());
  }
  string_permutation names_entity_end(Jali::Entity_kind entitykind) {
    const int ikind = static_cast<int>(entitykind);
    return boost::make_permutation_iterator(names_.begin(),
                                            entity_indexes_[ikind].end());
  }

  /// References to state vectors

  typedef std::shared_ptr<BaseStateVector> pointer;
  typedef const std::shared_ptr<BaseStateVector> const_pointer;

  /// Return pointer to i'th state vector
  pointer operator[](int i) { return state_vectors_[i]; }

  /// Return const pointer to the i'th state vector
  const_pointer operator[](int i) const { return state_vectors_[i]; }

  /// Number of state vectors
  int size() const {return state_vectors_.size();}


  /*! 
    @brief Find iterator to state vector by name
    @tparam T  Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name  String identifier for vector
    @param domain Shared pointer to the domain
    @param on_what What kind of entity data is defined on (CELL, NODE, etc.)
    @param parallel_type What type of parallel entity data is defined on (OWNED, GHOST, etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter on_what) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()
  */

  template<class T, class DomainType>
  iterator find(std::string const name,
                std::shared_ptr<DomainType> const domain,
                Entity_kind const on_what = Entity_kind::ANY_KIND,
                Parallel_type const parallel_type = Parallel_type::ALL) {

    iterator it = state_vectors_.begin();
    while (it != state_vectors_.end()) {
      std::shared_ptr<BaseStateVector> bv = *it;
      std::shared_ptr<StateVector<T, DomainType>> sv =
          std::dynamic_pointer_cast<StateVector<T, DomainType>>(bv);
      if (sv && (sv->name() == name) &&
          (sv->domain() == domain) &&
          ((on_what == Entity_kind::ANY_KIND) || (sv->on_what() == on_what)) &&
          ((parallel_type == Parallel_type::ALL) || (sv->parallel_type() == parallel_type)))
        break;
      else
        ++it;
    }
    return it;
  }


  /*! 
    @brief Find iterator to state vector by integer identifier
    @tparam T  Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name  String identifier for vector
    @param domain Shared pointer to the domain
    @param on_what What kind of entity data is defined on (CELL, NODE, etc.)
    @param parallel_type What type of parallel entity data is defined on (OWNED, GHOST, etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter on_what) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()
  */

  template <class T, class DomainType>
  iterator find(int const identifier,
                std::shared_ptr<DomainType> const domain,
                Jali::Entity_kind const on_what = Entity_kind::ANY_KIND,
                Parallel_type const parallel_type = Parallel_type::ALL) {
  
    return find<T>(BaseStateVector::int_to_string(identifier), domain, on_what,
                   parallel_type);  // DomainType deduced automatically
  }



  /*! 
    @brief Find iterator to state vector by integer identifier (const version)
    @tparam T  Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name  String identifier for vector
    @param domain Shared pointer to the domain
    @param on_what What kind of entity data is defined on (CELL, NODE, etc.)
    @param parallel_type What type of parallel entity data is defined on (OWNED, GHOST, etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter on_what) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()
  */

  template <class T, class DomainType>
  const_iterator find(std::string const name,
                      std::shared_ptr<DomainType> const domain,
                      Entity_kind const on_what = Entity_kind::ANY_KIND,
                      Parallel_type const parallel_type = Parallel_type::ALL)
      const {
    
    const_iterator it = state_vectors_.cbegin();
    while (it != state_vectors_.cend()) {
      std::shared_ptr<BaseStateVector> bv = *it;
      StateVector<T, DomainType> const& sv =
          dynamic_cast<StateVector<T, DomainType> *>(bv.get());
      if ((sv->name() == name) &&
          (sv->domain() == domain) &&
          ((on_what == Entity_kind::ANY_KIND) || (sv->on_what() == on_what)) &&
          ((parallel_type == Parallel_type::ALL) || (sv->parallel_type() == parallel_type)))
        break;
      else
        ++it;
    }
    return it;
  }
  
  /*! 
    @brief Find iterator to state vector by integer identifier (const version)
    @tparam T  Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name  String identifier for vector
    @param domain Shared pointer to the domain
    @param on_what What kind of entity data is defined on (CELL, NODE, etc.)
    @param parallel_type What type of parallel entity data is defined on (OWNED, GHOST, etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter on_what) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()
  */

  template <class T, class DomainType>
  const_iterator find(int const identifier,
                      std::shared_ptr<DomainType> const domain,
                      Entity_kind const on_what = Entity_kind::ANY_KIND,
                      Parallel_type const parallel_type = Parallel_type::ALL)
      const {
    return find<T>(BaseStateVector::int_to_string(identifier), domain, on_what,
                   parallel_type);  // DomainType deduced automatically
  }

  /*! 
    @brief Retrieve a state vector by name given the domain and entity type it
    is defined on
    @tparam T  Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name  String identifier for vector
    @param domain Shared pointer to the domain
    @param on_what What kind of entity data is defined on (CELL, NODE, etc.)
    @param parallel_type What type of parallel entity data is defined on (OWNED, GHOST, etc.)
    @param vector Raw pointer to the state vector

    @todo MUST GET RID OF THIS AND FORCE USAGE OF SHARED POINTER VERSION
  
    Get a statevector with the given name, particular entity kind
    and parallel type of entity it is defined on. The function returns
    true if such a vector was found; false otherwise. The caller must
    know the type of elements in the state vector, int, double,
    std::array<double, 3> or whatever else. The calling routine must
    declare the state vector as "T vector" where T is StateVector<int>
    or StateVector<double> or StateVector<some_other_type>. Even
    though this is a copy into *vector, its an inexpensive shallow
    copy of the meta data only
   */
  
  template <class T, class DomainType>
  bool get(std::string const name,
           std::shared_ptr<DomainType> const domain,
           Entity_kind const on_what,
           Parallel_type const parallel_type,
           StateVector<T, DomainType> *vector) {
    
    iterator it = find<T>(name, domain, on_what, parallel_type);
    if (it != state_vectors_.end()) {
      *vector = *(std::static_pointer_cast<StateVector<T, DomainType>>(*it));
      return true;
    } else {
      return false;
    }
  }

  /*! 
    @brief Retrieve a state vector by integer identifier given the domain and
    entity type it is defined on
    @tparam T  Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param identifier  Integer identifier for vector
    @param domain Shared pointer to the domain
    @param on_what What kind of entity data is defined on (CELL, NODE, etc.)
    @param parallel_type What type of parallel entity data is defined on (OWNED, GHOST, etc.)
    @param vector Raw pointer to the state vector

    @todo MUST GET RID OF THIS AND FORCE USAGE OF SHARED POINTER VERSION
  
    Get a statevector with the given integer identifier and particular
    entity kind and parallel type of entity it is defined on. The
    function returns true if such a vector was found; false
    otherwise. The caller must know the type of elements in the state
    vector, int, double, std::array<double, 3> or whatever else. The
    calling routine must declare the state vector as "T vector" where
    T is StateVector<int> or StateVector<double> or
    StateVector<some_other_type>. Even though this is a copy into
    *vector, its an inexpensive shallow copy of the meta data only
   */
  
  template <class T, class DomainType>
  bool get(int const identifier,
           std::shared_ptr<DomainType> const domain,
           Entity_kind const on_what,
           Parallel_type const parallel_type,
           StateVector<T, DomainType> *vector) {
    
    iterator it = find<T>(BaseStateVector::int_to_string(identifier), domain,
                          on_what, parallel_type);
    if (it != state_vectors_.end()) {
      *vector = *(std::static_pointer_cast<StateVector<T, DomainType>>(*it));
      return true;
    } else {
      return false;
    }
  }


  /*! 
    @brief Retrieve a state vector on the mesh by name (regardless of what
    type of entity it is on)
    @tparam T           Data type
    @tparam DomainType  Type of domain data is defined on (Mesh, MeshTile)
    @param name         String identifier for vector
    @param vector       Raw pointer to the state vector

    @todo MUST GET RID OF THIS AND FORCE USAGE OF SHARED POINTER VERSION
  
    Get a statevector on the mesh with the given name without regard
    to what kind of entity it is defined on. The function returns true
    if such a vector was found; false otherwise. The caller must know
    the type of elements in the state vector, int, double,
    std::array<double, 3> or whatever else. The calling routine must
    declare the state vector as "T vector" where T is StateVector<int>
    or StateVector<double> or StateVector<some_other_type>. Even
    though this is a copy into *vector, its an inexpensive shallow
    copy of the meta data only
   */

  template <class T, class DomainType>
  bool get(std::string const name, StateVector<T, DomainType> *vector_ptr) {
    return get(name, mymesh_, Entity_kind::ANY_KIND, Parallel_type::ALL,
               vector_ptr);
  }
  

  /*! 
    @brief Retrieve a shared pointer to a state vector by name given the domain
    and type of entity it is defined on
    @tparam T            Data type
    @tparam DomainType   Type of domain data is defined on (Mesh, MeshTile)
    @param name          String identifier for vector
    @param domain        Shared pointer to the domain
    @param on_what       What kind of entity data is defined on (CELL, NODE, etc.)
    @param parallel_type What type of parallel entity data is defined on (OWNED, GHOST, etc.)
    @param vector_ptr    Shared pointer to the state vector

    Get a statevector on given domain with the given name, particular
    entity kind and parallel type of entity it is defined on The
    function returns true if such a vector was found; false
    otherwise. The caller must know the type of elements in the state
    vector, int, double, std::array<double, 3> or whatever else. The
    calling routine must declare the state vector as "T vector" where
    T is StateVector<int> or StateVector<double> or
    StateVector<some_other_type>. Even though this is a copy into
    *vector, its an inexpensive shallow copy of the meta data only
   */

  template <class T, class DomainType>
  bool get(std::string const name,
           std::shared_ptr<DomainType> const domain,
           Entity_kind const on_what,
           Parallel_type const parallel_type,
           std::shared_ptr<StateVector<T, DomainType>> *vector_ptr) {
    
    iterator it = find<T>(name, domain, on_what, parallel_type);
    if (it != state_vectors_.end()) {
      *vector_ptr =
          std::static_pointer_cast<StateVector<T, DomainType>>(*it);
      return true;
    } else {
      return false;
    }
  }

  /*! 
    @brief Retrieve a shared pointer to a state vector by integer identifier
    given the domain and type of entity it is defined on
    @tparam T            Data type
    @tparam DomainType   Type of domain data is defined on (Mesh, MeshTile)
    @param identifier    Integer identifier for vector
    @param domain        Shared pointer to the domain
    @param on_what       What kind of entity data is defined on (CELL, NODE, etc.)
    @param parallel_type What type of parallel entity data is defined on (OWNED, GHOST, etc.)
    @param vector_ptr    Shared pointer to the state vector

    Get a statevector on given domain with the given integer
    identifier, particular entity kind and parallel type of entity it
    is defined on. The function returns true if such a vector was
    found; false otherwise. The caller must know the type of elements
    in the state vector, int, double, std::array<double, 3> or
    whatever else. The calling routine must declare the state vector
    as "T vector" where T is StateVector<int> or StateVector<double>
    or StateVector<some_other_type>. Even though this is a copy into
    *vector, its an inexpensive shallow copy of the meta data only
   */

  template <class T, class DomainType>
  bool get(int const identifier,
           std::shared_ptr<DomainType> const domain,
           Entity_kind const on_what,
           Parallel_type const parallel_type,
           std::shared_ptr<StateVector<T, DomainType>> *vector_ptr) {
 
    return get(BaseStateVector::int_to_string(identifier), domain,
               on_what, parallel_type, vector_ptr);
  }


  /*! 
    @brief Retrieve a shared pointer to a state vector by name given the 
    domain it is defined on (regardless of entity type it is on)
    @tparam T           Data type
    @tparam DomainType  Type of domain data is defined on (Mesh, MeshTile)
    @param name         String identifier for vector
    @param domain       Shared pointer to the domain
    @param vector_ptr   Shared pointer to the state vector

    Get a statevector with the given integer identifier without regard
    to what kind of entity it is defined on. The function returns true
    if such a vector was found; false otherwise. The caller must know
    the type of elements in the state vector, int, double,
    std::array<double, 3> or whatever else. The calling routine must
    declare the state vector as "T vector" where T is StateVector<int>
    or StateVector<double> or StateVector<some_other_type>. Even
    though this is a copy into *vector, its an inexpensive shallow
    copy of the meta data only
   */

  template <class T, class DomainType>
  bool get(std::string const name,
           std::shared_ptr<DomainType> const domain,
           std::shared_ptr<StateVector<T, DomainType>> *vector_ptr) {

    return get(name, domain, Entity_kind::ANY_KIND, Parallel_type::ALL,
               vector_ptr);
  }


  /*! 
    @brief Retrieve a shared pointer to a state vector by integer identifier
    given the domain it is defined on (regardless of entity type it is on)
    @tparam T           Data type
    @tparam DomainType  Type of domain data is defined on (Mesh, MeshTile)
    @param identifier   Integer identifier for vector
    @param domain       Shared pointer to the domain
    @param vector_ptr   Shared pointer to the state vector

    Get a statevector with the given integer identifier without regard
    to what kind of entity it is defined on. The function returns true
    if such a vector was found; false otherwise. The caller must know
    the type of elements in the state vector, int, double,
    std::array<double, 3> or whatever else. The calling routine must
    declare the state vector as "T vector" where T is StateVector<int>
    or StateVector<double> or StateVector<some_other_type>. Even
    though this is a copy into *vector, its an inexpensive shallow
    copy of the meta data only
   */

  template <class T, class DomainType>
  bool get(int const identifier,
           std::shared_ptr<DomainType> const domain,
           std::shared_ptr<StateVector<T, DomainType>> *vector_ptr) {

    return get(BaseStateVector::int_to_string(identifier), domain, vector_ptr);
  }

  /*! 
    @brief Add state vector using a string identifier
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param on_what     What kind of entity data is defined on (CELL, NODE, etc.)
    @param ptype       What type of parallel entity data is defined on (OWNED, GHOST, etc.)
    @param data        Raw pointer to data array

    Add state vector - returns reference to the added StateVector.
    The data is copied from the input memory location into an
    internal buffer.
  */

  template <class T, class DomainType>
  StateVector<T, DomainType>& add(std::string const name,
                                  std::shared_ptr<DomainType> domain,
                                  Entity_kind const on_what,
                                  Parallel_type const ptype,
                                  T* data) {
    
    iterator it = find<T>(name, domain, on_what, ptype);
    if (it == end()) {
      // a search of the state vectors by name and kind of entity turned up
      // empty, so add the vector to the list; if not, warn about duplicate
      // state data
      
      auto vector =
          std::make_shared<StateVector<T, DomainType>>(name, domain,
                                                       on_what, ptype,
                                                       data);
      state_vectors_.emplace_back(vector);

      // add the index of this vector in state_vectors_ to the vector of
      // indexes for this entity type, to allow iteration over state
      // vectors on this entity type with a permutation iterator

      int ikind = static_cast<int>(on_what);
      entity_indexes_[ikind].emplace_back(state_vectors_.size()-1);
      names_.emplace_back(name);

      // emplace back may cause reallocation of the vector so the iterator
      // may not be valid. Use [] operator to get reference to vector

      return (*vector);
    } else {
      // found a state vector by same name
      std::cerr <<
          "Attempted to add duplicate state vector. Ignoring\n" << std::endl;
      return
          (*(std::static_pointer_cast<StateVector<T, DomainType>>(*it)));
    }
  }



  /*! 
    @brief Add state vector using an integer identifier
    @tparam T            Data type
    @tparam DomainType   Type of domain data is defined on (Mesh, MeshTile)
    @param identifier    Integer identifier for vector
    @param domain        Shared pointer to the domain
    @param on_what       What kind of entity data is defined on (CELL, NODE, etc.)
    @param parallel_type What type of parallel entity data is defined on (OWNED, GHOST, etc.)
    @param data          Raw pointer to data array

    Add state vector - returns reference to the added StateVector.
    The data is copied from the input memory location into an
    internal buffer.
  */

  template <class T, class DomainType>
  StateVector<T, DomainType>& add(int const identifier,
                                  std::shared_ptr<DomainType> domain,
                                  Entity_kind const on_what,
                                  Parallel_type const parallel_type,
                                  T* data) {

    return add(BaseStateVector::int_to_string(identifier), domain, on_what,
               parallel_type, data);
  }


  /*!
    @brief Add state vectory by copying an input vector

    Add a new state vector to the state manager based on data from
    the input state vector. Meta data is copied from one vector to
    another and A DEEP COPY IS MADE of the input vector data. The
    need for a deep copy is why we cannot just assign in_vector to
    vector copy even if they are on the same domain
  */

  template <class T, class DomainType>
  StateVector<T, DomainType>&
  add(StateVector<T, DomainType> const& in_vec) {

    Entity_kind on_what = in_vec.on_what();
    Parallel_type ptype = in_vec.parallel_type();

    iterator it = find<T>(in_vec.name(), in_vec.domain(), on_what, ptype);
    if (it == end()) {

      std::shared_ptr<StateVector<T, DomainType>> vector_copy;

      // a search of the state vectors by name and kind of entity turned up
      // empty, so add the vector to the list

      if (mymesh_.get() != &(in_vec.mesh())) {

        // IT HAS BEEN SUGGESTED THAT THIS IS A BAD IDEA TO DO UNDER
        // THE THE HOOD SINCE SOMEONE MIGHT DO IT INADVERTENTLY AND
        // NOT KNOW BUT SINCE WE HAVE TO RETURN A REFERENCE WHAT WILL
        // WE RETURN?  WE CANNOT THROW AN EXCEPTION EITHER (SINCE GPUs
        // DON'T HANDLE THOSE)

        // the input vector is defined on a different mesh? copy the
        // vector data onto a vector defined on mesh and then add

        std::cerr << "Copying data from one mesh to another???\n";

        vector_copy =
            std::make_shared<StateVector<T, DomainType>>(in_vec.name(),
                                                         mymesh_,
                                                         on_what,
                                                         ptype,
                                                         &(in_vec[0]));
      } else {
        vector_copy = std::make_shared<StateVector<T, DomainType>>(in_vec);
      }
      state_vectors_.emplace_back(vector_copy);
      
      // add the index of this vector in state_vectors_ to the vector of
      // indexes for this entity type, to allow iteration over state
      // vectors on this entity type with a permutation iterator

      int nvec = state_vectors_.size();
      int ikind = static_cast<int>(on_what);
      entity_indexes_[ikind].emplace_back(nvec-1);
      names_.emplace_back(in_vec.name());

      // emplace back may cause reallocation of the vector so the iterator
      // may not be valid. Use [] operator to get reference to vector

      return (*(std::static_pointer_cast<StateVector<T, DomainType>>(state_vectors_[nvec-1])));
    } else {
      // found a state vector by same name
      std::cerr << "Attempted to add duplicate state vector. Ignoring\n" <<
          std::endl;
      return *(std::dynamic_pointer_cast<StateVector<T, DomainType>>(*it));
    }
  }


  /// @brief Import field data from mesh
  void init_from_mesh();


  /// @brief Export field data to mesh
  void export_to_mesh();


 private:

  // Constant pointer to the mesh associated with this state
  const std::shared_ptr<Mesh> mymesh_;

  // All the state vectors
  std::vector<std::shared_ptr<BaseStateVector>> state_vectors_;

  // Stores which indices of state_vectors_ correspond to data stored
  // on each entity kind
  std::vector<int> entity_indexes_[NUM_ENTITY_KINDS];

  // Names of the state vectors
  std::vector<std::string> names_;

};


std::ostream & operator<<(std::ostream & os, State const & s);

}  // namespace Jali

#endif  // JALI_STATE_H_
