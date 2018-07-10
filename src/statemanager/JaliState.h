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

  The State class (or statemanager) stores shared pointers to
  StateVectors.  The StateVectors store a weak pointer to the State
  object (hence, we inherit from std::enable_shared_from_this using a
  Curiously Recurring Template Pattern)
*/

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <boost/iterator/permutation_iterator.hpp>

#include "Mesh.hh"    // Jali mesh header
#include "JaliStateVector.h"

namespace Jali {

class State : public std::enable_shared_from_this<State> {
 public:

  // Constructor is protected - call the static method create instead as:
  // std::shared_ptr<State> mystate = Jali::State::create(mesh)

  static std::shared_ptr<State> create(std::shared_ptr<Mesh> mesh) {
    // Make a struct derived from State. Since its derived from State,
    // it can call the protected constructor

    struct shared_state_enabler : public State {
      shared_state_enabler(std::shared_ptr<Mesh> mesh) : State(mesh) {}
    };

    // Now call make_shared on shared_state_enabler whose constructor is public

    return std::make_shared<shared_state_enabler>(mesh);


    // See StackOverflow.com https://stackoverflow.com/questions/8147027/how-do-i-call-stdmake-shared-on-a-class-with-only-protected-or-private-const
    //
    // or search for "shared_ptr from protected or private constructor"
  }


  /// Copy constructor (disabled)

  State(const State &) = delete;

  /// Assignment operator (disabled)

  State & operator=(const State &) = delete;

  /// Destructor

  ~State() {}

  /// Mesh

  std::shared_ptr<Jali::Mesh> mesh() {return mymesh_;}

  /// Number of materials in simulation

  int num_materials() const {
    int nsets = material_cellsets_.size();
    return nsets;
  }

  /// Name of the i'th material

  std::string material_name(int m) const {
    int nsets = material_cellsets_.size();
    return (m < nsets) ? material_cellsets_[m]->name() : "UNDEFINED_MATERIAL";
  }

  /// material index by name (-1 if not found)

  int material_index_by_name(std::string const& name) const {
    int nsets = material_cellsets_.size();
    for (int m = 0; m < nsets; m++)
      if (name == material_cellsets_[m]->name())
        return m;
    return -1;
  }

  /// Get the mesh set associated with the i'th material
  //
  // if there is only one material in the problem, return a
  // nullptr. This means the code should work with the entire mesh

  std::shared_ptr<MeshSet> material_set(int m) const {
    return material_cellsets_[m];
  }


  /// Get the number of cells in the material

  int num_material_cells(int m) const {
    int nsets = material_cellsets_.size();
    if (nsets && m < nsets)
      return material_cellsets_[m]->num_entities();
  }


  /// Get the cells in the m'th material. If no materials have been defined
  /// return reference to a dummy vector thats empty

  std::vector<int> const& material_cells(int m) const {
    int nsets = material_cellsets_.size();
    if (nsets && m < nsets)
      return material_cellsets_[m]->entities();
    else
      return dummy_cellset_->entities();
  }


  /// Get the number of materials in a cell

  int num_cell_materials(int c) const {
    int nsets = material_cellsets_.size();
    return nsets ? cell_materials_[c].size() : 0;
  }

  /// Get the materials in a cell

  std::vector<int> const& cell_materials(int c) const {
    return material_cellsets_.size() ? cell_materials_[c] : dummy_cellmats_;
  }


  /// Get the mesh cell set associated with material with given name

  std::shared_ptr<MeshSet> material_set_by_name(std::string const& name) const {
    int nsets = material_cellsets_.size();
    for (int m = 0; m < nsets; m++)
      if (name == material_cellsets_[m]->name())
        return material_cellsets_[m];
    return nullptr;
  }

  /// Add new material to the state

  void add_material(std::string const& matname,
                    std::vector<int> const& matcells);

  /// Remove material from the state -- EXPENSIVE (CAN JUST MAKE SET EMPTY)

  void rem_material(int m);

  /// Remove material from the state -- EXPENSIVE

  void rem_material_by_name(std::string const& name) {
    int m = material_index_by_name(name);
    if (m >= 0)
      rem_material(m);
  }


  /// Cell index in material set

  int cell_index_in_material(int c, int m) {
    assert(m < material_cellsets_.size());
    return material_cellsets_[m]->index_in_set(c);
  }

  /// Add cells to a material

  void add_cells_to_material(int m, std::vector<int> const& cells);

  /// Remove cells from a material (EXPENSIVE - NOT IMPLEMENTED)

  void rem_cells_from_material(int m, std::vector<int> const& cells);

  //! Typedefs for iterators for going through all the state vectors

  typedef
  std::vector<std::shared_ptr<StateVectorBase>>::iterator
  iterator;

  typedef
  std::vector<std::shared_ptr<StateVectorBase>>::const_iterator
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
    std::vector<std::shared_ptr<StateVectorBase>>::iterator,
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

  typedef std::shared_ptr<StateVectorBase> pointer;
  typedef const std::shared_ptr<StateVectorBase> const_pointer;

  /// Return pointer to i'th state vector
  pointer operator[](int i) { return state_vectors_[i]; }

  /// Return const pointer to the i'th state vector
  const_pointer operator[](int i) const { return state_vectors_[i]; }

  /// Number of state vectors
  int size() const {return state_vectors_.size();}



  /*!
    @brief Find iterator to state vector by name (does not distinguish based
    on data type of vector, UNIVAL or MULTIVAL or what domain its defined on)
    @param name        String identifier for vector
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter kind) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()

    This is the coarsest form of find - it will scoop a vector of any
    type (UNIVAL or MULTIVAL) containing any data type and defined on
    any domain (mesh, meshset or meshtile)
  */

  iterator find(std::string name,
                Entity_kind kind = Entity_kind::ANY_KIND,
                Entity_type type = Entity_type::ALL) {

    iterator it = state_vectors_.begin();
    while (it != state_vectors_.end()) {
      std::shared_ptr<StateVectorBase> bv = *it;
      if (bv->name() == name &&
          ((kind == Entity_kind::ANY_KIND) || (bv->entity_kind() == kind)) &&
          ((type == Entity_type::ALL) || (bv->entity_type() == type)))
        break;
      else
        ++it;
    }

    return it;
  }


  /*!  
    @brief Find const iterator to state vector by name (does not
    distinguish based on data type of vector, UNIVAL or MULTIVAL or
    what domain its defined on)
    @param name        String identifier for vector
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter kind) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()

    This is the coarsest form of find - it will scoop a vector of any
    type (UNIVAL or MULTIVAL) containing any data type and defined on
    any domain (mesh, meshset or meshtile)
  */

  const_iterator find(std::string name,
                      Entity_kind kind = Entity_kind::ANY_KIND,
                      Entity_type type = Entity_type::ALL) const {

    const_iterator it = state_vectors_.cbegin();
    while (it != state_vectors_.cend()) {
      std::shared_ptr<StateVectorBase> bv = *it;
      if (bv->name() == name &&
          ((kind == Entity_kind::ANY_KIND) || (bv->entity_kind() == kind)) &&
          ((type == Entity_type::ALL) || (bv->entity_type() == type)))
        break;
      else
        ++it;
    }

    return it;
  }


  /*!
    @brief Find iterator to state vector on a given domain (regardless of vector type or data type)
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find iterator to state vector by NAME, the domain it is defined
    on, what kind of entity and what type of parallel entity it is
    defined on.  The type of entity (parameter kind) may be specified
    as ANY_KIND if caller does not care if the vector is on a
    particular entity type or knows that there is only one vector by
    this name defined on a specific entity kind. The function returns
    an iterator to a state vector in the state manager if found;
    otherwise, it returns State::end()
  */

  template <class DomainType>
  iterator find(std::string name,
                std::shared_ptr<DomainType> domain,
                Entity_kind kind = Entity_kind::ANY_KIND,
                Entity_type type = Entity_type::ALL) {

    iterator it = state_vectors_.begin();
    while (it != state_vectors_.end()) {
      std::shared_ptr<StateVectorBase> bv = *it;

      if (bv->name() == name &&
          ((kind == Entity_kind::ANY_KIND) || (bv->entity_kind() == kind)) &&
          ((type == Entity_type::ALL) || (bv->entity_type() == type))) {

        StateVector_type bvectype = bv->type();
        if (bvectype == StateVector_type::UNIVAL) {
          auto uvec =
              std::dynamic_pointer_cast<UniStateVectorBase<DomainType>>(bv);
          if (uvec && uvec->domain() == domain)
            break;
        } else if (bvectype == StateVector_type::MULTIVAL) {
          auto mvec =
              std::dynamic_pointer_cast<MultiStateVectorBase<DomainType>>(bv);
          if (mvec && mvec->domain() == domain)
            break;
        }
      }
      ++it;
    }

    return it;
  }




  /*!
    @brief Find const iterator to state vector on a given domain (regardless of vector type or data type)
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find const iterator to state vector by NAME, the domain it is defined
    on, what kind of entity and what type of parallel entity it is
    defined on.  The type of entity (parameter kind) may be specified
    as ANY_KIND if caller does not care if the vector is on a
    particular entity type or knows that there is only one vector by
    this name defined on a specific entity kind. The function returns
    an iterator to a state vector in the state manager if found;
    otherwise, it returns State::end()
  */

  template <class DomainType>
  const_iterator find(std::string name,
                      std::shared_ptr<DomainType> domain,
                      Entity_kind kind = Entity_kind::ANY_KIND,
                      Entity_type type = Entity_type::ALL) const {

    const_iterator it = state_vectors_.cbegin();
    while (it != state_vectors_.cend()) {
      std::shared_ptr<StateVectorBase> bv = *it;

      if (bv->name() == name &&
          ((kind == Entity_kind::ANY_KIND) || (bv->entity_kind() == kind)) &&
          ((type == Entity_type::ALL) || (bv->entity_type() == type))) {

        StateVector_type bvectype = bv->type();
        if (bvectype == StateVector_type::UNIVAL) {
          auto uvec =
              std::dynamic_pointer_cast<UniStateVectorBase<DomainType>>(bv);
          if (uvec && uvec->domain() == domain)
            break;
        } else if (bvectype == StateVector_type::MULTIVAL) {
          auto mvec =
              std::dynamic_pointer_cast<MultiStateVectorBase<DomainType>>(bv);
          if (mvec && mvec->domain() == domain)
            break;
        }
      }
      ++it;
    }

    return it;
  }


  /*!
    @brief Find iterator to state vector on a given domain (regardless of data type)
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param vectype     Enum type of state vector (UNIVAL or MULTIVAL)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find iterator to state vector by NAME, the domain it is defined
    on, what kind of entity and what type of parallel entity it is
    defined on.  The type of entity (parameter kind) may be specified
    as ANY_KIND if caller does not care if the vector is on a
    particular entity type or knows that there is only one vector by
    this name defined on a specific entity kind. The function returns
    an iterator to a state vector in the state manager if found;
    otherwise, it returns State::end()
  */

  template <class DomainType, StateVector_type vectype>
  iterator find(std::string name,
                std::shared_ptr<DomainType> domain,
                Entity_kind kind = Entity_kind::ANY_KIND,
                Entity_type type = Entity_type::ALL) {

    iterator it = state_vectors_.begin();
    while (it != state_vectors_.end()) {
      std::shared_ptr<StateVectorBase> bv = *it;

      if (bv->name() == name &&
          ((kind == Entity_kind::ANY_KIND) || (bv->entity_kind() == kind)) &&
          ((type == Entity_type::ALL) || (bv->entity_type() == type))) {

        StateVector_type bvectype = bv->type();
        if (vectype == StateVector_type::UNIVAL &&
            bvectype == StateVector_type::UNIVAL) {
          auto uvec =
              std::dynamic_pointer_cast<UniStateVectorBase<DomainType>>(bv);
          if (uvec && uvec->domain() == domain)
            break;
        } else if (vectype == StateVector_type::MULTIVAL &&
                   bvectype == StateVector_type::MULTIVAL) {
          auto mvec =
              std::dynamic_pointer_cast<MultiStateVectorBase<DomainType>>(bv);
          if (mvec && mvec->domain() == domain)
            break;
        }
      }
      ++it;
    }

    return it;
  }




  /*!
    @brief Find const iterator to state vector on a given domain (regardless of data type)
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param vectype     Enum type of state vector (UNIVAL or MULTIVAL)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find const iterator to state vector by NAME, the domain it is defined
    on, what kind of entity and what type of parallel entity it is
    defined on.  The type of entity (parameter kind) may be specified
    as ANY_KIND if caller does not care if the vector is on a
    particular entity type or knows that there is only one vector by
    this name defined on a specific entity kind. The function returns
    an iterator to a state vector in the state manager if found;
    otherwise, it returns State::end()
  */

  template <class DomainType, StateVector_type vectype>
  const_iterator find(std::string name,
                std::shared_ptr<DomainType> domain,
                Entity_kind kind = Entity_kind::ANY_KIND,
                Entity_type type = Entity_type::ALL) const {

    iterator it = state_vectors_.cbegin();
    while (it != state_vectors_.cend()) {
      std::shared_ptr<StateVectorBase> bv = *it;

      if (bv->name() == name &&
          ((kind == Entity_kind::ANY_KIND) || (bv->entity_kind() == kind)) &&
          ((type == Entity_type::ALL) || (bv->entity_type() == type))) {

        StateVector_type bvectype = bv->type();
        if (vectype == StateVector_type::UNIVAL &&
            bvectype == StateVector_type::UNIVAL) {
          auto uvec =
              std::dynamic_pointer_cast<UniStateVectorBase<DomainType>>(bv);
          if (uvec && uvec->domain() == domain)
            break;
        } else if (vectype == StateVector_type::MULTIVAL &&
                   bvectype == StateVector_type::MULTIVAL) {
          auto mvec =
              std::dynamic_pointer_cast<MultiStateVectorBase<DomainType>>(bv);
          if (mvec && mvec->domain() == domain)
            break;
        }
      }
      ++it;
    }

    return it;
  }


  /*!
    @brief Find iterator to state vector templated on data type T by name
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param vectype     Enum type of state vector (UNIVAL or MULTIVAL)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter kind) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()
  */

  template<class T, class DomainType,
           StateVector_type vectype = StateVector_type::UNIVAL>
  iterator find(std::string name,
                std::shared_ptr<DomainType> domain,
                Entity_kind kind = Entity_kind::ANY_KIND,
                Entity_type type = Entity_type::ALL) {

    if (vectype == StateVector_type::UNIVAL)
      return find<T, DomainType, UniStateVector>(name, domain, kind, type);
    else if (vectype == StateVector_type::MULTIVAL)
      return find<T, DomainType, MultiStateVector>(name, domain, kind, type);
  }


  /*!
    @brief Find const iterator to state vector of data type T by name
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param vectype     Enum type of state vector (UNIVAL or MULTIVAL)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter kind) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()
  */

  template<class T, class DomainType,
           StateVector_type vectype = StateVector_type::UNIVAL>
  const_iterator find(std::string name,
                std::shared_ptr<DomainType> domain,
                Entity_kind kind = Entity_kind::ANY_KIND,
                Entity_type type = Entity_type::ALL) const {

    if (vectype == StateVector_type::UNIVAL)
      return find<T, DomainType, UniStateVector>(name, domain, kind, type);
    else if (vectype == StateVector_type::MULTIVAL)
      return find<T, DomainType, MultiStateVector>(name, domain, kind, type);
  }




  /*!
    @brief Find iterator to state vector to data type T by name
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter kind) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()
  */

  template<class T, class DomainType,
           template<class /* T */, class /* DomainType */> class StateVecType>
  iterator find(std::string name,
                std::shared_ptr<DomainType> domain,
                Entity_kind kind = Entity_kind::ANY_KIND,
                Entity_type type = Entity_type::ALL) {

    iterator it = state_vectors_.begin();
    while (it != state_vectors_.end()) {
      std::shared_ptr<StateVecType<T, DomainType>> sv =
          std::dynamic_pointer_cast<StateVecType<T, DomainType>>(*it);

      // Check if we were able to cast the shared_ptr to BaseVector to
      // StateVecType and all the other characteristics match
      if (sv && (sv->name() == name) && (sv->domain() == domain) &&
          ((kind == Entity_kind::ANY_KIND) || (sv->entity_kind() == kind)) &&
          ((type == Entity_type::ALL) || (sv->entity_type() == type)))
        break;
      else
        ++it;
    }

    return it;
  }



  /*!
    @brief Find a const iterator to state vector of data type T by name
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, BOUNDARY_GHOST, ALL etc.)

    Find iterator to state vector by NAME and what kind of entity and
    what type of parallel entity it is defined on.  The type of entity
    (parameter kind) may be specified as ANY_KIND if caller does
    not care if the vector is on a particular entity type or knows
    that there is only one vector by this name defined on a specific
    entity kind. The function returns an iterator to a state vector in
    the state manager if found; otherwise, it returns State::end()
  */

  template<class T, class DomainType,
           template<class /* T */, class /* DomainType */> class StateVecType>
  const_iterator find(std::string name,
                      std::shared_ptr<DomainType> domain,
                      Entity_kind kind = Entity_kind::ANY_KIND,
                      Entity_type type = Entity_type::ALL) const {

    const_iterator it = state_vectors_.cbegin();
    while (it != state_vectors_.cend()) {
      std::shared_ptr<StateVecType<T, DomainType>> sv =
          std::dynamic_pointer_cast<StateVecType<T, DomainType>>(*it);

      // Check if we were able to cast the shared_ptr to BaseVector to
      // StateVecType and all the other characteristics match
      if (sv && (sv->name() == name) && (sv->domain() == domain) &&
          ((kind == Entity_kind::ANY_KIND) || (sv->entity_kind() == kind)) &&
          ((type == Entity_type::ALL) || (sv->entity_type() == type)))
        break;
      else
        ++it;
    }

    return it;
  }



  /*!
    @brief Retrieve a state vector by name given the domain and entity type it
    is defined on
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, etc.)
    @param vector      Pointer to the state vector

    Get a statevector with the given name, particular entity kind and
    type of entity it is defined on. The function returns true if such
    a vector was found; false otherwise. The caller must know the type
    of elements in the state vector, int, double, std::array<double,
    3> or whatever else. The calling routine must declare the state
    vector as "T vector" where T is UniStateVector<int> or
    UniStateVector<double> or StateVector<some_other_type>. Even though
    this is a copy into *vector, its an inexpensive shallow copy of
    the meta data only
  */

  template <class T, class DomainType,
            template<class /* T */, class /* DomainType */> class StateVecType>
  bool get(std::string name,
           std::shared_ptr<DomainType> domain,
           Entity_kind kind,
           Entity_type type,
           StateVecType<T, DomainType> *vector) {

    iterator it = find<T, DomainType, StateVecType>(name, domain, kind, type);
    if (it != state_vectors_.end()) {
      *vector = *(std::dynamic_pointer_cast<StateVecType<T, DomainType>>(*it));
      return true;
    } else {
      return false;
    }
  }


  /*!
    @brief Retrieve a const state vector by name given the domain and entity
    type it is defined on
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, etc.)
    @param vector      Pointer to the state vector

    Get a statevector with the given name, particular entity kind and
    type of entity it is defined on. The function returns true if such
    a vector was found; false otherwise. The caller must know the type
    of elements in the state vector, int, double, std::array<double,
    3> or whatever else. The calling routine must declare the state
    vector as "T vector" where T is UniStateVector<int> or
    UniStateVector<double> or UniStateVector<some_other_type>. Even though
    this is a copy into *vector, its an inexpensive shallow copy of
    the meta data only
  */

  template <class T, class DomainType,
            template<class /* T */, class /* DomainType */> class StateVecType>
  bool get(std::string name,
           std::shared_ptr<DomainType> domain,
           Entity_kind kind,
           Entity_type type,
           StateVecType<T, DomainType> const *vector) const {

    const_iterator it = find<T, DomainType, StateVecType>(name, domain, kind, type);
    if (it != state_vectors_.cend()) {
      *vector = *(std::dynamic_pointer_cast<StateVecType<T, DomainType> const>(*it));
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
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name         String identifier for vector
    @param vector       Pointer to the state vector

    Get a statevector on the mesh with the given name without regard
    to what kind of entity it is defined on. The function returns true
    if such a vector was found; false otherwise. The caller must know
    the type of elements in the state vector, int, double,
    std::array<double, 3> or whatever else. The calling routine must
    declare the state vector as "T vector" where T is UniStateVector<int>
    or UniStateVector<double> or UniStateVector<some_other_type>. Even
    though this is a copy into *vector, its an inexpensive shallow
    copy of the meta data only
  */

  template <class T, class DomainType,
            template <class /* T */, class /* DomainType */> class StateVecType>
  bool get(std::string name, StateVecType<T, DomainType> *vector) {
    return get(name, mymesh_, Entity_kind::ANY_KIND, Entity_type::ALL,
               vector);
  }





  /*!
    @brief Retrieve a const state vector on the mesh by name (regardless of
    what type of entity it is on)
    @tparam T           Data type
    @tparam DomainType  Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name         String identifier for vector
    @param vector       Pointer to the state vector

    Get a statevector on the mesh with the given name without regard
    to what kind of entity it is defined on. The function returns true
    if such a vector was found; false otherwise. The caller must know
    the type of elements in the state vector, int, double,
    std::array<double, 3> or whatever else. The calling routine must
    declare the state vector as "T vector" where T is UniStateVector<int>
    or UniStateVector<double> or UniStateVector<some_other_type>. Even
    though this is a copy into *vector, its an inexpensive shallow
    copy of the meta data only
  */

  template <class T, class DomainType,
            template <class /* T */, class /* DomainType */> class StateVecType>
  bool get(std::string name,
           StateVecType<T, DomainType> const *vector) const {
    return get(name, mymesh_, Entity_kind::ANY_KIND, Entity_type::ALL,
               vector);
  }





  /*!
    @brief Retrieve a shared pointer to a state vector by name given the domain
    and type of entity it is defined on
    @tparam T            Data type
    @tparam DomainType   Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name          String identifier for vector
    @param domain        Shared pointer to the domain
    @param kind          What kind of entity data is defined on (CELL, NODE, etc.)
    @param type          What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, etc.)
    @param vector_ptr    Shared pointer to the state vector

    Get a statevector on given domain with the given name, particular
    entity kind and type of entity it is defined on The function
    returns true if such a vector was found; false otherwise. The
    caller must know the type of elements in the state vector, int,
    double, std::array<double, 3> or whatever else. The calling
    routine must declare the state vector as "T vector" where T is
    UniStateVector<int> or UniStateVector<double> or
    UniStateVector<some_other_type>. Even though this is a copy into
    *vector, its an inexpensive shallow copy of the meta data only
  */

  template <class T, class DomainType,
            template <class /* T */, class /* DomainType */> class StateVecType>
  bool get(std::string name,
           std::shared_ptr<DomainType> domain,
           Entity_kind kind,
           Entity_type type,
           std::shared_ptr<StateVecType<T, DomainType>> *vector_ptr) {

    iterator it = find<T, DomainType, StateVecType>(name, domain, kind, type);
    if (it != state_vectors_.end()) {
      *vector_ptr = std::dynamic_pointer_cast<StateVecType<T, DomainType>>(*it);
      return true;
    } else {
      return false;
    }
  }


  /*!
    @brief Retrieve a shared pointer to a const state vector by name given
    the domain
    and type of entity it is defined on
    @tparam T            Data type
    @tparam DomainType   Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name          String identifier for vector
    @param domain        Shared pointer to the domain
    @param kind          What kind of entity data is defined on (CELL, NODE, etc.)
    @param type          What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, etc.)
    @param vector_ptr    Shared pointer to the state vector

    Get a statevector on given domain with the given name, particular
    entity kind and type of entity it is defined on The function
    returns true if such a vector was found; false otherwise. The
    caller must know the type of elements in the state vector, int,
    double, std::array<double, 3> or whatever else. The calling
    routine must declare the state vector as "T vector" where T is
    UniStateVector<int> or UniStateVector<double> or
    UniStateVector<some_other_type>. Even though this is a copy into
    *vector, its an inexpensive shallow copy of the meta data only
  */

  template <class T, class DomainType,
            template <class /* T */, class /* DomainType */> class StateVecType>
  bool get(std::string name,
           std::shared_ptr<DomainType> domain,
           Entity_kind kind,
           Entity_type type,
           std::shared_ptr<StateVecType<T, DomainType> const> *vector_ptr) const {

    const_iterator it = find<T, DomainType, StateVecType>(name, domain, kind, type);
    if (it != state_vectors_.cend()) {
      *vector_ptr = std::dynamic_pointer_cast<StateVecType<T, DomainType> const>(*it);
      return true;
    } else {
      return false;
    }
  }



  /*!
    @brief Retrieve a shared pointer to a state vector by name given the
    domain it is defined on (regardless of entity type it is on)
    @tparam T           Data type
    @tparam DomainType  Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name         String identifier for vector
    @param domain       Shared pointer to the domain
    @param vector_ptr   Shared pointer to the state vector

    Get a statevector with the given integer identifier without regard
    to what kind of entity it is defined on. The function returns true
    if such a vector was found; false otherwise. The caller must know
    the type of elements in the state vector, int, double,
    std::array<double, 3> or whatever else. The calling routine must
    declare the state vector as "T vector" where T is UniStateVector<int>
    or UniStateVector<double> or UniStateVector<some_other_type>. Even
    though this is a copy into *vector, its an inexpensive shallow
    copy of the meta data only

    NOT SURE THAT THERE IS ANY UTILITY IN RETURNING A SHARED POINTER
    HERE - THE DATA IS ALREADY STORED AS A SHARED POINTER, SO THERE IS
    NO RISK OF DELETING THE DATA WHEN AT LEAST ONE COPY OF THE STATE
    VECTOR POINTS TO THE DATA. THE STATE VECTOR ITSELF IS EASIER TO
    DEAL WITH DIRECTLY INSTEAD OF THROUGH A SHARED POINTER
  */

  template <class T, class DomainType,
            template <class /* T */, class /* DomainType */> class StateVecType>
  bool get(std::string name,
           std::shared_ptr<DomainType> domain,
           std::shared_ptr<StateVecType<T, DomainType>> *vector_ptr) {

    return get(name, domain, Entity_kind::ANY_KIND, Entity_type::ALL,
               vector_ptr);
  }




  /*!
    @brief Retrieve a shared pointer to a const state vector by name given
    the domain it is defined on (regardless of entity type it is on)
    @tparam T           Data type
    @tparam DomainType  Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name         String identifier for vector
    @param domain       Shared pointer to the domain
    @param vector_ptr   Shared pointer to the state vector

    Get a statevector with the given integer identifier without regard
    to what kind of entity it is defined on. The function returns true
    if such a vector was found; false otherwise. The caller must know
    the type of elements in the state vector, int, double,
    std::array<double, 3> or whatever else. The calling routine must
    declare the state vector as "T vector" where T is UniStateVector<int>
    or UniStateVector<double> or UniStateVector<some_other_type>. Even
    though this is a copy into *vector, its an inexpensive shallow
    copy of the meta data only

    NOT SURE THAT THERE IS ANY UTILITY IN RETURNING A SHARED POINTER
    HERE - THE DATA IS ALREADY STORED AS A SHARED POINTER, SO THERE IS
    NO RISK OF DELETING THE DATA WHEN AT LEAST ONE COPY OF THE STATE
    VECTOR POINTS TO THE DATA. THE STATE VECTOR ITSELF IS EASIER TO
    DEAL WITH DIRECTLY INSTEAD OF THROUGH A SHARED POINTER
  */

  template <class T, class DomainType,
            template <class /* T */, class /* DomainType */> class StateVecType>
  bool get(std::string name,
           std::shared_ptr<DomainType> domain,
           std::shared_ptr<StateVecType<T, DomainType> const> *vector_ptr) const {

    return get(name, domain, Entity_kind::ANY_KIND, Entity_type::ALL,
               vector_ptr);
  }




  /*!
    @brief Add an uninitialized state vector using a string identifier
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @tparam StateVecType  State vector class (UniStateVector or MultiStateVector)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, etc.)
    @param data        Raw pointer to data array (optional)

    Add state vector - returns reference to the added UniStateVector or
    MultiStateVector.

  */

  template <class T, class DomainType,
            template<class /* T */, class /* DomainType */> class StateVecType>
  StateVecType<T, DomainType>& add(std::string name,
                                   std::shared_ptr<DomainType> domain,
                                   Entity_kind kind,
                                   Entity_type type) {

    iterator it = find<T, DomainType, StateVecType>(name, domain, kind, type);
    if (it == end()) {
      // a search of the state vectors by name and kind of entity turned up
      // empty, so add the vector to the list; if not, warn about duplicate
      // state data

      // add the index of this vector in state_vectors_ to the vector of
      // indexes for this entity type, to allow iteration over state
      // vectors on this entity type with a permutation iterator

      int ikind = static_cast<int>(kind);
      entity_indexes_[ikind].emplace_back(state_vectors_.size()-1);
      names_.emplace_back(name);

      auto vector =
          std::make_shared<StateVecType<T, DomainType>>(name, domain,
                                                        shared_from_this(),
                                                        kind, type);
      state_vectors_.emplace_back(vector);
      return (*vector);
    } else {
      // found a state vector by same name
      std::cerr <<
          "Attempted to add duplicate state vector. Ignoring\n" << std::endl;
      return
          (*(std::dynamic_pointer_cast<StateVecType<T, DomainType>>(*it)));
    }
  }




  /*!
    @brief Add a single valued state vector (class UniStateVector) using a string identifier and optional array data
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, etc.)
    @param data        Raw pointer to data array

    Add state vector - returns reference to the added UniStateVector.
    The data is copied from the input memory location into an
    internal buffer.
  */

  template <class T, class DomainType>
  UniStateVector<T, DomainType>& add(std::string name,
                                  std::shared_ptr<DomainType> domain,
                                  Entity_kind kind,
                                  Entity_type type,
                                  T const * const data) {

    iterator it = find<T, DomainType, UniStateVector>(name, domain, kind, type);
    if (it == end()) {
      // a search of the state vectors by name and kind of entity turned up
      // empty, so add the vector to the list; if not, warn about duplicate
      // state data

      // add the index of this vector in state_vectors_ to the vector of
      // indexes for this entity type, to allow iteration over state
      // vectors on this entity type with a permutation iterator

      int ikind = static_cast<int>(kind);
      entity_indexes_[ikind].emplace_back(state_vectors_.size()-1);
      names_.emplace_back(name);

      auto vector =
          std::make_shared<UniStateVector<T, DomainType>>(name, domain,
                                                       shared_from_this(),
                                                       kind, type,
                                                       data);
      state_vectors_.emplace_back(vector);
      return (*vector);
    } else {  // found a state vector by same name
      std::cerr << "Attempted to add duplicate state vector. Ignoring\n";
      return (*(std::dynamic_pointer_cast<UniStateVector<T, DomainType>>(*it)));
    }
  }




  /*!
    @brief Add a multi-valued state vector (class MultiStateVector) using a string identifier and 2D array data
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, etc.)
    @param data        Raw pointer to data array

    Add state vector - returns reference to the added MultiStateVector.
    The data is copied from the input memory location into an
    internal buffer.
  */

  template <class T, class DomainType>
  MultiStateVector<T, DomainType>& add(std::string name,
                                    std::shared_ptr<DomainType> domain,
                                    Entity_kind kind,
                                    Entity_type type,
                                    Data_layout  layout,
                                    T const * const * data) {

    iterator it = find<T, DomainType, MultiStateVector>(name, domain, kind, type);
    if (it == end()) {
      // a search of the state vectors by name and kind of entity turned up
      // empty, so add the vector to the list; if not, warn about duplicate
      // state data

      // add the index of this vector in state_vectors_ to the vector of
      // indexes for this entity type, to allow iteration over state
      // vectors on this entity type with a permutation iterator

      int ikind = static_cast<int>(kind);
      entity_indexes_[ikind].emplace_back(state_vectors_.size()-1);
      names_.emplace_back(name);

      auto vector =
          std::make_shared<MultiStateVector<T, DomainType>>(name, domain,
                                                         shared_from_this(),
                                                         kind, type, layout,
                                                         data);
      state_vectors_.emplace_back(vector);
      return (*vector);
    } else {  // found a state vector by same name
      std::cerr << "Attempted to add duplicate state vector. Ignoring\n";
      return (*(std::dynamic_pointer_cast<MultiStateVector<T, DomainType>>(*it)));
    }
  }




  /*!
    @brief Add state vector using a string identifier and a single value
    @tparam T          Data type
    @tparam DomainType Type of domain data is defined on (Mesh, MeshTile)
    @param name        String identifier for vector
    @param domain      Shared pointer to the domain
    @param kind        What kind of entity data is defined on (CELL, NODE, etc.)
    @param type        What type of entity data is defined on (PARALLEL_OWNED, PARALLEL_GHOST, etc.)
    @param data        Value of all elements

    Add state vector - returns reference to the added StateVector.
    The data is copied from the input memory location into an
    internal buffer.

    ************************** NOTE *******************************
    This version of the overloaded operator is being DISABLED for
    pointer and array types (via the line 'typename
    std::enable_if....type') because complex template deduction rules
    are making the compiler invoke the reference version, when we call
    it with a non-const pointer.

    So, if we call

    double data[10];
    std::vector<double>& myvec = mystate.add(......, data)

    it thinks the template argument T is "double *" instead of "double"

    See, stackoverflow.com Q&A

    http://stackoverflow.com/questions/13665574/template-argument-deduction-and-pointers-to-constants


    We could make it work for some cases using

    template <class T, class DomainType,
              template<class, class> class StateVecType>
    auto add(........,
             T const& data) -> StateVecType<decltype(data+data), DomainType>&

    but this does not work if T is a double[3] or std::array<double, 3>
    as there is no + operator defined for these types

    ******************************************************************

    */

  template <class T, class DomainType,
            template<class /* T */, class /* DomainType */> class StateVecType>
  typename
  std::enable_if<(!std::is_pointer<T>::value && !std::is_array<T>::value),
      StateVecType<T, DomainType>>::type &
      add(std::string name, std::shared_ptr<DomainType> domain,
          Entity_kind kind, Entity_type type, T const& data) {

    iterator it = find<T, DomainType, StateVecType>(name, domain, kind, type);
    if (it == end()) {
      // a search of the state vectors by name and kind of entity turned up
      // empty, so add the vector to the list; if not, warn about duplicate
      // state data

      auto vector =
          std::make_shared<StateVecType<T, DomainType>>(name, domain,
                                                        shared_from_this(),
                                                        kind, type, data);
      state_vectors_.emplace_back(vector);

      // add the index of this vector in state_vectors_ to the vector of
      // indexes for this entity type, to allow iteration over state
      // vectors on this entity type with a permutation iterator

      int ikind = static_cast<int>(kind);
      entity_indexes_[ikind].emplace_back(state_vectors_.size()-1);
      names_.emplace_back(name);

      // emplace back may cause reallocation of the vector so the iterator
      // may not be valid. Use [] operator to get reference to vector

      return (*vector);
    } else {
      // found a state vector by same name
      std::cerr << "Attempted to add duplicate state vector. Ignoring\n" <<
          std::endl;
      return (*(std::dynamic_pointer_cast<StateVecType<T, DomainType>>(*it)));
    }
  }



  /*!
    @brief Add state vectory by copying an input vector

    Add a new state vector to the state manager based on data from
    the input state vector. Meta data is copied from one vector to
    another and A DEEP COPY IS MADE of the input vector data. The
    need for a deep copy is why we cannot just assign in_vector to
    vector copy even if they are on the same domain
  */

  template <class T, class DomainType,
            template<class /* T */, class /* DomainType */> class StateVecType>
    StateVecType<T, DomainType>&
  add(StateVecType<T, DomainType> const& in_vec) {

    Entity_kind kind = in_vec.entity_kind();
    Entity_type type = in_vec.entity_type();

    iterator it = find<T, DomainType, StateVecType>(in_vec.name(),
                                                    in_vec.domain(), kind,
                                                    type);
    if (it == end()) {

      std::shared_ptr<StateVecType<T, DomainType>> vector_copy;

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
            std::make_shared<StateVecType<T, DomainType>>(in_vec.name(),
                                                          mymesh_,
                                                          shared_from_this(),
                                                          kind, type,
                                                          &(in_vec[0]));
      } else {
        vector_copy = std::make_shared<StateVecType<T, DomainType>>(in_vec);
      }
      state_vectors_.emplace_back(vector_copy);

      // add the index of this vector in state_vectors_ to the vector of
      // indexes for this entity type, to allow iteration over state
      // vectors on this entity type with a permutation iterator

      int nvec = state_vectors_.size();
      int ikind = static_cast<int>(kind);
      entity_indexes_[ikind].emplace_back(nvec-1);
      names_.emplace_back(in_vec.name());

      // emplace back may cause reallocation of the vector so the iterator
      // may not be valid. Use [] operator to get reference to vector

      return (*(std::dynamic_pointer_cast<StateVecType<T, DomainType>>(state_vectors_[nvec-1])));
    } else {
      // found a state vector by same name
      std::cerr << "Attempted to add duplicate state vector. Ignoring\n" <<
          std::endl;
      return *(std::dynamic_pointer_cast<StateVecType<T, DomainType>>(*it));
    }
  }


  /// @brief Import field data from mesh
  void init_from_mesh();


  /// @brief Export field data to mesh
  void export_to_mesh();

 protected:

  /// Constructor (Private - Use create_state)
  //
  // Private because we will have to call shared_from_this() within
  //  the add functions to send a shared_ptr to state vector
  //  constructors in which case State cannot be created on the stack,
  //  only on the heap

  explicit State(std::shared_ptr<Jali::Mesh> mesh) : mymesh_(mesh) {
    Entity_ID_List dummy_owned_cells, dummy_ghost_cells;
    dummy_cellset_ = std::make_shared<MeshSet>("dummy_cellset_",
                                               *mesh, Entity_kind::CELL,
                                               dummy_owned_cells,
                                               dummy_ghost_cells, false);
  }


 private:

  // Constant pointer to the mesh associated with this state
  const std::shared_ptr<Mesh> mymesh_;

  // Meshsets associated with materials. If there is only one
  // material, there will be no meshset stored because its the whole
  // mesh
  std::vector<std::shared_ptr<MeshSet>> material_cellsets_;

  // One dummy material set (to return when no materials are defined)
  std::shared_ptr<MeshSet> dummy_cellset_;

  // Dummy material vector (to return in cell_materials when no mats
  // are defined)
  std::vector<int> dummy_cellmats_;

  // Lists of materials in cells
  std::vector<std::vector<int>> cell_materials_;

  // All the state vectors
  std::vector<std::shared_ptr<StateVectorBase>> state_vectors_;

  // Stores which indices of state_vectors_ correspond to data stored
  // on each entity kind
  std::vector<int> entity_indexes_[NUM_ENTITY_KINDS];

  // Names of the state vectors
  std::vector<std::string> names_;

};

std::ostream & operator<<(std::ostream & os, State const & s);

}  // namespace Jali

#endif  // JALI_STATE_H_
