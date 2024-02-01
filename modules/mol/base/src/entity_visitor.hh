//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------
#ifndef OST_ENTITY_VISITOR_HH
#define OST_ENTITY_VISITOR_HH

#include "module_config.hh"
#include "entity_visitor_fw.hh"

#include <ost/mol/view_type_fw.hh>
#include <ost/mol/handle_type_fw.hh>

namespace ost { namespace mol {
/// \section The Entity Visitor interface
///
/// For hierarchical traversal of the entity-chain-residue-atom hierarchy
/// OpenStructure offers the concept of (simplified) visitors. By passing a
/// subclass of EntityVisitor to the the EntityHandle::Apply(EntityVisitor&)
/// member function, the whole subtree rooted at the object where Apply was
/// invoked is then traversed in a recursive manner.
///  
/// Visitors may be applied to both Handles and Views, so they are a good way
/// to generalise your algorithms for EntityView and EntitHandle the like.
///
/// Algorithms are implemented by subclassing EntityVisitor and overriding one
/// or more of the callbacks.
/// 
/// EntityVisitor provides stub implementations for all of the callbacks. You
/// only have to overload the callbacks you are interested in. Every callback
/// returns a bool to indicate whether the hierarchical traversal should
/// continue. Returning false indicates that recursive traversal should stop
/// and continue with the next sibling.
///  
/// \subsection Example
///
/// The following example prints residues and their atoms to stdout 
/// \code
/// class Printer : public EntityVisitor { public:
///   virtual bool VisitResidue(const ResidueHandle& residue) {
///     std::cout << residue.GetKey() << "." << residue.GetNumber() << std::endl;
///     return true;
///   }
///   virtual bool VisitAtom(const AtomHandle& atom) {
///     std::cout << " " << atom.GetName() << std::endl;
///     return true;
///   }
/// };
/// \endcode
///
/// \note
/// The order of traversal for chains and atoms is implementation specific and
/// your algorithms should not rely on a particular order. However, residues
/// are visited in ascending order of their residue numbers.
///


/// \brief EntityVisitor interface
class DLLEXPORT_OST_MOL EntityVisitor {
public:
  virtual ~EntityVisitor();

  virtual void OnEntry();
  virtual void OnExit();
  /// \brief Chain callback
  virtual bool VisitChain(const ChainHandle& c);
  
  /// \brief Residue callback
  virtual bool VisitResidue(const ResidueHandle& r);
  
  /// \brief Atom callback
  virtual bool VisitAtom(const AtomHandle& a);
  
  /// \brief Bond callback
  /// \todo not implemented for views
  virtual bool VisitBond(const BondHandle& b);
  
  /// \brief Torsion callback.
  /// \todo Not implemented for views
  virtual bool VisitTorsion(const TorsionHandle& t);
};

/// \brief EntityViewVisitor interface
class DLLEXPORT_OST_MOL EntityViewVisitor {
public:
  virtual ~EntityViewVisitor();

  virtual void OnEntry();
  virtual void OnExit();

  virtual bool VisitChain(const ChainView& c);
  
  virtual bool VisitResidue(const ResidueView& r);
  
  virtual bool VisitAtom(const AtomView& a);

  // is handle until bond view is implemented
  virtual bool VisitBond(const BondHandle& b);
};

}} //ns

#endif
