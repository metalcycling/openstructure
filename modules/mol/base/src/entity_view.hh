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
#ifndef OST_ENTITY_VIEW_HN
#define OST_ENTITY_VIEW_HN
/*
   Author: Marco Biasini
 */
#include <ost/mol/module_config.hh>

#include <ost/mol/residue_prop.hh>
#include <ost/mol/entity_base.hh>
#include <ost/mol/view_type_fw.hh>
#include <ost/mol/impl/entity_impl_fw.hh>
#include <ost/mol/entity_visitor_fw.hh>
#include <ost/mol/handle_type_fw.hh>
#include <ost/mol/query.hh>
#include <ost/mol/property_id.hh>
#include <ost/geom/geom.hh>

#ifdef _MSC_VER
# ifdef FindAtom
#  undef FindAtom
# endif
# ifdef AddAtom
#  undef AddAtom
# endif
#endif

namespace ost { namespace mol {

/// \brief definition of EntityView
/// 
/// \section views Views
///
/// Views are very widely used as an input to algorithms and the link between
/// graphical objects and entities.
/// 
/// Views point to the same data as their corresponding handle classes. When
/// modifying and entity, such as applying a coordinate transform to the atoms,
/// these changes will also be reflected in the views.
///
/// The entity - chain -residue - atom hierarchy of the entity handles is
/// replicated in the views, the most important difference being that not all
/// atoms, residues and chains present in the entity must be present in the 
/// view. One can think of the views as a \em selection.
/// 
/// \subsection creating_views Creating Views
/// 
/// Views can be assembled in several different ways: To create a full view
/// containing all the atoms, residues, chains and bonds of the entity, use
/// EntityHandle::CreateFullView(), an empty view containing basically nothing
/// is created with EntityHandle::CreateEmptyView().
/// 
/// Elements can be added with EntitView::AddAtom(), 
/// EntityView::AddResidue(), EntityView::AddChain(). The behaviour of these
/// methods is influenced by a combination of ViewAddFlag. Note that none of
/// the \em add methods insert bonds. Bonds can be added with 
/// EntityView::AddAllInclusiveBonds() or manually with EntityView::AddBond().
/// The elements can be removed with EntityView::RemoveChain(), 
/// EntityView::RemoveResidue() and EntityView::RemoveAtom().
///
/// A very convenient way to create views is with so-called \ref Query query
/// statements. The queries are executed via EntityHandle::Select() and
/// EntityView::Select().
///
/// \subsection bool_view_ops Boolean Operations on Views
///
/// Two views can be combined by applying boolean operators such as Union(), 
/// Intersection() and Difference().
class DLLEXPORT_OST_MOL EntityView: public EntityBase {
  friend class AtomView;
  friend class ResidueView;
public:
  EntityView();
  /// \brief create new entity view from entity handle
  ///
  /// \param entity The handle the view is created from.
  /// \param flags An ORed together list of ViewAddFlag
  explicit EntityView(const EntityHandle& entity,
                      ViewAddFlags flags=0);

  /// \name View validity
  //@{
  /// \brief check validity of handle
  /// 
  /// check, whether the entity view points to a valid entity.
  /// \note It is an error to use any method other than #IsValid, #Impl and 
  ///       #operator bool() when the handle is invalid. An InvalidHandle
  ///       exception will be thrown.
  operator bool() const { return this->IsValid(); }
  /// \brief check validity of handle
  /// \sa #operator bool()
  bool IsValid() const { return data_.get()!=0; }
  //@}
  /// \brief Get entity's mass
  Real GetMass() const;
             
  /// \brief Get entity's center of mass (mass weighted)
  geom::Vec3 GetCenterOfMass() const;
  
  /// \brief Get entity's center of atoms (not mass weighted)
  ///
  /// Returns the center of all the atoms in an entity. This is
  /// similar to GetCenterOfMass(), but the atoms are not mass weighted
  geom::Vec3 GetCenterOfAtoms() const;
  
  geom::AlignedCuboid GetBounds() const;
  
public:
  /// \name Internal
  //@{
  /// \brief construct entity view from data and entity impl
  EntityView(const EntityViewDataPtr& data,
             const impl::EntityImplPtr impl);

  /// \brief Get entity view-specific data
  const EntityViewDataPtr& ViewData() const; 
  /// \brief Get entity view-specific data
  EntityViewDataPtr& ViewData();
  //@}
public:
  /// \brief create copy of view
  /// 
  /// The copy has the semantics of a deep copy.
  EntityView Copy() const;
  /// \brief Add chain to view
  ///
  /// If ViewAddFlag::CHECK_DUPLICATES is set, the method will ensure that the
  /// handle is not already included in the view.
  ///
  /// \param chain_handle The chain handle to add
  /// \param flags Can be a ORed together combination of ViewAddFlag
  /// \return  The added chain view, or a reference to the existing chain view
  ///     that points to the very same handle.
  ChainView AddChain(const ChainHandle& chain_handle,
                     ViewAddFlags flags=0);
                     
  /// \brief Add chain to view
  ///
  /// Same as AddChain(const ChainHandle&) but with the difference that only
  /// the residues or atoms included in the chain view are included when one of
  /// ViewAddFlag::INCLUDE_ATOMS or ViewAddFlag::INCLUDE_RESIDUES is set.
  /// \param chain_view the handle that the view points to will be added to the view.
  /// \param flags can be a ORed together combination of ViewAddFlag
  /// \return  The added chain view, or a reference to the existing chain view
  ///      that points to the very same handle.
  ChainView AddChain(const ChainView& chain_view,
                     ViewAddFlags flags=0);

  /// \brief Add residue to view
  ///
  /// If the residues's chain has not been added to the chain already, it will
  /// be added.
  /// \param residue_handle is the residue to add
  /// \param flags can be a ORed together combination of ViewAddFlag
  /// \return The added residue view, or a reference to the existing residue
  ///     view that points to the very same handle.
  ResidueView AddResidue(const ResidueHandle& residue_handle,
                         ViewAddFlags flags=0);
  /// \brief Add residue to view
  /// 
  /// Same as AddReside(const ResidueHandle&, ViewAddFlags) but with the 
  /// difference that only atoms included in the residue view are added
  /// when ViewAddFlag::INCLUDE_ATOMS is set.
  /// 
  /// \sa AddResidue(const ResidueHandle&, ViewAddFlags)
  ResidueView AddResidue(const ResidueView& residue_view,
                         ViewAddFlags flags=0);
  /// \brief Add atom to view
  ///
  /// If the atom's chain  or residue have not been added to the view, they
  /// will be added
  ///
  /// \param atom_handle is the atom to add
  /// \param flags can be a ORed together combination of ViewAddFlag
  /// \return The added atom view, or a reference to the existing atom view
  ///     that points to the very same handle.
  AtomView AddAtom(const AtomHandle& atom_handle,
                   ViewAddFlags flags=0);

  AtomView AddAtom(const AtomView& atom_view,
                  ViewAddFlags flags=0);
  /// \brief Select a structural subset of the entity view
  /// \sa Query
  EntityView Select(const String& query_string, QueryFlags flags=0) const;

  /// \brief Select a structural subset of the entity view
  /// \sa Query
  EntityView Select(const Query& query, QueryFlags flags=0) const;

  /// \brief Add a bond
  bool AddBond(const BondHandle& bond_handle);

  /// \brief Get total number of atoms in this entity
  int GetAtomCount() const;

  /// \brief Get total number of residues in this entity
  int GetResidueCount() const;

  /// \brief Get number of chains in this entity
  int GetChainCount() const;

  /// \brief Find chain by name
  ChainView FindChain(const String& chain_name) const;

  /// \brief Find atom for given residue handle
  ///
  /// deprecated. Use ViewForHandle instead
  ResidueView FindResidue(const ResidueHandle& residue) const;

  ResidueView FindResidue(const String& chain, const ResNum& num) const;

  /// \brief locate residue view for given residue handle
  ///
  /// returns an invalid ResidueView if the residue could not be
  /// located in the structure
  ResidueView ViewForHandle(const ResidueHandle& residue) const;
  /// \brief Find view for given atom handle
  ///
  /// deprecated. Use ViewForHandle instead
  AtomView FindAtom(const AtomHandle& atom) const;
  AtomView ViewForHandle(const AtomHandle& atom) const;

  /// \brief locate atom view for given atom handle
  ///
  /// returns an invalid AtomView if the residue could not be
  /// located in the structure
  AtomView FindAtom(const String& chain_name, const ResNum& num,
                    const String& atom_name) const;
  /// \brief Find chain for given chain handle
  ///
  /// deprecated: Use ViewForHandle instead

  ChainView FindChain(const ChainHandle& chain) const;
  /// \brief locate chain view for given chain handle
  ///
  /// returns an invalid AtomView if the residue could not be
  /// located in the structure
  ChainView ViewForHandle(const ChainHandle& chain) const;
  /// \brief  Check whether the view includes the given chain
  /// \return True, if the chain is included in the view, false if not.
  bool IsChainIncluded(const ChainHandle& chain) const;

  /// \brief  remove atom from view
  void RemoveAtom(AtomView view);

  /// \brief  remove residue from view
  void RemoveResidue(ResidueView view);

  /// \brief remove chain from view
  void RemoveChain(ChainView view);

  /// \brief clear bonds and add all inclusive ones (ie where both atoms are 
  ///     present)
  void AddAllInclusiveBonds();

  const BondHandleList& GetBondList() const;
  
  /// \brief  find atoms close to point
  AtomViewList FindWithin(const geom::Vec3& center, Real radius) const;
  
  /// \brief  Get number of bonds
  int GetBondCount() const;

  const ChainViewList& GetChainList() const;

  /// \brief Get list of all residues included in the view
  ResidueViewList GetResidueList() const;

  /// \brief Get list of all atoms included in the view
  AtomViewList GetAtomList() const;
  
  /// \brief Apply entity visitor to whole view
  void Apply(EntityVisitor& visitor);

  /// \brief Apply entity view visitor to whole view
  void Apply(EntityViewVisitor& visitor);

  /// \brief Get handle this view points to
  EntityHandle GetHandle() const;

  /// \brief Return the min/max pair for the given property
  std::pair<Real,Real> GetMinMax(const String& prop, 
                                     Prop::Level hint=Prop::UNSPECIFIED) const;

  /// \brief  Get angle in radians between bonds a1-a2 and a2-a3
  Real GetAngle(const AtomHandle& a1, const AtomHandle& a2,
                const AtomHandle& a3) const;

  /// \brief  Get angle in radians between bonds a1-a2 and a2-a3
  Real GetAngle(const AtomView& a1, const AtomView& a2,
                const AtomView& a3) const;

  /// \brief  Extend current view to include all atoms of each residue where
  ///         at least one atom is selected currently
  EntityView ExtendViewToResidues() const;

  /// \brief  Extend current view to include all atoms that are within the sum
  ///         of their vdw radii + gap
  ///
  /// Include all atoms within: at1.GetRadius() + at2.GetRadius() + gap
  EntityView ExtendViewToSurrounding(Real gap) const;

  /// \brief returns a string containing a human-readable summary of the
  ///     entity view
  String Dump() const;
  /// \brief create a new full view
  /// 
  /// This is an alias of EntityView::Copy()
  EntityView CreateFullView() const;
  /// \brief create empty view
  EntityView CreateEmptyView() const;
  
  /// \brief true if view contains at least one atom
  bool HasAtoms() const;
#ifdef _MSC_VER
  // alternative routines for msvc
  AtomView AddXAtom(const AtomHandle& ah, ViewAddFlags flags=0);
  AtomView FindXAtom(const AtomHandle& ah);
#endif
  void RemoveBond(const BondHandle& bond);
protected:
  /// \internal
  void AddAtomInternal(const AtomView& av);
  void RemoveAtomInternal(const AtomView& av);

private:


  EntityViewDataPtr  data_;
};

typedef std::pair<EntityView, EntityView> EntityViewPair;

}} // ns

#endif // OST_ENTITY_VIEW_HH

