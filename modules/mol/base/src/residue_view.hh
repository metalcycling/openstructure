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
#ifndef OST_RESIDUE_VIEW_HH
#define OST_RESIDUE_VIEW_HH

#include <ost/mol/query.hh>
#include <ost/mol/module_config.hh>
#include <ost/geom/geom.hh>
#include <ost/mol/residue_base.hh>
#include <ost/mol/view_type_fw.hh>
#include <ost/mol/impl/residue_impl_fw.hh>
#include <ost/mol/entity_visitor_fw.hh>
#include <ost/mol/handle_type_fw.hh>

namespace ost { namespace mol {

/// \brief residue view
///
/// Atoms are added with ResidueView::AddAtom() and removed with 
/// ResidueView::RemoveAtom(). Like other view operations, these operations only
/// affect the view and do not alter the structure and topology of the underlying
/// \ref ResidueHandle "residue".
class DLLEXPORT_OST_MOL ResidueView : public ResidueBase {

  friend class ChainView;

public:
  /// \brief Create invalid ResidueView
  ResidueView();

  /// \brief Create new residue view
  ///
  /// \param chain is the parent chain view to insert into
  /// \param residue is residue handle this view points to.
  /// \note For better readability consider using ChainView::AddResidue()
  ResidueView(const ChainView&      chain,
              const ResidueHandle&  residue);

public:
  
  /// \name Handle validity
  //@{
  /// \brief check validity of view
  ///
  /// check, whether the residue view points to a valid residue.
  /// \note It is an error to use any method other than #IsValid, #Impl and
  ///       #operator bool() when the handle is invalid. An InvalidHandle
  ///       exception will be thrown.
  operator bool() const { return this->IsValid(); }
  /// \brief check validity of handle
  /// \sa #operator bool()
  bool IsValid() const;
  //@}
  
  /// \name internal
  //@{
  /// \brief construct new view from view data and impl.
  ///
  /// This has the same effect as using the default copy constructor of the view.
  ResidueView(const ResidueViewDataPtr& data,
              const impl::ResidueImplPtr& impl);

  /// \brief Get internal view data
  const ResidueViewDataPtr& ViewData() const {
    return data_;
  }
  /// \brief Get internal view data
  ResidueViewDataPtr& ViewData() {
    return data_;
  }
  //@}
public:
  /// \name Adding and removing atoms
  //@{
  /// \brief Add atom to view
  ///
  /// If ViewAddFlag::CHECK_DUPLICATES is set, the method will ensure that the
  /// view is not already included in this view. This method is an alias for
  /// `residue_view.AddAtom(atom_view.GetHandle())`
  ///
  /// \param atom_view is the atom view to add
  /// \param flags can be a ORed together combination of ViewAddFlag
  /// \return  The added atom view, or a reference to the existing atom
  ///     view that points to the very same view.
  AtomView AddAtom(const AtomView& atom_view,
                   ViewAddFlags flags=0);

  /// \brief Add atom to view
  ///
  /// If ViewAddFlag::CHECK_DUPLICATES is set, the method will ensure that the
  /// handle is not already included in the view.
  ///
  /// \param atom_handle the atom handle to add
  /// \param flags can be a ORed together combination of ViewAddFlag
  /// \return  The added atom view, or a reference to the existing atom
  ///      view that points to the very same handle.
  AtomView AddAtom(const AtomHandle& atom_handle,
                   ViewAddFlags flags=0);
  //@}

  /// \brief Find atom by name
  AtomView FindAtom(const String& atom_name) const;

  /// \brief Find residue by residue handle
  AtomView ViewForHandle(const AtomHandle& handle) const;
  AtomView FindAtom(const AtomHandle& handle) const;

  /// \brief Check whether the view includes the the given atom.
  bool IsAtomIncluded(const AtomHandle& handle) const;

  /// \brief remove given atom from view
  ///
  /// all bonds involving this atom will be removed from the view as well.
  void RemoveAtom(AtomView view);

  /// \brief remove all atoms
  ///
  /// all bonds involving one of the atoms will be removed from the view as
  /// well.
  void RemoveAtoms();

  /// \brief Apply entity visitor to whole chain
  void Apply(EntityVisitor& visitor);
  void Apply(EntityViewVisitor& visitor);

  /// \brief     return number of atoms in this residue view.
  int GetAtomCount() const;

  /// \brief     get list of atoms in this view
  const AtomViewList& GetAtomList() const;

  /// \brief get handle this view points to
  ResidueHandle GetHandle() const;

  /// \brief get entity
  EntityView GetEntity() const;
  
  /// \brief get index of residue view in chain
  int GetIndex() const;
  /// \brief    get parent chain view.
  ChainView GetChain() const;
  
  /// \brief Get residue's mass
  double GetMass() const;
  
  /// \brief Get residue's center of mass (mass weighted)
  geom::Vec3 GetCenterOfMass() const;
  
  /// \brief Get residue's center of atoms (not mass weighted)
  ///
  /// Returns the center of all the atoms in this residue. This is
  /// similar to GetCenterOfMass(), but the atoms are not mass weighted
  geom::Vec3 GetCenterOfAtoms() const;
  
  /// \brief Get residue's axis-aligned bounding box
  geom::AlignedCuboid GetBounds() const;

  /// \brief return view based on a query object
  /// \sa Query
  EntityView Select(const Query& q, QueryFlags flags=0) const;

  /// \brief return view based on query String.
  /// \sa Query
  EntityView Select(const String& query_string, QueryFlags flags=0) const;

  /// \brief get unique id
  ///
  /// The unique id is the same for all ResidueViews pointing to the same residue
  /// view data.
  unsigned long GetHashCode() const;
  
  bool operator==(const ResidueView& rhs) const;
  
  bool operator!=(const ResidueView& rhs) const;

  bool HasAtoms() const;
protected:
  /// \brief set the index of residiue view in chain
  /// should be called from chainview whenever indexes change
  void SetIndex(int index);

private:
  ResidueViewDataPtr  data_;
};

}} // ns

#endif // OST_RESIDUE_VIEW_HH

