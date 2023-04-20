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
#ifndef OST_CHAIN_VIEW_HH
#define OST_CHAIN_VIEW_HH

#include <ost/mol/module_config.hh>
#include <ost/geom/geom.hh>

#include <ost/mol/query.hh>
#include <ost/mol/chain_base.hh>
#include <ost/mol/impl/chain_impl_fw.hh>
#include <ost/mol/impl/entity_impl_fw.hh>
#include <ost/mol/view_type_fw.hh>
#include <ost/mol/entity_visitor_fw.hh>
#include <ost/mol/handle_type_fw.hh>
#include <ost/mol/residue_prop.hh>

namespace ost { namespace mol {

/// \brief definition of ChainView
class DLLEXPORT_OST_MOL ChainView : public ChainBase {
public:
  ChainView();

  ChainView(ChainViewDataPtr data, 
            impl::ChainImplPtr impl);
  
  ChainView(const EntityView& entity,
            const ChainHandle& chain);
public:
  
  /// \name Handle validity
  //@{
  /// \brief check validity of handle
  /// 
  /// check, whether the chain handle points to a valid handle.
  /// \note It is an error to use any method other than #IsValid, #Impl() and 
  ///       #operator bool() when the handle is invalid. An InvalidHandle
  ///       exception will be thrown.
  operator bool() const { return this->IsValid(); }
  /// \brief check validity of handle
  /// \sa #operator bool()
  bool IsValid() const { return data_.get()!=0; }
  //@}
  
  /// \brief Get parent entity view
  EntityView GetEntity() const;
  
  void Apply(EntityVisitor& visitor);
  void Apply(EntityViewVisitor& visitor);
  
  /// \brief Get number of residues in this chain view
  int GetResidueCount() const;
  
  /// \brief Get number of atoms in this chain view
  int GetAtomCount() const;
  
  /// \brief   Get number of all inclusive bonds of this chain.
  ///
  /// \warning For a bond between two different chains, it is not clearly
  ///          defined to which chain it belongs. Therefore, this method will
  ///          only return the number of bonds where both bonded atoms are
  ///          within the same chain (all inclusive bonds).
  int GetBondCount() const;

public:
  /// \brief Add residue to view
  ///
  ///
  /// If ViewAddFlag::CHECK_DUPLICATES is set, the method will ensure that the 
  /// handle is not already included in the view.
  ///
  /// \param residue_handle The residue handle to add
  /// \param flags Can be a ORed together combination of ViewAddFlag. Flags that
  ///      do not affect the operation are silently ignored.
  /// \return the added residue view, or a reference to the existing residue 
  ///      view that points to the very same handle.
  ResidueView AddResidue(const ResidueHandle& residue_handle, 
                         ViewAddFlags flags=0);      
                         
  /// \brief Add residue to view
  ///
  /// This method is identical to AddResidue(const ResidueHandle&) except that 
  /// in case ViewAddFlag::INCLUDE_ATOMS is set only the atoms included in the
  /// residue view are considered.
  /// 
  /// If ViewAddFlag::CHECK_DUPLICATES is set, the method will ensure that the 
  /// view is not already included in the view.
  ///
  /// \param residue_view The residue view to add
  /// \param flags Can be a ORed together combination of ViewAddFlag. Flags that
  ///     do not affect the operation are silently ignored.
  /// \return The added residue view, or a reference to the existing residue 
  ///     view that points to the very same handle.                         
  ResidueView AddResidue(const ResidueView& residue_view, 
                         ViewAddFlags flags=0);                                 
  /// \brief Add atom to view
  /// 
  /// If the atoms's residue has not been added to the view it will be added.
  /// 
  /// \param atom_handle The atom to add
  /// \param flags Can be a ORed together combination of ViewAddFlag
  /// \return The added atom view, or a reference to the existing atom view 
  ///      that points to the very same handle.
  AtomView AddAtom(const AtomHandle& atom_handle, 
                   ViewAddFlags flags=0);  

  /// \brief Find residue by number                       
  ResidueView FindResidue(const ResNum& number) const;
    
  /// \brief Find view for given atom handle
  /// Deprecated: Use ViewForHandle instead
  AtomView FindAtom(const AtomHandle& atom) const;
  
  AtomView ViewForHandle(const AtomHandle& atom) const;
  AtomView FindAtom(const ResNum& num, const String& name) const;
  /// \brief Find residue by residue handle
  ResidueView ViewForHandle(const ResidueHandle& handle) const;
  
  ResidueView FindResidue(const ResidueHandle& handle) const;
  /// \brief Check whether the view includes the given residue.
  bool IsResidueIncluded(const ResidueHandle& handle) const;
  
  /// \brief Get handle this view points to
  ChainHandle GetHandle() const;

  /// \brief Get residues of this chain view
  /// 
  /// \todo Residue list may be unordered. Is it the responsibility of the 
  ///     method building the view to have them ordered or is it our job?
  const ResidueViewList& GetResidueList() const;
  
  /// \brief remove given residue from the view
  void RemoveResidue(ResidueView view);

  int GetResidueIndex(const ResNum& number) const;
  
  /// \brief get residue by index
  /// 
  /// returns an invalid residue handle if the index is out of bounds
  ResidueView GetResidueByIndex(int index) const;

  
  /// \brief Get entity's mass
  Real GetMass() const;
  
  /// \brief Get entity's center of mass (mass weighted)
  geom::Vec3 GetCenterOfMass() const;
  
  /// \brief Get entity's center of atoms (not mass weighted)
  ///
  /// Returns the center of all the atoms in an entity. This is
  /// similar to GetCenterOfMass(), but the atoms are not mass weighted
  geom::Vec3 GetCenterOfAtoms() const;
  
  /// \brief Get entity's axis aligned bounding box
  geom::AlignedCuboid GetBounds() const;
  
  /// /// \brief remove all residues from the view
  void RemoveResidues();
  /// \name Internal
  //@{
  /// \brief Get view data
  ChainViewDataPtr& ViewData() {
    return data_;
  }
  
  /// \brief Get view data
  const ChainViewDataPtr& ViewData() const {
    return data_;
  }
  /// \brief whether the residues form an ordered sequence with respect to
  ///     their residue numbers.
  bool InSequence() const;

  /// \brief return view based on a query object
  /// \sa Query
  EntityView Select(const Query& q, QueryFlags flags=0) const;

  /// \brief return view based on query String.
  /// \sa Query
  EntityView Select(const String& query_string, QueryFlags flags=0) const;

  /// \brief get unique id
  ///
  /// The unique id is the same for all ChainViews pointing to the same chain
  /// view data.
  unsigned long GetHashCode() const;
    
  bool operator==(const ChainView& rhs) const;
  bool operator!=(const ChainView& rhs) const;
  bool HasAtoms() const;
  //@}
private:
  ChainViewDataPtr  data_;
};

}} // ns
#endif // OST_CHAIN_VIEW_HH

