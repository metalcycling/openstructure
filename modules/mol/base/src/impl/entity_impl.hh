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
#ifndef OST_ENTITY_IMPL_HH
#define OST_ENTITY_IMPL_HH

#include <map>
#include <vector>

#include <boost/enable_shared_from_this.hpp>

#include <ost/mol/module_config.hh>
#include <ost/geom/geom.hh>

#include <ost/mol/impl/entity_impl_fw.hh>

#include <ost/mol/entity_view.hh>

#include <ost/mol/residue_prop.hh>
#include <ost/mol/impl/atom_impl_fw.hh>
#include <ost/mol/impl/residue_impl_fw.hh>
#include <ost/mol/impl/chain_impl_fw.hh>
#include <ost/mol/impl/connector_impl_fw.hh>
#include <ost/mol/impl/torsion_impl_fw.hh>
#include <ost/mol/impl/fragment_impl_fw.hh>
#include <ost/mol/impl/pointer_iterator.hh>
#include <ost/mol/entity_visitor_fw.hh>
#include <ost/mol/entity_observer_fw.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/spatial_organizer.hh>


#include <ost/generic_property.hh>

namespace ost { namespace mol { namespace impl {

/// \internal
typedef std::map<AtomImpl*,AtomImplPtr> AtomImplMap;
/// \internal
typedef std::map<ResidueImpl*,ResidueImplPtr> ResidueImplMap;
/// \internal
typedef std::vector<ChainImplPtr> ChainImplList;
/// \internal
typedef std::map<ConnectorImpl*,ConnectorImplP> ConnectorImplMap;
/// \internal
typedef std::map<TorsionImpl*,TorsionImplP> TorsionImplMap;
/// \internal
typedef std::vector<FragmentImplP> FragmentImplList;
/// \internal
typedef std::map<EntityObserver*,EntityObserverPtr> EntityObserverMap;
/// \internal
typedef SpatialOrganizer<AtomImplPtr> SpatialAtomOrganizer;

/// \internal
typedef enum {
  DirtyICS=1,
  DirtyXCS=2,
  // dirty trace (implies dirty ICS)
  DirtyTrace=DirtyICS+4,
  DirtyOrganizer=8,
  DisableICS=16
} EntityDirtyFlags;


/// \internal
class EntityImpl: public GenericPropContainerImpl,
                  public boost::enable_shared_from_this<EntityImpl>                  
{
public:
  EntityImpl();
  ~EntityImpl();

  /// \brief copy entity
  EntityImplPtr Copy();
  
  Real GetMass() const;
  geom::Vec3 GetCenterOfMass() const;
  geom::Vec3 GetCenterOfAtoms() const;
  /// \brief returns the axis-aligned bounding box of the entity
  geom::AlignedCuboid GetBounds() const;
  


  // default copy ctor and assignment op should work for now
  AtomImplPtr CreateAtom(const ResidueImplPtr& rp, const String& name,
                         const geom::Vec3& pos, const String& ele);

  ResidueImplPtr CreateResidue(const ChainImplPtr& cp,
                               const ResNum& n,
                               const ResidueKey& k);

  ChainImplPtr InsertChain(const String& cname);
  /// \brief insert a new chain based on parameters of the given chain
  /// 
  /// The chain will have no residues and atoms
  // force deep to be set explicitely, because it is better than implicit
  // (and since we are on the impl level interface consistency isn't that critical)
  ChainImplPtr InsertChain(const String& chain_name, const ChainImplPtr& chain, bool deep);
  ConnectorImplP Connect(const AtomImplPtr& first, const AtomImplPtr& second,
                         Real len, Real theta, Real phi,
                         unsigned char bond_order);

  TorsionImplP AddTorsion(const String& name, const AtomImplPtr& a1,
                          const AtomImplPtr& a2, const AtomImplPtr& a3,
                          const AtomImplPtr& a4);

  void TraceDirectionality();

  //! Search for torsion.
  TorsionImplP FindTorsion(const AtomImplPtr& a1,
                           const AtomImplPtr& a2,
                           const AtomImplPtr& a3,
                           const AtomImplPtr& a4) const;

  void RenameChain(ChainImplPtr chain, const String& new_name);
  bool SetAngle(const AtomImplPtr& a1, const AtomImplPtr& a2,
                const AtomImplPtr& a3, Real angle);

  Real GetAngle(const AtomImplPtr& a1, const AtomImplPtr& a2,
                  const AtomImplPtr& a3) const;
  Real GetAngleXCS(const AtomImplPtr& a1, const AtomImplPtr& a2,
                     const AtomImplPtr& a3) const;
  Real GetAngleICS(const AtomImplPtr& a1, const AtomImplPtr& a2,
                     const AtomImplPtr& a3) const;
  // update positions from internal coordinate system
  void UpdateFromICS();

  // update internal coordinate system from positions
  void UpdateFromXCS();

  void Apply(EntityVisitor& v);
  void ApplyTransform(const geom::Transform& t);

  void SetTransform(const geom::Transform& t);
  const geom::Transform& GetTransform() const {return transform_;}
  bool HasTransform() const {return has_transform_;}
  void ClearTransform();
  void FixTransform();

  void AttachObserver(const EntityObserverPtr& o);
  void DetachObserver(const EntityObserverPtr& o);
  void NotifyObserver();

  
  void UpdateOrganizer();
  
  AtomImplList FindWithin(const geom::Vec3& pos, Real radius);
  // use query flag defaults
  EntityView Select(const EntityHandle& h, const Query& q) const;
  // override query flag defaults with given flags
  EntityView Select(const EntityHandle& h, const Query& q, 
                    QueryFlags flags) const;
  EntityView CreateFullView(const EntityHandle& h) const;
  void SetDefaultQueryFlags(QueryFlags f) {default_query_flags_=f;}
  QueryFlags GetDefaultQueryFlags() const {return default_query_flags_;}


  /// Get chain by name. Returns an invalid ChainImplPtr if no chain with the
  /// given name exists.
  ChainImplPtr FindChain(const String& name) const;

  ResidueImplPtr FindResidue(const String& chain_name,
                            const ResNum& residue) const;

  AtomImplPtr FindAtom(const String& chain_name,
                       const ResNum& residue,
                       const String& atom_name) const;
  /*
    this is a deep-swap which will exchange the internal structure
    of two implementations, which by design affects _all_ entity handles
    that point two either of these two implementations
  */
  void Swap(EntityImpl& impl);

  //! Get total number of atoms of all chains belonging to this entity
  int GetAtomCount() const;

  //! Get total number of bonds belonging to this entity
  int GetBondCount() const;

  //! Get list of bonds
  mol::BondHandleList GetBondList() const;


  //! Get total number of residues of all chains belonging to this entity
  int GetResidueCount() const;

  //! Get number of chains
  int GetChainCount() const;

  const ChainImplList& GetChainList() const { return chain_list_; }

  ChainImplList& GetChainList() { return chain_list_; }

  void DeleteFromConnMap(const ConnectorImplP& conn);

  void DeleteChain(const ChainImplPtr& chain);

  void DeleteAtom(const AtomImplPtr& atom);

  void IncXCSEditorCount();
  void DecXCSEditorCount();
    
  void IncICSEditorCount();
  void DecICSEditorCount();

  TorsionImplMap& GetTorsionMap();
  
  void UpdateXCSIfNeeded();
  
  void UpdateICSIfNeeded();
  
  void UpdateOrganizerIfNeeded();
  
  void EnableICS();
  
  bool HasICS() const;
  
  bool IsXCSDirty() const;
  
  void MarkXCSDirty();
  void MarkICSDirty();
  void MarkTraceDirty();
  void MarkOrganizerDirty();

  void UpdateTransformedPos();

  const String& GetName() const;

  impl::ChainImplList::iterator GetChain(const String& name);

  pointer_it<ChainImplPtr> GetChainIter(const String& name);
  void SetName(const String& ent_name);

  void ReorderAllResidues();

  void RenumberAllResidues(int start, bool keep_spacing);

  void RenumberChain(const String& name, int start, bool keep_spacing);

private:
  void DoCopy(EntityImplPtr dest);
  
  void ReplicateHierarchy(EntityImplPtr dest);
  
  void DoCopyBondsAndTorsions(EntityImplPtr dest);

  AtomImplMap atom_map_;
  ChainImplList chain_list_;
  ConnectorImplMap connector_map_;
  TorsionImplMap torsion_map_;

  geom::Transform transform_;
  bool has_transform_;

  SpatialAtomOrganizer atom_organizer_;
  FragmentImplList fragment_list_;
  EntityObserverMap observer_map_;

  int xcs_editor_count_;
  int ics_editor_count_;
  int dirty_flags_;
  String name_;

  unsigned long next_index_;

  QueryFlags default_query_flags_;

  template <bool always_true>
  EntityView do_selection(const EntityHandle&, const Query&, QueryFlags) const;
};

}}} // ns

#endif
