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
#ifndef OST_ATOM_IMPL_HH
#define OST_ATOM_IMPL_HH

#include <boost/enable_shared_from_this.hpp>

#include <ost/mol/module_config.hh>
#include <ost/geom/geom.hh>
#include <ost/mol/entity_visitor_fw.hh>
#include <ost/mol/impl/atom_prop.hh>

#include <ost/mol/impl/atom_impl_fw.hh>
#include <ost/mol/impl/residue_impl_fw.hh>
#include <ost/mol/impl/connector_impl_fw.hh>
#include <ost/mol/impl/fragment_impl_fw.hh>
#include <ost/mol/impl/entity_impl_fw.hh>

#include <ost/generic_property.hh>
#include <ost/mol/property_id.hh>

namespace ost { namespace mol { namespace impl {

/// \brief Atom implementation.
/// 
/// \section Connectors
/// each atom distinguishes between two types of connectors: primary and
/// secondary. each atoms has at most one primary connector: this connects the
/// atom with the atom that serves as the frame of reference for the internal
/// coordinates. the secondary connectors connect this atom with the atoms that
/// themselves depend on the coordinate system of this atom.
/// \internal
class AtomImpl: public GenericPropContainerImpl,
                public boost::enable_shared_from_this<AtomImpl> {
public:
   AtomImpl(const EntityImplPtr& ent, const ResidueImplPtr& res,
            const String& name, const geom::Vec3& pos, const String& ele,
            unsigned long index);

  ~AtomImpl();
  void Apply(EntityVisitor& h);

  // for efficiency reasons, the simple setter/getter methods are
  // replaced by direct access - this is the impl layer after all

  const String& Name() const {return name_;}
  String& Name() {return name_;}

  // DEPRECATED
  const String& GetName() const {return name_;}

  const geom::Vec3& TransformedPos() const {return tf_pos_;}
  geom::Vec3& TransformedPos() {return tf_pos_;}

  const geom::Vec3& OriginalPos() const {return pos_;}
  geom::Vec3& OriginalPos() {return pos_;}

  ResidueImplPtr GetResidue() const;

  void SetPrimaryConnector(const ConnectorImplP& bp) {
    prim_connector_=bp;
  }
  
  const ConnectorImplP& GetPrimaryConnector() const {
    return prim_connector_;
  }

  const ConnectorImplList& GetSecondaryConnectors() const {
    return connector_list_;
  }

  void AddSecondaryConnector(const ConnectorImplP& bp);
  
  // updates position and then follows secondary connectors
  void UpdateFromICS();

  // derive internal coordinates from position and linkage
  void UpdateFromXCS();
  
  String GetQualifiedName() const;
  
  EntityImplPtr GetEntity() const;
  
  int GetConnectorCount() const {
    return connector_list_.size()+(prim_connector_ ? 1 : 0);
  }

  void SetVisited(bool f) {set_state_bit(0,f);}
  bool IsVisited() const {return get_state_bit(0);}

  void SetTraced(bool f) {set_state_bit(1,f);}
  bool IsTraced() const {return get_state_bit(1);}


  Real GetBFactor() const { return b_factor_; }
  
  const String& GetElement() const { return element_; }  
  
  void SetElement(const String& ele) 
  { 
    if (element_!=ele) {
      element_=ele;
      AtomProp* old_prop = prop_;
      prop_=impl::AtomProp::GetDefaultProps(element_);
      if (old_prop && !old_prop->is_default) {
        if(old_prop->has_anisou) this->SetAnisou(old_prop->anisou);
        if(old_prop->charge != 0.0) this->SetCharge(old_prop->charge);
        delete old_prop;
      }
    }
  }
  bool HasDefaultProps() const { return prop_->is_default; }
  void SetAnisou(const geom::Mat3& anisou)
  {
    if (prop_->is_default && prop_->anisou!=anisou) {
      prop_=new AtomProp(*prop_);
      prop_->is_default=false;
      prop_->has_anisou=true;
    }

    prop_->anisou=anisou;
  }
  
  const geom::Mat3& GetAnisou() const
  {
    return prop_->anisou;
  }
  void SetBFactor(Real factor) 
  { 
    b_factor_=factor;
  }  
  
  void SetOccupancy(Real occ)
  {
    occupancy_=occ;
  }
  
  Real GetOccupancy() const 
  { 
    return occupancy_; 
  }
  
  Real GetRadius() const { return prop_->radius; }
  
  Real GetMass() const { return prop_->mass; }
  Real GetCharge() const { return prop_->charge; }
  
  bool IsHetAtom() { return is_hetatm_; }
  
  void SetHetAtom(bool het) { is_hetatm_=het; }
  
  void SetMass(Real mass)
  {
    if (prop_->is_default && prop_->mass!=mass) {
      prop_=new AtomProp(*prop_);
      prop_->is_default=false;
    }
    prop_->mass=mass;
  }
  
  void SetRadius(Real radius)
  {
    if (prop_->is_default && prop_->radius!=radius) {
      prop_=new AtomProp(*prop_);
      prop_->is_default=false;
    }
    prop_->radius=radius;
  }
  
  void SetCharge(Real charge)
  {
    if (prop_->is_default && prop_->charge!=charge) {
      prop_=new AtomProp(*prop_);
      prop_->is_default=false;
    }
    prop_->charge=charge;
  }
  
  
  unsigned int GetState() const 
  {
    return state_;
  }
  
  void SetState(int state) 
  {
    state_=state;
  }
  
  void ClearDirectionality();
  /// \brief trace directionality
  void TraceDirectionality(FragmentImplP frag, ConnectorImplP conn,  
                           int n, unsigned int& c);

  bool HasPrevious() const {return prim_connector_.get()!=NULL;}

  void DeleteAllConnectors();
  
  void DeleteConnector(const ConnectorImplP& conn, 
                       bool delete_other=true);

  void DeleteAllTorsions();
  
  String GetStringProperty(Prop::ID prop_id) const;
                                                  
  Real GetFloatProperty(Prop::ID prop_id) const;
  
  int GetIntProperty(Prop::ID prop_id) const;                     

  unsigned long GetIndex() const {return index_;}
  void SetIndex(unsigned long index) {index_=index;}
                     
private:
  ResidueImplW res_;
  String       name_;
  geom::Vec3   pos_;
  geom::Vec3   tf_pos_;
  Real         occupancy_;
  Real         b_factor_;
  AtomProp*    prop_;
  bool         is_hetatm_;
  String       element_;  
  ConnectorImplP prim_connector_;
  ConnectorImplList connector_list_;
  FragmentImplP fragment_;
  
  /// total rotation matrix. used to optimise recursive traversal. Do not rely
  /// on its value!
  geom::Mat3 total_rot_;
  
  /// state used by many algorithms that need to associate data with an
  /// atom
  unsigned int state_;

  void set_state_bit(unsigned int bit, bool state) 
  {
    unsigned int mask = 0x1<<bit;
    if(state) {
      state_ |= mask;
    } else {
      state_ &= ~mask;
    }
  }

  bool get_state_bit(unsigned int bit) const {
    unsigned int mask = 0x1<<bit;
    return (state_ & mask)!=0;
  }

  unsigned long index_;
};

/// \internal
std::ostream& operator<<(std::ostream& o, const AtomImplPtr ap);

/// \internal
bool ConnectorExists(const AtomImplPtr& a, const AtomImplPtr& b);

/// \internal
ConnectorImplP GetConnector(const AtomImplPtr& a, const AtomImplPtr& b);

}}} // ns

#endif

