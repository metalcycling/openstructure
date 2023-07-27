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
#include <algorithm>
#include <limits>

#include <boost/bind/bind.hpp>
#include <ost/log.hh>
#include <ost/mol/bond_handle.hh>
#include <ost/mol/chain_view.hh>
#include <ost/mol/atom_view.hh>
#include <ost/mol/entity_visitor.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/chain_impl.hh>
#include <ost/mol/entity_view.hh>

#include "residue_view.hh"
#include "residue_handle.hh"


namespace ost { namespace mol {

class DLLEXPORT_OST_MOL ResidueViewData {
public:
  ResidueViewData(int residue_index): index(residue_index){}
  ResidueViewData(const ChainViewDataPtr& chain_data, int residue_index)
    : chain(chain_data), index(residue_index) { }
  ChainViewDataWeakPtr  chain;
  AtomViewList          atoms;
  int                   index;
};

ResidueView::ResidueView(const ChainView&      chain,
                         const ResidueHandle&  residue)
  : ResidueBase(residue.Impl()), data_(new ResidueViewData(chain.ViewData(), chain.GetResidueCount())) {
}

ResidueView::ResidueView(const ResidueViewDataPtr& data,
                         const impl::ResidueImplPtr& impl)
 : ResidueBase(impl), data_(data) {
}

ResidueView::ResidueView() {
}

ResidueHandle ResidueView::GetHandle() const {
  return ResidueHandle(Impl());
}

void ResidueView::Apply(EntityVisitor& visitor)
{
  this->CheckValidity();
  if (visitor.VisitResidue(this->GetHandle())) {
    AtomViewList::iterator i;
    for (i=data_->atoms.begin(); i!=data_->atoms.end(); ++i) {
      (*i).Apply(visitor);
    }
  }
}

void ResidueView::Apply(EntityViewVisitor& visitor)
{
  this->CheckValidity();
  if (visitor.VisitResidue(*this)) {
    AtomViewList::iterator i;
    for (i=data_->atoms.begin(); i!=data_->atoms.end(); ++i) {
      (*i).Apply(visitor);
    }
  }
}


int ResidueView::GetAtomCount() const {
  this->CheckValidity();
  return static_cast<int>(data_->atoms.size());
}

AtomView ResidueView::AddAtom(const AtomHandle& atom_handle,
                              ViewAddFlags flags) {
  this->CheckValidity();
  AtomView v;
  if ((flags & ViewAddFlag::CHECK_DUPLICATES) &&
      (v=this->ViewForHandle(atom_handle))) {
    return v;
  }
  v=AtomView(*this, atom_handle);
  data_->atoms.push_back(v);
  return v;
}

AtomView ResidueView::AddAtom(const AtomView& atom_view,
                              ViewAddFlags flags) {
  return this->AddAtom(atom_view.GetHandle(), flags);
}


AtomView ResidueView::FindAtom(const String& atom_name) const {
  this->CheckValidity();
  const AtomViewList& l=data_->atoms;
  AtomViewList::const_iterator i;
  i=std::find_if(l.begin(), l.end(), boost::bind(&AtomView::GetName, boost::placeholders::_1)==atom_name);
  return i==data_->atoms.end() ? AtomView() : *i;
}

AtomView ResidueView::ViewForHandle(const AtomHandle& handle) const {
  this->CheckValidity();
  const AtomViewList& l=data_->atoms;
  AtomViewList::const_iterator i;
  i=std::find_if(l.begin(), l.end(), bind(&AtomView::GetHandle, boost::placeholders::_1)==handle);
  return i==data_->atoms.end() ? AtomView() : *i;
}

bool ResidueView::IsAtomIncluded(const AtomHandle& handle) const {
  // validity is checked by FindAtom
  return this->ViewForHandle(handle).IsValid();
}


const AtomViewList& ResidueView::GetAtomList() const {
  this->CheckValidity();
  return data_->atoms;
}

ChainView ResidueView::GetChain() const {
  this->CheckValidity();
  if (!data_->chain.expired()) {
    return ChainView(data_->chain.lock(), Impl()->GetChain());
  }
  throw InvalidHandle();
}


void ResidueView::RemoveAtom(AtomView view) {
  this->CheckValidity();
  if (!view.IsValid())
    return;
  view.RemoveBonds();
  data_->atoms.erase(std::find(data_->atoms.begin(), 
                               data_->atoms.end(), view));
  this->GetEntity().RemoveAtomInternal(view);
}

EntityView ResidueView::GetEntity() const 
{
  this->CheckValidity();
  return this->GetChain().GetEntity();
}

void ResidueView::RemoveAtoms() 
{
  this->CheckValidity();
  std::for_each(data_->atoms.begin(), data_->atoms.end(),
                bind(&AtomView::RemoveBonds, boost::placeholders::_1));
  AtomViewList::iterator i;
  EntityView ev=this->GetEntity();
  for (i=data_->atoms.begin(); i!=data_->atoms.end(); ++i) {
    ev.RemoveAtomInternal(*i);
  }
  data_->atoms.clear();
}

bool ResidueView::operator==(const ResidueView& rhs) const 
{
  return rhs.Impl()==Impl() && (rhs.data_==data_);
}

bool ResidueView::operator!=(const ResidueView& rhs) const 
{
  return !this->operator==(rhs);
}


int ResidueView::GetIndex() const
{
  this->CheckValidity();
  return data_->index;
}

void ResidueView::SetIndex(int index)
{
  this->CheckValidity();
  data_->index = index;
}

double ResidueView::GetMass() const
{
  this->CheckValidity();
  double mass = 0;
  AtomViewList::const_iterator i;
  for (i=data_->atoms.begin(); i!=data_->atoms.end(); ++i) {
    mass+=(*i).GetMass();
  }
  return mass;
}

geom::AlignedCuboid ResidueView::GetBounds() const 
{
  this->CheckValidity();
  geom::Vec3 mmin( std::numeric_limits<Real>::max());
  geom::Vec3 mmax(-std::numeric_limits<Real>::max());  
  AtomViewList::const_iterator i;
  for (i=data_->atoms.begin(); i!=data_->atoms.end(); ++i) {
    mmin=geom::Min(mmin, (*i).GetPos());
    mmax=geom::Max(mmax, (*i).GetPos());
  }
  return geom::AlignedCuboid(mmin, mmax);
}

geom::Vec3 ResidueView::GetCenterOfMass() const 
{
  this->CheckValidity();  
  geom::Vec3 center;
  double mass = this->GetMass();
  if (!data_->atoms.empty() && mass > 0) {
    AtomViewList::const_iterator i;
    for (i=data_->atoms.begin(); i!=data_->atoms.end(); ++i) {
      center+=(*i).GetPos()*(*i).GetMass();
    }
      center/=mass;
  }
  return center;
}

geom::Vec3 ResidueView::GetCenterOfAtoms() const 
{
  this->CheckValidity();  
  geom::Vec3 sum;
  if (!data_->atoms.empty()) {
    AtomViewList::const_iterator i;
    for (i=data_->atoms.begin(); i!=data_->atoms.end(); ++i) {
      sum+=(*i).GetPos();
    }
    sum/=data_->atoms.size();
  }
  return sum;
}

EntityView ResidueView::Select(const Query& q, QueryFlags flags) const
{
  this->CheckValidity();
  std::stringstream s;
  if (q.GetQueryString() != "") {
     s << "cname='" << Impl()->GetChain()->GetName() << "' and rindex=" << Impl()->GetIndex()
      << " and (" << q.GetQueryString() << ")";
  }
  else s << "cname='" << Impl()->GetChain()->GetName() << "' and rindex=" << Impl()->GetIndex();
  return this->GetEntity().Select(Query(s.str()), flags);
}

EntityView ResidueView::Select(const String& q, QueryFlags flags) const {
  this->CheckValidity();
  std::stringstream s;
  if (q != "") {
     s << "cname='" << Impl()->GetChain()->GetName() << "' and rindex=" << Impl()->GetIndex()
       << " and (" << q << ")";
  }
  else s << "cname='" << Impl()->GetChain()->GetName() << "' and rindex=" << Impl()->GetIndex();
  return this->GetEntity().Select(Query(s.str()), flags);
}

AtomView ResidueView::FindAtom(const AtomHandle& handle) const
{
  LOG_WARNING("ResidueView::FindAtom(handle) is deprecated. "
              "Use ResidueView::ViewForHandle instead.");
  return this->ViewForHandle(handle);
}

bool ResidueView::HasAtoms() const {
  this->CheckValidity();
  return data_->atoms.size()>0;
}

unsigned long ResidueView::GetHashCode() const
{
  this->CheckValidity();
  return reinterpret_cast<unsigned long>(data_.get());
}

bool ResidueView::IsValid() const
{
  return Impl().get()!=0 && Impl()->GetEntity();
}

}} //ns
