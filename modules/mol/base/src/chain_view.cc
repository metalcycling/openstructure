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
#include <ost/mol/residue_handle.hh>
#include <ost/mol/chain_view.hh>
#include <ost/mol/residue_view.hh>
#include <ost/mol/entity_visitor.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/atom_view.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/chain_impl.hh>
#include <ost/mol/entity_handle.hh>


namespace ost { namespace mol {

namespace {

// Ideally, the ResNumComp could be written with two call operators, one of
// the form (ResNum, ResidueView) and one (ResidueView, ResNum). Unfortunately,
// the STL on windows wants two argumentso of identical type. That's why we
// work around it by storing a ResNum internally and the ResidueView that is
// valid against it. But because this version is slightly slower, we use the
// fast version where possible.
struct ResNumComp {
  ResNumComp(const ResNum& num): number_(num) {}
  ResNumComp(): number_(1) {};
  bool operator()(const ResNum& a,
                  const ResidueView& b) const {
    return a<b.GetNumber();
  }
  bool operator()(const ResidueView& a,
                  const ResNum& b) const {
    return a.GetNumber()<b;
  }   
  bool operator()(const ResidueView& a, const ResidueView& b)
  {
    assert(!(a.IsValid() && b.IsValid()));
    return a.IsValid() ? a.GetNumber()<number_ : number_<b.GetNumber();
  }
private:
  ResNum number_;
};

}

class DLLEXPORT_OST_MOL ChainViewData {
public:
  ChainViewData(const EntityView& view)
    : entity(view.ViewData()), in_sequence(true) {  
  }
  EntityViewDataWeakPtr entity;
  ResidueViewList residues;
  std::map<unsigned long, ResidueView> handle_to_view;

  ResidueView ViewForHandle(const ResidueHandle& r)
  {
    std::map<unsigned long, ResidueView>::iterator i=handle_to_view.find(r.GetHashCode());
    if (i!=handle_to_view.end()) {
      return i->second;
    }
    return ResidueView();
  }

  bool            in_sequence;
};

ChainView::ChainView() {

}

EntityView ChainView::GetEntity() const {
  this->CheckValidity();
  if (!data_->entity.expired()) {
    return EntityView(data_->entity.lock(), Impl()->GetEntity());    
  }
  throw InvalidHandle();
}

ChainView::ChainView(const EntityView& entity,
                     const ChainHandle& chain)
 : ChainBase(chain.Impl()), data_(new ChainViewData(entity)) {
}

ChainView::ChainView(ChainViewDataPtr data, 
                     impl::ChainImplPtr impl) 
 : ChainBase(impl), data_(data) {
}

ChainHandle ChainView::GetHandle() const {
  return ChainHandle(Impl());
}

void ChainView::Apply(EntityVisitor& visitor) 
{
  this->CheckValidity();  
  ResidueViewList::iterator i;
  if (visitor.VisitChain(this->GetHandle())) {
    for (i=data_->residues.begin(); i!=data_->residues.end(); ++i) {
      (*i).Apply(visitor);
    }    
  }
}

void ChainView::Apply(EntityViewVisitor& visitor) 
{
  this->CheckValidity();  
  if (visitor.VisitChain(*this)) {
    ResidueViewList::iterator i;
    for (i=data_->residues.begin(); i!=data_->residues.end(); ++i) {
      (*i).Apply(visitor);
    }    
  }
}

int ChainView::GetResidueCount() const {
  this->CheckValidity();  
  return static_cast<int>(data_->residues.size());
}

int ChainView::GetAtomCount() const {
  this->CheckValidity();  
  int count=0;
  ResidueViewList::const_iterator it=data_->residues.begin();
  while(it!=data_->residues.end()) {
    count+=(*it).GetAtomCount();
    ++it;
  }
  return count;
}

int ChainView::GetBondCount() const {
  this->CheckValidity();
  int count=0;
  EntityView ev = this->GetEntity();
  const BondHandleList& bhl = ev.GetBondList();
  for (BondHandleList::const_iterator i=bhl.begin(); i!=bhl.end(); ++i) {
    if (i->GetFirst().GetResidue().GetChain().GetName()==this->GetName() &&
        i->GetSecond().GetResidue().GetChain().GetName()==this->GetName()) {
      count++;
    }
  }
  return count;
}

ResidueView ChainView::FindResidue(const ResNum& number) const {
  this->CheckValidity();  
  const ResidueViewList& l=data_->residues;
  ResidueViewList::const_iterator i;
  if (data_->in_sequence) {
    std::pair<ResidueViewList::const_iterator, 
              ResidueViewList::const_iterator> p;
#if defined(_MSC_VER)            
    p=std::equal_range(l.begin(), l.end(), ResidueView(), ResNumComp(number));    
#else
    p=std::equal_range(l.begin(), l.end(), number, ResNumComp());    
#endif
    if (p.first!=p.second) {
     return *p.first;
    }
    return ResidueView();
  } else {
    i=std::find_if(l.begin(), l.end(), 
                   bind(&ResidueView::GetNumber, boost::placeholders::_1)==number);
    return i==data_->residues.end() ? ResidueView() : *i;    
  }
}

ResidueView ChainView::FindResidue(const ResidueHandle& residue) const {
  LOG_WARNING("ChainView::FindResidue(handle) is deprecated. "
              "Use ChainView::ViewForHandle instead.");
  return this->ViewForHandle(residue);
}

ResidueView ChainView::ViewForHandle(const ResidueHandle& handle) const {
  this->CheckValidity();  
  return data_->ViewForHandle(handle);
}


bool ChainView::IsResidueIncluded(const ResidueHandle& handle) const {
  this->CheckValidity();  
  return this->ViewForHandle(handle).IsValid();
}

ResidueView ChainView::AddResidue(const ResidueHandle& residue_handle, 
                                  ViewAddFlags flags) {
  this->CheckValidity();                                    
  ResidueView rv;
  if ((flags & ViewAddFlag::CHECK_DUPLICATES) && 
      (rv=this->ViewForHandle(residue_handle)))
    return rv;
  rv=ResidueView(*this, residue_handle);
  if (!data_->residues.empty()) {
    if (data_->residues.back().GetNumber()>=rv.GetNumber())
      data_->in_sequence=false;
  }  
  data_->residues.push_back(rv);
  data_->handle_to_view[rv.GetHandle().GetHashCode()] = rv;
  if (flags & ViewAddFlag::INCLUDE_ATOMS) {
    const impl::AtomImplList& l=residue_handle.Impl()->GetAtomList();
    for (impl::AtomImplList::const_iterator i=l.begin(); i!=l.end(); ++i) {
      rv.AddAtom(AtomHandle(*i), flags);
    }
  }  
  return rv;
}

AtomView ChainView::FindAtom(const AtomHandle& atom) const {
  LOG_WARNING("ChainView::FindAtom(handle) is deprecated. "
              "Use ChainView::ViewForHandle instead.");
  return this->ViewForHandle(atom);
}

AtomView ChainView::ViewForHandle(const AtomHandle& atom) const {
  ResidueHandle residue=atom.GetResidue();
  if (atom.GetEntity()!=this->GetEntity().GetHandle())
    return AtomView();

  ResidueView v=this->FindResidue(residue.GetNumber());
  return v.IsValid() ? v.ViewForHandle(atom) : AtomView();
}

AtomView ChainView::AddAtom(const AtomHandle& atom_handle, 
                            ViewAddFlags flags) {
  ResidueView rv=this->AddResidue(atom_handle.GetResidue(), 
                                  ViewAddFlag::CHECK_DUPLICATES);
  return rv.AddAtom(atom_handle, flags);
}

const ResidueViewList& ChainView::GetResidueList() const {
  this->CheckValidity();
  return data_->residues;
}



void ChainView::RemoveResidue(ResidueView view) {
  this->CheckValidity();
  if (!view.IsValid())
    return;
  view.RemoveAtoms();
  ResidueViewList::iterator i=data_->residues.begin();
  for (; i!=data_->residues.end(); ++i) {
    if (*i==view) {
      break;
    }
  }
  ResidueViewList::iterator to_del = i;
  int index = view.GetIndex();
  for(; i!=data_->residues.end(); ++i) {
    int res_index = (*i).GetIndex();
    if(index < res_index){
      (*i).SetIndex(res_index-1);
    }
  }
  data_->residues.erase(to_del);
  data_->handle_to_view.erase(view.GetHandle().GetHashCode());
}

ResidueView ChainView::AddResidue(const ResidueView& residue_view, 
                                  ViewAddFlags flags) {
  this->CheckValidity();
  ResidueView rv;
  if ((flags & ViewAddFlag::CHECK_DUPLICATES) && 
      (rv=this->ViewForHandle(residue_view.GetHandle()))) {
    if (!(flags & ViewAddFlag::INCLUDE_ATOMS)) {
      return rv;
    }
  } else {
    rv=ResidueView(*this, residue_view.GetHandle());    
  }

  if (!data_->residues.empty()) {
    if (data_->residues.back().GetNumber()>=rv.GetNumber())
      data_->in_sequence=false;
  }
  data_->residues.push_back(rv);
  data_->handle_to_view[rv.GetHandle().GetHashCode()] = rv;
  if (flags & ViewAddFlag::INCLUDE_ATOMS) {
    AtomViewList::const_iterator i=residue_view.GetAtomList().begin();
    for (; i!=residue_view.GetAtomList().end(); ++i) {
      rv.AddAtom(*i);
    }
  }
  return rv;  
}
                       
void ChainView::RemoveResidues() {
  this->CheckValidity();
  std::for_each(data_->residues.begin(), data_->residues.end(),
                bind(&ResidueView::RemoveAtoms, boost::placeholders::_1));
  data_->residues.clear();
  data_->handle_to_view.clear();
}



ResidueView ChainView::GetResidueByIndex(int i) const
{
  this->CheckValidity();
  bool bounded=(i>=0 && i<static_cast<int>(data_->residues.size()));
  return bounded ? data_->residues[i] : ResidueView();
}

int ChainView::GetResidueIndex(const ResNum& number) const
{
  this->CheckValidity();
  ResidueViewList::iterator i;
  if (data_->in_sequence) {
    std::pair<ResidueViewList::iterator, ResidueViewList::iterator> p;
#if defined(_MSC_VER)  
    p=std::equal_range(data_->residues.begin(), data_->residues.end(), 
                       ResidueView(), ResNumComp(number));    
#else
    p=std::equal_range(data_->residues.begin(), data_->residues.end(), 
                       number, ResNumComp());    
#endif
    i=p.first;
  } else {
    i=std::find_if(data_->residues.begin(), data_->residues.end(),
                   bind(&ResidueView::GetNumber, boost::placeholders::_1)==number);    
  }

  return i==data_->residues.end() ? -1 : i-data_->residues.begin();
}

bool ChainView::operator==(const ChainView& rhs) const
{
  return data_==rhs.data_;
}

bool ChainView::operator!=(const ChainView& rhs) const
{
  return !this->operator==(rhs);
}

Real ChainView::GetMass() const {
  this->CheckValidity();
  double mass = 0;
  ResidueViewList::const_iterator i;
  for (i=data_->residues.begin(); i!=data_->residues.end(); ++i) {
    mass+=i->GetMass();
  }
  return mass;
}

geom::AlignedCuboid ChainView::GetBounds() const
{
  this->CheckValidity();
  geom::Vec3 mmin( std::numeric_limits<Real>::max());
  geom::Vec3 mmax(-std::numeric_limits<Real>::max());

  
  bool atoms=false;
  for (ResidueViewList::const_iterator
       i=data_->residues.begin(); i!=data_->residues.end(); ++i) {
    ResidueView r=*i;
    for (AtomViewList::const_iterator j=r.GetAtomList().begin(),
         e2=r.GetAtomList().end(); j!=e2; ++j) {
      mmin=geom::Min(mmin, j->GetPos());
      mmax=geom::Max(mmax, j->GetPos());
      atoms=true;
    }
  }
  if (atoms) {
    return geom::AlignedCuboid(mmin, mmax);
  } else {
    return geom::AlignedCuboid(geom::Vec3(), geom::Vec3());
  }
}

geom::Vec3 ChainView::GetCenterOfAtoms() const
{
  this->CheckValidity();
  geom::Vec3 center;
  if(this->HasAtoms()) {
    int atom_count = 0;
    ResidueViewList::const_iterator i;
    for (i=data_->residues.begin(); i!=data_->residues.end(); ++i) {
      ResidueView r=*i;
      for (AtomViewList::const_iterator j=r.GetAtomList().begin(),
           e2=r.GetAtomList().end(); j!=e2; ++j) {
        center+=j->GetPos();
        atom_count+=1;
      }
    }
    center/=atom_count;
  }
  return center;
}

geom::Vec3 ChainView::GetCenterOfMass() const
{
  this->CheckValidity();
  geom::Vec3 center;
  Real mass = this->GetMass();
  if(this->HasAtoms() && mass > 0) {
    ResidueViewList::const_iterator i;
    for (i=data_->residues.begin(); i!=data_->residues.end(); ++i) {
      ResidueView r=*i;
      for (AtomViewList::const_iterator j=r.GetAtomList().begin(),
          e2=r.GetAtomList().end(); j!=e2; ++j) {
        center+=j->GetPos() * j->GetMass();
      }
    }
    center/=mass;
  }
  return center;
}

AtomView ChainView::FindAtom(const ResNum& num, 
                             const String& name) const
{
  this->CheckValidity();
  ResidueView res=this->FindResidue(num);
  return res.IsValid() ? res.FindAtom(name) : AtomView();
}


bool ChainView::InSequence() const
{
  this->CheckValidity();
  return data_->in_sequence;
}

EntityView ChainView::Select(const Query& q, QueryFlags flags) const
{
  this->CheckValidity();
  if (q.GetQueryString() != "") {
    return this->GetEntity().Select(Query("cname='"+Impl()->GetName()+"' and ("+
                                        q.GetQueryString()+")"), flags);
  }
  else {
  return this->GetEntity().Select(Query("cname='"+Impl()->GetName()+"'"), flags);
}
}

EntityView ChainView::Select(const String& q, QueryFlags flags) const {
  this->CheckValidity();
  if (q != "") {
    return this->GetEntity().Select(Query("cname='"+Impl()->GetName()+"' and ("+
                                        q+")"), flags);
  }
  else return this->GetEntity().Select(Query("cname='"+Impl()->GetName()+"'"), flags);
}


bool ChainView::HasAtoms() const {
  this->CheckValidity();  
  for (ResidueViewList::const_iterator it=data_->residues.begin(), 
       e=data_->residues.end(); it!=e; ++it) {
    if ((*it).HasAtoms()) {
      return true;
    }
  }
  return false;
}

unsigned long ChainView::GetHashCode() const
{
  this->CheckValidity();
  return reinterpret_cast<unsigned long>(data_.get());
}

bool ChainView::IsValid() const
{
  return Impl().get()!=0 && Impl()->GetEntity();
}


}} // ns

