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
#include <ost/log.hh>
#include "atom_impl.hh"
#include "chain_impl.hh"
#include "residue_impl.hh"
#include "entity_impl.hh"
#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/entity_visitor.hh>
#include <boost/bind/bind.hpp>
#include <algorithm>
#include <numeric>

using std::accumulate;

namespace ost { namespace mol { namespace impl {

ChainImpl::ChainImpl(const EntityImplPtr& e, const String& name):
  ent_(e), 
  name_(name),
  residue_list_(),
  in_sequence_(true),
  type_(CHAINTYPE_UNKNOWN),
  description_()
{}

String ChainImpl::GetName() const
{
  return name_;
}

struct ResNumComp {
  bool operator()(const ResNum& a, 
                  const ResidueImplPtr& b) {
    return a<b->GetNumber();
  }
  bool operator()(const ResidueImplPtr& a, 
                  const ResNum& b) {
    return a->GetNumber()<b;
  }  
};

bool ChainImpl::InSequence() const {
  return in_sequence_;
}

ResidueImplPtr ChainImpl::AppendResidue(const ResidueImplPtr& res, bool deep)
{
  ResidueImplPtr dst_res=this->AppendResidue(res->GetKey(), res->GetNumber());
  dst_res->Assign(*res.get());                                                
  dst_res->SetOneLetterCode(res->GetOneLetterCode());
  dst_res->SetSecStructure(res->GetSecStructure());
  dst_res->SetChemClass(res->GetChemClass());  
  dst_res->SetProtein(res->IsProtein());
  dst_res->SetChemType(res->GetChemType());
  dst_res->SetIsLigand(res->IsLigand());

  if(deep){
    AtomImplList::iterator it=res->GetAtomList().begin(),
                           it_end=res->GetAtomList().end();

    if(res->HasAltAtoms()){

      for(;it!=it_end;++it){

        AtomImplPtr atom=*it;
        std::vector<String> alt_gn = res->GetAltAtomGroupNames(atom);

        if(alt_gn.empty()){
          //there are no alternative locations for this atoms =>just add it
          dst_res->InsertAtom(atom);
        }

        else{
          AtomImplPtr dst_atom;
          
          dst_atom = dst_res->InsertAltAtom(atom->GetName(),alt_gn[0],
                                            res->GetAltAtomPos(atom,alt_gn[0]),
                                            atom->GetElement(),
                                            res->GetAltAtomOcc(atom,alt_gn[0]),
                                            res->GetAltAtomBFactor(atom,alt_gn[0]));
          
          //in the InsertAltAtom function only name, position and 
          //element get passed when inserting the actual atom into the residue. 
          //All the rest we have to map by ourself as we would call the 
          //dst_res->InsertAtom(atom) function
          dst_atom->Assign(*atom.get());
          dst_atom->SetState(atom->GetState());
          dst_atom->SetBFactor(atom->GetBFactor());
          dst_atom->SetOccupancy(atom->GetOccupancy());
          dst_atom->SetHetAtom(atom->IsHetAtom());

          if (!atom->HasDefaultProps()) {
            dst_atom->SetRadius(atom->GetRadius());
            dst_atom->SetCharge(atom->GetCharge());
            dst_atom->SetMass(atom->GetMass());
            dst_atom->SetAnisou(atom->GetAnisou());
          }

          //the pos for the first alt group name is already added in the
          //InsertAltAtom function, all the others get added here
          for(std::vector<String>::iterator alt_gn_it = alt_gn.begin() + 1,
              alt_gn_it_e = alt_gn.end(); alt_gn_it != alt_gn_it_e; ++alt_gn_it){

            dst_res->AddAltAtomPos(*alt_gn_it,dst_atom,
                                   res->GetAltAtomPos(atom,*alt_gn_it),
                                   res->GetAltAtomOcc(atom,*alt_gn_it),
                                   res->GetAltAtomBFactor(atom,*alt_gn_it));
          }
        }
      }
    }
    else{
      for(;it!=it_end;++it){
        dst_res->InsertAtom(*it);
      }
    }
  }
  return dst_res;
}

ResidueImplPtr ChainImpl::InsertResidueBefore(int index, 
                                              const ResNum& num2, 
                                              const ResidueKey& key) {
  ResNum num=num2;
  if (index>=static_cast<int>(residue_list_.size()) || index<0) {
    return ResidueImplPtr();
  } 
  // check if we break the InSequence property
  if (in_sequence_) {
    ResNum next=residue_list_[index]->GetNumber();
    if (next<=num2 || (index>0 && residue_list_[index-1]->GetNumber()>=num))
      in_sequence_=false;
  }
  ResidueImplPtr rp = ent_.lock()->CreateResidue(shared_from_this(), num, key);  
  residue_list_.insert(residue_list_.begin()+index, rp);
  if (in_sequence_) {
    this->UpdateShifts();
  }
  return rp;
}

ResidueImplPtr ChainImpl::InsertResidueAfter(int index, 
                                             const ResNum& num2, 
                                             const ResidueKey& key) {
  ResNum num=num2;
  if (index>=static_cast<int>(residue_list_.size()) || index<0) {
    return ResidueImplPtr();
  } 
  // check if we break the InSequence property
  if (in_sequence_) {
    ResNum prev=residue_list_[index]->GetNumber();
    if (prev>=num || (index<static_cast<int>(residue_list_.size())-1 && 
                      residue_list_[index+1]->GetNumber()<=num)) {
      in_sequence_=false;
    }
  }
  ResidueImplPtr rp = ent_.lock()->CreateResidue(shared_from_this(), num, key);
  residue_list_.insert(residue_list_.begin()+index+1, rp);
  if (in_sequence_) {
    this->UpdateShifts();
  }
  return rp;
}

void ChainImpl::UpdateShifts() 
{
  shifts_.clear();
  ResidueImplList::iterator i=residue_list_.begin(); 
  ResNum last(0);
  for (int count=0; i!=residue_list_.end(); ++i, ++count) {
    const ResNum& next=(*i)->GetNumber();
    if (last.GetNum()==next.GetNum()) {
      // introduce negative shift, due to inserted residue. We should do 
      // something more clever here, since it is very likely that there are
      // more than one consecutively inserted residues.
      Shift s;
      s.start=count;
      s.shift=-1;
      last=next;      
      shifts_.push_back(s);      
      continue;
    } else if (last.GetNum()+1==next.GetNum()) {
      last=next;
      continue;
    }
    Shift s;
    s.start=count;
    s.shift=next.GetNum()-last.GetNum()-1;
    last=next;
    shifts_.push_back(s);
  }  
}

void ChainImpl::DeleteAllResidues() {
  ResidueImplList::iterator i=residue_list_.begin();
  for (; i!=residue_list_.end(); ++i) {
    (*i)->DeleteAllAtoms();
  }
  residue_list_.clear();
  this->UpdateShifts();
}

void ChainImpl::DeleteResidue(const ResNum& number) {
  int index=this->GetIndexForResNum(number);
  if (index>=0) {
    ResidueImplPtr r=residue_list_[index];
    r->DeleteAllAtoms();
    residue_list_.erase(residue_list_.begin()+index);
    this->UpdateShifts();    
  }
}

void ChainImpl::DeleteResidue(const ResidueImplPtr& residue) {
  if (residue->GetChain().get()!=this)
    return;
  int index=this->GetIndex(residue);
  if (index>=0) {
    ResidueImplPtr r=residue_list_[index];
    r->DeleteAllAtoms();
    residue_list_.erase(residue_list_.begin()+index);
    this->UpdateShifts();    
  }
}

ResidueImplPtr ChainImpl::AppendResidue(const ResidueKey& key, 
                                        const ResNum& num)
{
  ResidueImplPtr rp = ent_.lock()->CreateResidue(shared_from_this(), num, key);  
  if (residue_list_.empty()) {
    if (num!=ResNum(1)) {
    // insert shift here...      
      Shift s;
      s.start=std::min(rp->GetNumber().GetNum()-1, 0);
      s.shift=rp->GetNumber().GetNum()-1;
      shifts_.push_back(s);    
    }
  } else {
    if (residue_list_.back()->GetNumber()>=num) {
      in_sequence_=false;
    }
    if (in_sequence_) {
      LOG_DEBUG("appending residue " << num);    
      // update is only needed if we introduce a new gap.
      if (residue_list_.back()->GetNumber().GetNum()+1<num.GetNum()) {
        Shift s;
        s.start=residue_list_.size();
        s.shift=num.GetNum()-residue_list_.back()->GetNumber().GetNum()-1;
        shifts_.push_back(s);
      } else if (residue_list_.back()->GetNumber().GetNum()==num.GetNum()) {
        Shift s;
        s.start=residue_list_.size();
        s.shift=-1;
        shifts_.push_back(s);
      }      
    }
  }
  residue_list_.push_back(rp);
  return rp;
}

ResidueImplPtr ChainImpl::AppendResidue(const ResidueKey& key)
{
  ResNum num(1);
  if(!residue_list_.empty()) {
    num=residue_list_.back()->GetNumber()+1;
  }
  ResidueImplPtr rp = ent_.lock()->CreateResidue(shared_from_this(),num,key);
  residue_list_.push_back(rp);
  return rp;
}

ResidueImplPtr ChainImpl::GetPrev(const ResidueImplPtr& r) const
{
  if (!r)
    return ResidueImplPtr();
  int index=this->GetIndex(r)-1;
  if (index>-1 && index<static_cast<int>(residue_list_.size())-1) {
    return residue_list_[index];
  }
  return ResidueImplPtr();
}

int ChainImpl::GetResidueCount() const {
  return static_cast<int>(residue_list_.size());
}

int ChainImpl::GetAtomCount() const {
  int count=0;
  ResidueImplList::const_iterator it=residue_list_.begin();
  while(it!=residue_list_.end()) {
    count+=(*it)->GetAtomCount();
    ++it;
  }
  return count;
}

int ChainImpl::GetBondCount() const {
  int count=0;
  ResidueImplList::const_iterator it=residue_list_.begin();
  while(it!=residue_list_.end()) {
    count+=(*it)->GetBondCount();
    ++it;
  }
  return count;
}


ResidueImplPtr ChainImpl::GetNext(const ResidueImplPtr& r) const
{
  if (!r)
    return ResidueImplPtr();
  int index=this->GetIndex(r)+1;
  if (index>0 && index<=static_cast<int>(residue_list_.size())-1) {
    return residue_list_[index];
  }
  return ResidueImplPtr();
}


const ResidueImplList& ChainImpl::GetResidueList() const
{
  return residue_list_;
}

void ChainImpl::Apply(EntityVisitor& v)
{
  LOG_TRACE("visitor @" << &v << " visiting chain impl @" << this);
  if (v.VisitChain(ChainHandle(shared_from_this()))) {
    for (ResidueImplList::iterator it=residue_list_.begin();
         it!=residue_list_.end();++it) {
      (*it)->Apply(v);
    }
  }
}

ResidueImplPtr ChainImpl::FindResidue(const ResNum& number) const 
{
  int index=this->GetIndexForResNum(number);
  bool invalid=index<0 || index>static_cast<int>(residue_list_.size())-1;
  return   invalid ? ResidueImplPtr() : residue_list_[index];
}

AtomImplPtr ChainImpl::FindAtom(const ResNum& num, 
                                const String& atom_name) const 
{
  ResidueImplPtr r=this->FindResidue(num);                      
  if (r) {
    return r->FindAtom(atom_name);
  }
  return AtomImplPtr();
}

EntityImplPtr ChainImpl::GetEntity() const 
{
  return ent_.lock();
}
int ChainImpl::GetIndexForResNum(const ResNum& number) const
{
  if (in_sequence_) {
    return this->GetIndexForResNumInSequence(number);
  } else {
      ResidueImplList::const_iterator k;    
      k=std::find_if(residue_list_.begin(), 
                     residue_list_.end(), 
                     bind(&ResidueImpl::GetNumber, boost::placeholders::_1)==number);

      if (k==residue_list_.end())
        return -1;
      int pos=std::distance(residue_list_.begin(), k);
      assert(residue_list_[pos]->GetNumber()==number);
      return pos;
  }
}

int ChainImpl::GetIndexForResNumInSequence(const ResNum& number) const
{
  int pos=number.GetNum()-1;    
  std::list<Shift>::const_iterator i;
  for (i=shifts_.begin(); i!=shifts_.end(); ++i) {
    const Shift& s=*i;
    if (pos<s.start) {
      break;
    } else if (pos<s.start+s.shift) {
      return -1;
    }
    pos-=s.shift;
  }
  while (pos>=0 && pos<static_cast<int>(residue_list_.size()) && 
         residue_list_[pos]->GetNumber()<number) {
    pos++;
  }
  if (pos<0 || pos>=static_cast<int>(residue_list_.size())) {
    return -1;
  }
  assert(residue_list_[pos]->GetNumber()==number);
  return pos;
}

int ChainImpl::GetIndex(const ResidueImplPtr& res) const 
{
  if (!res) {
    return -1;
  }
  ResNum number=res->GetNumber();
  if (in_sequence_) {
    return this->GetIndexForResNumInSequence(number);
  } else {
      ResidueImplList::const_iterator k=residue_list_.begin()-1;
      do {
        k=std::find_if(k+1, residue_list_.end(), 
                       bind(&ResidueImpl::GetNumber, boost::placeholders::_1)==number);
      } while(k!=residue_list_.end() && (*k)!=res);

      if (k==residue_list_.end())
        return -1;
      int pos=std::distance(residue_list_.begin(), k);
      assert(residue_list_[pos]->GetNumber()==number);
      return pos;
  }
}

void ChainImpl::AssignSecondaryStructure(SecStructure ss, 
                                         const ResNum& start, 
                                         const ResNum& end)
{
  int start_index=this->GetIndexForResNum(start);
  int i=start_index;
  bool found_end=false;
  if (i>=0) {
    while (i<static_cast<int>(residue_list_.size())) {
      if (residue_list_[i]->GetNumber()==end) {
        found_end=true;
        break;
      }
      ++i;
    }
  }
  if (!found_end) {
    return;
  }
  for (int end=i, i=start_index; i<=end; ++i) {
    residue_list_[i]->SetSecStructure(ss);
  }
}

Real ChainImpl::GetMass() const
{
  double mass = 0;
  for (ResidueImplList::const_iterator i=residue_list_.begin(); 
        i!=residue_list_.end(); ++i) {
    ResidueImplPtr r=*i;
    mass+=r->GetMass();
  }
  return mass;
}

geom::AlignedCuboid ChainImpl::GetBounds() const 
{
  geom::Vec3 mmin( std::numeric_limits<Real>::infinity());
  geom::Vec3 mmax(-std::numeric_limits<Real>::infinity());
  bool atoms=false;
  for (ResidueImplList::const_iterator i=residue_list_.begin(); 
        i!=residue_list_.end(); ++i) {
    ResidueImplPtr r=*i;
    for (AtomImplList::iterator j=r->GetAtomList().begin(); 
          j!=r->GetAtomList().end(); ++j) {
      mmin=geom::Min(mmin, (*j)->TransformedPos());
      mmax=geom::Max(mmax, (*j)->TransformedPos());
      atoms=true;
    }
  }
  if (!atoms) {
    return geom::AlignedCuboid(geom::Vec3(), geom::Vec3());
  }
  return geom::AlignedCuboid(mmin, mmax);
}


geom::Vec3 ChainImpl::GetCenterOfAtoms() const
{
  geom::Vec3 sum;
  if(this->GetAtomCount() > 0) {
    for (ResidueImplList::const_iterator i=residue_list_.begin(); 
        i!=residue_list_.end(); ++i) {
      ResidueImplPtr r=*i;
      for (AtomImplList::iterator j=r->GetAtomList().begin(); 
          j!=r->GetAtomList().end(); ++j) {
          sum+=(*j)->TransformedPos();
      }
    }
    sum/=this->GetAtomCount();
  }
  return sum;  
}

geom::Vec3 ChainImpl::GetCenterOfMass() const
{
  geom::Vec3 center;
  Real mass = this->GetMass();
  if(this->GetAtomCount() > 0 && mass > 0) {
    for (ResidueImplList::const_iterator i=residue_list_.begin(); 
         i!=residue_list_.end(); ++i) {
      ResidueImplPtr r=*i;
      for (AtomImplList::iterator j=r->GetAtomList().begin(); 
           j!=r->GetAtomList().end(); ++j) {
        center+=(*j)->TransformedPos() * (*j)->GetMass();
      }
    }
    center/=mass;
  }
  return center;  
}

void ChainImpl::SetName(const String& new_name)
{
  name_=new_name;
}

namespace {
  bool rnum_cmp(const ResidueImplPtr& r1, const ResidueImplPtr& r2)
  {
    return r1->GetNumber()<r2->GetNumber();
  }
}

void ChainImpl::ReorderResidues()
{
  std::sort(residue_list_.begin(),residue_list_.end(),rnum_cmp);
  UpdateShifts();
}

void ChainImpl::RenumberAllResidues(int start, bool keep_spacing)
{
  ResNum actual_num=ResNum(start);
  ResNum original_start_num=residue_list_[0]->GetNumber();
  ResNum start_diff=ResNum(start)-original_start_num;


  for (ResidueImplList::const_iterator i=residue_list_.begin(); 
         i!=residue_list_.end(); ++i) {

      if(keep_spacing){
         ResNum temp=(*i)->GetNumber();
         (*i)->SetNumber(temp+start_diff);
        continue;
      }

      (*i)->SetNumber(actual_num);
      actual_num++;
  }
  UpdateShifts();
}

void ChainImpl::RenumberAllResidues(const ResNumList& new_numbers)
{
  if (new_numbers.size() != residue_list_.size()) {
    throw Error("number of residues and residue numbers must match");
  }
  for (size_t i = 0; i<new_numbers.size(); ++i) {
      residue_list_[i]->SetNumber(new_numbers[i]);
  }
  this->UpdateShifts();
}

void ChainImpl::SetInSequence(const int index)
{
  ResNum num=residue_list_[index]->GetNumber();
  //Check if rp is in sequence
  if (in_sequence_) {
    if (index>0 && residue_list_[index-1]->GetNumber()>=num)
      in_sequence_=false;
    if (index<static_cast<int>(residue_list_.size())-1 && residue_list_[index+1]->GetNumber()<=num)
      in_sequence_=false;
  }
  if (in_sequence_) {
    this->UpdateShifts();
  }
}

void ChainImpl::UpdateTransformedPos()
{
  for (ResidueImplList::iterator rit=residue_list_.begin(); rit!=residue_list_.end(); ++rit) {
    (*rit)->UpdateTransformedPos();
  }
}
  
}}} // ns
