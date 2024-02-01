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

/*
  Author: Marco Biasini
 */

#include <ost/invalid_handle.hh>
#include <ost/seq/sequence_handle.hh>
#include <ost/seq/impl/sequence_impl.hh>


namespace ost { namespace seq {

ConstSequenceHandle::ConstSequenceHandle()
{ }

SequenceHandle::SequenceHandle(const impl::SequenceImplPtr& impl):
  impl_(impl)
{ }

ConstSequenceHandle::ConstSequenceHandle(const impl::SequenceImplPtr& impl):
  impl_(impl)
{
  
}

bool ConstSequenceHandle::operator==(const ConstSequenceHandle& rhs) const
{
  return impl_==rhs.impl_;
}

bool ConstSequenceHandle::operator!=(const ConstSequenceHandle& rhs) const
{
  return impl_!=rhs.impl_;
}

char ConstSequenceHandle::operator[](int index) const
{
  this->CheckValidity();
  return this->GetOneLetterCode(index);
}



const String& ConstSequenceHandle::GetRole() const
{
  this->CheckValidity();
  return Impl()->GetRole();
}


const String& SequenceHandle::GetRole() const
{
  this->CheckValidity();
  return Impl()->GetRole();
}
  
void SequenceHandle::SetRole(const String& role) const
{
  this->CheckValidity();
  Impl()->SetRole(role);
}


void ConstSequenceHandle::CheckValidity() const
{
  if (!impl_) {
    throw InvalidHandle();
  }
}

bool ConstSequenceHandle::IsValid() const
{
  return impl_.get()!=0;
}

impl::SequenceImplPtr& ConstSequenceHandle::Impl() const 
{
  return impl_;
}


SequenceHandle CreateSequence(const String& name, const String& seq, 
                              const String& role)
{
  return SequenceHandle(impl::SequenceImpl::FromString(name, seq, role));
}


String ConstSequenceHandle::GetGaplessString() const
{
  this->CheckValidity();
  return Impl()->GetGaplessString();
}

int ConstSequenceHandle::GetOffset() const
{
  this->CheckValidity();
  return Impl()->GetOffset();  
}


int ConstSequenceHandle::GetLength() const
{
  this->CheckValidity();
  return Impl()->GetLength();
}


char ConstSequenceHandle::GetOneLetterCode(int position) const
{
  this->CheckValidity();
  return Impl()->GetOneLetterCode(position);
}

mol::ResidueView ConstSequenceHandle::GetResidue(int position) const
{
  this->CheckValidity();
  return Impl()->GetResidue(position);  
}

mol::EntityView ConstSequenceHandle::GetAttachedView() const
{
  this->CheckValidity();
  return Impl()->GetAttachedView(); 
}

SequenceHandle ConstSequenceHandle::Copy() const
{
  this->CheckValidity();
  return Impl()->Copy();     
}

bool ConstSequenceHandle::HasAttachedView() const
{
  this->CheckValidity();
  return Impl()->HasAttachedView();   
}


int ConstSequenceHandle::GetResidueIndex(int pos) const
{
  this->CheckValidity();
  return Impl()->GetResidueIndex(pos);   
}

int ConstSequenceHandle::GetPos(int residue_index) const
{
  this->CheckValidity();
  return Impl()->GetPos(residue_index);   
}


int ConstSequenceHandle::GetFirstNonGap() const
{
  this->CheckValidity();
  return Impl()->GetFirstNonGap(); 
}


int ConstSequenceHandle::GetLastNonGap() const
{
  this->CheckValidity();
  return Impl()->GetLastNonGap(); 
}

const String& ConstSequenceHandle::GetName() const
{
  this->CheckValidity();
  return Impl()->GetName();
}

const String& ConstSequenceHandle::GetString() const
{
  this->CheckValidity();
  return Impl()->GetString(); 
}

void SequenceHandle::SetName(const String& name)
{
  this->CheckValidity();
  return Impl()->SetName(name);   
}

void SequenceHandle::SetString(const String& seq)
{
  this->CheckValidity();
  return Impl()->SetString(seq);   
}

void SequenceHandle::SetOffset(int offset)
{
  this->CheckValidity();
  return Impl()->SetOffset(offset);   
}


void SequenceHandle::AttachView(const mol::EntityView& view)
{
  this->CheckValidity();
  Impl()->AttachView(view);   
}


void SequenceHandle::AttachView(const mol::EntityView& view, 
                                const String& chain_name)
{
  this->CheckValidity();
  Impl()->AttachView(view, chain_name);   
}

#if(OST_INFO_ENABLED)

/// \brief export sequence to info
void  SequenceToInfo(const ConstSequenceHandle& sequence,
                     info::InfoGroup& group)
{
  impl::SequenceImplToInfo(sequence.Impl(), group);
}

/// \brief create sequence from info
SequenceHandle SequenceFromInfo(info::InfoGroup& group)
{
  return SequenceHandle(impl::SequenceImplFromInfo(group));
}
#endif



std::ostream& operator<<(std::ostream& os, const ConstSequenceHandle& sequence)
{
  if (sequence.IsValid()) {
    os << sequence.GetName() << ": " << sequence.GetString();    
  } else {
    os << "invalid";
  }

  return os;
}

bool SequenceHandle::operator==(const SequenceHandle& rhs) const
{
  return impl_==rhs.impl_;
}

bool SequenceHandle::operator!=(const SequenceHandle& rhs) const
{
  return impl_!=rhs.impl_;
}


void SequenceHandle::CheckValidity() const
{
  if (!impl_) {
    throw InvalidHandle();
  }
}

bool SequenceHandle::IsValid() const
{
  return impl_.get()!=0;
}

impl::SequenceImplPtr& SequenceHandle::Impl() const 
{
  return impl_;
}


SequenceHandle::SequenceHandle():
  impl_()
{ }


String SequenceHandle::GetGaplessString() const
{
  this->CheckValidity();
  return Impl()->GetGaplessString();
}

void SequenceHandle::Normalise() {
  this->CheckValidity();
  Impl()->Normalise();
}

int SequenceHandle::GetOffset() const
{
  this->CheckValidity();
  return Impl()->GetOffset();  
}

SequenceHandle::operator ConstSequenceHandle() const
{
  return ConstSequenceHandle(impl_);
}

int SequenceHandle::GetLength() const
{
  this->CheckValidity();
  return Impl()->GetLength();
}


char SequenceHandle::GetOneLetterCode(int position) const
{
  this->CheckValidity();
  return Impl()->GetOneLetterCode(position);
}

void SequenceHandle::SetOneLetterCode(int position, char new_char)
{
  this->CheckValidity();
  Impl()->SetOneLetterCode(position, new_char);
}

mol::ResidueView SequenceHandle::GetResidue(int position) const
{
  this->CheckValidity();
  return Impl()->GetResidue(position);  
}

mol::EntityView SequenceHandle::GetAttachedView() const
{
  this->CheckValidity();
  return Impl()->GetAttachedView(); 
}

SequenceHandle SequenceHandle::Copy() const
{
  this->CheckValidity();
  return Impl()->Copy();     
}

bool SequenceHandle::HasAttachedView() const
{
  this->CheckValidity();
  return Impl()->HasAttachedView();   
}


int SequenceHandle::GetResidueIndex(int pos) const
{
  this->CheckValidity();
  return Impl()->GetResidueIndex(pos);   
}

int SequenceHandle::GetPos(int residue_index) const
{
  this->CheckValidity();
  return Impl()->GetPos(residue_index);   
}


int SequenceHandle::GetFirstNonGap() const
{
  this->CheckValidity();
  return Impl()->GetFirstNonGap(); 
}


int SequenceHandle::GetLastNonGap() const
{
  this->CheckValidity();
  return Impl()->GetLastNonGap(); 
}

const String& SequenceHandle::GetName() const
{
  this->CheckValidity();
  return Impl()->GetName();
}
void SequenceHandle::Append(char olc)
{
  this->CheckValidity();
  Impl()->Append(olc);
}

int ConstSequenceHandle::GetIndex(const String& substr) const
{
  this->CheckValidity();
  return Impl()->GetIndex(substr);
}


const String& SequenceHandle::GetString() const
{
  this->CheckValidity();
  return Impl()->GetString(); 
}

GenericPropContainerImpl* ConstSequenceHandle::GpImpl()
{
  return Impl().get();
}

const GenericPropContainerImpl* ConstSequenceHandle::GpImpl() const
{
  return Impl().get();
}

GenericPropContainerImpl* SequenceHandle::GpImpl()
{
  return Impl().get();
}

const GenericPropContainerImpl* SequenceHandle::GpImpl() const
{
  return Impl().get();
}

int SequenceHandle::GetIndex(const String& substr) const
{
  this->CheckValidity();
  return Impl()->GetIndex(substr);
}

char SequenceHandle::operator[](size_t index) const
{
  this->CheckValidity();
  return this->GetString()[index];
}

bool Match(const ConstSequenceHandle& s1, const ConstSequenceHandle& s2)
{
  return Match(s1.GetString(), s2.GetString());
}

bool Match(const String& s1, const String& s2)
{
  if (s1.size()!=s2.size()) {
    return false;
  }
  for (size_t i=0; i<s1.size(); ++i) {
    char c1=s1[i];
    char c2=s2[i];
    if (toupper(c1)!=toupper(c2) && toupper(c1)!='X' && toupper(c2)!='X') {
      return false;
    }
  }
  return true;
}
}}
