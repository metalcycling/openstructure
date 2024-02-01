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
#include <ost/config.hh>
#if(OST_INFO_ENABLED)
#include <ost/info/info.hh>
#endif

#include <ost/seq/impl/sequence_list_impl.hh>
#include <ost/seq/sequence_list.hh>

namespace ost { namespace seq {


SequenceList CreateSequenceList()
{
  return SequenceList(impl::SequenceListImplPtr(new impl::SequenceListImpl));
}

SequenceList::SequenceList() {}

void SequenceList::AddSequence(const SequenceHandle& sequence)
{
  this->CheckValidity();
  impl_->AddSequence(sequence.Impl());
}

// \brief get number of sequences
int SequenceList::GetCount() const
{
  this->CheckValidity();
  return impl_->GetCount();
}

SequenceList::SequenceList(const impl::SequenceListImplPtr& impl):
  impl_(impl)
{ }

SequenceHandle SequenceList::operator[](int index) const
{
  this->CheckValidity();
  return SequenceHandle(impl_->GetSequence(index));
}

SequenceList::Iterator SequenceList::Begin() const
{
  this->CheckValidity();
  return SequenceList::Iterator(impl_->Begin(), impl_->End());
}

SequenceList::Iterator SequenceList::End() const
{
  this->CheckValidity();
  return SequenceList::Iterator(impl_->End(), impl_->End());
}

SequenceHandle SequenceList::FindSequence(const String& name) const
{
  this->CheckValidity();
  return SequenceHandle(impl_->FindSequence(name));
}

bool SequenceList::IsValid() const
{
  return impl_.get()!=NULL;
}

void SequenceList::CheckValidity() const 
{
  if (!impl_) {
    throw InvalidHandle();
  }
}

bool ConstSequenceList::IsValid() const
{
  return impl_.get()!=NULL;
}

void ConstSequenceList::CheckValidity() const 
{
  if (!impl_) {
    throw InvalidHandle();
  }
}

ConstSequenceList::ConstSequenceList() {}

void ConstSequenceList::AddSequence(const ConstSequenceHandle& sequence)
{
  this->CheckValidity();
  impl_->AddSequence(sequence.Impl());
}

/// \brief get number of sequences in list
int ConstSequenceList::GetCount() const
{
  this->CheckValidity();
  return impl_->GetCount();
}

ConstSequenceHandle ConstSequenceList::operator[](int index) const
{
  return ConstSequenceHandle(impl_->GetSequence(index));
}

ConstSequenceList::Iterator ConstSequenceList::Begin() const
{
  return ConstSequenceList::Iterator(impl_->Begin(), impl_->End());
}

ConstSequenceList::Iterator ConstSequenceList::End() const
{
  return ConstSequenceList::Iterator(impl_->End(), impl_->End());
}

ConstSequenceList::ConstSequenceList(const impl::SequenceListImplPtr& impl):
  impl_(impl)
{ }

bool ConstSequenceList::SequencesHaveEqualLength() const
{
  return impl_->SequencesHaveEqualLength();
}

bool SequenceList::SequencesHaveEqualLength() const
{
  return impl_->SequencesHaveEqualLength();
}

impl::SequenceListImplPtr& SequenceList::Impl() const
{
  return impl_;
}

impl::SequenceListImplPtr& ConstSequenceList::Impl() const
{
  return impl_;
}

int ConstSequenceList::GetMinLength() const
{
  this->CheckValidity();
  return impl_->GetMinLength();
}

int ConstSequenceList::GetMaxLength() const
{
  this->CheckValidity();
  return impl_->GetMaxLength();
}

int SequenceList::GetMinLength() const
{
  this->CheckValidity();
  return impl_->GetMinLength();
}

int SequenceList::GetMaxLength() const
{
  this->CheckValidity();
  return impl_->GetMaxLength();
}

#if(OST_INFO_ENABLED)
/// \brief export sequence list to info
void SequenceListToInfo(const ConstSequenceList& seq_list, 
                        info::InfoGroup& group)
{
  impl::SequenceListImplToInfo(seq_list.Impl(), group);
}

/// \brief create sequence list object from info
SequenceList SequenceListFromInfo(info::InfoGroup& group)
{
  return SequenceList(impl::SequenceListImplFromInfo(group));
}
#endif

SequenceList::operator ConstSequenceList() const
{
  return ConstSequenceList(impl_);
}

SequenceList SequenceList::Take(int n) const
{
  this->CheckValidity();
  if (n<0) {
    return SequenceList(impl_->Slice(impl_->GetCount()+n, -n));
  } else {
    return SequenceList(impl_->Slice(0, n));
  }
}

SequenceList SequenceList::Slice(int first, int n) const
{
  this->CheckValidity();
  return SequenceList(impl_->Slice(first, n));
}

ConstSequenceList ConstSequenceList::Take(int n) const
{
  this->CheckValidity();
  if (n<0) {
    return ConstSequenceList(impl_->Slice(impl_->GetCount()+n, -n));
  } else {
    return ConstSequenceList(impl_->Slice(0, n));
  }
}

ConstSequenceList ConstSequenceList::Slice(int first, int n) const
{
  this->CheckValidity();
  return ConstSequenceList(impl_->Slice(first, n));
}

ConstSequenceHandle ConstSequenceList::FindSequence(const String& name) const
{
  this->CheckValidity();
  return ConstSequenceHandle(impl_->FindSequence(name));
}

ConstSequenceList DLLEXPORT_OST_SEQ CreateConstSequenceList()
{
  return ConstSequenceList(impl::SequenceListImplPtr(new impl::SequenceListImpl));  
}

}}
