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
#include <ost/mol/residue_view.hh>
#include <ost/mol/atom_view.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/mol/residue_view.hh>
#include <ost/seq/impl/sequence_list_impl.hh>
#include <ost/seq/impl/sequence_impl.hh>
#include <ost/seq/sequence_list.hh>
#include <ost/seq/aligned_region.hh>
#include <ost/seq/aligned_column_iterator.hh>
#include <ost/integrity_error.hh>

namespace ost { namespace seq {

AlignmentHandle CreateAlignment()
{
  return AlignmentHandle(impl::SequenceListImplPtr(new impl::SequenceListImpl));
}

AlignmentHandle::AlignmentHandle():
  impl_()
{
}

AlignmentHandle::AlignmentHandle(const impl::SequenceListImplPtr& impl):
  impl_(impl)
{ }


int AlignmentHandle::GetPos(int seq_index, int residue_index) const
{
  this->CheckValidity();
  return impl_->GetPos(seq_index, residue_index);
}

int AlignmentHandle::GetResidueIndex(int seq_index, int pos) const
{
  this->CheckValidity();
  return impl_->GetResidueIndex(seq_index, pos);
}


void AlignmentHandle::AddSequence(const ConstSequenceHandle& sequence)
{
  this->CheckValidity();
  if (!sequence.IsValid()) {
    throw InvalidSequence();
  }
  if (impl_->GetCount()>0 &&
      impl_->GetSequence(0)->GetLength()!=sequence.GetLength()) {
    throw Error("sequence doesn't have the same length as the "
                             "alignment");
  }
  return impl_->AddSequence(sequence.Impl());
}

ConstSequenceHandle AlignmentHandle::GetSequence(int seq_id) const
{
  this->CheckValidity();
  return impl_->GetSequence(seq_id);
}


String AlignmentHandle::ToString(int width) const
{
  this->CheckValidity();
  return impl_->ToString(width);
}

void AlignmentHandle::RemoveSequence(int seq_index)
{
  this->CheckValidity();
  impl_->RemoveSequence(seq_index);
}

int AlignmentHandle::GetLength() const
{
  this->CheckValidity();
  if (impl_->GetCount()==0) {
    return 0;
  } else {
    return impl_->GetSequence(0)->GetLength();
  }
}


AlignmentHandle AlignmentHandle::Copy() const
{
  this->CheckValidity();
  return impl_->Copy();
}

bool AlignmentHandle::operator==(const AlignmentHandle& rhs) const
{
  return rhs.impl_==impl_;
}

bool AlignmentHandle::operator!=(const AlignmentHandle& rhs) const
{
    return rhs.impl_!=impl_;
}

void AlignmentHandle::CheckValidity() const
{
  if (!impl_) {
    throw InvalidHandle();
  }
}

int AlignmentHandle::GetCount() const
{
  this->CheckValidity();
  return impl_->GetCount();
}

AlignmentHandle AlignmentFromSequenceList(const SequenceList& seq_list)
{
  if (seq_list.IsValid() && seq_list.SequencesHaveEqualLength()) {
    return AlignmentHandle(seq_list.Impl());
  }
  throw Error("sequences have different lengths");
}

ConstSequenceList AlignmentHandle::GetSequences() const
{
  this->CheckValidity();
  return ConstSequenceList(impl_);
}


void AlignmentHandle::AttachView(int seq_index, const mol::EntityView& view)
{
  this->CheckValidity();
  impl_->GetSequence(seq_index)->AttachView(view);
}

void AlignmentHandle::AttachView(int seq_index, const mol::EntityView& view,
                                 const String& chain_name)
{
  this->CheckValidity();
  impl_->GetSequence(seq_index)->AttachView(view, chain_name);
}

ConstSequenceHandle AlignmentHandle::FindSequence(const String& name) const
{
  this->CheckValidity();
  return ConstSequenceHandle(impl_->FindSequence(name));
}

int AlignmentHandle::FindSequenceIndex(const String& name) const
{
  this->CheckValidity();
  return impl_->FindSequenceIndex(name);
}  

void AlignmentHandle::Cut(int start, int end)
{
  this->CheckValidity();
  for (impl::SequenceListImpl::Iterator i=impl_->Begin(),
       e=impl_->End(); i!=e; ++i) {
    (*i)->Cut(start, end-start);
  }
}

void AlignmentHandle::Replace(const AlignedRegion& aln_r, int start, int end){
  this->CheckValidity();
  //check that alignment handle and aligned region contain same number of sequences
  if (impl_->GetCount() != aln_r.GetAlignmentHandle().GetCount()) {
    throw IntegrityError("alignment handle and aligned region are required "
                         "to have the same number of sequences");
  }
  int aln_r_start=aln_r.GetStart();
  int aln_r_length=aln_r.GetLength();
  AlignmentHandle aln=aln_r.GetAlignmentHandle();
  //iterate over sequences and replace part of sequences with the substrings
  //from aligned region
  for (int i=0;i<impl_->GetCount() ;++i) {
    String s=aln.GetSequence(i).GetString().substr(aln_r_start, aln_r_length);
    this->GetSequence(i).Impl()->Replace(s, start, end);
  }
}

void AlignmentHandle::ShiftRegion(int start, int end, int amount,
                                  int master)
{
  this->CheckValidity();
  int cnt=0;
  if(master!=-1){
    impl::SequenceImplPtr handle = this->GetSequence(master).Impl();
    handle->ShiftRegion(start, end, amount);
  }
  else{
    for (impl::SequenceListImpl::Iterator i=impl_->Begin(),
         e=impl_->End(); i!=e; ++i, ++cnt) {
        (*i)->ShiftRegion(start, end, amount);
    }
  }
}

AlignedRegion AlignmentHandle::MakeRegion(int start, int n, int master) const
{
  this->CheckValidity();
  if((start<0 )||( n < 0 )|| ((start + n) > this->GetLength())){
    throw std::out_of_range("invalid region");
  }
  return AlignedRegion(*this, start, start+n, master);
}

mol::ResidueView AlignmentHandle::GetResidue(int seq_index, int pos) const
{
  this->CheckValidity();
  return impl_->GetSequence(seq_index)->GetResidue(pos);
}

char AlignmentHandle::GetOneLetterCode(int seq_index, int pos) const
{
  this->CheckValidity();
  return impl_->GetSequence(seq_index)->GetOneLetterCode(pos);
}

AlignedColumnIterator AlignmentHandle::begin() const
{
  // no need to check for validity here
  return AlignedColumnIterator(*this, 0, this->GetLength());
}

AlignedColumnIterator AlignmentHandle::end() const
{
  // no need to check for validity here
  return AlignedColumnIterator(*this, this->GetLength(), this->GetLength());
}

AlignedColumn AlignmentHandle::operator[](int index) const
{
  return AlignedColumn(*this, index);
}

void AlignmentHandle::SetSequenceName(int seq_index, const String& name)

{
  this->CheckValidity();
  impl_->GetSequence(seq_index)->SetName(name);
}

void AlignmentHandle::SetSequenceOffset(int seq_index, int offset)
{
  this->CheckValidity();
  impl_->GetSequence(seq_index)->SetOffset(offset);
}

int AlignmentHandle::GetSequenceOffset(int seq_index)
{
  this->CheckValidity();
  return impl_->GetSequence(seq_index)->GetOffset();
}

Real AlignmentHandle::GetCoverage(int seq_index) const
{
  this->CheckValidity();
  return impl_->GetCoverage(seq_index);
}

mol::EntityViewPair AlignmentHandle::GetMatchingBackboneViews(int idx0, int idx1) const
{
  this->CheckValidity();
  const impl::SequenceImpl& s1=*impl_->GetSequence(idx0).get();
  const impl::SequenceImpl& s2=*impl_->GetSequence(idx1).get();
  if (!s1.HasAttachedView() || !s2.HasAttachedView()) {
    throw Error("both sequences must have a view attached");
  }
  mol::EntityView v1=s1.GetAttachedView().CreateEmptyView();
  mol::EntityView v2=s2.GetAttachedView().CreateEmptyView();
  for (int i=0; i<s1.GetLength(); ++i) {
    if (s1[i]=='-' && s2[i]=='-') {
      continue;
    }
    mol::ResidueView r1=s1.GetResidue(i);
    mol::ResidueView r2=s2.GetResidue(i);
    if (!r1.IsValid() || !r2.IsValid()) {
      continue;
    }
    const char* bb_anames[]={"N", "CA", "C", "O"};
    //for (size_t j=0; )
    for (size_t j=0; j<4; ++j) {
      mol::AtomView a1=r1.FindAtom(bb_anames[j]);
      mol::AtomView a2=r2.FindAtom(bb_anames[j]);
      if (!a1.IsValid() || !a2.IsValid()) {
        continue;
      }
      v1.AddAtom(a1);
      v2.AddAtom(a2);
    }
  }
  return mol::EntityViewPair(v1, v2);
}


const String& AlignmentHandle::GetSequenceRole(int seq_index)
{
  this->CheckValidity();
  return impl_->GetSequence(seq_index)->GetRole();
  
}
  
void AlignmentHandle::SetSequenceRole(int seq_index, const String& role)
{
  this->CheckValidity();
  impl_->GetSequence(seq_index)->SetRole(role);
  
}


}}
