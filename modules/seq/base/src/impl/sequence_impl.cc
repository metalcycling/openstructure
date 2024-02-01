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

#include <boost/format.hpp>
#include <ost/config.hh>
#if(OST_INFO_ENABLED)
#include <ost/info/info.hh>
#endif

#include <ost/integrity_error.hh>
#include <ost/mol/chain_view.hh>

#include "sequence_impl.hh"

namespace ost { namespace seq { namespace impl {

using boost::format;

namespace {

bool is_sane_code(char code)
{
  return (isalpha(code) || code=='?' || code=='.');
}

}

bool SequenceImpl::IsSequenceStringSane(const String& seq_string)
{
  for (size_t i=0; i<seq_string.length(); ++i) {
    if (is_sane_code(seq_string[i]) || seq_string[i]=='-')
      continue;
    return false;
  }
  return true;
}

SequenceImplPtr SequenceImpl::FromString(const String& seq_name,
                                         const String& seq_string, 
                                         const String& role)
{

  if (SequenceImpl::IsSequenceStringSane(seq_string)) {
    return SequenceImplPtr(new SequenceImpl(seq_name, seq_string, role));
  }
  throw InvalidSequence();
}

void SequenceImpl::Append(char olc)
{
  seq_string_.push_back(olc);
  if (olc=='-') {
    if (seq_string_.size()>1 && 
        seq_string_[seq_string_.size()-2]=='-') {
      shifts_.back().shift++;
    } else {
      Shift new_shift;
      new_shift.start=seq_string_.size()-1;
      new_shift.shift=1;
      shifts_.push_back(new_shift);
    }
  }

}

int SequenceImpl::GetIndex(const String& substr) const
{
  size_t pos=seq_string_.find(substr);
  if (pos==String::npos) {
    return -1;
  }
  return pos;
}

void SequenceImpl::SetString(const String& seq)
{
  if (SequenceImpl::IsSequenceStringSane(seq)) {
    seq_string_=seq;
    this->ShiftsFromSequence();
  }
  else {
    throw InvalidSequence();
  }
}

SequenceImpl::SequenceImpl(const String& seq_name,
                   const String& seq_string, const String& role)
  : seq_name_(seq_name), seq_string_(seq_string), seq_role_(role), offset_(0)
{
  this->ShiftsFromSequence();
}

SequenceImplPtr SequenceImpl::Copy() const
{
  SequenceImplPtr new_seq(new SequenceImpl(seq_name_, seq_string_, seq_role_));
  new_seq->offset_=offset_;
  new_seq->shifts_=shifts_;
  new_seq->attached_view_=attached_view_;
  return new_seq;
}

void SequenceImpl::ShiftsFromSequence()
{
  shifts_.clear();
  for (size_t i=0; i<seq_string_.length(); ++i) {
    Shift new_shift;
    new_shift.shift=0;
    new_shift.start=i;
    while (i<seq_string_.length() && seq_string_[i]=='-') {
      ++new_shift.shift;
      ++i;
    }
    if (new_shift.shift>0) {
      shifts_.push_back(new_shift);
    }
  }
}

void SequenceImpl::SetOneLetterCode(int position, char new_char)
{
  if (position<0 || position>=static_cast<int>(seq_string_.length()))
    throw std::out_of_range("Position is not covered in sequence");
  if (new_char=='-') {
    if (seq_string_[position]!='-') {
      seq_string_[position]=new_char;
      this->ShiftsFromSequence();
      return;
    }
  } else  {
    if (seq_string_[position]=='-') {
      seq_string_[position]=new_char;
      this->ShiftsFromSequence();
      return;
    } else {
      seq_string_[position]=new_char;
    }
  }
}

int SequenceImpl::GetResidueIndex(int pos) const
{
  if (pos<0 || pos>=static_cast<int>(seq_string_.length()))
    throw std::out_of_range("Position is not covered in sequence");
  if (seq_string_[pos]=='-')
    throw Error("Requested position contains a gap");
  std::list<Shift>::const_iterator i;
  int shift=0;
  for (i=shifts_.begin(); i!=shifts_.end(); ++i) {
    const Shift& s=*i;
    if (s.start<=pos)
      shift+=s.shift;
  }
  return pos-shift+offset_;
}

void SequenceImpl::SetName(const String& seq_name)
{
  seq_name_=seq_name;
}

const String& SequenceImpl::GetName() const
{
  return seq_name_;
}

int SequenceImpl::GetPos(int index) const
{
  int shifted_index=index-offset_;
  int pos=this->GetPosNoBounds(shifted_index);
  if (pos<0 || pos>=static_cast<int>(seq_string_.length()))
    throw Error("number not covered in sequence");
  return pos;
}

int SequenceImpl::GetPosNoBounds(int index) const
{
  int shift=0;
  std::list<Shift>::const_iterator i;
  for (i=shifts_.begin(); i!=shifts_.end(); ++i) {
    const Shift& s=*i;
    if (s.start-shift<=index)
      shift+=s.shift;
  }
  return index+shift;
}

/// \brief  Set offset for sequence alignment.
void SequenceImpl::SetOffset(int offset)
{
  offset_=offset;
}

/// \brief Get sequence offset
int SequenceImpl::GetOffset() const
{
  return offset_;
}

int SequenceImpl::GetLength() const {
  return static_cast<int>(seq_string_.length());
}

char SequenceImpl::GetOneLetterCode(int position) const
{
  if (position<0 || position>=static_cast<int>(seq_string_.length()))
    throw std::out_of_range("Position is not covered in sequence");
  return seq_string_[position];
}

String SequenceImpl::GetGaplessString() const {
  String seq_string_no_gaps = "";
  for (size_t i=0; i<seq_string_.length(); ++i) {
    if ( seq_string_[i]!='-') {
      seq_string_no_gaps = seq_string_no_gaps + seq_string_[i];
    }
  }
  return seq_string_no_gaps;
}

mol::ResidueView SequenceImpl::GetResidue(int position) const
{
  if (!attached_view_.IsValid())
    return mol::ResidueView();

  if (static_cast<int>(seq_string_.size())<=position || position <0 ||
      seq_string_[position]=='-') {
    return mol::ResidueView();
  }
  int p=this->GetResidueIndex(position);
  if (p!=-1) {
    mol::ChainView chain=(attached_view_.GetChainList())[0];
    return chain.GetResidueByIndex(p);
  } else {
    return mol::ResidueView();
  }
}


mol::EntityView SequenceImpl::GetAttachedView() const
{
  return attached_view_;
}


void SequenceImpl::AttachView(const mol::EntityView& view)
{
  static const char* msg="Expected 1 chain, but %d chains found";
  attached_view_=view;
  if (view.IsValid() && attached_view_.GetChainCount()!=1) {
    throw IntegrityError(str(format(msg) % attached_view_.GetChainCount()));
  }
}

void SequenceImpl::AttachView(const mol::EntityView& view, 
                              const String& chain_name)
{
  attached_view_=view.Select(chain_name);
  if (attached_view_.GetChainCount()==0) {
    throw IntegrityError("No chain with given name in view");
  }
}

bool SequenceImpl::HasAttachedView() const
{
  return attached_view_.IsValid();
}

#if(OST_INFO_ENABLED)
void SequenceImplToInfo(const SequenceImplPtr& sequence, info::InfoGroup& group)
{
  group.SetTextData(sequence->GetString());
  group.SetAttribute("name", sequence->GetName());
  std::ostringstream ss;
  ss << sequence->GetOffset();
  group.SetAttribute("offset", ss.str());
}

SequenceImplPtr SequenceImplFromInfo(const info::InfoGroup& group)
{
  String text=group.GetTextData();
  std::istringstream ss(group.GetAttribute("offset"));
  int offset=0;
  ss >> offset;
  SequenceImplPtr sequence=SequenceImpl::FromString(group.GetAttribute("name"), 
                                                    text);
  sequence->SetOffset(offset);
  return sequence;
}
#endif

int SequenceImpl::GetFirstNonGap() const
{
  size_t pos=seq_string_.find_first_not_of("-");
  return pos==String::npos? -1 : pos;  
}

int SequenceImpl::GetLastNonGap() const
{
  size_t pos=seq_string_.find_last_not_of("-");
  return pos==String::npos? -1 : pos;
}

void SequenceImpl::Cut(int start, int n)
{
  seq_string_.replace(start, n, "");
  this->ShiftsFromSequence();
}

void SequenceImpl::Replace(const String& str,int start, int end)
{
  seq_string_.replace(start, end-start, str);
  this->ShiftsFromSequence();
}

void SequenceImpl::Normalise() {
  const char gaps[]={'.', '-', '\0'};
  bool is_gap;
  String new_seq_string = "";
  int n_gaps = sizeof(gaps) / sizeof(gaps[0]);
  for (size_t i=0; i<seq_string_.length(); ++i) {
    is_gap = false;
    for (int j=0; j<n_gaps; ++j) {
      if (seq_string_[i]==gaps[j]) {
        is_gap = true;
        break;
      }
    }
    if (is_gap == false) {
      new_seq_string = new_seq_string + char(toupper(seq_string_[i]));
    }
  }
  seq_string_ = new_seq_string;
  this->ShiftsFromSequence();
}

void SequenceImpl::ShiftRegion(int start, int end, int amount)
{
  if(start > end || start + amount < 0 || end + amount > this->GetLength()){
    throw std::out_of_range("ShiftRegion: invalid region");
  }
  String str1=seq_string_.substr(start, end-start);
  if (amount<0) {
    String str2=seq_string_.substr(start+amount, -amount);
    seq_string_.replace(start+amount, end-start, str1);
    seq_string_.replace(end+amount, -amount, str2);
  } else {
    String str2=seq_string_.substr(end, amount);
    seq_string_.replace(start+amount, end-start, str1);
    seq_string_.replace(start, amount, str2);
  }
}

}}} //ns

