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
#include <limits>

#include <ost/config.hh>
#if(OST_INFO_ENABLED)
#include <ost/info/info.hh>
#endif
#include <ost/seq/impl/sequence_list_impl.hh>

namespace ost { namespace seq { namespace impl {

void SequenceListImpl::AddSequence(const SequenceImplPtr& sequence)
{
  list_.push_back(sequence);
}

int SequenceListImpl::GetPos(int seq_index, int number) const
{
  SequenceImplPtr s=this->GetSequence(seq_index);
  if (s) {
    return s->GetPos(number);
  }
  return -1;
}

int SequenceListImpl::GetResidueIndex(int seq_index, int pos) const
{
  SequenceImplPtr s=this->GetSequence(seq_index);
  if (s) {
    return s->GetResidueIndex(pos);
  }

  return -1;
}

SequenceImplPtr SequenceListImpl::FindSequence(const String& name) const
{
  for (size_t i=0; i<list_.size(); ++i) {
    if (list_[i]->GetName()==name) {
      return list_[i];
    }

  }
  return SequenceImplPtr();
}


int SequenceListImpl::FindSequenceIndex(const String& name) const
{
  for (size_t i=0; i<list_.size(); ++i) {
    if (list_[i]->GetName()==name) {
      return int(i);
    }
  }
  return -1;
}
  
  
String SequenceListImpl::ToString(int width) const
{
  std::stringstream buffer;
  // determine label width
  int label_size=0;
  for (std::vector<SequenceImplPtr>::const_iterator i=list_.begin(),
       e=list_.end(); i!=e; ++i) {
    label_size=std::max(label_size, static_cast<int>((*i)->GetName().size()));
  }
  label_size+=2;
  if (label_size>width/4) {
    label_size=width/4;
  }  
  int offset=0;
  bool done=false;
  int text_len=width-label_size;
  while (!done) {
    done=true;
    for (std::vector<SequenceImplPtr>::const_iterator i=list_.begin(),
         e=list_.end(); i!=e; ++i) {
      SequenceImplPtr s=*i;
      if (static_cast<int>(s->GetName().size())>label_size-2) {
        buffer << s->GetName().substr(0, label_size-4) << "... ";
      } else {
        buffer << s->GetName() << String(label_size-s->GetName().size(), ' ');
      }

      if (offset<s->GetLength()) {
        int rem=s->GetLength()-offset;
        int actual=text_len>rem ? rem : text_len;
        buffer << s->GetString().substr(offset, actual) << std::endl;
        if (s->GetLength()>offset+actual) {
          done=false;
        }
      } else {
        buffer << std::endl;
      }
    }
    if (!done) {
      buffer << std::endl;
    }
    offset+=text_len;
  }
  return buffer.str();
}

void SequenceListImpl::RemoveSequence(int index)
{
  if (index<0 || index>=static_cast<int>(list_.size())) {
    return;
  }
  list_.erase(list_.begin()+index);
}
int SequenceListImpl::GetMinLength() const
{
  int min_length=std::numeric_limits<int>::max();
  for (std::vector<SequenceImplPtr>::const_iterator i=list_.begin(),
       e=list_.end(); i!=e; ++i) {
    min_length=std::min(min_length, (*i)->GetLength());
  }
  return min_length;
}


int SequenceListImpl::GetMaxLength() const
{
  int max_length=0;
  for (std::vector<SequenceImplPtr>::const_iterator i=list_.begin(),
       e=list_.end(); i!=e; ++i) {
    max_length=std::max(max_length, (*i)->GetLength());
  }
  return max_length;
}

SequenceListImplPtr SequenceListImpl::Copy() const
{
  SequenceListImplPtr new_ali(new SequenceListImpl);
  for (std::vector<SequenceImplPtr>::const_iterator i=list_.begin(),
       e=list_.end(); i!=e; ++i) {
    new_ali->AddSequence((*i)->Copy());
  }
  return new_ali;
}

#if(OST_INFO_ENABLED)
void SequenceListImplToInfo(const SequenceListImplPtr& seq_list,
                            info::InfoGroup& group)
{
  for (SequenceListImpl::ConstIterator i=seq_list->Begin(),
       e=seq_list->End(); i!=e; ++i) {
    info::InfoGroup seq_group=group.CreateGroup("sequence");
    SequenceImplToInfo(*i, seq_group);
  }
} 


SequenceListImplPtr SequenceListImplFromInfo(info::InfoGroup& group)
{
  SequenceListImplPtr seq_list(new SequenceListImpl);
  info::InfoGroupList sequences=group.GetGroups("sequence");
  for (info::InfoGroupList::const_iterator i=sequences.begin(),
       e=sequences.end(); i!=e; ++i) {
    seq_list->AddSequence(SequenceImplFromInfo(*i));
  }
  return seq_list;
}
#endif

SequenceListImplPtr SequenceListImpl::Slice(int first, int n) const
{
  int last=std::min(first+n, this->GetCount());
  if (first>=0) {
    SequenceListImplPtr seq_list(new SequenceListImpl);
    for (int i=first; i<last; ++i) {
      seq_list->AddSequence(this->GetSequence(i));
    }
    return seq_list;
  } else {
    return SequenceListImplPtr(new SequenceListImpl);
  }
}

Real SequenceListImpl::GetCoverage(int seq_index) const
{
  int a=0, b=0;
  String seq_string_0 = GetSequence(0)->GetString();
  String seq_string_index = GetSequence(seq_index)->GetString();
  for (int i=0; i<GetSequence(0)->GetLength(); ++i) {
    if (seq_string_0[i]!='-') {
      a+=1;
      if (seq_string_index[i]!='-') {
        b+=1;
      }
    }
  }
  Real coverage=Real(b)/Real(a);
  return coverage;
}


}}}
