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
#ifndef OST_SEQ_SEQUENCE_LIST_IMPL_HH
#define OST_SEQ_SEQUENCE_LIST_IMPL_HH

/*
  Author: Marco Biasini
 */
#include <vector>
#include <boost/shared_ptr.hpp>
#include <ost/seq/impl/sequence_list_impl_fw.hh>
#include <ost/seq/impl/sequence_impl.hh>

namespace ost { namespace seq { namespace impl {

/// \internal
class DLLEXPORT_OST_SEQ SequenceListImpl {
public:
  typedef SequenceImplListConstIterator  ConstIterator;
  typedef SequenceImplListIterator  Iterator;  
  SequenceListImpl() { } 
  
  /// \brief add sequence
  void AddSequence(const SequenceImplPtr& sequence);
  /// \brief get number of sequences in list
  int GetCount() const { return list_.size(); }
  
  SequenceImplPtr& GetSequence(int i) {
    unsigned int index = static_cast<unsigned int>(i);
    if (index<list_.size()) {
      return list_[index];
    }
    throw std::out_of_range("Index not covered SequenceList");
  }
  
  const SequenceImplPtr& GetSequence(unsigned int i) const {
    unsigned int index = static_cast<unsigned int>(i);
    if (index<list_.size()) {
      return list_[index];
    }
    throw std::out_of_range("Index not covered SequenceList");
  }
  
  
  void RemoveSequence(int index);
  int GetPos(int seq_index, int residue_index) const;

  int GetResidueIndex(int seq_index, int pos) const;

  SequenceImplPtr FindSequence(const String& name) const;

  int FindSequenceIndex(const String& name) const;
  
  String ToString(int width=80) const;

  SequenceListImplPtr Slice(int first, int n) const;
  int GetMinLength() const;


  int GetMaxLength() const;
  
  bool SequencesHaveEqualLength() const 
  { 
    return this->GetMaxLength()==this->GetMinLength();
  }

  Real GetCoverage(int seq_index) const;

  SequenceListImplPtr Copy() const;

  Iterator Begin() { return list_.begin(); }
  Iterator End() { return list_.end(); }
  ConstIterator Begin() const { return list_.begin(); }
  ConstIterator End() const { return list_.end(); }  
  
  
private:
  std::vector<SequenceImplPtr> list_;
};

#if(OST_INFO_ENABLED)
/// \brief export sequence list impl to info
/// \internal
void DLLEXPORT_OST_SEQ 
SequenceListImplToInfo(const SequenceListImplPtr& seq_list,
                        info::InfoGroup& group);

/// \brief create sequence impl list object from info
/// \internal
SequenceListImplPtr DLLEXPORT_OST_SEQ 
SequenceListImplFromInfo(info::InfoGroup& group);
#endif

}}}

#endif
