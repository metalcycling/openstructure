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
#ifndef OST_SEQ_SEQUENCE_LIST_HH
#define OST_SEQ_SEQUENCE_LIST_HH

#include <ost/seq/impl/sequence_list_impl_fw.hh>
#include <ost/seq/impl/sequence_list_iterator.hh>
#include <ost/seq/sequence_handle.hh>

namespace ost { namespace seq {

class SequenceList;

/// \brief list of immutable sequences. 
/// 
/// The sequence list is reference counted and thus very cheap to pass around.
class DLLEXPORT_OST_SEQ ConstSequenceList {
public:
  friend class SequenceList;
  typedef impl::SequenceListIterator<ConstSequenceHandle, 
                                     impl::SequenceImplListIterator> Iterator;
  typedef Iterator iterator;
  ConstSequenceList();
  
  void AddSequence(const ConstSequenceHandle& sequence);
  
  /// \brief get number of sequences in list
  int GetCount() const;
  
  ConstSequenceHandle operator[](int index) const;
  Iterator Begin() const;
  Iterator End() const;
  
  // for stl compat
  Iterator begin() const { return Begin(); }
  Iterator end() const { return End(); }
  
  bool SequencesHaveEqualLength() const;
  bool IsValid() const;  
  
  /// \brief create a sequence list consisting of n first sequences
  /// 
  /// if n is negative, the sequences taken from the end of the list
  ConstSequenceList Take(int n) const;
  
  /// \brief create a sequence list from the given splice interval
  ConstSequenceList Slice(int first, int n) const;
  
  ConstSequenceHandle FindSequence(const String& name) const;
  int GetMinLength() const;
  int GetMaxLength() const;
  /// \internal
  impl::SequenceListImplPtr& Impl() const;  
  ConstSequenceList(const impl::SequenceListImplPtr& impl);    
protected:

  void CheckValidity() const;
  mutable impl::SequenceListImplPtr impl_;
};

/// \brief list of sequences. 
/// 
/// The sequence list is reference counted and thus very cheap to pass around.
class DLLEXPORT_OST_SEQ SequenceList {
public:
  
  typedef impl::SequenceListIterator<SequenceHandle, 
                                     impl::SequenceImplListIterator> Iterator;
  typedef Iterator iterator;
  SequenceList();
  void AddSequence(const SequenceHandle& sequence);
  
  bool SequencesHaveEqualLength() const;  
  // \brief get number of sequences
  int GetCount() const;
  SequenceHandle operator[](int index) const;  
  Iterator Begin() const;
  Iterator End() const;
  
  /// \brief create a sequence list consisting of n first sequences
  /// 
  /// if n is negative, the sequences taken from the end of the list
  SequenceList Take(int n) const;
  
  /// \brief create a sequence list from the given splice interval
  SequenceList Slice(int first, int n) const;
  // for stl compat
  Iterator begin() const { return Begin(); }
  Iterator end() const { return End(); }  
  
  bool IsValid() const;
  
  int GetMinLength() const;
  int GetMaxLength() const;  
  
  SequenceHandle FindSequence(const String& name) const;
  // \internal
  impl::SequenceListImplPtr& Impl() const;
  SequenceList(const impl::SequenceListImplPtr& impl);  
  
  operator ConstSequenceList() const;
protected:
  void CheckValidity() const;

  mutable impl::SequenceListImplPtr impl_;
};

SequenceList DLLEXPORT_OST_SEQ CreateSequenceList();
ConstSequenceList DLLEXPORT_OST_SEQ CreateConstSequenceList();

#if(OST_INFO_ENABLED)
/// \brief export sequence list to info
void DLLEXPORT_OST_SEQ SequenceListToInfo(const ConstSequenceList& seq_list,
                                          info::InfoGroup& group);

/// \brief create sequence list object from info
SequenceList DLLEXPORT_OST_SEQ SequenceListFromInfo(info::InfoGroup& group);
#endif
}}

#endif
