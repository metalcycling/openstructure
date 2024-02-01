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
#ifndef OST_SEQ_SEQUENCE_HANDLE_HH
#define OST_SEQ_SEQUENCE_HANDLE_HH

/*
  Author: Marco Biasini
 */

#include <ost/base.hh>
#include <ost/generic_property.hh>
#include <ost/config.hh>
#if(OST_INFO_ENABLED)
#include <ost/info/info_fw.hh>
#endif
#include <ost/seq/module_config.hh>

#include <ost/mol/residue_view.hh>
#include <ost/mol/entity_view.hh>

#include <ost/seq/impl/sequence_impl_fw.hh>

namespace ost { namespace seq {

class SequenceHandle;
class AlignmentHandle;

/// \brief immutable sequence
/// 
/// The ConstSequenceHandle provides all read-only methods of the 
/// \ref SequenceHandle "sequence handle".
class DLLEXPORT_OST_SEQ ConstSequenceHandle : 
  public ConstGenericPropContainer<ConstSequenceHandle> {
public:
  friend class ConstGenericPropContainer<ConstSequenceHandle>;
  friend class AlignmentHandle;
  friend class ConstSequenceList;
  friend class SequenceList;
  
  typedef String::const_iterator iterator;
  
  /// \brief create invalid sequence handle
  /// 
  /// \sa IsValid()
  ConstSequenceHandle();
  
  /// \brief       Get residue index corresponding to given sequence position
  /// \param pos   zero-based index
  /// \throws      out_of_range if pos is not in [0, length-1].
  int GetResidueIndex(int pos) const;

  /// \brief Get zero-based index for given residue number.
  /// \throws out_of_range, if number is not included in sequence.
  int GetPos(int residue_index) const;

  /// \brief Get position of first non-gap character in sequence
  int GetFirstNonGap() const;
  
  /// \brief Get position of last non-gap character in sequence
  int GetLastNonGap() const;  
  
  /// \brief Get name of sequence
  /// 
  /// \sa SequenceHandle::GetName()
  const String& GetName() const;  
  
  /// \brief get sequence as a string, including all the gaps
  /// 
  /// \sa GetGaplessString()
  const String& GetString() const;
  
  /// \brief Get sequence as string ignoring gaps
  String GetGaplessString() const;

  /// \brief Get sequence offset from N-terminus
  ///
  /// \sa SequenceHandle::SetOffset
  int GetOffset() const;  
  
  /// \brief Get lenght of sequence, including gaps.
  int GetLength() const;
  
  /// \brief get index of substring
  int GetIndex(const String& substr) const;

  /// \brief get one letter code of residue at position
  char GetOneLetterCode(int position) const;

  /// \brief get residue at position
  ///
  /// will return the residue view at the given sequence position or an invalid
  /// residue view when no view is attached, the index is out of bounds or the
  /// position contains a gap.
  mol::ResidueView GetResidue(int position) const;

  /// \brief get attached view. may be an invalid entity view
  /// \sa SequenceHandle::AttachView(const mol::EntityView&, const String&)
  mol::EntityView GetAttachedView() const;
  
  /// \brief create copy sequence
  /// The newly created sequence has the same attached view.
  SequenceHandle Copy() const;
  
  /// \brief whether the sequence has an attached view
  /// 
  /// \sa SequenceHandle::AttachView(const mol::EntityView&, const String&)
  bool HasAttachedView() const;
  
  const String& GetRole() const;

  bool operator==(const ConstSequenceHandle& rhs) const;
  bool operator!=(const ConstSequenceHandle& rhs) const;  
  
  char operator[](int index) const;
  
  iterator begin() const { return this->GetString().begin(); }
  iterator end() const { return this->GetString().end(); }
  
  /// \brief whether the sequence is valid
  bool IsValid() const;
  /// \internal
protected:
  GenericPropContainerImpl* GpImpl();

  const GenericPropContainerImpl* GpImpl() const;
public:
  ConstSequenceHandle(const impl::SequenceImplPtr& impl);
  impl::SequenceImplPtr& Impl() const;  
private:
  void CheckValidity() const;
  mutable impl::SequenceImplPtr impl_;  
};

/// \brief mutable sequence handle.
/// 
/// Sequences are mostly used as part of a \ref AlignmentHandle
/// "multiple sequence alignment". The class allows for fast mapping between
/// residue index and position in the sequence. The GetResidueIndex() method
/// maps from position in the sequence to residue index, taking the
/// \ref GetOffset() "sequence offset" into account. The reverse
/// mapping is done with GetPos().
/// 
/// Optionally, an entity view may be attached to the sequence with
/// #AttachView(). If #HasAttachedView() is true, the #GetResidue() method
/// can be used directly to map from sequence position to residue. Note that
/// it is required that the entity view has only one chain.
//
/// Sequences are created from string as
/// \code
/// String seq_str="ALGW--ALN"
/// SequenceHandle seq=CreateSequence("NAME", seq_str);
/// \endcode
///
/// Sequences IO is dealt with in the sequence module. For more information, 
/// consult \ref module_seq "this page".
class DLLEXPORT_OST_SEQ SequenceHandle : 
  public GenericPropContainer<SequenceHandle> {
public:
  typedef String::const_iterator iterator;
  friend class GenericPropContainer<SequenceHandle>;  
  friend class SequenceList;

  friend class AlignmentHandle;
  friend class ConstSequenceList;
  
  /// \brief       Get residue index corresponding to given sequence position
  /// \param pos   zero-based index
  /// \throws      out_of_range if pos is not in [0, length-1].
  int GetResidueIndex(int pos) const;

  /// \brief Get zero-based index for given residue number.
  /// \throws out_of_range, if number is not included in sequence.
  int GetPos(int residue_index) const;

  /// \brief Get position of first non-gap character in sequence
  int GetFirstNonGap() const;
  
  /// \brief Get position of last non-gap character in sequence
  int GetLastNonGap() const;  
  
  /// \brief Get name of sequence
  /// 
  /// \sa SequenceHandle::GetName()
  const String& GetName() const;  
  
  /// \brief get sequence as a string, including all the gaps
  /// 
  /// \sa GetGaplessString()
  const String& GetString() const;
  
  /// \brief Get sequence as string ignoring gaps
  String GetGaplessString() const;

  /// \brief Get sequence offset from N-terminus
  ///
  /// \sa SequenceHandle::GetOffset
  int GetOffset() const;  
  
  /// \brief Get lenght of sequence, including gaps.
  int GetLength() const;

  /// \brief get one letter code of residue at position
  char GetOneLetterCode(int position) const;

  /// \brief get index of substring
  int GetIndex(const String& substr) const;
  
  /// \brief get residue at position
  ///
  /// will return the residue view at the given sequence position or an invalid
  /// residue view when no view is attached, the index is out of bounds or the
  /// position contains a gap.
  mol::ResidueView GetResidue(int position) const;

  /// \brief get attached view. may be an invalid entity view
  /// \sa SequenceHandle::AttachView(const mol::EntityView&, const String&)
  mol::EntityView GetAttachedView() const;

  /// \brief remove '.', '-' as gaps and make sequence all-uppercase
  /// Changes happen in place.
  void Normalise();

  /// \brief create copy sequence
  /// The newly created sequence has the same attached view.
  SequenceHandle Copy() const;
  
  /// \brief whether the sequence has an attached view
  /// 
  /// \sa SequenceHandle::AttachView(const mol::EntityView&, const String&)
  bool HasAttachedView() const;
  
  bool operator==(const SequenceHandle& rhs) const;
  bool operator!=(const SequenceHandle& rhs) const;  
  
  /// \brief whether the sequence is valid
  bool IsValid() const;
    
  /// \brief create invalid sequence handle
  /// 
  /// \sa IsValid()
  SequenceHandle();
  
  /// \brief set name of sequence
  void SetName(const String& name);

  /// \brief Set sequence String
  void SetString(const String& seq);

  /// \brief Set sequence offset
  ///
  /// By default the sequence offset is zero, i.e. the beginning of the sequence
  /// lies exactly at the N-terminus. Setting the sequence offset to a positive
  /// number will shift the sequence towards the C-terminus.
  void SetOffset(int offset);

  /// \brief attach entity view to sequence
  ///
  /// \throw IntegrityError when the view contains more than one chain
  void AttachView(const mol::EntityView& view);
  
  void SetOneLetterCode(int position, char new_char);
  
  char operator[](size_t index) const;
    
  iterator begin() const { return this->GetString().begin(); }
  iterator end() const { return this->GetString().end(); }
  
  void Append(char olc);

  operator ConstSequenceHandle() const;
  /// \brief attach entity view to sequence
  ///
  /// The sequence is mapped onto the chain with given name
  void AttachView(const mol::EntityView& view, const String& chain_name);
   
  /// \internal
  SequenceHandle(const impl::SequenceImplPtr& impl);  
  
  const String& GetRole() const;
  
  void SetRole(const String& role) const;
  
  impl::SequenceImplPtr& Impl() const;  

  GenericPropContainerImpl* GpImpl();

  const GenericPropContainerImpl* GpImpl() const;
private:
  void CheckValidity() const;  
  mutable impl::SequenceImplPtr impl_;
};

SequenceHandle DLLEXPORT_OST_SEQ CreateSequence(const String& name, 
                                                const String& seq, 
                                                const String& role="UNKNOWN");
#if(OST_INFO_ENABLED)                                                
/// \brief export sequence to info
void DLLEXPORT_OST_SEQ SequenceToInfo(const ConstSequenceHandle& sequence,
                                      info::InfoGroup& group);

/// \brief create sequence from info
SequenceHandle DLLEXPORT_OST_SEQ SequenceFromInfo(info::InfoGroup& group);
#endif


DLLEXPORT_OST_SEQ std::ostream& operator<<(std::ostream& os, 
                                           const ConstSequenceHandle& sequence);

bool DLLEXPORT_OST_SEQ Match(const ConstSequenceHandle& s1,
                             const ConstSequenceHandle& s2);

bool DLLEXPORT_OST_SEQ Match(const String& s1,
                             const String& s2);
}}

#endif
