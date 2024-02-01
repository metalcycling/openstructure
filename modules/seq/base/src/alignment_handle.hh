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
#ifndef OST_SEQ_ALIGNMENT_HANDLE_HH
#define OST_SEQ_ALIGNMENT_HANDLE_HH

/*
  Author: Marco Biasini
 */
#include <ost/seq/module_config.hh>

#include <ost/seq/impl/sequence_list_impl_fw.hh>
#include <ost/seq/sequence_handle.hh>
#include <ost/seq/sequence_list.hh>

namespace ost { namespace seq {

class SequenceList;
class AlignedRegion;
class AlignedColumn;
class AlignedColumnIterator;

/// \brief representation of a multiple sequence alignemnt consisting of two or
///     more sequences
///
/// A sequence alignment consists of two or more sequences. The number of
/// sequences in the alignment can be obtained by #GetCount().
/// All sequences are of length #GetLength().
///
/// Typically sequence alignments are used column-based, i.e by looking at an
/// aligned columns in the sequence alignment. To get a row-based (sequence)
/// view on the sequence list, use AlignmentHandle::GetSequences(). For an
/// overview of how to use the sequence module, see \ref module_seq "here"
///
/// All operators that operate on an alignment will again  produce a valid
/// alignment. This mean that it is not possible to change the length of one
/// sequence, without  adjusting the other sequences, too.
class DLLEXPORT_OST_SEQ AlignmentHandle {
public:

  typedef AlignedColumnIterator iterator;
  AlignmentHandle();

  /// \brief  Get position in sequence with index seq_index that corresponds to
  ///         the given residue index.
  ///
  ///  Identical to \c msa->GetSequence(seq_id)->GetPos(index)
  int GetPos(int seq_index, int residue_index) const;

  /// \brief  Get residue index for seq_index and alignment position.
  ///
  ///  Identitical to \c msa->GetSequence(seq_id)->GetResidueIndex(pos)
  int GetResidueIndex(int seq_index, int pos) const;

  mol::ResidueView GetResidue(int seq_index, int pos) const;

  char GetOneLetterCode(int seq_index, int pos) const;

  /// \brief     Add new sequence to multiple sequence alignment.
  ///
  /// If the sequence length does not match with the length of the other
  /// sequences in the alignment, an InvalidSequence() exception is thrown.
  void AddSequence(const ConstSequenceHandle& sequence);

  /// \brief     Get sequence with given index.
  /// \return    sequence or invalid handle if the index is out of bounds
  ConstSequenceHandle GetSequence(int seq_id) const;
  
  /// \brief remove sequence at seq_index from alignment
  void RemoveSequence(int seq_index);
  /// \brief   Convert multiple sequence alignment to string.
  String ToString(int width=80) const;

  /// \brief   Get length of sequences in alignment
  int GetLength() const;

  /// \brief deep-copy multi sequence alignment
  AlignmentHandle Copy() const;

  /// \brief find sequence by name.
  ///
  /// If several sequences have the same name, the first matching sequence will
  /// be returned.
  ConstSequenceHandle FindSequence(const String& name) const;
  
  int FindSequenceIndex(const String& name) const;
  
  /// \brief attach view to given sequence
  /// \sa SequenceHandle::AttachView(const mol::EntityView&)
  void AttachView(int seq_index, const mol::EntityView& view);

  /// \brief attach view to given sequence
  /// \sa SequenceHandle::AttachView(const mol::EntityView&, const String&)
  void AttachView(int seq_index, const mol::EntityView& view,
                  const String& chain_name);


  mol::EntityViewPair GetMatchingBackboneViews(int idx0=0, int idx1=1) const;
  /// \brief set name of sequence
  void SetSequenceName(int seq_index, const String& name);

  void SetSequenceOffset(int seq_index, int offset); 
  int  GetSequenceOffset(int seq_index);
  /// \brief Get list of sequences (read-only)
  ConstSequenceList GetSequences() const;

  /// \brief create an aligned region.
  ///
  /// \param start is the index of the first column
  /// \param n is the length of the sequence
  /// \param master is the reference system for operations such as shifting.
  ///      If set to -1, no master sequence is defined and the operations will
  ///      affect all sequences
  ///
  /// If the aligned region is out of bounds, a std::out_of_bounds exeception will be thrown.
  AlignedRegion MakeRegion(int start, int n, int master=-1) const;

  /// \brief get number of sequences in alignment
  int GetCount() const;

  bool operator==(const AlignmentHandle& rhs) const;
  bool operator!=(const AlignmentHandle& rhs) const;

  /// \brief get aligned column at index
  ///
  /// This method does not throw any exception. Upon accessing methods of the
  /// aligned column, exceptions might be thrown when the index is out of
  /// bounds.
  AlignedColumn operator[](int index) const;

  AlignmentHandle(const impl::SequenceListImplPtr& impl);

  /// \brief cut out half-closed interval start, end
  void Cut(int start, int end);

  ///\brief Replace part of an alignment
  void Replace(const AlignedRegion& rhs, int start, int end);
  /// \brief shift half-closed interval by amount
  ///
  /// if master is -1, all sequences of the alignment will be shifted. Otherwise
  /// only the sequence with given index is affected.
  void ShiftRegion(int start, int end, int amount, int master=-1);

  /// \brief Column iterator start-point
  iterator begin() const;
  /// \brief Column iterator end-point
  iterator end() const;

  bool IsValid() const { return impl_.get()!=0; }
    
  ///\brief get coverage of a specifi sequence 
  ///
  /// returns a value representing how extensively the specified sequence
  /// covers the first sequence (sequence 0). The function return a value
  /// between 0 (no coverage) and 1 (full coverage)
  Real GetCoverage(int seq_index) const;
  
  
  const String& GetSequenceRole(int seq_index);
  
  void SetSequenceRole(int seq_index, const String& role);
private:
  void CheckValidity() const;
  impl::SequenceListImplPtr impl_;
};

AlignmentHandle DLLEXPORT_OST_SEQ CreateAlignment();

/// \brief convert alignment from sequence list
///
/// If the sequences in the SequenceList have different lengths, an
/// InvalidAlignment exception is thrown.
///
/// \return alignment consisting of the sequences in seq_list.
AlignmentHandle DLLEXPORT_OST_SEQ
AlignmentFromSequenceList(const SequenceList& seq_list);

typedef std::vector<AlignmentHandle> AlignmentList;

}}

#endif
