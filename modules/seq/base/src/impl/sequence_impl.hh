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
#ifndef OST_SEQUENCE_HH
#define OST_SEQUENCE_HH

/*
  Author: Marco Biasini
 */
#include <ost/message.hh> 
#include <list>

#include <boost/shared_ptr.hpp>
#include <ost/generic_property.hh>
#include <ost/config.hh>
#if(OST_INFO_ENABLED)
#include <ost/info/info_fw.hh>
#endif
#include <ost/mol/residue_prop.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/residue_view.hh>

#include <ost/seq/module_config.hh>
#include <ost/seq/invalid_sequence.hh>
namespace ost { namespace seq { namespace impl {

/// \internal
class SequenceImpl;

/// \internal
typedef boost::shared_ptr<SequenceImpl> SequenceImplPtr;

/// \internal
class DLLEXPORT_OST_SEQ SequenceImpl : public GenericPropContainerImpl {
public:
  /// \brief       Construct new sequence object from sequence_string.
  static SequenceImplPtr FromString(const String& seq_name,
                                    const String& sequence_string,
                                    const String& role="UNKNOWN");

  /// \brief       Get residue index corresponding to given sequence position
  /// \param pos   zero-based index
  /// \throws      out_of_range if pos is not in [0, length-1].
  int GetResidueIndex(int pos) const;

  /// \brief       Get zero-based index for given residue number.
  /// \throws      out_of_range, if number is not included in sequence.
  int GetPos(int residue_index) const;

  /// \brief      Get position of first non-gap character in sequence
  int GetFirstNonGap() const;
  
  /// \brief      Get position of last non-gap character in sequence
  int GetLastNonGap() const;  
  
  const String& GetName() const;

  void SetName(const String& name);
  
  void Cut(int start, int n);

  void ShiftRegion(int start, int end, int amount);
  const String& GetString() const {
    return seq_string_;
  }
  /// \brief Set sequence String
  void SetString(const String& seq);

  /// \brief Remove everything 'unusual' from sequence
  void Normalise();

  /// \brief replace substring starting from start to end
  void Replace(const String& str,int start, int end);

  /// \brief Get sequence as String ignoring gaps
  String GetGaplessString() const;

  /// \brief Get sequence offset from N-terminus
  ///
  /// \sa #SetOffset
  int GetOffset() const;
  
  int GetIndex(const String& substr) const;
  /// \brief Set sequence offset
  ///
  /// By default the sequence offset is zero, i.e. the beginning of the sequence
  /// lies exactly at the N-terminus. Setting the sequence offset to a positive
  /// number will shift the sequence towards the C-terminus.
  void SetOffset(int offset);

  /// \brief Get lenght of sequence, including gaps.
  int GetLength() const;
  /// \brief       Construct new sequence object from sequence_string
  ///
  /// If you want to check whether the sequence String does only contain
  /// valid characters use \c CreateSequence instead.
  SequenceImpl(const String& seq_name, const String& sequence_string, 
               const String& role);

  /// \brief get one letter code of residue at position
  char GetOneLetterCode(int position) const;

  void SetOneLetterCode(int position, char new_char);
  /// \brief get residue at position
  ///
  /// will return the residue view at the given sequence position or an invalid
  /// residue view when no view is attached, the index is out of bounds or the
  /// position contains a gap.
  mol::ResidueView GetResidue(int position) const;

  /// \brief get attached view. may be an invalid entity view
  mol::EntityView GetAttachedView() const;

  /// \brief attach entity view to sequence
  ///
  /// \throw IntegrityError when the view contains more than one chain
  void AttachView(const mol::EntityView& view);

  /// \brief create copy sequence
  ///
  /// The newly created sequence has the same attached view.
  SequenceImplPtr Copy() const;
  /// \brief attach entity view to sequence
  ///
  /// The sequence is mapped onto the chain with given name
  void AttachView(const mol::EntityView& view, const String& chain_name);

  /// \brief whether the sequence has an attached view
  bool HasAttachedView() const;
  
  void Append(char olc);
  
  char& operator[](size_t index)
  {
    return seq_string_[index];
  }
  char operator[](size_t index) const
  {
    return seq_string_[index];
  }
  const String& GetRole() const
  {
    return seq_role_;
  }
  
  void SetRole(const String& role)
  {
    seq_role_=role;
  }
private:

  /// \brief       Recalculates gap shifts from sequence.
  ///
  /// Each hyphen in  sequence is interpreted as a gap. Lower-and-uppercase
  /// symbols are interpreted as amino acid residues. If other symbols are
  /// found in the String, an InvalidSequence exception is thrown.
  void ShiftsFromSequence();


  int GetPosNoBounds(int index) const;

  static bool IsSequenceStringSane(const String& seq_string);

  typedef struct {
    int start;
    int shift;
  } Shift;
  String              seq_name_;
  String              seq_string_;
  String              seq_role_;
  std::list<Shift>    shifts_;
  bool                editing_;
  int                 offset_;
  mol::EntityView     attached_view_;
};

/// \internal
typedef std::vector<SequenceImplPtr>  SequenceList;

#if(OST_INFO_ENABLED)
/// \internal
SequenceImplPtr DLLEXPORT_OST_SEQ SequenceImplFromInfo(const info::InfoGroup& group);

/// \internal
void DLLEXPORT_OST_SEQ SequenceImplToInfo(const SequenceImplPtr& sequence,
                                     info::InfoGroup& group);
#endif


}}} //ns
#endif

