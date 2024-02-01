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
#ifndef OST_SEQ_ALIGNED_REGION_HH
#define OST_SEQ_ALIGNED_REGION_HH

/*
  Author: Marco Biasini
 */
#include <ost/seq/module_config.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/aligned_column_iterator.hh>

namespace ost { namespace seq {


/// \brief represents a region in a \ref AlignmentHandle "sequence alignment"
/// 
/// The region covers the half-closed interval start/stop.
/// 
/// You should never need to construct an aligned region manually. Use one of 
/// the helper methods of AlignmentHandle.
class DLLEXPORT_OST_SEQ AlignedRegion {
public:
  typedef AlignedColumnIterator iterator;
  AlignedRegion(const AlignmentHandle& aln, 
                int start, int end, int master=-1);
  
  /// \brief get end of aligned region
  int GetEnd() const;
  /// \brief get start of aligned region
  int GetStart() const;
  int GetMaster() const;
  void SetMaster(int master);
  /// \brief delete interval and update length of AlignedRegion
  void Delete();
  /// \brief replace region with content of AlignedRegion and set length to 
  ///       length of \p aln_r
  void Replace(const AlignedRegion& aln_r);
  
  /// \brief shift the aligned region of the master sequence to the left by
  ///        n characters.
  void ShiftLeft(int n);
  
  /// \brief shift the aligned region of the master sequence to the right by
  ///        n characters.  
  /// \sa
  void ShiftRight(int n);
  
  /// \brief length of the aligned region
  /// 
  /// \sa AlignmentHandle::GetLength()
  int GetLength() const;
  
  ///\brief retrieve alignment handle for aligned region
  AlignmentHandle GetAlignmentHandle() const;

  /// \brief get aligned column at given index
  /// 
  /// The indices range from 0 to GetLength()-1.
  AlignedColumn operator[](int index) const;
  
  /// \brief iterator pointing at start of aligned region
  AlignedColumnIterator begin();
  /// \brief iterator pointing at end of aligned region.
  AlignedColumnIterator end();
  
  bool operator==(const AlignedRegion& rhs) const;
  bool operator!=(const AlignedRegion& rhs) const;  
private:
  AlignmentHandle aln_;
  int             start_;
  int             end_ ;
  int             master_;
};

typedef std::vector<AlignedRegion> AlignedRegionList;

}}

#endif
