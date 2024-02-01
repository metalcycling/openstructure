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
#ifndef OST_SEQ_ALIGNED_COLUMN_HH
#define OST_SEQ_ALIGNED_COLUMN_HH

/*
  Author: Marco Biasini
 */
#include <ost/seq/module_config.hh>
#include <ost/seq/alignment_handle.hh>
namespace ost { namespace seq {

/// \brief Provides access to a column in a \ref AlignedRegion "aligned region"
///     or a \ref AlignmentHandle "sequence alignment"
class DLLEXPORT_OST_SEQ AlignedColumn {
public:
  AlignedColumn(const AlignmentHandle& aln, int index);
  /// \brief construct empty (invalid) aligned column
  AlignedColumn();
  /// \brief get index of column
  /// 
  /// \sa AlignmentHandle::GetIndex()
  int GetIndex() const;
  
  char operator[](int row) const;
  
  /// \brief get residue index at row
  /// 
  /// \sa AlignmentHandle::GetResidueIndex()
  int GetResidueIndex(int row) const;
  
  /// \brief get number of rows
  /// 
  /// \sa AlignmentHandle::GetCount()
  int GetRowCount() const;
  
  /// \brief get residue at row
  /// 
  /// \sa AlignmentHandle::GetResidue()
  mol::ResidueView GetResidue(int row) const;
private:
  void CheckRowValidity(int row) const;

  AlignmentHandle aln_;
  int             index_;
};

DLLEXPORT_OST_SEQ std::ostream& operator<<(std::ostream& os, 
                                           const AlignedColumn& col);

}}

#endif
