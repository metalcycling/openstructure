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
#include <ost/seq/aligned_column.hh>
namespace ost { namespace seq {
  
AlignedColumn::AlignedColumn(const AlignmentHandle& aln, int index):
  aln_(aln), index_(index)
{
  
}

AlignedColumn::AlignedColumn():
  aln_(), index_(0)
{
}

/// \brief get index of column
int AlignedColumn::GetIndex() const
{
  return index_;
}

char AlignedColumn::operator[](int row) const
{
  this->CheckRowValidity(row);
  return aln_.GetOneLetterCode(row, index_);
}

int AlignedColumn::GetResidueIndex(int row) const
{
  this->CheckRowValidity(row);
  return aln_.GetResidueIndex(row, index_);
}

mol::ResidueView AlignedColumn::GetResidue(int row) const
{
  this->CheckRowValidity(row);
  return aln_.GetResidue(row, index_);
}

std::ostream& operator<<(std::ostream& os, const AlignedColumn& col)
{
  for (int i=0; i<col.GetRowCount(); ++i) {
    os << col[i];
  }
  return os;
}

int AlignedColumn::GetRowCount() const
{
  return aln_.GetCount();
}

void AlignedColumn::CheckRowValidity(int row) const
{
  if(row<0 || row>=this->GetRowCount()){
    throw std::out_of_range("Row out of bounds");
  }
}


}}
