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
#ifndef OST_SEQ_ALIGNED_COLUMN_ITERATOR_HH
#define OST_SEQ_ALIGNED_COLUMN_ITERATOR_HH

#include <ost/seq/aligned_column.hh>

/*
  Author: Marco Biasini
 */
namespace ost { namespace seq {

class AlignedColumnIterator : public std::iterator<std::forward_iterator_tag,
                                                   AlignedColumn> {
private:
  void UpdateVal()
  {
    if (curr_<end_) {
      val_=AlignedColumn(aln_, curr_);
    }
  }
public:
  AlignedColumnIterator(const AlignmentHandle& aln, int start, int end):
     aln_(aln), curr_(start), end_(end)
  {
    this->UpdateVal();
  }

  AlignedColumnIterator& operator++()
  {
    ++curr_;
    this->UpdateVal();
    return *this;
  }

  AlignedColumnIterator& operator--()
  {
    --curr_;
    this->UpdateVal();
    return *this;
  }

  AlignedColumnIterator operator++(int)
  {
    AlignedColumnIterator ans(*this);
    ++(*this);
    return ans;
  }

  AlignedColumn* operator->()
  {
    return &val_;
  }

  AlignedColumn& operator*()
  {
    return val_;
  }

  bool operator==(const AlignedColumnIterator& rhs) const
  {
    return aln_==rhs.aln_ && rhs.curr_==curr_;
  }

  bool operator!=(const AlignedColumnIterator& rhs) const
  {
    return !(*this==rhs);
  }
private:
  AlignmentHandle aln_;
  int             curr_;
  int             end_;
  AlignedColumn   val_;
};
}}

#endif
