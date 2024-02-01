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

/*
  Author: Marco Biasini
 */
#include <ost/integrity_error.hh>
#include <ost/seq/aligned_region.hh>
#include <ost/seq/alignment_handle.hh>
namespace ost { namespace seq {

AlignedRegion::AlignedRegion(const AlignmentHandle& aln, 
                             int start, int end, int master):
   aln_(aln), start_(start), end_(end), master_(master)
{
}


int AlignedRegion::GetEnd() const
{
 return end_; 
}


int AlignedRegion::GetStart() const
{
  return start_;
}


void AlignedRegion::Delete()
{
  aln_.Cut(start_, end_);
  end_=start_;
}

void AlignedRegion::Replace(const AlignedRegion& aln_r)
{
  aln_.Replace(aln_r,start_, end_);
  end_=start_+end_;
}

void AlignedRegion::ShiftLeft(int n)
{
  if (start_-n<0) {
    std::stringstream s;
    s << "Can't shift aligned region" << this->GetStart() << ":" 
      << this->GetEnd() << " by " << n << " characters to the left";
    throw IntegrityError(s.str());
  }
  aln_.ShiftRegion(start_, end_, -n, master_);
  start_-=n;
  end_-=n;
}

void AlignedRegion::ShiftRight(int n)
{
  if (end_+n>aln_.GetLength()) {
    std::stringstream s;
    s << "Can't shift aligned region" << this->GetStart() << ":" 
      << this->GetEnd() << " by " << n << " characters to the right";
    throw IntegrityError(s.str());
  }
  aln_.ShiftRegion(start_, end_, n, master_);
  start_+=n;
  end_+=n;
}

AlignedColumnIterator AlignedRegion::begin()
{
  return AlignedColumnIterator(aln_, start_, end_);
}

AlignedColumnIterator AlignedRegion::end()
{
  return AlignedColumnIterator(aln_, end_, end_);
}

AlignedColumn AlignedRegion::operator[](int index) const
{
  int col_index = start_+index;
  if(col_index<start_ || col_index>=end_){
    throw std::out_of_range("Index out of region");
  }
  return AlignedColumn(aln_, col_index);
}

bool AlignedRegion::operator==(const AlignedRegion& rhs) const
{
  return aln_==rhs.aln_ && rhs.start_==start_ && 
         end_==rhs.end_ && master_==rhs.master_;
}

bool AlignedRegion::operator!=(const AlignedRegion& rhs) const
{
  return !(*this==rhs);
}

int AlignedRegion::GetMaster() const
{
  return master_;
}

void AlignedRegion::SetMaster(int master)
{
  if(master < -1 || master >= aln_.GetCount()){
    throw IntegrityError("Master out of bounds");
  }
  master_=master;
}

int AlignedRegion::GetLength() const
{
  return end_-start_;
}

AlignmentHandle AlignedRegion::GetAlignmentHandle() const
{
  return aln_;
}

}}
