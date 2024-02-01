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
#include <boost/python.hpp>
#include <boost/python/slice.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


#include <ost/export_helper/pair_to_tuple_conv.hh>
#include <ost/generic_property.hh>
#include <ost/export_helper/generic_property_def.hh>
#include <ost/config.hh>
#if(OST_INFO_ENABLED)
#include <ost/info/info.hh>
#endif
#include <ost/mol/mol.hh>
#include <ost/seq/sequence_handle.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/sequence_op.hh>
#include <ost/seq/sequence_list.hh>
#include <ost/seq/aligned_region.hh>
#include <ost/seq/views_from_sequences.hh>
#include "const_seq_list_export_def.hh"

using namespace ost;
using namespace ost::seq;
using namespace boost::python;

void (SequenceHandle::*attach_one)(const mol::EntityView&)=&SequenceHandle::AttachView;
void (SequenceHandle::*attach_two)(const mol::EntityView&,
                             const String&)=&SequenceHandle::AttachView;

void (AlignmentHandle::*attach_view_a)(int, const mol::EntityView&)=&AlignmentHandle::AttachView;
void (AlignmentHandle::*attach_view_b)(int, const mol::EntityView&,
                             const String&)=&AlignmentHandle::AttachView;
SequenceHandle (*seq_from_chain_a)(const String&,const mol::ChainView&)=&SequenceFromChain;
SequenceHandle (*seq_from_chain_b)(const String&,const mol::ChainHandle&)=&SequenceFromChain;

bool (*m1)(const String&, const String&)=&Match;
bool (*m2)(const ConstSequenceHandle&, const ConstSequenceHandle&)=&Match;
template <typename T>
T do_slice(const T& t, slice& sl) {
  int start=0, end=t.GetCount();
  try {
    start=extract<int>(sl.start());
    if (start<0) {
      start=t.GetCount()+start;
    }    
  } catch(error_already_set&) {
    PyErr_Clear();
  }
  try {
    end=extract<int>(sl.stop());
    if (end<0) {
      end=t.GetCount()+end;
    }    
  } catch(error_already_set&) {
    PyErr_Clear();
  }
  return t.Slice(start, end-start);
}

ConstSequenceList do_slice_a(ConstSequenceList& t, slice sl)
{
  return do_slice<ConstSequenceList>(t, sl);
}
String aln_to_str(const AlignmentHandle& aln) {
  return aln.ToString();
}

SequenceList do_slice_b(SequenceList& t, slice sl)
{
  return do_slice<SequenceList>(t, sl);
}

String slice_seq(const ConstSequenceHandle& sh, slice& sl) {
  int start=0, end=sh.GetLength();
  try {
    start=extract<int>(sl.start());
    if (start<0) {
      start=sh.GetLength()+start;
    }    
  } catch(error_already_set&) {
    PyErr_Clear();
  }
  try {
    end=extract<int>(sl.stop());
    if (end<0) {
      end=sh.GetLength()+end;
    }
  } catch(error_already_set&) {
    PyErr_Clear();
  }
  String s=sh.GetString();
  return s.substr(start, end-start);
}

AlignedRegion slice_aln(const AlignmentHandle& aln, slice sl) {
  int start=0, end=aln.GetLength();
  try {
    start=extract<int>(sl.start());
    if (start<0) {
      start=aln.GetLength()+start;
    }    
  } catch(error_already_set&) {
    PyErr_Clear();
  }
  try {
    end=extract<int>(sl.stop());
    if (end<0) {
      end=aln.GetLength()+end;
    }
  } catch(error_already_set&) {
    PyErr_Clear();
  }
  return aln.MakeRegion(start, end-start);
}


struct RevRegionRangeIter {
  RevRegionRangeIter(AlignedColumnIterator b,
                     AlignedColumnIterator e):
    b_(b), e_(e)
  { }

  AlignedColumn next()
  {
    if (b_==e_) {
      boost::python::objects::stop_iteration_error();
    }
    --e_;
    AlignedColumn col=*e_;

    return col;
  }
private:
  AlignedColumnIterator b_;
  AlignedColumnIterator e_;
};

struct RegionRangeIter {
  RegionRangeIter(AlignedColumnIterator b,
                  AlignedColumnIterator e):
    b_(b), e_(e)
  { }

  AlignedColumn next()
  {
    if (b_==e_) {
      boost::python::objects::stop_iteration_error();
    }
    AlignedColumn col=*b_;
    ++b_;
    return col;
  }
private:
  AlignedColumnIterator b_;
  AlignedColumnIterator e_;
};

struct ConstSeqListIter {
  ConstSeqListIter(ConstSequenceList& list):
    l_(list), b_(l_.begin()), e_(l_.end())
  { }

  ConstSequenceHandle next()
  {
    if (b_==e_) {
      boost::python::objects::stop_iteration_error();
    }
    ConstSequenceHandle s=*b_;  
    ++b_;
    return s;
  }
private:
  ConstSequenceList           l_;
  ConstSequenceList::iterator b_;
  ConstSequenceList::iterator e_;
};

struct SeqListIter {
  SeqListIter(SequenceList& list):
    l_(list), b_(l_.begin()), e_(l_.end())
  { }

  SequenceHandle next()
  {
    if (b_==e_) {
      boost::python::objects::stop_iteration_error();
    }
    SequenceHandle s=*b_;
    ++b_;
    return s;
  }
private:
  SequenceList           l_;
  SequenceList::iterator b_;
  SequenceList::iterator e_;
};

ConstSeqListIter iter_cs(ConstSequenceList& sl) 
{
  return ConstSeqListIter(sl);
}

SeqListIter iter_sl(SequenceList& sl) 
{
  return SeqListIter(sl);
}

RegionRangeIter iter_range1(AlignmentHandle& aln)
{
  return RegionRangeIter(aln.begin(), aln.end());
}

RegionRangeIter iter_range2(AlignedRegion& aln_region)
{
  return RegionRangeIter(aln_region.begin(), aln_region.end());
}

RevRegionRangeIter iter_range_rev(AlignedRegion& aln_region)
{
  return RevRegionRangeIter(aln_region.begin(), aln_region.end());
}

template <typename C, typename O>
void const_seq_handle_def(O& bp_class)
{
  bp_class
    .def("GetResidueIndex", &C::GetResidueIndex)
    .def("GetPos", &C::GetPos)
    .def("GetLength", &C::GetLength)
    .def("GetResidue", &C::GetResidue)
    .def("GetOneLetterCode", &C::GetOneLetterCode)
    .def("__iter__", iterator<C>())
    .def("__getitem__", &C::GetOneLetterCode)
    .def("__getitem__", slice_seq)
    .def("GetOffset", &C::GetOffset)
    .def("Copy", &C::Copy)
    .def("IsValid", &C::IsValid)
    .def("GetFirstNonGap", &C::GetFirstNonGap)
    .def("GetLastNonGap", &C::GetLastNonGap)
    .add_property("first_non_gap", &C::GetFirstNonGap)
    .add_property("last_non_gap", &C::GetLastNonGap)
    .def("GetAttachedView", &C::GetAttachedView)
    .def("GetGaplessString", &C::GetGaplessString)
    .add_property("role", make_function(&C::GetRole, 
                                        return_value_policy<copy_const_reference>()))
    .def("GetString", &C::GetString,
         return_value_policy<copy_const_reference>())
         .def("GetName", &C::GetName,
              return_value_policy<copy_const_reference>())
    .def("HasAttachedView", &C::HasAttachedView)
    .def("__len__", &C::GetLength)
    .add_property("length", &C::GetLength)
    .add_property("attached_view", &C::GetAttachedView)
    .add_property("name",
                  make_function(&C::GetName,
                                return_value_policy<copy_const_reference>()))
    .add_property("offset", &C::GetOffset)
    .add_property("gapless_string", &C::GetGaplessString)
    .add_property("string",
                  make_function(&C::GetString,
                                return_value_policy<copy_const_reference>()))
    .def("__str__", &C::GetString,
         return_value_policy<copy_const_reference>())
  ;
}


void export_sequence()
{
  class_<ConstSequenceHandle> const_seq("ConstSequenceHandle", init<>());
  const_seq_handle_def<ConstSequenceHandle>(const_seq);
  const_generic_prop_def<ConstSequenceHandle>(const_seq);
  class_<SequenceHandle> seq_handle("SequenceHandle", init<>());
  const_seq_handle_def<SequenceHandle>(seq_handle);
  generic_prop_def<SequenceHandle>(seq_handle);  
  seq_handle
    .def("__setitem__", &SequenceHandle::SetOneLetterCode)
    .def("SetOneLetterCode", &SequenceHandle::SetOneLetterCode)  
    .def("SetOffset", &SequenceHandle::SetOffset)
    .def("AttachView", attach_one)
    .def("AttachView", attach_two)
    .def("Append", &SequenceHandle::Append)
    .def("SetString", &SequenceHandle::SetString)
    .def("Normalise", &SequenceHandle::Normalise)
    .add_property("string",
                  make_function(&SequenceHandle::GetString,
                                return_value_policy<copy_const_reference>()),
                  &SequenceHandle::SetString)
    .add_property("role", make_function(&SequenceHandle::GetRole, 
                                        return_value_policy<copy_const_reference>()),
                  &SequenceHandle::SetRole)
    .def("SetName", &SequenceHandle::SetName)
    .add_property("name",
                  make_function(&SequenceHandle::GetName,
                                return_value_policy<copy_const_reference>()),
                  &SequenceHandle::SetName)
    .add_property("offset", &SequenceHandle::GetOffset,
                  &SequenceHandle::SetOffset)
  ;
  def("Match", m1);
  def("Match", m2);

  implicitly_convertible<SequenceHandle, ConstSequenceHandle>();
  
  def("CreateSequence", &CreateSequence, 
      (arg("name"), arg("seq"), arg("role")="UNKNOWN"));
  /*class_<SequenceHandleList>("SequenceHandleList", init<>())
    .def(vector_indexing_suite<SequenceHandleList>())
  ;*/
  class_<RegionRangeIter>("RegionRangeIter", no_init)
    .def("next", &RegionRangeIter::next)
    .def("__next__", &RegionRangeIter::next)
  ;
    class_<RevRegionRangeIter>("RevRegionRangeIter", no_init)
    .def("next", &RevRegionRangeIter::next)
    .def("__next__", &RevRegionRangeIter::next)
  ;
  class_<ConstSeqListIter>("ConstSeqListIter", no_init)
    .def("next", &ConstSeqListIter::next)
    .def("__next__", &ConstSeqListIter::next)
  ;
  class_<SeqListIter>("SeqListIter", no_init)
    .def("next", &SeqListIter::next)
    .def("__next__", &SeqListIter::next)
  ;
  class_<AlignmentHandle>("AlignmentHandle", init<>())
    .def("GetCount", &AlignmentHandle::GetCount)
    .add_property("sequence_count", &AlignmentHandle::GetCount)
    .def("GetSequence", &AlignmentHandle::GetSequence)
    .def("GetPos", &AlignmentHandle::GetPos)
    .def("GetResidueIndex", &AlignmentHandle::GetResidueIndex)
    .def("GetResidue", &AlignmentHandle::GetResidue)
    .def("AddSequence", &AlignmentHandle::AddSequence)
    .def("GetMatchingBackboneViews", &AlignmentHandle::GetMatchingBackboneViews,
         (arg("idx_a")=0, arg("idx_b")=1))
    .def("FindSequence", &AlignmentHandle::FindSequence)
    .def("FindSequenceIndex", &AlignmentHandle::FindSequenceIndex)
    .def("Copy", &AlignmentHandle::Copy)
    .def("ToString", &AlignmentHandle::ToString, (arg("width")=80))
    .def("__str__", aln_to_str)
    .def("GetLength", &AlignmentHandle::GetLength)
    .def("__len__", &AlignmentHandle::GetLength)
    .def("GetSequences", &AlignmentHandle::GetSequences)
    .def("GetCoverage", &AlignmentHandle::GetCoverage)
    .def("AttachView", attach_view_a)
    .def("SetSequenceRole", &AlignmentHandle::SetSequenceRole)
    .def("GetSequenceRole", &AlignmentHandle::GetSequenceRole, 
         return_value_policy<copy_const_reference>())
    .def("AttachView", attach_view_b)
    .def("Cut", &AlignmentHandle::Cut)
    .def("MakeRegion", &AlignmentHandle::MakeRegion)
    .def("RemoveSequence", &AlignmentHandle::RemoveSequence)
    .def("Replace",&AlignmentHandle::Replace)
    .def("__getitem__", &slice_aln)
    .def("__getitem__", &AlignmentHandle::operator[])
    .def("__iter__", iter_range1)
    .add_property("sequences", &AlignmentHandle::GetSequences)
    .def("SetSequenceName",  &AlignmentHandle::SetSequenceName)
    .def("SetSequenceOffset", &AlignmentHandle::SetSequenceOffset)
    .def("GetSequenceOffset", &AlignmentHandle::GetSequenceOffset)
  ;
  class_<AlignedColumn>("AlignedColumn", no_init)
    .def("GetIndex", &AlignedColumn::GetIndex)
    .def("__getitem__", &AlignedColumn::operator[])
    .def("GetRowCount", &AlignedColumn::GetRowCount)
    .def("GetResidue", &AlignedColumn::GetResidue)
    .def(self_ns::str(self))
  ;
  class_<AlignedRegion>("AlignedRegion", no_init)
    .def("Delete", &AlignedRegion::Delete)
    .def("Replace", &AlignedRegion::Replace)
    .def("ShiftLeft", &AlignedRegion::ShiftLeft)
    .def("ShiftRight", &AlignedRegion::ShiftRight)
    .def("GetStart", &AlignedRegion::GetStart)
    .def("GetEnd", &AlignedRegion::GetEnd)
    .def("GetLength", &AlignedRegion::GetLength)
    .def("SetMaster", &AlignedRegion::SetMaster)
    .def("GetMaster", &AlignedRegion::GetMaster)
    .def("GetAlignmentHandle",&AlignedRegion::GetAlignmentHandle)
    .add_property("master", &AlignedRegion::GetMaster,
                  &AlignedRegion::SetMaster)
    .def("__getitem__", &AlignedRegion::operator[])
    .def("__len__", &AlignedRegion::GetLength)
    .def("__iter__", iter_range2)
    .add_property("start", &AlignedRegion::GetStart)
    .add_property("end", &AlignedRegion::GetEnd)
  ;
  class_<ConstSequenceList>("ConstSequenceList", init<>())
    CONST_SEQ_LIST_DEF(ConstSequenceList)
    .def("__getitem__", &do_slice_a)
    .def("__iter__", &iter_cs)    
  ;
  class_<SequenceList>("SequenceList", init<>())
    CONST_SEQ_LIST_DEF(SequenceList)
    .def("__getitem__", &do_slice_b)
    .def("__iter__", &iter_sl)
  ;
  class_<AlignmentList>("AlignmentList", init<>())
    .def(vector_indexing_suite<AlignmentList>())
    .def("__getitem__", &do_slice_b)
  ;
  implicitly_convertible<SequenceList, ConstSequenceList>();

  def("CreateSequenceList", &CreateSequenceList);
  def("SequenceFromChain", seq_from_chain_a);
  def("SequenceFromChain", seq_from_chain_b);
#if(OST_INFO_ENABLED)
  def("SequenceToInfo", &SequenceToInfo);
  def("SequenceListToInfo", &SequenceListToInfo);
  def("SequenceFromInfo", &SequenceFromInfo);  
  def("SequenceListFromInfo", &SequenceListFromInfo);
#endif
  def("ViewsFromSequences", &ViewsFromSequences, (arg("seq1"), arg("seq2")));
  def("ViewsFromAlignment", &ViewsFromAlignment, 
      (arg("aln"), arg("index1")=0, arg("index2")=1));
  def("CreateAlignment", &CreateAlignment);
  def("AlignmentFromSequenceList", &AlignmentFromSequenceList);
}
