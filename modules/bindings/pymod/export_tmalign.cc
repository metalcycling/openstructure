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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <ost/bindings/wrap_tmalign.hh>
using namespace boost::python;

ost::bindings::TMAlignResult WrapTMAlignPos(const geom::Vec3List& pos_one, 
                                            const geom::Vec3List& pos_two, 
                                            const ost::seq::SequenceHandle& seq1,
                                            const ost::seq::SequenceHandle& seq2,
                                            bool fast, bool rna) {

  return ost::bindings::WrappedTMAlign(pos_one, pos_two, seq1, seq2, fast, rna);
}

ost::bindings::TMAlignResult WrapTMAlignView(const ost::mol::ChainView& chain1,
                                             const ost::mol::ChainView& chain2, 
                                             bool fast) {

  return ost::bindings::WrappedTMAlign(chain1, chain2, fast);
}

void export_TMAlign() {
  class_<ost::bindings::TMAlignResult>("TMAlignResult", init<Real, Real, int, const geom::Mat4&, 
                                                             const ost::seq::AlignmentHandle&>())
    .add_property("rmsd", make_function(&ost::bindings::TMAlignResult::GetRMSD))
    .add_property("tm_score", make_function(&ost::bindings::TMAlignResult::GetTMScore))
    .add_property("aligned_length", make_function(&ost::bindings::TMAlignResult::GetAlignedLength))
    .add_property("transform", make_function(&ost::bindings::TMAlignResult::GetTransform,
                               return_value_policy<reference_existing_object>()))
    .add_property("alignment", make_function(&ost::bindings::TMAlignResult::GetAlignment,
                               return_value_policy<reference_existing_object>()))
  ;

  def("WrappedTMAlign", &WrapTMAlignPos, (arg("pos1"), arg("pos2"), arg("seq1"), arg("seq2"),
                                          arg("fast")=false, arg("rna")=false));

  def("WrappedTMAlign", &WrapTMAlignView, (arg("chain1"), arg("chain2"),
                                           arg("fast")=false));
}
