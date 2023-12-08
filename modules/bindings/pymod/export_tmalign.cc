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

ost::bindings::MMAlignResult WrapMMAlignView(const ost::mol::EntityView& ent1,
                                             const ost::mol::EntityView& ent2, 
                                             bool fast,
                                             boost::python::dict& mapping) {
  boost::python::list keys(mapping.keys());
  boost::python::list values(mapping.values());
  std::map<String, String> m_mapping;
  for(uint i = 0; i < boost::python::len(keys); ++i) {
    m_mapping[boost::python::extract<String>(keys[i])] =
    boost::python::extract<String>(values[i]);
  }

  return ost::bindings::WrappedMMAlign(ent1, ent2, fast, m_mapping);
}

void export_TMAlign() {
  class_<ost::bindings::TMAlignResult>("TMAlignResult", init<Real, Real, Real, int, const geom::Mat4&, 
                                                             const ost::seq::AlignmentHandle&>())
    .add_property("rmsd", make_function(&ost::bindings::TMAlignResult::GetRMSD))
    .add_property("tm_score", make_function(&ost::bindings::TMAlignResult::GetTMScore))
    .add_property("tm_score_swapped", make_function(&ost::bindings::TMAlignResult::GetTMScoreSwapped))
    .add_property("aligned_length", make_function(&ost::bindings::TMAlignResult::GetAlignedLength))
    .add_property("transform", make_function(&ost::bindings::TMAlignResult::GetTransform,
                               return_value_policy<reference_existing_object>()))
    .add_property("alignment", make_function(&ost::bindings::TMAlignResult::GetAlignment,
                               return_value_policy<reference_existing_object>()))
  ;

  class_<ost::bindings::MMAlignResult>("MMAlignResult", init<Real, Real, Real, int, const geom::Mat4&,
                                                             const ost::seq::AlignmentList&,
                                                             const std::vector<String>&,
                                                             const std::vector<String>&>())
    .add_property("rmsd", make_function(&ost::bindings::MMAlignResult::GetRMSD))
    .add_property("tm_score", make_function(&ost::bindings::MMAlignResult::GetTMScore))
    .add_property("tm_score_swapped", make_function(&ost::bindings::MMAlignResult::GetTMScoreSwapped))
    .add_property("transform", make_function(&ost::bindings::MMAlignResult::GetTransform,
                               return_value_policy<reference_existing_object>()))
    .add_property("aligned_length", make_function(&ost::bindings::MMAlignResult::GetAlignedLength))
    .add_property("alignments", make_function(&ost::bindings::MMAlignResult::GetAlignments,
                               return_value_policy<reference_existing_object>()))
    .add_property("ent1_mapped_chains", make_function(&ost::bindings::MMAlignResult::GetEnt1MappedChains,
                               return_value_policy<reference_existing_object>()))
    .add_property("ent2_mapped_chains", make_function(&ost::bindings::MMAlignResult::GetEnt2MappedChains,
                               return_value_policy<reference_existing_object>()))
  ;


  def("WrappedTMAlign", &WrapTMAlignPos, (arg("pos1"), arg("pos2"), arg("seq1"), arg("seq2"),
                                          arg("fast")=false, arg("rna")=false));

  def("WrappedTMAlign", &WrapTMAlignView, (arg("chain1"), arg("chain2"),
                                           arg("fast")=false));

  def("WrappedMMAlign", &WrapMMAlignView, (arg("ent1"), arg("ent2"),
                                           arg("fast")=false,
                                           arg("mapping")=boost::python::dict()));
}
