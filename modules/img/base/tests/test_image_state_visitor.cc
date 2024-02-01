//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  Author: Ansgar Philippsen
*/

#include "test_image_state_visitor.hh"

#include <ost/img/image_state.hh>

using namespace ost::img;
using namespace ost::img::image_state;

namespace test_image_state_visitor {

struct MyConstVisitorFnc
{
  template <typename V, class DP>
  void VisitState(const ImageStateImpl<V,DP>& isi) {
    BOOST_CHECK(true);
  }
};

typedef ImageStateNonModVisitor<MyConstVisitorFnc> MyConstVisitor;

struct MyIPVisitorFnc
{
  template <typename V, class DP>
  void VisitState(ImageStateImpl<V,DP>& isi) {
    BOOST_CHECK(true);
  }
};
typedef ImageStateModIPVisitor<MyIPVisitorFnc> MyIPVisitor;


struct MyOPVisitorFnc
{
  template <typename V, class DP>
  ImageStateBasePtr VisitState(const ImageStateImpl<V,DP>& isi) {
    BOOST_CHECK(true);
    return ImageStateBasePtr();
  }
};
typedef ImageStateModOPVisitor<MyOPVisitorFnc> MyOPVisitor;


void test_all()
{
  RealSpatialImageState rs_is;
  RealFrequencyImageState rf_is;
  ComplexSpatialImageState cs_is;
  ComplexFrequencyImageState cf_is;
  ComplexHalfFrequencyImageState hf_is;
  
  MyConstVisitor const_visitor;
  MyIPVisitor ip_visitor;
  MyOPVisitor op_visitor;

  rs_is.ApplyIP(const_visitor);
  rf_is.ApplyIP(const_visitor);
  rs_is.ApplyIP(const_visitor);
  rf_is.ApplyIP(const_visitor);
  hf_is.ApplyIP(const_visitor);

  rs_is.ApplyIP(ip_visitor);
  rf_is.ApplyIP(ip_visitor);
  rs_is.ApplyIP(ip_visitor);
  rf_is.ApplyIP(ip_visitor);
  hf_is.ApplyIP(ip_visitor);

  rs_is.Apply(op_visitor);
  rf_is.Apply(op_visitor);
  rs_is.Apply(op_visitor);
  rf_is.Apply(op_visitor);
  hf_is.Apply(op_visitor);
}

} // namespace

test_suite* CreateImageStateVisitorTest()
{
  using namespace test_image_state_visitor;
  test_suite* ts=BOOST_TEST_SUITE("ImageStateVisitor Test");

  ts->add(BOOST_TEST_CASE(&test_all));

  return ts;
}
