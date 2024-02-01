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

#include <sstream>

#include <boost/mpl/vector.hpp>


#include "test_image_state.hh"

#include <ost/img/extent.hh>
#include <ost/img/extent_iterator.hh>
#include <ost/img/size.hh>
#include <ost/img/image_state/image_state_impl.hh>
#include <ost/img/image_state/image_state_factory.hh>
#include <ost/img/image_state/image_state_def.hh>

#include <ost/img/alg/randomize.hh>

using namespace ost::img;
using namespace ost::img::image_state;

namespace test_image_state {

alg::Randomize randomizer;

bool is_clear(const ImageStateBase& b)
{
  for(ExtentIterator it(b.GetExtent()); !it.AtEnd(); ++it) {
    if(b.GetReal(it) != 0.0) return false;
  }
  return true;
}

void TestPointToCoord() 
{
  RealSpatialImageState spatial(Extent(Size(10,10,10)),
                                PixelSampling(Vec3(2.0,2.0,2.0), SPATIAL, 
                                              Extent(Size(10,10, 10))));
  ComplexFrequencyImageState frequency(Extent(Size(10,10,10)),
                                       PixelSampling(Vec3(2.0,2.0,2.0), FREQUENCY, 
                                                     Extent(Size(10,10,10))));
  spatial.SetAbsoluteOrigin(Vec3(5.0,5.0,5.0));
  frequency.SetAbsoluteOrigin(Vec3(5.0,5.0,5.0));
  BOOST_CHECK(spatial.IndexToCoord(Point(2.0,2.0,2.0)) == Vec3(14.0,14.0,14.0));
  BOOST_CHECK(frequency.IndexToCoord(Point(2.0,2.0,2.0)) == Vec3(0.1,0.1,0.1));
  BOOST_CHECK(frequency.CoordToIndex(Vec3(0.2,0.2,0.2)) == Vec3(4.0,4.0,4.0));
  BOOST_CHECK(spatial.CoordToIndex(Vec3(16.0,16.0,16.0)) == Vec3(3.0,3.0,3.0));
}

class TestRealSpatial { 
public:

  void test_Construction() {
    Extent ext(Size(8,10,13));
    RealSpatialImageState st(ext,PixelSampling());

    BOOST_CHECK(st.GetSize()==Size(8,10,13));

    BOOST_CHECK(is_clear(st));
  }

  static test_suite* Create()
  {
    test_suite* ts=BOOST_TEST_SUITE("image state test: real/spatial");
    
    boost::shared_ptr<TestRealSpatial> p(new TestRealSpatial());
    ts->add( BOOST_CLASS_TEST_CASE(&TestRealSpatial::test_Construction,p));
    return ts;
  }
};

class TestComplexHalfFrequency { 
public:

  void test_Construction1D() {
    std::ostringstream ostr;
    ComplexHalfFrequencyImageState st1( Extent(Size(15)),PixelSampling() );
    ostr.str("");
    ostr << "expected " << Extent(Point(0),Point(7)) << " but got " << st1.GetExtent();
    BOOST_CHECK_MESSAGE(st1.GetExtent()==Extent(Point(0),Point(7)),ostr.str());
    ostr.str("");
    ostr << "expected " << Extent(Point(-7),Point(7)) << " but got " << st1.GetLogicalExtent();
    BOOST_CHECK_MESSAGE(st1.GetLogicalExtent()==Extent(Point(-7),Point(7)),ostr.str());
    BOOST_CHECK(is_clear(st1));
  }

  void test_Construction2D() {
    std::ostringstream ostr;
    ComplexHalfFrequencyImageState st2( Extent(Size(8,10)),PixelSampling() );
    ostr.str("");
    ostr << "expected " << Extent(Point(-3,0),Point(4,5)) << " but got " << st2.GetExtent();
    BOOST_CHECK_MESSAGE(st2.GetExtent()==Extent(Point(-3,0),Point(4,5)),ostr.str());
    ostr.str("");
    ostr << "expected " << Extent(Point(-3,-4),Point(4,5)) << " but got " << st2.GetLogicalExtent();
    BOOST_CHECK_MESSAGE(st2.GetLogicalExtent()==Extent(Point(-3,-4),Point(4,5)),ostr.str());
    BOOST_CHECK(is_clear(st2));
  }

  void test_Construction3D() {
    std::ostringstream ostr;
    ComplexHalfFrequencyImageState st3( Extent(Size(8,10,13)),PixelSampling() );
    ostr.str("");
    ostr << "expected " << Extent(Point(-3,-4,0),Point(4,5,6)) << " but got " << st3.GetExtent();
    BOOST_CHECK_MESSAGE(st3.GetExtent()==Extent(Point(-3,-4,0),Point(4,5,6)),ostr.str());
    ostr.str("");
    ostr << "expected " << Extent(Point(-3,-4,-6),Point(4,5,6)) << " but got " << st3.GetLogicalExtent();
    BOOST_CHECK_MESSAGE(st3.GetLogicalExtent()==Extent(Point(-3,-4,-6),Point(4,5,6)),ostr.str());
    BOOST_CHECK(is_clear(st3));
  }

  void test_ConjugateAccess1D() {
    ComplexHalfFrequencyImageState st1( Extent(Size(15)),PixelSampling() );
    st1.ApplyIP(randomizer);
    st1.SetComplex(Point(0),Complex(0,0));
    for(int u=1;u<=7;++u) {
      BOOST_REQUIRE(st1.GetComplex(Point(u))==std::conj(st1.GetComplex(Point(-u))));
    }
  }

  void test_ConjugateAccess2D() {
    std::ostringstream ostr;
    ComplexHalfFrequencyImageState st2( Extent(Size(4,7)),PixelSampling() );
    st2.ApplyIP(randomizer);
    st2.SetComplex(Point(0,0),Complex(0,0));
    for(int v=1;v<=3;++v) {
      for(int u=-1;u<=1;++u) {
	ostr.str("");
	ostr << "@" << Point(u,v) << " expected conjugate value of " << st2.GetComplex(Point(u,v));
	BOOST_REQUIRE_MESSAGE(st2.GetComplex(Point(u,v))==std::conj(st2.GetComplex(Point(-u,-v))),ostr.str());
      }
    }
  }

  void test_ConjugateAccess3D() {
    std::ostringstream ostr;
    ComplexHalfFrequencyImageState st3( Extent(Size(3,4,5)),PixelSampling() );
    st3.ApplyIP(randomizer);
    st3.SetComplex(Point(0,0),Complex(0,0));
    for(int w=1;w<=2;++w) {
      for(int v=-1;v<=1;++v) {
        for(int u=-1;u<=1;++u) {
          ostr.str("");
          ostr << "@" << Point(u,v,w) << " expected conjugate value of " 
               << st3.GetComplex(Point(u,v,w));
          BOOST_REQUIRE_MESSAGE(st3.GetComplex(Point(u,v,w))==std::conj(st3.GetComplex(Point(-u,-v,-w))),ostr.str());
        }
      }
    }
  }

  static test_suite* Create()
  {
    test_suite* ts=BOOST_TEST_SUITE("image state test: half/complex/frequency");
    
    boost::shared_ptr<TestComplexHalfFrequency> p(new TestComplexHalfFrequency());
    ts->add(BOOST_CLASS_TEST_CASE(&TestComplexHalfFrequency::test_Construction1D,p));
    ts->add(BOOST_CLASS_TEST_CASE(&TestComplexHalfFrequency::test_Construction2D,p));
    ts->add(BOOST_CLASS_TEST_CASE(&TestComplexHalfFrequency::test_Construction3D,p));
    ts->add(BOOST_CLASS_TEST_CASE(&TestComplexHalfFrequency::test_ConjugateAccess1D,p));
    ts->add(BOOST_CLASS_TEST_CASE(&TestComplexHalfFrequency::test_ConjugateAccess2D,p));
    ts->add(BOOST_CLASS_TEST_CASE(&TestComplexHalfFrequency::test_ConjugateAccess3D,p));
    //ts->add(&TestPointToCoord);
    return ts;
  }
};

} // namespace

test_suite* CreateImageStateTest()
{
  using namespace test_image_state;
  test_suite* ts=BOOST_TEST_SUITE("image state test");

  ts->add( TestRealSpatial::Create() );
  ts->add( TestComplexHalfFrequency::Create() );
  
  return ts;
}
