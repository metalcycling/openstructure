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
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <ost/img/image.hh>
#include <ost/img/alg/filter.hh>
#include <ost/img/alg/randomize.hh>
#include <ost/img/normalizer_impl.hh>

using namespace ost::img;
using namespace ost::img::alg;

void LineIterator1DTest() 
{
  ImageHandle image = CreateImage(Size(3));
  LineIterator li1(image.GetExtent(),Axis::X);
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li1;
  BOOST_CHECK(li1.AtEnd());

  LineIterator li2(image.GetExtent(),Axis::Y);
  BOOST_CHECK(li2.AtEnd());

  LineIterator li3(image.GetExtent(),Axis::Z);
  BOOST_CHECK(li3.AtEnd());

  LineIterator li4(image.GetExtent(),(Axis::X | Axis::Y));
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li4;
  BOOST_CHECK(li4.AtEnd());

  LineIterator li5(image.GetExtent(),(Axis::X | Axis::Z));
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li5;
  BOOST_CHECK(li5.AtEnd());

  LineIterator li6(image.GetExtent(),(Axis::Y | Axis::Z));
  BOOST_CHECK(li6.AtEnd());

  LineIterator li7(image.GetExtent(),(Axis::X | Axis::Y | Axis::Z ));
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li7;
  BOOST_CHECK(li7.AtEnd());
}

void LineIterator2DTest() 
{
  ImageHandle image = CreateImage(Size(3,3)); 
  LineIterator li1(image.GetExtent(),Axis::X);
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,1,0),Point(3,2,1)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,2,0),Point(3,3,1)));
  ++li1;
  BOOST_CHECK(li1.AtEnd());

  LineIterator li2(image.GetExtent(),Axis::Y);
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,3,1)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,3,1)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,3,1)));
  ++li2;
  BOOST_CHECK(li2.AtEnd());

  LineIterator li3(image.GetExtent(),Axis::Z);
  BOOST_CHECK(li3.AtEnd());

  LineIterator li4(image.GetExtent(),(Axis::X | Axis::Y));
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,1,0),Point(3,2,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,2,0),Point(3,3,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,3,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,3,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,3,1)));
  ++li4;
  BOOST_CHECK(li4.AtEnd());

  LineIterator li5(image.GetExtent(),Axis::X | Axis::Z);
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,1,0),Point(3,2,1)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,2,0),Point(3,3,1)));
  ++li5;
  BOOST_CHECK(li5.AtEnd());

  LineIterator li6(image.GetExtent(),(Axis::Y | Axis::Z));
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,3,1)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,3,1)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,3,1)));
  ++li6;
  BOOST_CHECK(li6.AtEnd());

  LineIterator li7(image.GetExtent(),(Axis::X | Axis::Y | Axis::Z));
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,1,0),Point(3,2,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,2,0),Point(3,3,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,3,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,3,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,3,1)));
  ++li7;
  BOOST_CHECK(li7.AtEnd());

}

void LineIterator3DTest() 
{
  ImageHandle image = CreateImage(Size(3,3,3));  
  LineIterator li1(image.GetExtent(),Axis::X);
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,1,0),Point(3,2,1)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,2,0),Point(3,3,1)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,0,1),Point(3,1,2)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,1,1),Point(3,2,2)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,2,1),Point(3,3,2)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,0,2),Point(3,1,3)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,1,2),Point(3,2,3)));
  ++li1;
  BOOST_CHECK(li1.GetLine().GetExtent()==Extent(Point(0,2,2),Point(3,3,3)));
  ++li1;
  BOOST_CHECK(li1.AtEnd());

  LineIterator li2(image.GetExtent(),Axis::Y);
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,3,1)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,3,1)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,3,1)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(0,0,1),Point(1,3,2)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(1,0,1),Point(2,3,2)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(2,0,1),Point(3,3,2)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(0,0,2),Point(1,3,3)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(1,0,2),Point(2,3,3)));
  ++li2;
  BOOST_CHECK(li2.GetLine().GetExtent()==Extent(Point(2,0,2),Point(3,3,3)));
  ++li2;
  BOOST_CHECK(li2.AtEnd());

  LineIterator li3(image.GetExtent(),Axis::Z);
  BOOST_CHECK(li3.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,1,3)));
  ++li3;
  BOOST_CHECK(li3.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,1,3)));
  ++li3;
  BOOST_CHECK(li3.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,1,3)));
  ++li3;
  BOOST_CHECK(li3.GetLine().GetExtent()==Extent(Point(0,1,0),Point(1,2,3)));
  ++li3;
  BOOST_CHECK(li3.GetLine().GetExtent()==Extent(Point(1,1,0),Point(2,2,3)));
  ++li3;
  BOOST_CHECK(li3.GetLine().GetExtent()==Extent(Point(2,1,0),Point(3,2,3)));
  ++li3;
  BOOST_CHECK(li3.GetLine().GetExtent()==Extent(Point(0,2,0),Point(1,3,3)));
  ++li3;
  BOOST_CHECK(li3.GetLine().GetExtent()==Extent(Point(1,2,0),Point(2,3,3)));
  ++li3;
  BOOST_CHECK(li3.GetLine().GetExtent()==Extent(Point(2,2,0),Point(3,3,3)));
  ++li3;
  BOOST_CHECK(li3.AtEnd());

  LineIterator li4(image.GetExtent(),(Axis::X | Axis::Y));
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,1,0),Point(3,2,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,2,0),Point(3,3,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,0,1),Point(3,1,2)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,1,1),Point(3,2,2)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,2,1),Point(3,3,2)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,0,2),Point(3,1,3)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,1,2),Point(3,2,3)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,2,2),Point(3,3,3)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,3,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,3,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,3,1)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,0,1),Point(1,3,2)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(1,0,1),Point(2,3,2)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(2,0,1),Point(3,3,2)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(0,0,2),Point(1,3,3)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(1,0,2),Point(2,3,3)));
  ++li4;
  BOOST_CHECK(li4.GetLine().GetExtent()==Extent(Point(2,0,2),Point(3,3,3)));
  ++li4;
  BOOST_CHECK(li4.AtEnd());

  LineIterator li5(image.GetExtent(),(Axis::X | Axis::Z));
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,1,0),Point(3,2,1)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,2,0),Point(3,3,1)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,0,1),Point(3,1,2)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,1,1),Point(3,2,2)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,2,1),Point(3,3,2)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,0,2),Point(3,1,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,1,2),Point(3,2,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,2,2),Point(3,3,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,1,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,1,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,1,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,1,0),Point(1,2,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(1,1,0),Point(2,2,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(2,1,0),Point(3,2,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(0,2,0),Point(1,3,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(1,2,0),Point(2,3,3)));
  ++li5;
  BOOST_CHECK(li5.GetLine().GetExtent()==Extent(Point(2,2,0),Point(3,3,3)));
  ++li5;
  BOOST_CHECK(li5.AtEnd());

  LineIterator li6(image.GetExtent(),(Axis::Y | Axis::Z));
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,3,1)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,3,1)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,3,1)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(0,0,1),Point(1,3,2)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(1,0,1),Point(2,3,2)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(2,0,1),Point(3,3,2)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(0,0,2),Point(1,3,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(1,0,2),Point(2,3,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(2,0,2),Point(3,3,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,1,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,1,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,1,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(0,1,0),Point(1,2,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(1,1,0),Point(2,2,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(2,1,0),Point(3,2,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(0,2,0),Point(1,3,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(1,2,0),Point(2,3,3)));
  ++li6;
  BOOST_CHECK(li6.GetLine().GetExtent()==Extent(Point(2,2,0),Point(3,3,3)));
  ++li6;
  BOOST_CHECK(li6.AtEnd());

  LineIterator li7(image.GetExtent(),(Axis::X | Axis::Y | Axis::Z));
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,0),Point(3,1,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,1,0),Point(3,2,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,2,0),Point(3,3,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,1),Point(3,1,2)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,1,1),Point(3,2,2)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,2,1),Point(3,3,2)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,2),Point(3,1,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,1,2),Point(3,2,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,2,2),Point(3,3,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,3,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,3,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,3,1)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,1),Point(1,3,2)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(1,0,1),Point(2,3,2)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(2,0,1),Point(3,3,2)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,2),Point(1,3,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(1,0,2),Point(2,3,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(2,0,2),Point(3,3,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,0,0),Point(1,1,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(1,0,0),Point(2,1,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(2,0,0),Point(3,1,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,1,0),Point(1,2,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(1,1,0),Point(2,2,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(2,1,0),Point(3,2,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(0,2,0),Point(1,3,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(1,2,0),Point(2,3,3)));
  ++li7;
  BOOST_CHECK(li7.GetLine().GetExtent()==Extent(Point(2,2,0),Point(3,3,3)));
  ++li7;
  BOOST_CHECK(li7.AtEnd());

}

void FilterTest() 
{

  ImageHandle image = CreateImage(Size(7,7,7));
  Extent filled_ext(Size(10,10,10),Point(32,32,32));
  for (ExtentIterator it(filled_ext);!it.AtEnd();++it)
  {
    image.SetReal(Point(it),10.0);
  }

  GaussianFilter gf(1);

  ImageHandle image_gauss=image.Apply(gf);

  GaussianGradientMagnitudeFilter gfdf(1);

  ImageHandle image_gauss_first=image.Apply(gfdf);

  GaussianLaplacianFilter gsdf(1);

  ImageHandle image_gauss_second=image.Apply(gsdf);

  AnisotropicFilter af(1.0);

  ImageHandle image_anisotropic=image.Apply(af);
}

void FourierFiltersTest() 
{

  using namespace ost::img;
  using namespace ost::img::alg;
 
  ImageHandle image =  CreateImage(Size(40,40,40),REAL,SPATIAL);
  Vec3 vox_sampl(1.0,1.0,1.0); 
  image.SetPixelSampling(vox_sampl);

  Randomize rand;
  image.ApplyIP(rand);

  image = image+1.0;  

  LowPassFilter lp(4.0);
  ImageHandle image2 = image.Apply(lp);

  HighPassFilter hp(4.0);
  ImageHandle image3 = image.Apply(hp);

  GaussianLowPassFilter gauss_lp(4.0);
  ImageHandle image4 = image.Apply(gauss_lp);

  GaussianHighPassFilter gauss_hp(4.0);
  ImageHandle image5 = image.Apply(gauss_hp);

  FermiLowPassFilter fermi_lp(4.0,0.05);
  ImageHandle image6 = image.Apply(fermi_lp);

  FermiHighPassFilter fermi_hp(4.0,0.05);
  ImageHandle image7 = image.Apply(fermi_hp);

  ButterworthLowPassFilter butterworth_lp(4.0,3.5);
  ImageHandle image8 = image.Apply(butterworth_lp);

  ButterworthHighPassFilter butterworth_hp(4.0,4.5);
  ImageHandle image9 = image.Apply(butterworth_hp);

}


boost::unit_test::test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    boost::unit_test::test_suite* test = BOOST_TEST_SUITE( "FiltersTest" );

    test->add( BOOST_TEST_CASE( &LineIterator1DTest ) );
    test->add( BOOST_TEST_CASE( &LineIterator2DTest ) );
    test->add( BOOST_TEST_CASE( &LineIterator3DTest ) );
    test->add( BOOST_TEST_CASE( &FourierFiltersTest ) );
    test->add( BOOST_TEST_CASE( &FilterTest ) );

    return test;
}
