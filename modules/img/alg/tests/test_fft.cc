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
#include <ost/base.hh>
#include <ost/img/image_state.hh>
#include <ost/img/alg/fft.hh>
#include <ost/img/alg/dft.hh>
#include <ost/img/alg/randomize.hh>
#include <ost/img/alg/alg_shift.hh>

#include "tests.hh"
#include "test_utils.hh"

namespace {

template <typename V, class DP>
void DumpStates(const image_state::ImageStateImpl<V,DP>& one,
		const image_state::ImageStateImpl<V,DP>& two)
{
  if(one.GetSize()!=two.GetSize()) {
    return;
  }

  for(image_state::IndexIterator it(one.GetSize()); !it.AtEnd(); ++it) {
    std::cerr << it << " " << one.Value(it) << "  " << two.Value(it) << std::endl;
  }
}


void Test_C2C_1D_calc(int N)
{
  try {
    ComplexSpatialImageState in=ComplexSpatialImageState( Extent(Size(N)),PixelSampling() );
    ComplexFrequencyImageState ref=ComplexFrequencyImageState( Extent(Size(N)),PixelSampling() );

    int step;

    if (N==64) {
      step=4;
    } else {
      step=3;
    }

    int sum =0;

    for (ExtentIterator ex_it(in.GetExtent());!ex_it.AtEnd();++ex_it) {
      if (ex_it[0] % step == 0) {
        in.Value(ex_it)=Complex(4.0,0.0);
        sum=sum+4;
      }
    }

    ImageStateBasePtr out_base;

    out_base = in.Apply(alg::FFT());

    ComplexFrequencyImageState* out_state = dynamic_cast<ComplexFrequencyImageState*>(out_base.get());

    BOOST_REQUIRE(out_state!=0);

    //std::cout << "FFT" << std::endl;
    //for (ExtentIterator ex_it(out_state->GetExtent());!ex_it.AtEnd();++ex_it) {
    //  std::cout << ex_it[0] << " " << out_state->Value(ex_it) << std::endl;
    //}

    BOOST_CHECK(out_state->GetExtent()==Extent(Size(N),Point(0)));

    Point peak(N/step,0);
    Complex peak_value = out_state->Value(peak);
    Complex calc_peak_value = Complex(sum,0.0);
    //std::cout << std::abs(peak_value-calc_peak_value);
    BOOST_REQUIRE(std::abs(peak_value-calc_peak_value)<1e-06);

  } catch (alg::FFTException& e) {
    BOOST_ERROR("FFT Exception caught!");
  }
}

void Test_C2C_1D_even_calc()
{
  Test_C2C_1D_calc(64);
}

void Test_C2C_1D_odd_calc()
{
  Test_C2C_1D_calc(63);
}

void Test_C2C_3D_calc()
{
  int N=16;
  try {
    ComplexSpatialImageState in( Extent(Size(N,N,N)),PixelSampling() );
    ComplexFrequencyImageState ref( Extent(Size(N,N,N)),PixelSampling() );

    int step = 4;
    int sum = 0;
    for (ExtentIterator ex_it(in.GetExtent());!ex_it.AtEnd();++ex_it) {
      if (ex_it[0] % step == 0 && ex_it[1] % step == 0 && ex_it[2] % step == 0) {
         in.Value(ex_it)=Complex(4.0,0.0);
         sum=sum+4;
      }
    }

    ImageStateBasePtr out_base;

    out_base = in.Apply(alg::FFT());

    ComplexFrequencyImageState* out_state = dynamic_cast<ComplexFrequencyImageState*>(out_base.get());

    BOOST_REQUIRE(out_state!=0);

    BOOST_CHECK(out_state->GetExtent()==Extent(Size(N,N,N),Point(0,0,0)));

    Point peak(N/step,0,0);
    Complex peak_value = out_state->Value(peak);
    Complex calc_peak_value = Complex(sum,0.0);
    BOOST_REQUIRE(peak_value == calc_peak_value);
    Point peak2(0,N/step,0);
    Complex peak_value2 = out_state->Value(peak2);
    Complex calc_peak_value2 = Complex(sum,0.0);
    BOOST_REQUIRE(peak_value2 == calc_peak_value2);
    Point peak3(0,0,N/step);
    Complex peak_value3 = out_state->Value(peak3);
    Complex calc_peak_value3 = Complex(sum,0.0);
    BOOST_REQUIRE(peak_value3 == calc_peak_value3);

  } catch (alg::FFTException& e) {
    BOOST_ERROR("FFT Exception caught!");
  }

}

void Test_R2H_1D_calc(int N)
{
  try {
    RealSpatialImageState in=RealSpatialImageState( Extent(Size(N)),PixelSampling() );
    ComplexHalfFrequencyImageState ref=ComplexHalfFrequencyImageState( Extent(Size(N)),PixelSampling() );

    int step;

    if (N==64) {
      step=4;
    } else {
      step=3;
    }

    int sum =0;

    for (ExtentIterator ex_it(in.GetExtent());!ex_it.AtEnd();++ex_it) {
      if (ex_it[0] % step == 0) {
        in.Value(ex_it)=4.0;
        sum=sum+4;
      }
    }

    ImageStateBasePtr out_base;

    out_base = in.Apply(alg::FFT());

    ComplexHalfFrequencyImageState* out_state = dynamic_cast<ComplexHalfFrequencyImageState*>(out_base.get());

    BOOST_REQUIRE(out_state!=0);

    BOOST_CHECK(out_state->GetExtent()==Extent(Point(0),Size(N/2+1)));
    BOOST_CHECK(out_state->GetLogicalExtent()==Extent(Size(N),Point(0)));

    //std::cout << "FFT" << std::endl;
    //for (ExtentIterator ex_it(out_state->GetExtent());!ex_it.AtEnd();++ex_it) {
    //  std::cout << ex_it[0] << " " << out_state->Value(ex_it) << std::endl;
    //}

    BOOST_CHECK(out_state->GetExtent()==Extent(Point(0),Size(N/2+1)));
    BOOST_CHECK(out_state->GetLogicalExtent()==Extent(Size(N),Point(0)));

    Point peak(N/step);
    Complex peak_value = out_state->Value(peak);
    Complex calc_peak_value = Complex(sum,0.0);

    BOOST_REQUIRE(peak_value == calc_peak_value);

  } catch (alg::FFTException& e) {
    BOOST_ERROR("FFT Exception caught in Test_R2H_1D_explicit!");
  }
}

void Test_R2H_1D_even_calc()
{
  Test_R2H_1D_calc(64);
}

void Test_R2H_1D_odd_calc()
{
  Test_R2H_1D_calc(63);
}

void Test_R2H_2D_explicit()
{
  try {
    static Real din[] = {
      0.732276, 0.514038, 0.614191, 0.812722, 0.029531, 0.544793,
      0.541531, 0.365433, 0.446773, 0.126797, 0.759645, 0.199166,
      0.609194, 0.310389, 0.491145, 0.360568, 0.916562, 0.54668,
      0.236379, 0.131892, 0.0490496, 0.0249416, 0.254728, 0.316022,
      0.908825, 0.0584649, 0.300042, 0.252541, 0.476233, 0.824687,
      0.902688, 0.468777, 0.042249, 0.590042, 0.0424999, 0.214492
    };
    
    RealSpatialImageState in_state( Extent(Size(6,6)),PixelSampling() );
    for(int v=0;v<6;++v) {
      for(int u=0;u<6;++u) {
	in_state.Value(Index(u,v,0)) = din[u*6+v];
      }
    }

    static Real tout[][2] = {
      {15.016,0}, {1.79937,1.15406}, {1.63976,0.226115}, 
      {1.6911,0}, {1.63976,-0.226115}, {1.79937,-1.15406},
      {1.55692,-0.512984}, {-0.0550974,-0.471976}, {0.416683,1.63666},
      {-0.851851,-1.37895}, {1.92232,1.51598}, {-0.234998,2.64483},
      {-1.11715,0.203644}, {0.227241,-0.979282}, {0.46483,1.3786}, 
      {-1.48793,-0.946194}, {-0.367004,0.0618105}, {-0.794775,0.601116},
      {3.58978,0}, {-0.179098,0.664966}, {0.0857343,1.5342}, 
      {0.0151358,0}, {0.0857343,-1.5342}, {-0.179098,-0.664966},
      {-1.11715,-0.203644}, {-0.794775,-0.601116}, {-0.367004,-0.0618105}, 
      {-1.48793,0.946194}, {0.46483,-1.3786}, {0.227241,0.979282},
      {1.55692,0.512984}, {-0.234998,-2.64483}, {1.92232,-1.51598}, 
      {-0.851851,1.37895}, {0.416683,-1.63666}, {-0.0550974,0.471976}
    };
    
    ComplexHalfFrequencyImageState ref_state(Extent(Size(6,6)),PixelSampling());
  
    for(int v=0;v<4;++v) {
      for(int u=0;u<6;++u) {
	int i=u*6+v;
	ref_state.Value(Index(u,v,0)) = Complex(tout[i][0],tout[i][1]);
      }
    }

    ImageStateBasePtr out_base;

    out_base = in_state.Apply(alg::FFT());

    ComplexHalfFrequencyImageState* out_state = dynamic_cast<ComplexHalfFrequencyImageState*>(out_base.get());

    BOOST_REQUIRE(out_state!=0);

    BOOST_REQUIRE(CheckEqualStates(ref_state,*out_state,1e-4));

    // bw tf
    ImageStateBasePtr out_rev_base;

    out_rev_base = out_state->Apply(alg::FFT());

    RealSpatialImageState* out_rev = dynamic_cast<RealSpatialImageState*>(out_rev_base.get());

    BOOST_REQUIRE(out_rev!=0);

    BOOST_REQUIRE(CheckEqualStates(in_state, *out_rev, 1e-4));

    BOOST_CHECK(in_state.GetExtent()==out_rev->GetExtent());

  } catch (alg::FFTException& e) {
    BOOST_ERROR("FFT Exception caught in Test_R2H_2D_explicit!");
  }
}

void Test_R2H_2D_calc()
{
  int N=16;
  try {
    RealSpatialImageState in( (Extent(Point(-2,-3),Size(N,N))),PixelSampling() );
    ComplexHalfFrequencyImageState ref( Extent(Size(N,N)),PixelSampling() );

    int step = 4;
    int sum = 0;
    for (ExtentIterator ex_it(in.GetExtent());!ex_it.AtEnd();++ex_it) {
      if (ex_it[0] % step == 0 && ex_it[1] % step == 0) {
         in.Value(ex_it)=4.0;
         sum=sum+4;
      }
    }

    ImageStateBasePtr out_base;

    out_base = in.Apply(alg::FFT());

    ComplexHalfFrequencyImageState* out_state = dynamic_cast<ComplexHalfFrequencyImageState*>(out_base.get());

    BOOST_REQUIRE(out_state!=0);

    BOOST_CHECK(out_state->GetLogicalExtent()==Extent(Size(N,N),Point(0,0)));
    Point st=Extent(Size(N,N),Point(0,0)).GetStart();
    Point pst(st[0],0);
    BOOST_CHECK(out_state->GetExtent()==Extent(pst,Size(N,N/2+1)));

    Point peak(N/step,0);
    Complex peak_value = out_state->Value(peak);
    Complex calc_peak_value = Complex(-sum,0.0);
    BOOST_REQUIRE(peak_value == calc_peak_value);
    Point peak2(0,N/step);
    Complex peak_value2 = out_state->Value(peak2);
    Complex calc_peak_value2 = Complex(0,sum);

    BOOST_REQUIRE(peak_value2 == calc_peak_value2);

  } catch (alg::FFTException& e) {
    BOOST_ERROR("FFT Exception caught in Test_R2H_2D_explicit!");
  }
}

void Test_R2H_3D_calc()
{
  int N=16;
  try {
    RealSpatialImageState in=RealSpatialImageState( Extent(Size(N,N,N),Point(0,0,0)),PixelSampling() );
    ComplexHalfFrequencyImageState ref=ComplexHalfFrequencyImageState( Extent(Size(N,N,N)),PixelSampling() );

    int step = 4;
    int sum = 0;
    for (ExtentIterator ex_it(in.GetExtent());!ex_it.AtEnd();++ex_it) {
      if (ex_it[0] % step == 0 && ex_it[1] % step == 0 && ex_it[2] % step == 0) {
         in.Value(ex_it)=4.0;
         sum=sum+4;
      }
    }

    ImageStateBasePtr out_base;

    out_base = in.Apply(alg::FFT());

    ComplexHalfFrequencyImageState* out_state = dynamic_cast<ComplexHalfFrequencyImageState*>(out_base.get());

    BOOST_REQUIRE(out_state!=0);

    BOOST_CHECK(out_state->GetLogicalExtent()==Extent(Size(N,N,N),Point(0,0,0)));
    Point st=Extent(Size(N,N,N),Point(0,0,0)).GetStart();
    Point pst(st[0],st[1],0);
    BOOST_CHECK(out_state->GetExtent()==Extent(pst,Size(N,N,N/2+1)));

    Point peak(N/step,0,0);
    Complex peak_value = out_state->Value(peak);
    Complex calc_peak_value = Complex(0.0,sum);
    BOOST_REQUIRE(peak_value == calc_peak_value);
    Point peak2(0,N/step,0);
    Complex peak_value2 = out_state->Value(peak2);
    Complex calc_peak_value2 = Complex(0.0,sum);
    BOOST_REQUIRE(peak_value2 == calc_peak_value2);
    Point peak3(0,0,N/step);
    Complex peak_value3 = out_state->Value(peak3);
    Complex calc_peak_value3 = Complex(0.0,sum);
    BOOST_REQUIRE(peak_value3 == calc_peak_value3);
    
  } catch (alg::FFTException& e) {
    BOOST_ERROR("FFT Exception caught in Test_R2H_3D_explicit!");
  }
}

void Test_DFT(DataType TYPE)
{
  Point st(2,3,-1);
  ImageHandle ri1=CreateImage(Extent(st,Size(4,3,2)), TYPE ,SPATIAL);
  ri1.ApplyIP(alg::Randomize());
  ImageHandle fi1=ri1.Apply(alg::DFT());

  ImageHandle ri2=ri1.Apply(alg::Shift(-st));
  ImageHandle fi2=ri2.Apply(alg::FFT());

  std::ostringstream msg;

  msg << "fi1 fi2 Extent mismatch: " << fi1.GetExtent() << "!=" << fi2.GetExtent();
  BOOST_REQUIRE_MESSAGE(fi1.GetExtent()==fi2.GetExtent(),msg.str());

  for(ExtentIterator it(fi1.GetExtent());!it.AtEnd();++it) {
    Real absdiff=std::abs(fi1.GetComplex(it)-fi2.GetComplex(it));
    msg.str("");
    msg << "@" << Point(it) << ": abs(" << fi1.GetComplex(it) << "-" <<fi2.GetComplex(it) << ")="<<absdiff;
    BOOST_REQUIRE_MESSAGE(absdiff<1e-4,msg.str());
  }

  ImageHandle ri3 = fi1.Apply(alg::DFT());

  msg.str("");
  msg << "ri1 ri3 Extent mismatch: " << ri3.GetExtent() << "!=" << ri1.GetExtent();
  BOOST_REQUIRE_MESSAGE(ri3.GetExtent()==ri1.GetExtent(),msg.str());

  for(ExtentIterator it(ri3.GetExtent());!it.AtEnd();++it) {
    Real absdiff=std::abs(ri3.GetReal(it)-ri1.GetReal(it));
    msg.str("");
    msg << "@" << Point(it) << ": abs(" << ri3.GetReal(it) << "-" <<ri1.GetReal(it) << ")="<<absdiff;
    BOOST_REQUIRE_MESSAGE(absdiff<1e-5,msg.str());
  }
}

void Test_DFT_REAL()
{
  Test_DFT(REAL);
}

void Test_DFT_COMPLEX()
{
  Test_DFT(COMPLEX);
}

void Test_Invalid()
{
  RealFrequencyImageState rfs;
  BOOST_CHECK_THROW(rfs.Apply(alg::FFT()), alg::FFTException);
}

void Test_Sampling()
{
  std::ostringstream msg;
  ImageHandle ih1=CreateImage(Size(4,10,20));
  ih1.SetPixelSampling(geom::Vec3(2.0,4.0,1.0));
  ih1.ApplyIP(alg::Randomize());
  BOOST_CHECK(geom::Length(ih1.GetPixelSampling()-geom::Vec3(2.0,4.0,1.0))<1e-5);
  BOOST_CHECK(geom::Length(ih1.GetSpatialSampling()-geom::Vec3(2.0,4.0,1.0))<1e-5);
  ih1.ApplyIP(alg::FFT());
  msg << "expected (0.125,0.025,0.05) but got " << ih1.GetPixelSampling();
  BOOST_CHECK_MESSAGE(geom::Length(ih1.GetPixelSampling()-geom::Vec3(0.125,0.025,0.05))<1e-5,msg.str());
  BOOST_CHECK(geom::Length(ih1.GetSpatialSampling()-geom::Vec3(2.0,4.0,1.0))<1e-5);
  ih1.ApplyIP(alg::FFT());
  BOOST_CHECK(geom::Length(ih1.GetPixelSampling()-geom::Vec3(2.0,4.0,1.0))<1e-5);
  BOOST_CHECK(geom::Length(ih1.GetSpatialSampling()-geom::Vec3(2.0,4.0,1.0))<1e-5);
}

struct test_memalloc_fnc {
  test_memalloc_fnc() {}

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& is1) {
    ImageStateImpl<T,D> is2 = is1;

    for(ExtentIterator it(is1.GetExtent());!it.AtEnd();++it) {
      T value = is1.Value(it);
      value+=1.0;
      is1.Value(it)=value;
    }

    for(ExtentIterator it(is2.GetExtent());!it.AtEnd();++it) {
      T value = is2.Value(it);
      value+=1.0;
      is2.Value(it)=value;
    }
  }

  String GetAlgorithmName() const {return "";}
};

typedef ImageStateModIPAlgorithm<test_memalloc_fnc> test_memalloc_alg;

void Test_Memalloc()
{
  test_memalloc_alg talg;
  ImageHandle i1=CreateImage(Size(6,5,4));
  i1.ApplyIP(alg::Randomize());
  i1.ApplyIP(talg);
  ImageHandle i2=i1.Apply(alg::FFT());
  i2.ApplyIP(talg);
  ImageHandle i3=i2.Apply(alg::FFT());
  i3.ApplyIP(talg);

  i2=CreateImage(Size(6,5,4),COMPLEX,HALF_FREQUENCY);
  i2.ApplyIP(talg);
  i3=i2.Apply(alg::FFT());
  i3.ApplyIP(talg);

  ImageHandle i4 =  CreateImage(Size(40,40,40),REAL,SPATIAL);
  Vec3 vox_sampl(1.0,1.0,1.0); 
  i4.SetPixelSampling(vox_sampl);

  alg::Randomize rand;
  alg::FFT fft;  

  i4.ApplyIP(rand);
  i4.ApplyIP(fft);
  i4.ApplyIP(fft);
  i4.ApplyIP(talg);
}


} // namespace 

test_suite* CreateFFTTest()
{
  test_suite* ts=BOOST_TEST_SUITE("FFT Test");

  ts->add(BOOST_TEST_CASE(&Test_C2C_1D_even_calc));
  ts->add(BOOST_TEST_CASE(&Test_C2C_1D_odd_calc));
  ts->add(BOOST_TEST_CASE(&Test_C2C_3D_calc));
  ts->add(BOOST_TEST_CASE(&Test_R2H_1D_even_calc));
  ts->add(BOOST_TEST_CASE(&Test_R2H_1D_odd_calc));
  ts->add(BOOST_TEST_CASE(&Test_R2H_2D_explicit));
  ts->add(BOOST_TEST_CASE(&Test_R2H_2D_calc));
  ts->add(BOOST_TEST_CASE(&Test_R2H_3D_calc));
  ts->add(BOOST_TEST_CASE(&Test_DFT_REAL));
  ts->add(BOOST_TEST_CASE(&Test_DFT_COMPLEX));
  ts->add(BOOST_TEST_CASE(&Test_Sampling));
  ts->add(BOOST_TEST_CASE(&Test_Invalid));
  ts->add(BOOST_TEST_CASE(&Test_Memalloc));

  return ts;
}

