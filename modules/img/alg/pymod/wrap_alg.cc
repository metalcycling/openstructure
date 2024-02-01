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
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include <boost/python.hpp>
#include <ost/img/alg/auto_correlate.hh>
#include <ost/img/alg/clear.hh>
#include <ost/img/alg/conjugate.hh>
#include <ost/img/alg/correlate.hh>
#include <ost/img/alg/cross_correlate.hh>
#include <ost/img/alg/convolute.hh>
#include <ost/img/alg/fft.hh>
#include <ost/img/alg/dft.hh>
#include <ost/img/alg/fill.hh>
#include <ost/img/alg/threshold.hh>
#include <ost/img/alg/density_slice.hh>
#include <ost/img/alg/histogram.hh>
#include <ost/img/alg/negate.hh>
#include <ost/img/alg/stat.hh>
#include <ost/img/alg/stat_min_max.hh>
#include <ost/img/alg/randomize.hh>
#include <ost/img/alg/power_spectrum.hh>
#include <ost/img/alg/mask_image.hh>
#include <ost/img/alg/smooth_mask_image.hh>
#include <ost/img/alg/clip_min_max.hh>
#include <ost/img/alg/local_sigma_threshold.hh>
#include <ost/img/extent.hh>
#include <ost/img/alg/transform.hh>
#include <ost/img/alg/discrete_shrink.hh>
#include <ost/img/alg/fractional_shift.hh>
#include <ost/img/alg/line_average.hh>
#include <ost/img/alg/rscrosscorr.hh>
#include <ost/img/alg/highest_peak_search_3d.hh>

using namespace ost;
using namespace ost::img;
//using namespace ost::img::alg;
using namespace boost::python;

// fw decl
void export_Filter();
void export_Normalizer();
void export_Transcendentals();
void export_Polar();

namespace {
list get_histo_bins(const alg::Histogram& hist) {
  list result;
  const std::vector<int>& bins = hist.GetBins();
  for(std::vector<int>::const_iterator it=bins.begin();it!=bins.end();++it) {
    result.append(*it);
  }
  return result;
}

void frac_shift0(alg::FractionalShift* s) 
{
  s->SetShift();
}
  
void frac_shift1(alg::FractionalShift* s, Real sx) 
{
  s->SetShift(sx);
}

void frac_shift2(alg::FractionalShift* s, Real sx, Real sy) 
{
  s->SetShift(sx,sy);
}
  
void frac_shift3(alg::FractionalShift* s, Real sx, Real sy, Real sz) 
{
  s->SetShift(sx,sy,sz);
}
  
void frac_shift4(alg::FractionalShift* s, const Vec2& v) 
{
  s->SetShift(v);
}
  
void frac_shift5(alg::FractionalShift* s, const Vec3& v) 
{
  s->SetShift(v);
}
  
Real real_spatial_cross_correlation_1(const ConstImageHandle& image1, 
                                      const ConstImageHandle& image2)
{
  return alg::RealSpatialCrossCorrelation(image1,image2);
}


Real real_spatial_cross_correlation_2(const ConstImageHandle& image1, 
                                      const ConstImageHandle& image2, 
                                      const Extent& extent)
{
  return alg::RealSpatialCrossCorrelation(image1,image2,extent);
}

} //anon ns

BOOST_PYTHON_MODULE(_ost_img_alg)
{
  class_<alg::AutoCorrelate, bases<ConstModOPAlgorithm> >("AutoCorrelate",init<>());

  class_<alg::Clear,bases<ConstModIPAlgorithm> >("Clear", init<>());

  class_<alg::Conj,bases<ConstModIPAlgorithm> >("Conj", init<>());

  class_<alg::ExplicitConvolute,bases<ConstModOPAlgorithm> >("ExplicitConvolute", init<const ConstImageHandle&,optional<bool> >());

  class_<alg::ExplicitCorrelate,bases<ConstModOPAlgorithm> >("ExplicitCorrelate", init<const ConstImageHandle&>());

  class_<alg::CrossCorrelate,bases<ConstModIPAlgorithm> >("CrossCorrelate", init<const ConstImageHandle&>());

  class_<alg::DFT,bases<ConstModOPAlgorithm> >("DFT", init<>());

  class_<alg::FFT,bases<ConstModOPAlgorithm> >("FFT", init<>());

  class_<alg::PowerSpectrum,bases<ConstModOPAlgorithm> >("PowerSpectrum",init<>());

  class_<alg::Fill,bases<ConstModIPAlgorithm> >("Fill", init<const Complex&>())
    .def(init<Real>())
    ;

  class_<alg::Threshold,bases<ConstModIPAlgorithm> >("Threshold", init<const Complex&>())
    .def(init<Real>())
    ;

  class_<alg::DensitySlice,bases<ConstModIPAlgorithm> >("DensitySlice",
						   init<const Complex&,const Complex&>())
    .def(init<Real,Real>())
    ;

  class_<alg::ClipMinMax,bases<ConstModIPAlgorithm> >("ClipMinMax", 
						 init<const Complex&,const Complex&>())
    .def(init<Real,Real>())
    ;

  class_<alg::Histogram, bases<NonModAlgorithm> >("Histogram", init<int,Real,Real>() )
    .def("GetBins",get_histo_bins)
    ;

  class_<alg::Negate,bases<ConstModIPAlgorithm> >("Negate", init<>());

  class_<alg::Stat, bases<NonModAlgorithm> >("Stat", init<>() )
    .def("GetMean",&alg::Stat::GetMean)
    .def("SetMean",&alg::Stat::SetMean)
    .def("GetMinimum",&alg::Stat::GetMinimum)
    .def("GetMinimumPosition",&alg::Stat::GetMinimumPosition)
    .def("SetMinimum",&alg::Stat::SetMinimum)
    .def("GetMaximum",&alg::Stat::GetMaximum)
    .def("GetMaximumPosition",&alg::Stat::GetMaximumPosition)
    .def("SetMaximum",&alg::Stat::SetMaximum)
    .def("GetSum",&alg::Stat::GetSum)
    .def("GetRootMeanSquare",&alg::Stat::GetRootMeanSquare)
    .def("GetSkewness",&alg::Stat::GetSkewness)
    .def("GetKurtosis",&alg::Stat::GetKurtosis)
    .def("GetVariance",&alg::Stat::GetVariance)
    .def("SetVariance",&alg::Stat::SetVariance)
    .def("GetStandardDeviation",&alg::Stat::GetStandardDeviation)
    .def("SetStandardDeviation",&alg::Stat::SetStandardDeviation)
    .def("GetCenterOfMass",&alg::Stat::GetCenterOfMass)
    .def(self_ns::str(self))
  ;
  class_<alg::StatMinMax, bases<NonModAlgorithm> >("StatMinMax", init<>() )
    .def("GetMinimum",&alg::StatMinMax::GetMinimum)
    .def("GetMinimumPosition",&alg::StatMinMax::GetMinimumPosition)
    .def("SetMinimum",&alg::StatMinMax::SetMinimum)
    .def("GetMaximum",&alg::StatMinMax::GetMaximum)
    .def("GetMaximumPosition",&alg::StatMinMax::GetMaximumPosition)
    .def("SetMaximum",&alg::StatMinMax::SetMaximum)
    .def(self_ns::str(self))
  ;

  class_<alg::Randomize,bases<ConstModIPAlgorithm> >("Randomize", init<>());
  
  class_<alg::MaskImage, bases<ConstModIPAlgorithm> >("MaskImage",init<const MaskPtr&,optional<bool> >())
    ;

  class_<alg::SmoothMaskImage, bases<ModIPAlgorithm> >("SmoothMaskImage",init<const MaskPtr&, Real, optional<bool> >())
    ;

  class_<alg::LocalSigmaThreshold, bases<ConstModOPAlgorithm> >("LocalSigmaThreshold", init<int,Real>() );

  export_Filter();
  export_Normalizer();
  export_Polar();

  export_Transcendentals();

  class_<alg::Transformation>("Transformation", init<>())
    .def(init<const Mat4&>())
    .def(init<const Mat4&, const Mat4&>())
    ;

  class_<alg::Rotate2D, bases<alg::Transformation> >("Rotate2D", init<Real>() );
  class_<alg::Rotate2D90, bases<alg::Transformation> >("Rotate2D90", init<int>() );
  class_<alg::Scale, bases<alg::Transformation> >("Scale", init<Real>() );
  class_<alg::Translate, bases<alg::Transformation> >("Translate", init<const Vec3&>() );

  class_<alg::Transform,bases<ConstModOPAlgorithm> >("Transform",init<>())
    .def(init<const alg::Transformation&, optional<const Vec3&> >())
    ;

  class_<alg::Shift,bases<ConstModOPAlgorithm> >("Shift",init<>())
    .def(init<const Point&>())
    ;

  class_<alg::Mirror,bases<ConstModOPAlgorithm> >("Mirror",init<>())
    .def(init<int>())
    ;

  class_<alg::DiscreteShrink, bases<ConstModOPAlgorithm> >("DiscreteShrink",init<optional<const Size&> >())
    .def("SetBlocksize",&alg::DiscreteShrink::SetBlocksize)
    .def("GetBlocksize",&alg::DiscreteShrink::GetBlocksize)
    ;

  class_<alg::FractionalShift, bases<ConstModIPAlgorithm> >("FractionalShift", init<optional <Real,Real,Real> >() )
    .def(init<const Vec3&>())
    .def("SetShift",frac_shift0)
    .def("SetShift",frac_shift1)
    .def("SetShift",frac_shift2)
    .def("SetShift",frac_shift3)
    .def("SetShift",frac_shift4)
    .def("SetShift",frac_shift5)
    .def("GetShift",&alg::FractionalShift::GetShift)
    ;
  class_<alg::LineAverage, bases<ConstModOPAlgorithm> >("LineAverage", init<>() )
    .def(init<unsigned int>())
    .def("GetGeometry",&alg::LineAverage::GetGeometry)
    .def("SetGeometry",&alg::LineAverage::SetGeometry)
  ;
  def("RealSpatialCrossCorrelation",real_spatial_cross_correlation_1);
  def("RealSpatialCrossCorrelation",real_spatial_cross_correlation_2);

  class_<alg::HighestPeakSearch3D, bases<NonModAlgorithm> >("HighestPeakSearch3D" ,init<int,int,Real>())
    .def(init<int,int>())
    .def(init<int>())
    .def(init<>())
    .def("GetPeakList",&alg::HighestPeakSearch3D::GetPeakList)
    .def("ClearPeakList",&alg::HighestPeakSearch3D::ClearPeakList)
  ;

}
