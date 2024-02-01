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
  boost::python wrapper for alg/filter

  Author: Ansgar Philippsen
*/

#include <boost/python.hpp>
using namespace boost::python;

#include <ost/img/alg/filter.hh>

using namespace ost::img;
using namespace ost::img::alg;

void export_Filter()
{
  class_<GaussianFilter, bases<ModIPAlgorithm> >("GaussianFilter",
						 init<optional<Real> >() )
    .def("SetSigma",&GaussianFilter::SetSigma)
    .def("SetQ",&GaussianFilter::SetQ)
  ;

  class_<GaussianGradientMagnitudeFilter, bases<ModIPAlgorithm> >("GaussianGradientMagnitudeFilter",
						 init<optional<Real> >() )
    .def("SetSigma",&GaussianGradientMagnitudeFilter::SetSigma)
    .def("SetQ",&GaussianGradientMagnitudeFilter::SetQ)
  ;

  class_<GaussianLaplacianFilter, bases<ModIPAlgorithm> >("GaussianLaplacianFilter",
						 init<optional<Real> >() )
    .def("SetSigma",&GaussianLaplacianFilter::SetSigma)
    .def("SetQ",&GaussianLaplacianFilter::SetQ)
  ;





  class_<AnisotropicFilter, bases<GaussianFilter> >("AnisotropicFilter",
						    init<optional<Real> >() )
  ;

  class_<LowPassFilter, bases<ConstModIPAlgorithm> >("LowPassFilter",init<Real>())
    .def("GetLimit", &LowPassFilter::GetLimit)
    .def("SetLimit", &LowPassFilter::SetLimit)
  ;  

  class_<HighPassFilter, bases<ConstModIPAlgorithm> >("HighPassFilter",init<Real>())
    .def("GetLimit", &HighPassFilter::GetLimit)
    .def("SetLimit", &HighPassFilter::SetLimit)
  ;  

  class_<GaussianLowPassFilter, bases<ConstModIPAlgorithm> >("GaussianLowPassFilter",init<Real>())
    .def("GetLimit", &GaussianLowPassFilter::GetLimit)
    .def("SetLimit", &GaussianLowPassFilter::SetLimit)
  ;  
  
  class_<GaussianHighPassFilter, bases<ConstModIPAlgorithm> >("GaussianHighPassFilter",init<Real>())
    .def("GetLimit", &GaussianHighPassFilter::GetLimit)
    .def("SetLimit", &GaussianHighPassFilter::SetLimit)
  ;  

  class_<FermiLowPassFilter, bases<ConstModIPAlgorithm> >("FermiLowPassFilter",init<Real,Real>())
    .def("GetLimit", &FermiLowPassFilter::GetLimit)
    .def("SetLimit", &FermiLowPassFilter::SetLimit)
    .def("GetT", &FermiLowPassFilter::GetT)
    .def("SetT", &FermiLowPassFilter::SetT)
  ;  

  class_<FermiHighPassFilter, bases<ConstModIPAlgorithm> >("FermiHighPassFilter",init<Real,Real>())

    .def("GetLimit", &FermiHighPassFilter::GetLimit)
    .def("SetLimit", &FermiHighPassFilter::SetLimit)
    .def("GetT", &FermiHighPassFilter::GetT)
    .def("SetT", &FermiHighPassFilter::SetT)
  ;  

  class_<ButterworthLowPassFilter, bases<ConstModIPAlgorithm> >("ButterworthLowPassFilter",init<Real,Real>())
    .def("GetLimit", &ButterworthLowPassFilter::GetLimit)
    .def("SetLimit", &ButterworthLowPassFilter::SetLimit)
    .def("GetStop", &ButterworthLowPassFilter::GetStop)
    .def("SetStop", &ButterworthLowPassFilter::SetStop)
    .def("GetEps", &ButterworthLowPassFilter::GetEps)
    .def("SetEps", &ButterworthLowPassFilter::SetEps)
    .def("GetA", &ButterworthLowPassFilter::GetA)
    .def("SetA", &ButterworthLowPassFilter::SetA)
  ;
  class_<ButterworthHighPassFilter, bases<ConstModIPAlgorithm> >("ButterworthHighPassFilter",init<Real,Real>())
    .def("GetLimit", &ButterworthHighPassFilter::GetLimit)
    .def("SetLimit", &ButterworthHighPassFilter::SetLimit)
    .def("GetStop", &ButterworthHighPassFilter::GetStop)
    .def("SetStop", &ButterworthHighPassFilter::SetStop)
    .def("GetEps", &ButterworthHighPassFilter::GetEps)
    .def("SetEps", &ButterworthHighPassFilter::SetEps)
    .def("GetA", &ButterworthHighPassFilter::GetA)
    .def("SetA", &ButterworthHighPassFilter::SetA)
  ;

}
