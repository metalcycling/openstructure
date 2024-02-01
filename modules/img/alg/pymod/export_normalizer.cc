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
  boost::python wrapper for normalizer

  Author: Andreas Schenk
*/

#include <boost/python.hpp>
using namespace boost::python;

#include <ost/img/alg/norm.hh>

using namespace ost;
using namespace ost::img;
using namespace ost::img::alg;
namespace {

Normalizer (*linrange1)(const ConstImageHandle&,Real,Real) = CreateLinearRangeNormalizer;
Normalizer (*linrange2)(Real,Real,Real,Real) = CreateLinearRangeNormalizer;
Normalizer (*logrange1)(const ConstImageHandle&,Real,Real) = CreateLogRangeNormalizer;
Normalizer (*logrange2)(Real,Real,Real,Real) = CreateLogRangeNormalizer;
BOOST_PYTHON_FUNCTION_OVERLOADS(od_overloads, CreateTransmissionToODNormalizer, 2, 3)

} //anon ns

void export_Normalizer()
{
  Complex (Normalizer::*convert_complex)(Complex) const = &Normalizer::Convert;
  Real (Normalizer::*convert_double)(Real) const = &Normalizer::Convert;
  Complex (Normalizer::*backconvert_complex)(Complex) const = &Normalizer::BackConvert;
  Real (Normalizer::*backconvert_double)(Real) const = &Normalizer::BackConvert;
  
  class_<Normalizer,bases<ConstModIPAlgorithm> >("Normalizer",no_init)
   //complex defs first (otherwise Real will be converted to complex)
    .def("Convert",convert_complex)
    .def("Convert",convert_double)
    .def("BackConvert",backconvert_complex)
    .def("BackConvert",backconvert_double)
    ;
  
  def("CreateNoOpNormalizer",CreateNoOpNormalizer);
  def("CreateLinearNormalizer",CreateLinearNormalizer);
  def("CreateLinearRangeNormalizer",linrange1);
  def("CreateLinearRangeNormalizer",linrange2);
  def("CreateLogNormalizer",CreateLogNormalizer);
  def("CreateLogRangeNormalizer",logrange1);
  def("CreateLogRangeNormalizer",logrange2);
  def("CreateTransmissionToODNormalizer",CreateTransmissionToODNormalizer,od_overloads());
}
