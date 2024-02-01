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

#include <boost/python.hpp>
using namespace boost::python;

#include <ost/img/extent.hh>
#include <ost/img/alg/transform.hh>
#include <ost/img/alg/discrete_shrink.hh>
#include <ost/img/alg/fractional_shift.hh>

using namespace ost::img;
namespace tf = img::alg;

#if 0
Extent (*transform_extent)(const Extent&,const Transformer&) = tf::Transform;
Point (*transform_point)(const Point&,const Transformer&) = tf::Transform;

struct TransformerWrap : Transformer, wrapper<Transformer>
{
  Mat4 GetMatrix() const
  {
    return this->get_override("GetMatrix")();
  }
  Mat4 GetInverseMatrix() const
  {
    return this->get_override("GetInverseMatrix")();
  }
};
#endif

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(setshift_overloads, SetShift, 0,3)

void    (tf::FractionalShift::*setshift1)(Real,Real,Real)    = &tf::FractionalShift::SetShift;
void    (tf::FractionalShift::*setshift2)(const Vec3&)                    = &tf::FractionalShift::SetShift;


BOOST_PYTHON_MODULE(_ost_tf)
{
  class_<tf::Transformation>("Transformation", init<>())
    .def(init<const Mat4&>())
    .def(init<const Mat4&, const Mat4&>())
    ;

  class_<tf::Rotate2D, bases<tf::Transformation> >("Rotate2D", init<Real>() );
  class_<tf::Rotate2D90, bases<tf::Transformation> >("Rotate2D90", init<int>() );
  class_<tf::Scale, bases<tf::Transformation> >("Scale", init<Real>() );
  class_<tf::Translate, bases<tf::Transformation> >("Translate", init<const Vec3&>() );
 
  class_<tf::Transform,bases<ConstModOPAlgorithm> >("Transform",init<>())
    .def(init<const tf::Transformation&, optional<const Vec3&> >())
    ;

  class_<tf::Shift,bases<ConstModOPAlgorithm> >("Shift",init<>())
    .def(init<const Point&>())
    ;

  class_<tf::Mirror,bases<ConstModOPAlgorithm> >("Mirror",init<>())
    .def(init<int>())
    ;

  class_<tf::DiscreteShrink, bases<ConstModOPAlgorithm> >("DiscreteShrink",init<optional<const Size&> >())
    .def("SetBlocksize",&tf::DiscreteShrink::SetBlocksize)
    .def("GetBlocksize",&tf::DiscreteShrink::GetBlocksize)
    ;
    
  class_<tf::FractionalShift, bases<ModIPAlgorithm> >("FractionalShift", init<optional <Real,Real,Real> >() )
    .def(init<const Vec3&>())
    .def("SetShift",setshift1,setshift_overloads())
    .def("SetShift",setshift2)
    .def("GetShift",&tf::FractionalShift::GetShift)
    ;


#if 0
  // transformer wrapper
  class_<TransformerWrap, boost::noncopyable>("Transformer", no_init)
    .def("GetMatrix",pure_virtual(&Transformer::GetMatrix))
    .def("GetInverseMatrix",pure_virtual(&Transformer::GetInverseMatrix))
    ;
  
  class_<tf::Rotate2D, bases<Transformer> >("Rotate2D", init<Real>() );
  class_<tf::Rotate2D90, bases<Transformer> >("Rotate2D90", init<int>() );
  class_<tf::Scale, bases<Transformer> >("Scale", init<Real>() );
  class_<tf::Translate, bases<Transformer> >("Translate", init<const Vec3&>() );

  class_<tf::CustomTransformer, bases<Transformer> >("CustomTransformer", init<const Mat4&, optional<const Mat4&> >());

  def("Transform",transform_extent, args("Extent","Transformer"));
  def("Transform",transform_point, args("Point","Transformer"));
#endif

}
