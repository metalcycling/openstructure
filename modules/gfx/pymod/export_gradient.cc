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

/*
  Authors: Marco Biasini, Ansgar Philippsen
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include <ost/message.hh>
#include <ost/gfx/gradient.hh>
using namespace ost;
using namespace ost::gfx;

namespace {
  Gradient* make_gradient_d(const dict& d)
  {
    std::unique_ptr<Gradient> grad(new Gradient);
    list keys = d.keys();
    for(int i=0;i<len(keys);++i) {
      extract<float> fex(keys[i]);
      if(!fex.check()) {
        throw Error("expected floats as keys");
      }
      float mark = fex();
      Color col;
      object val = d[keys[i]];
      extract<Color> cex(val);
      if(cex.check()) {
        // use gfx.Color
        col=cex();
      } else {
        // try simple sequence
        if(len(val)!=3) {
          throw Error("expected values of gfx.Color or float triplets");
        }
        try {
          col=gfx::Color(extract<float>(val[0]),extract<float>(val[1]),extract<float>(val[2]));
        } catch (...) {
          throw Error("expected values of gfx.Color or float triplets");
        }
      }
      grad->SetColorAt(mark,col);
    }
    return grad.release();
  }

  Gradient* make_gradient_l(const list& l)
  {
    std::unique_ptr<Gradient> grad(new Gradient);
    float mf = len(l)<2 ? 0.0 : 1.0/static_cast<float>(len(l)-1); 
    for(int i=0;i<len(l);++i) {
      float mark = static_cast<float>(i)*mf;
      Color col;
      object val = l[i];
      extract<Color> cex(val);
      if(cex.check()) {
        // use gfx.Color
        col=cex();
      } else {
        // try simple sequence
        if(len(val)!=3) {
          throw Error("expected values of gfx.Color or float triplets");
        }
        try {
          col=gfx::RGB(extract<float>(val[0]),extract<float>(val[1]),extract<float>(val[2]));
        } catch (...) {
          throw Error("expected values of gfx.Color float triplets");
        }
      }
      grad->SetColorAt(mark,col);
    }
    return grad.release();
  }

  std::string sl_repr(const Gradient::StopList& sl) {
    std::ostringstream m;
    m << "[";
    for(size_t i=0;i<sl.size();++i) {
      Color c = sl[i].color;
      m << "(" << sl[i].t << "," << "gfx.RGB(" << c[0] << "," << c[1] << "," << c[2] << "))";
      if(i<sl.size()-1) m << ",";
    }
    m << "]";
    return m.str();
  }
}

void export_gradient()
{
  class_<Gradient>("Gradient", init<>())
    .def(init<const String&>())
    .def("__init__", make_constructor(make_gradient_d))
    .def("__init__", make_constructor(make_gradient_l))
    .def("SetColorAt", &Gradient::SetColorAt)
    .def("GetColorAt", &Gradient::GetColorAt)
    .def("GetStops", &Gradient::GetStops)
    .add_property("stops", &Gradient::GetStops)
    .def("GradientToInfo", &Gradient::GradientToInfo)
    .def("GradientFromInfo", &Gradient::GradientFromInfo).staticmethod("GradientFromInfo")
    .add_property("hsv_mode",&Gradient::GetHSVMode,&Gradient::SetHSVMode)
  ;
  implicitly_convertible<String, Gradient>();

  class_<Gradient::StopList>("GradientStopList", init<>())
    .def(vector_indexing_suite<Gradient::StopList>())
    .def("__repr__",sl_repr)
  ;

  class_<Gradient::Stop>("GradientStop", init<>())
    .def("GetColor", &Gradient::Stop::GetColor)
    .add_property("color", &Gradient::Stop::GetColor)
    .def("GetRel", &Gradient::Stop::GetRel)
    .add_property("rel", &Gradient::Stop::GetRel)
  ;

}
