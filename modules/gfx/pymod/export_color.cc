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
using namespace boost::python;

#include <ost/message.hh>
#include <ost/gfx/color.hh>
using namespace ost;
using namespace ost::gfx;

namespace {
  float get(const Color& c, int i) {
    if(i<0 || i>3) {
      throw Error("Color: index out of bounds");
    }
    return c[i];
  }

  void set(Color& c, int i, float v) {
    if(i<0 || i>3) {
      throw Error("Color: index out of bounds");
    }
    c[i]=v;
  }

  std::string repr(const Color& c) {
    std::ostringstream m;
    m << "gfx.RGBA(" << c.GetRed() << "," << c.GetGreen() << "," << c.GetBlue() << "," << c.GetAlpha() << ")";
    m << " gfx.HSV(" << c.GetHue() << "," << c.GetSat() << "," << c.GetVal() << ")";
    return m.str();
  }

  Color rgbh(uint rgb) {return RGB(static_cast<uchar>((rgb>>16)&0xff),
                                    static_cast<uchar>((rgb>>8)&0xff),
                                    static_cast<uchar>((rgb)&0xff));}
  Color rgbah(uint rgba) {return RGBA(static_cast<uchar>((rgba>>24)&0xff),
                                       static_cast<uchar>((rgba>>16)&0xff),
                                       static_cast<uchar>((rgba>>8)&0xff),
                                       static_cast<uchar>((rgba)&0xff));}

  tuple get_rgb(const Color& c) {return make_tuple(c.GetRed(),c.GetGreen(),c.GetBlue());}
  void set_rgb(Color& c, object rgb) {
    extract<geom::Vec3&> vec3(rgb);
    if(vec3.check()) {
      c.SetRGB(vec3()[0],vec3()[1],vec3()[2]);
    } else {
      // assume sequence
      c.SetRGB(extract<float>(rgb[0]),
               extract<float>(rgb[1]),
               extract<float>(rgb[2]));

    }
  }

  tuple get_rgba(const Color& c) {return make_tuple(c.GetRed(),c.GetGreen(),c.GetBlue(),c.GetAlpha());}
  void set_rgba(Color& c, object rgba) {
    extract<geom::Vec4&> vec4(rgba);
    if(vec4.check()) {
      c.SetRGB(vec4()[0],vec4()[1],vec4()[2]);
      c.SetAlpha(vec4()[3]);
    } else {
      // assume sequence
      c.SetRGB(extract<float>(rgba[0]),
               extract<float>(rgba[1]),
               extract<float>(rgba[2]));
      c.SetAlpha(extract<float>(rgba[3]));
    }
  }

  tuple get_hsv(const Color& c) {return make_tuple(c.GetHue(),c.GetSat(),c.GetVal());}
  void set_hsv(Color& c, object hsv) {
    extract<geom::Vec3&> vec3(hsv);
    if(vec3.check()) {
      c.SetHSV(vec3()[0],vec3()[1],vec3()[2]);
    } else {
      // assume sequence
      c.SetHSV(extract<float>(hsv[0]),
               extract<float>(hsv[1]),
               extract<float>(hsv[2]));
    }
  }

  tuple get_hsva(const Color& c) {return make_tuple(c.GetHue(),c.GetSat(),c.GetVal(),c.GetAlpha());}
  void set_hsva(Color& c, object hsva) {
    extract<geom::Vec4&> vec4(hsva);
    if(vec4.check()) {
      c.SetHSV(vec4()[0],vec4()[1],vec4()[2]);
      c.SetAlpha(vec4()[3]);
    } else {
      // assume sequence
      c.SetHSV(extract<float>(hsva[0]),
               extract<float>(hsva[1]),
               extract<float>(hsva[2]));
      c.SetAlpha(extract<float>(hsva[3]));
    }
  }

}

void export_color()
{
  class_<Color>("Color",init<>())
    .def(self_ns::str(self))
    .def("__repr__",repr)
    .add_property("r",&Color::GetRed,&Color::SetRed)
    .add_property("g",&Color::GetGreen,&Color::SetGreen)
    .add_property("b",&Color::GetBlue,&Color::SetBlue)
    .add_property("a",&Color::GetAlpha,&Color::SetAlpha)
    .add_property("h",&Color::GetHue,&Color::SetHue)
    .add_property("s",&Color::GetSat,&Color::SetSat)
    .add_property("v",&Color::GetVal,&Color::SetVal)
    .add_property("red",&Color::GetRed,&Color::SetRed)
    .add_property("green",&Color::GetGreen,&Color::SetGreen)
    .add_property("blue",&Color::GetBlue,&Color::SetBlue)
    .add_property("alpha",&Color::GetAlpha,&Color::SetAlpha)
    .add_property("hue",&Color::GetHue,&Color::SetHue)
    .add_property("sat",&Color::GetSat,&Color::SetSat)
    .add_property("val",&Color::GetVal,&Color::SetVal)
    .add_property("rgb",get_rgb,set_rgb)
    .add_property("rgba",get_rgba,set_rgba)
    .add_property("hsv",get_hsv,set_hsv)
    .add_property("hsva",get_hsva,set_hsva)
    .def("__getitem__",get)
    .def("__setitem__",set)
    // DEPRECATED
    .def(init<float, float, float, optional<float> >())
    .def("Red",&Color::GetRed)
    .def("Green",&Color::GetGreen)
    .def("Blue",&Color::GetBlue)
    .def("Alpha",&Color::GetAlpha)
    .def("ToHSV",&Color::ToHSV)
    .def("FromRGBA",&Color::FromRGB)
    ;

  def("RGB",RGB);
  def("RGBh",rgbh);
  def("RGBb",RGBb);
  def("RGBi",RGBi);
  def("RGBA",RGBA);
  def("RGBAh",rgbah);
  def("RGBAb",RGBAb);
  def("RGBAi",RGBAi);
  def("HSV",HSV);
  def("HSVi",HSVi);
  def("HSVA",HSVA);
  def("HSVAi",HSVAi);
}
