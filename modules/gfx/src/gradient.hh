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
#ifndef OST_GFX_GRADIENT_HH
#define OST_GFX_GRADIENT_HH

/*
  Authors: Marco Biasini, Ansgar Philippsen
*/

#include <vector>

#include "color.hh"

#include <ost/info/info.hh>
#include <ost/info/info_fw.hh>

namespace ost { namespace gfx {

/*!
 \brief color gradient
 
 Gradients map a scalar value in the range of 0 to 1 to a color. The 
 gradient consists of zero or more color stops. These stop points define the 
 output of GetColorAt(). New stops can be added with SetColorAt().

 Usually, the gradient operates in RGB colorspace, which means a linear
 interpolation of each color component between stops. As an alternative,
 interpolation can also be done in HSV colorspace, by setting HSV mode
 to true with SetHSVMode(). However, carefully using HSV colorspace will
 yield better results, if some potential caveats are taken into consideration.
 The best example is a simple gradient from red to blue. The half-point in RGB
 colorspace is (0.5,0.0,0.5), so a darkish purple. Converting to HSV one can see
 that the saturation is at 50%, which is unfortunate. Using hsv mode can circument
 this issue, but in the above example, naively using hsv mode will give green
 at the half-point and not purple. Why? Because red translates to HSV(0,1,1) and blue to
 HSV(2/3,1,1), which at the half-point is HSV(1/3,1,1), or green. To interpolate
 on the hue in the other direction, use HSV and unwrapped hue values, e.g. 
 [HSV(0.0,1,1),HSV(-1/3,1,1)] or [HSV(1,1,1),HSV(2/3,1,1)] - both of these will
 give a half-point or a nice purple with full saturation
 
 \sa \ref gradient.py "Gradient Example"
*/
class DLLEXPORT_OST_GFX Gradient {
public:
  struct Stop {
    Stop(): t(0.0) {}
    
    Stop(float tt, const Color& c): t(tt), color(c) {}
    float t;
    Color color;
    
    Color GetColor() const {
      return color;
    }
    
    float GetRel() const {
      return t;
    }

    bool operator==(const Stop& ref) const
    {
      return t==ref.t && color==ref.color;
    }
    
  };
  
  typedef std::vector<Stop> StopList;
  
public:
  /*!
    In python, the gradient can also be initialize with a dictionary, mapping
    stops to either float triplets or gfx.Color objects
  */
  Gradient();

  /// \brief initialize with a pre-define gradient name
  Gradient(const String& name);

  /// \brief get color
  Color GetColorAt(float t) const;
  /// \brief set color
  /// 
  /// \param t is the position where the new \ref Stop "stop" should be inserted
  /// \param color is the color of the stop
  void SetColorAt(float t, const Color& color);  
  
  /// \brief convert gradient to info
  // TODO: rename to ToInfo
  // TODO:: add stand-alone function named GradientToInfo(Gradient, InfoGroup)
  void GradientToInfo(info::InfoGroup& group) const;

  /// \brief get list of color stops
  StopList GetStops() const;
  
  /// \brief create gradient from info
  // TODO: move outside
  // TODO: add FromInfo member function that overwrites this
  static gfx::Gradient GradientFromInfo(const info::InfoGroup& group);

  /// \brief turn HSV interpolation mode on or off
  void SetHSVMode(bool m) {hsv_mode_=m;}
  /// \brief return the HSV interpolation mode flag
  bool GetHSVMode() const {return hsv_mode_;}
    
private:
  std::vector<Stop>  stops_;
  bool hsv_mode_;
};


///  \example gradient.py
///
/// Draw a set of lines and color them using a gradient
}}

#endif
