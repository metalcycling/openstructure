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
#ifndef OST_COLOR_HH
#define OST_COLOR_HH

/*
  Author: Ansgar Philippsen, Marco Biasini
*/

#include <iosfwd>

#include <boost/operators.hpp>

#include <ost/gfx/module_config.hh>
#include <ost/geom/geom.hh>

namespace ost { namespace gfx {

/*!
  A color is defined in both RGB as well as HSV space, with an 
  additional A component. To initialize a color use the factory
  functions RGB(r,g,b), RGBA(r,g,b,a), HSV(h,s,v) and HSVA(h,s,v,a)

  All values are defined in the range of 0.0 to 1.0, but are not clamped
  on reading or writing, except when updating cross-domain, i.e. changes 
  to R,G or B cause HSV to be updated accordingly, and vice versa.
  This may cause previously out-of-bounds values to become clamped or wrapped
  to the defined domain. An out-of-bounds H will be wrapped to be between 0 and 1, 
  all others will be clamped.

  For interaction with OpenGL, an operator float* method exists that
  return a pointer to a float[4], in the order RGBA

  Operator overloading works on the RGB components, i.e. C1+C2, C+s, C*s
*/
class DLLEXPORT_OST_GFX Color :   
  private boost::additive<Color>,
  private boost::additive<Color, float>,
  private boost::multiplicative<Color, float>{

public:
  //! \brief initialize to white
  Color();

  //! convenience to set RGB triplet (ranges 0-1)
  void SetRGB(float r, float g, float b);
  //! convenience to get RGB triplet as Vec3
  geom::Vec3 GetRGB() const;
  //! convenience to get RGBA as Vec4
  geom::Vec4 GetRGBA() const;
  //! retrieve red
  float GetRed() const;
  //! set red
  void SetRed(float);
  //! retrieve green
  float GetGreen() const;
  //! set green
  void SetGreen(float);
  //! retrieve blue
  float GetBlue() const;
  //! set blue
  void SetBlue(float);
  //! convenience to set HSV triplet (ranges 0-1)
  void SetHSV(float h, float s, float v);
  //! convenience to get HSV triplet as Vec3
  geom::Vec3 GetHSV() const;
  //! convenience to get HSVA as Vec4
  geom::Vec4 GetHSVA() const;
  //! retrieve hue
  float GetHue() const;
  //! set hue
  void SetHue(float);
  //! retrieve saturatuin
  float GetSat() const;
  //! set saturation
  void SetSat(float);
  //! retrieve value
  float GetVal() const;
  //! set value
  void SetVal(float);
  //! retrieve alpha
  float GetAlpha() const;
  //! set alpha
  void SetAlpha(float);

  
  /*!
    \brief direct access to RGBA components

    In the context of a OpenGL call that requires a pointer
    to an RGB or RGBA float triplet or quadruplet, this will
    do the automatic casting from a Color object, i.e. 
    glColor3fv(my_color);

    This will also allow access to the RGBA values
    via the array notation, i.e. my_color[0]
  */
  operator const float* () const;
  operator float* ();
  
  Color& operator*=(float rhs);
  Color& operator+=(float rhs);
    
  Color& operator+=(const Color& rhs);  
  Color& operator-=(const Color& rhs);    
  Color& operator-=(float rhs);
  Color& operator/=(float rhs);  

  ////////////////////////////////////////////////////////////////

  //! DEPRECATED
  geom::Vec3 ToHSV() const {return GetHSV();}

  //! DEPRECATED
  Color(float r, float g, float b, float a=1.0);

  //! DEPRECATED
  float Red() const {return GetRed();}
  //! DEPRECATED
  float Green() const {return GetGreen();}
  //! DEPRECATED
  float Blue() const {return GetBlue();}
  //! DEPRECATED
  float Alpha() const {return GetAlpha();}

  //! DEPRECATED
  static Color FromRGB(unsigned char r, unsigned char g, 
                       unsigned char b, unsigned char a = 0xff) {
    static float f=1.0/255.0;
    return Color(f*static_cast<float>(r),f*static_cast<float>(g),
                 f*static_cast<float>(b),f*static_cast<float>(a));
  }

  ////////////////////////////////////////////////////////////////

private:
  void to_hsv() const;
  void to_rgb() const;

  mutable float rgba_[4];
  mutable float hsv_[3];
  mutable bool rgb_dirty_;
  mutable bool hsv_dirty_;
};

#undef RGB
/// \brief RGB color spec from floats (0.0-1.0)
Color DLLEXPORT_OST_GFX RGB(float r, float g, float b);

/// \brief RGB color spec from bytes (0-255)
Color DLLEXPORT_OST_GFX RGBb(uchar r, uchar g, uchar b);

/// \brief RGB color spec from integers (0-65535)
Color DLLEXPORT_OST_GFX RGBi(unsigned int r, unsigned int g, unsigned int b);

/// \brief RGBA color spec from floats (0.0-1.0)
Color DLLEXPORT_OST_GFX RGBA(float r, float g, float b, float a);

/// \brief RGBA color spec from bytes (0-255)
Color DLLEXPORT_OST_GFX RGBAb(uchar r, uchar g, uchar b, uchar a);

/// \brief RGBA color spec from integers (0-65535)
Color DLLEXPORT_OST_GFX RGBAi(unsigned int r, unsigned int g, unsigned int b, unsigned int a);


/*!
  \brief HSV color spec from floats

  h: Hue from 0 to 1 (0=red, 2/6=green, 4/6=blue)
     outside values are modulus 1, not clamped
  s: Saturation from 0 (no color) to 1 (full color)
     outside values are clamped
  v: Value from 0 (no light, black) to 1 (full light)
     outside values are clamped

  out-of-bounds values are perfectly fine, as they are
  sometimes needed in hsv gradients, e.g.

  g1=gfx.Gradient([gfx.HSV(-60.0/360.0,1,1),gfx.HSV(60.0/360.0,1,1)])
  g1.hsv_mode=True
  g2=gfx.Gradient([gfx.HSV(300.0/360.0,1,1),gfx.HSV(420.0/360.0,1,1)])
  g2.hsv_mode=True

  both will blend from purple to yellow via red
*/
Color DLLEXPORT_OST_GFX HSV(float h, float s, float v);

/*!
  \brief HSV color spec from integers

  h: Hue from 0 to 359 (0=red, 120=green, 240=blue)
     outside values are modulus 360, not clamped
  s: Saturation from 0 (no color) to 100 (full color)
     outside values are clamped
  v: Value from 0 (no light, black) to 100 (full light)
     outside values are clamped

  signed integers on purpose, allows negative values to be
  used for hsv gradients, see comments for HSV(float,float,float)
*/
  Color DLLEXPORT_OST_GFX HSVi(int h, int s, int v);

/*!
  \brief HSVA color spec from floats

  see HSV(float,float,float); alpha from 0 to 1
*/
Color DLLEXPORT_OST_GFX HSVA(float h, float s, float v, float a);

/*!
  \brief HSVA color spec from integers

  see HSVi(int, int, int); alpha from 0 to 100
*/
Color DLLEXPORT_OST_GFX HSVAi(int h, int s, int v, int a);


//! \brief string form
DLLEXPORT_OST_GFX std::ostream& operator<<(std::ostream&, const Color& c);

}}

#endif

