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
  Author: Ansgar Philippsen, Marco Biasini
*/

#include <ostream>

#include "color.hh"

namespace ost { namespace gfx {

Color::Color()
{
  rgba_[0]=1.0;
  rgba_[1]=1.0;
  rgba_[2]=1.0;
  rgba_[3]=1.0;
  rgb_dirty_=false;
  hsv_dirty_=true;
}

Color::Color(float r, float g, float b, float a) 
{
  rgba_[0]=r;
  rgba_[1]=g;
  rgba_[2]=b;
  rgba_[3]=a;
  rgb_dirty_=false;
  hsv_dirty_=true;
}
  
void Color::SetRGB(float r, float g, float b)
{
  rgba_[0]=r;
  rgba_[1]=g;
  rgba_[2]=b;
  rgb_dirty_=false;
  hsv_dirty_=true;
}

geom::Vec3 Color::GetRGB() const
{
  if(rgb_dirty_) to_rgb();
  return geom::Vec3(rgba_);
}

geom::Vec4 Color::GetRGBA() const
{
  if(rgb_dirty_) to_rgb();
  return geom::Vec3(rgba_);
}

float Color::GetRed() const
{
  if(rgb_dirty_) to_rgb();
  return rgba_[0];
}

void Color::SetRed(float x)
{
  if(rgb_dirty_) to_rgb();
  rgba_[0]=x;
  hsv_dirty_=true;
}

float Color::GetGreen() const
{
  if(rgb_dirty_) to_rgb();
  return rgba_[1];
}

void Color::SetGreen(float x)
{
  if(rgb_dirty_) to_rgb();
  rgba_[1]=x;
  hsv_dirty_=true;
}

float Color::GetBlue() const
{
  if(rgb_dirty_) to_rgb();
  return rgba_[2];
}

void Color::SetBlue(float x)
{
  if(rgb_dirty_) to_rgb();
  rgba_[2]=x;
  hsv_dirty_=true;
}

void Color::SetHSV(float h, float s, float v)
{
  hsv_[0]=h;
  hsv_[1]=s;
  hsv_[2]=v;
  hsv_dirty_=false;
  rgb_dirty_=true;
}

geom::Vec3 Color::GetHSV() const
{
  if(hsv_dirty_) to_hsv();
  return geom::Vec3(hsv_);
}

geom::Vec4 Color::GetHSVA() const
{
  if(hsv_dirty_) to_hsv();
  return geom::Vec4(hsv_[0],hsv_[1],hsv_[2],rgba_[3]);
}

float Color::GetHue() const
{
  if(hsv_dirty_) to_hsv();
  return hsv_[0];
}

void Color::SetHue(float x)
{
  if(hsv_dirty_) to_hsv();
  hsv_[0]=x;
  rgb_dirty_=true;
}

float Color::GetSat() const
{
  if(hsv_dirty_) to_hsv();
  return hsv_[1];
}

void Color::SetSat(float x)
{
  if(hsv_dirty_) to_hsv();
  hsv_[1]=x;
  rgb_dirty_=true;
}

float Color::GetVal() const
{
  if(hsv_dirty_) to_hsv();
  return hsv_[2];
}

void Color::SetVal(float x)
{
  if(hsv_dirty_) to_hsv();
  hsv_[2]=x;
  rgb_dirty_=true;
}

float Color::GetAlpha() const
{
  return rgba_[3];
}

void Color::SetAlpha(float x)
{
  rgba_[3]=x;
}

////////////////////////////////////
// operators

Color::operator const float* () const 
{
  if(rgb_dirty_) to_rgb();
  return rgba_;
}

Color::operator float* () 
{
  if(rgb_dirty_) to_rgb();
  hsv_dirty_=true; 
  return rgba_;
}

Color& Color::operator*=(float rhs)
{
  if(rgb_dirty_) to_rgb();
  rgba_[0]*=rhs;
  rgba_[1]*=rhs;
  rgba_[2]*=rhs;
  rgba_[3]*=rhs; 
  hsv_dirty_=true; 
  return *this;
}

Color& Color::operator+=(float rhs)
{
  if(rgb_dirty_) to_rgb();
  rgba_[0]+=rhs;
  rgba_[1]+=rhs;
  rgba_[2]+=rhs;
  rgba_[3]+=rhs;
  hsv_dirty_=true; 
  return *this;
}
  
Color& Color::operator+=(const Color& rhs)
{
  if(rgb_dirty_) to_rgb();
  rgba_[0]+=rhs[0];
  rgba_[1]+=rhs[1];
  rgba_[2]+=rhs[2];
  rgba_[3]+=rhs[3];  
  hsv_dirty_=true; 
  return *this;  
}

Color& Color::operator-=(const Color& rhs)
{
  if(rgb_dirty_) to_rgb();
  rgba_[0]-=rhs[0];
  rgba_[1]-=rhs[1];
  rgba_[2]-=rhs[2];
  rgba_[3]-=rhs[3];  
  hsv_dirty_=true; 
  return *this;
}

Color& Color::operator-=(float rhs)
{
  if(rgb_dirty_) to_rgb();
  rgba_[0]-=rhs;
  rgba_[1]-=rhs;
  rgba_[2]-=rhs;
  rgba_[3]-=rhs;
  hsv_dirty_=true; 
  return *this;
}

Color& Color::operator/=(float rhs)
{
  if(rgb_dirty_) to_rgb();
  rgba_[0]/=rhs;
  rgba_[1]/=rhs;
  rgba_[2]/=rhs;
  rgba_[3]/=rhs;  
  hsv_dirty_=true; 
  return *this;
}

//////////////////////////////////////
// private methods

void Color::to_rgb() const
{
  float hh=fmod(hsv_[0],1.0f);
  if(hh<0.0) hh+=1.0;
  float ss=std::min(1.0f,std::max(0.0f,hsv_[1]));
  float vv=std::min(1.0f,std::max(0.0f,hsv_[2]));
  if (ss<1e-9){
    rgba_[0]=vv;
    rgba_[1]=vv;
    rgba_[2]=vv;
  } else {
    float var_h=hh*6.0;
    int var_i =static_cast<int>(var_h);
    float var_1 = vv*(1-ss);
    float var_2 = vv*(1-ss*( var_h - static_cast<float>(var_i)));
    float var_3 = vv*(1-ss*(1-(var_h-static_cast<float>(var_i))));
    switch(var_i){
    case 0:
      rgba_[0]=vv;
      rgba_[1]=var_3;
      rgba_[2]=var_1;
      break;
    case 1:
      rgba_[0] = var_2;
      rgba_[1] = vv;
      rgba_[2] = var_1;
      break;
    case 2:
      rgba_[0] = var_1;
      rgba_[1] = vv;
      rgba_[2] = var_3;
      break;
    case 3:
      rgba_[0] = var_1 ; 
      rgba_[1] = var_2 ; 
      rgba_[2] = vv;
      break;
    case 4:
      rgba_[0] = var_3 ;
      rgba_[1] = var_1 ;
      rgba_[2] = vv;
      break;
    case 5:
      rgba_[0] = vv;
      rgba_[1] = var_1 ;
      rgba_[2] = var_2;
      break;
    }
  }
  rgb_dirty_=false;
}

void Color::to_hsv() const
{
  float rr = std::min(1.0f,std::max(0.0f,rgba_[0]));
  float gg = std::min(1.0f,std::max(0.0f,rgba_[1]));
  float bb = std::min(1.0f,std::max(0.0f,rgba_[2]));

  float vmin = std::min(std::min(rr, gg), bb);
  float vmax = std::max(std::max(rr, gg), bb);
  float vdelta = vmax - vmin;

  hsv_[2] = vmax;

  if (vdelta < 1.0e-9) {
    hsv_[0] = 0.0;
    hsv_[1] = 0.0;
  } else {
    hsv_[1] = vdelta / vmax;
    float dr = (((vmax - rr) / 6.0) + (vdelta / 2.0)) / vdelta;
    float dg = (((vmax - gg) / 6.0) + (vdelta / 2.0)) / vdelta;
    float db = (((vmax - bb) / 6.0) + (vdelta / 2.0)) / vdelta;
    
    if (rr == vmax) {
      hsv_[0] = db - dg;
    } else if (gg == vmax) {
      hsv_[0] = (1.0 / 3.0) + dr - db;
    } else if (bb == vmax) {
      hsv_[0] = (2.0 / 3.0) + dg - dr;
    }
    if (hsv_[0] < 0.0) {
      hsv_[0] += 1.0;
    }
    if (hsv_[0] > 1.0) {
      hsv_[0] -= 1.0;
    }
  }
  hsv_dirty_=false;
}

////////////////////////////////////
// stand alone functions

Color RGB(float r, float g, float b)
{
  Color nrvo;
  nrvo.SetRGB(r,g,b);
  return nrvo;
}

Color RGBb(uchar r, uchar g, uchar b)
{
  static float f=1.0/255.0;
  Color nrvo;
  nrvo.SetRGB(static_cast<float>(r)*f,static_cast<float>(g)*f,static_cast<float>(b)*f);
  return nrvo;
}

Color RGBi(unsigned int r, unsigned int g, unsigned int b)
{
  static float f=1.0/65535.0;
  Color nrvo;
  nrvo.SetRGB(static_cast<float>(r)*f,static_cast<float>(g)*f,static_cast<float>(b)*f);
  return nrvo;
}

Color RGBA(float r, float g, float b, float a)
{
  Color nrvo;
  nrvo.SetRGB(r,g,b);
  nrvo.SetAlpha(a);
  return nrvo;
}

Color RGBAb(uchar r, uchar g, uchar b, uchar a)
{
  static float f=1.0/255.0;
  Color nrvo;
  nrvo.SetRGB(static_cast<float>(r)*f,static_cast<float>(g)*f,static_cast<float>(b)*f);
  nrvo.SetAlpha(static_cast<float>(a)*f);
  return nrvo;
}

Color RGBAi(unsigned int r, unsigned int g, unsigned int b, unsigned a)
{
  static float f=1.0/65535.0;
  Color nrvo;
  nrvo.SetRGB(static_cast<float>(r)*f,static_cast<float>(g)*f,static_cast<float>(b)*f);
  nrvo.SetAlpha(static_cast<float>(a)*f);
  return nrvo;
}

Color HSV(float h, float s, float v)
{
  Color nrvo;
  nrvo.SetHSV(h,s,v);
  return nrvo;
}

Color HSVi(int r, int g, int b)
{
  static float f=1.0/65535.0;
  Color nrvo;
  nrvo.SetHSV(static_cast<float>(r)*f,static_cast<float>(g)*f,static_cast<float>(b)*f);
  return nrvo;
}

Color HSVA(float h, float s, float v, float a)
{
  Color nrvo;
  nrvo.SetHSV(h,s,v);
  nrvo.SetAlpha(a);
  return nrvo;
}

Color HSVAi(int r, int g, int b, int a)
{
  static float f=1.0/65535.0;
  Color nrvo;
  nrvo.SetHSV(static_cast<float>(r)*f,static_cast<float>(g)*f,static_cast<float>(b)*f);
  return nrvo;
}

std::ostream& operator<<(std::ostream& s, const Color& c)
{
  s << "Color(r=" << c.GetRed() << ",g=" <<c.GetGreen() << ",b=" << c.GetBlue() << ",h=" << c.GetHue() << ",s=" << c.GetSat() << ",v=" << c.GetVal() << ",a=" << c.GetAlpha() << ")";
  return s;
}

}} // ns
