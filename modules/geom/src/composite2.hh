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
#ifndef GEOM_COMPOSITE2_HH
#define GEOM_COMPOSITE2_HH

#include "def.hh"

#include <vector>
#include <algorithm>

#include "circular_iterator.hh"
#include "vec2.hh"

namespace geom {

//! Line2
class DLLEXPORT_OST_GEOM Line2 {
public:
  Line2();
  Line2(const Vec2& from, const Vec2& to);

  Vec2 At(Real r) const;
  Vec2 GetOrigin() const;
  Vec2 GetDirection() const;
private:
  Vec2 ori_,dir_;
};

class DLLEXPORT_OST_GEOM Rectangle2{
public:
  Rectangle2();
  Rectangle2(Vec2 topleft, Vec2 bottomright);
  Real GetWidth() const;
  Real GetHeight() const;
  Vec2 GetStart() const;
  Vec2 GetEnd() const;
  void SetStart(const Vec2& v);
  void SetEnd(const Vec2& v);
  Real GetArea() const;
  Vec2 operator[] (unsigned index) const;
  Vec2& operator[] (unsigned index);
private:
  Vec2 topleft_;
  Vec2 bottomright_;
};

//! Polygon2
///\sa \ref create_split_image.py "Create Split Image Example"
class DLLEXPORT_OST_GEOM Polygon2: public std::vector<Vec2>
{
public:
  typedef const_circular_iter<Polygon2> const_circular_iterator;
  typedef circular_iter<Polygon2> circular_iterator;

    Polygon2();
    Polygon2(const Polygon2& p);
    Polygon2(const std::vector<Vec2>& v);
    circular_iterator cbegin(){return circular_iterator(begin(),end());}
    circular_iterator cend(){return circular_iterator(begin(),end(),--end());}
    const_circular_iterator cbegin() const {return const_circular_iterator(begin(),end());}
    const_circular_iterator cend() const {return const_circular_iterator(begin(),end(),--end());}
    void AddNode(const Vec2& n){push_back(n);};
    unsigned int GetNodeCount() const {return size();};
    Vec2 GetNode(unsigned int i) const {return operator[](i);};
    void SetNode(unsigned int i,const Vec2& v);
    circular_iterator FindSegment(std::vector<Vec2>::const_iterator start,std::vector<Vec2>::const_iterator end);
    circular_iterator FindNode(const Vec2& v) {return circular_iterator(begin(),end(),find(begin(),end(),v));}
    void Erase(circular_iterator first,circular_iterator last);
    Polygon2 operator+(const Vec2& v) const;
    Polygon2 operator*(Real d) const;
    Real GetArea() const;
    Vec2 GetCentroid() const;
    void Expand(Real val);
    Rectangle2 GetBoundingBox() const;
private:
        ;
};

class DLLEXPORT_OST_GEOM Ellipse2
{
public:
  Ellipse2();
  Ellipse2(Vec2 ori,Real a,Real  b,Real gamma);
  Vec2 At(Real t) const;
  Vec2 AtAngle(Real angle_) const;
  Rectangle2 GetBoundingBox() const;
  Real GetA() const {return a_;}
  Real GetB() const {return b_;}
  Real GetGamma() const {return gamma_;}
  Real GetArea() const {return a_*b_*M_PI;}
  Vec2 GetOrigin() const {return origin_;}
  void SetA(Real a) {a_=a;}
  void SetB(Real b) {b_=b;}
  void SetGamma(Real gamma) {gamma_=gamma;}
  void SetOrigin(Vec2 ori) {origin_=ori;}
private:
  Vec2 calc_(Real t) const;
  Vec2 origin_;
  Real a_;
  Real b_;
  Real gamma_;
};

class DLLEXPORT_OST_GEOM Hyperbola2
{
public:
  Hyperbola2();
  Hyperbola2(Vec2 ori,Real a,Real  b,Real gamma);
  Vec2 At(Real t, bool righthalf=true) const;
  Real GetA() const {return a_;}
  Real GetB() const {return b_;}
  Real GetGamma() const {return gamma_;}
  Vec2 GetOrigin() const {return origin_;}
  void SetA(Real a) {a_=a;}
  void SetB(Real b) {b_=b;}
  void SetGamma(Real gamma) {gamma_=gamma;}
  void SetOrigin(Vec2 ori) {origin_=ori;}
private:
  Vec2 origin_;
  Real a_;
  Real b_;
  Real gamma_;
};

class DLLEXPORT_OST_GEOM Circle2
{
public:
  Circle2();
  Circle2(const Circle2& c);
  Circle2(const Vec2& center, Real radius);
  void SetCenter(const Vec2& center);
  void SetRadius(Real r);
  Vec2 GetCenter() const;
  Real GetRadius() const;
  Real GetCircumference() const;
  Real GetArea() const;
protected:
  Vec2 center_;
  Real radius_;
};

} // ns

#endif

/// \example create_split_image.py
///
///This scripts loads one or more images and creates a split image where each of the images is displayed in a cone of equal angle.
///
///Usage:
///
/// \code dng create_split_image.py <input image 1> [<input image 2> <input image 3> .... ] <output image> \endcode
///
/// <BR>
/// <BR>


