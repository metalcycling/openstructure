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
  unsigned integer triplet based size

  Author: Ansgar Philippsen
*/

#ifndef IMG_SIZE_H
#define IMG_SIZE_H

#include <iosfwd>

#include <ost/img/module_config.hh>

namespace ost { namespace img {

class Point;

//! class encapsulating 1D to 3D size
class DLLEXPORT_OST_IMG_BASE Size {
 public:
  Size();
  Size(const Size &s);

  //! 1D constructor
  explicit Size(unsigned int w);

  //! 2D constructor
  explicit Size(unsigned int w, unsigned int h);

  //! 3D constructor
  explicit Size(unsigned int w, unsigned int h, unsigned int d);

  //! initialize from two points
  Size(const Point& p1, const Point& p2);

  //! return the point at half the size
  /*!
    for even dimension, its N/2-1, for odd dimension its (N-1)/2
  */
  Point GetHalf() const;
  
  //! get width component (convenience to index 0)
  unsigned int GetWidth() const;
  //! get height component (convenience to index 1)
  unsigned int GetHeight() const;
  //! get depth component (convenience to index 2)
  unsigned int GetDepth() const;
  //! get area (width * height)
  unsigned int GetSlab() const;
  //! get volume (depcrecated, use GetVolume() instead)
  unsigned int GetVol() const;
  //! get volume
  unsigned int GetVolume() const;

  //! return dimensionality (1-3)
  unsigned int GetDim() const;

  // some operators
  Size& operator=(const Size& s);
  unsigned int operator[](unsigned int index) const;
  unsigned int& operator[](unsigned int index);
  bool operator==(const Size& s) const;
  bool operator!=(const Size& s) const;

  Size& operator+=(const Size& s);
  Size& operator-=(const Size& s);

 private:
  unsigned int w_,h_,d_;

  bool equal(const Size& s) const;
};

DLLEXPORT_OST_IMG_BASE Size operator+(const Size& s1, const Size& s2);
DLLEXPORT_OST_IMG_BASE Size operator-(const Size& s1, const Size& s2);

DLLEXPORT_OST_IMG_BASE std::ostream& operator<<(std::ostream& os, const Size& s);


}} // namespace img


#endif
