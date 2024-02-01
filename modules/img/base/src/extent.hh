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
  Extent

  Author: Ansgar Philippsen
*/

#ifndef IMG_EXTENT_H
#define IMG_EXTENT_H

#include <ost/message.hh>
#include <iosfwd>
#include <ost/message.hh>

#include "point.hh"
#include "size.hh"

/*
  TODO:

  add CenterOrigin() -> return extent who's (0,0,0) is in the center
  add TLFOrigin()    -> return extent with start=(0,0,0)
  add operators      -> add/sub with point and extent
*/

namespace ost { namespace img {

struct DLLEXPORT InvalidExtentException: public Error
{
  InvalidExtentException(const String& m= String("") ):
    Error(String("an invalid extent occurred ("+m+String(")")))
  {}
};

//! Defines lower and upper valid indices
/*!
  Is part of the basic characteristic of each Data instance.
  Encapsulates a start and end index. Size and Dimension are 
  directly extracted from Extent.
*/
class DLLEXPORT_OST_IMG_BASE Extent {
 public:
  //! default constructor
  Extent();

  //! copy constructor
  Extent(const Extent &r);

  //! initialize with start and end point
  Extent(const Point& p1, const Point& p2);

  //! initialize with Size
  /*!
    The extent will go from -(size-1)/2 to (size-1)/2 for odd sizes,
    and from -(size/2-1) to size/2 for even images

    not explicit on purpose, facilitates easier creation code
  */
  Extent(const Size& s);

  //! initialize with center point and size
  /*
    The given point will be the in the center of the
    extent, with the size according to the even/odd rule
    as described for Extent(const Size& s)
  */
  Extent(const Size& size, const Point& center);

  //! initialize with Size and start point
  /*
    The extent will go from start to (size-1)+start
  */
  Extent(const Point& start, const Size& size);

  //! Return lower/left/front corner
  const Point& GetStart() const;
  //! Return upper/right/back corner
  const Point& GetEnd() const;

  //! Set start point, changing size
  void SetStart(const Point& o);
  //! Set end point, changing size
  void SetEnd(const Point& e);

  //! Returns true if the given point is within extent
  bool Contains(const Point& p) const;
  //! Returns true if the given exten is within extent
  bool Contains(const Extent& e) const;

  //! return center
  Point GetCenter() const;

  //! Return size of extent
  Size GetSize() const;
  int GetWidth() const;
  int GetHeight() const;
  int GetDepth() const;
  int GetVolume() const;

  //! Return dimension
  int GetDim() const;

  //! Wrap point around
  /* 
     Ensures that the given point lies within Extent, 
     using wrap around if necessary
   */
  Point WrapAround(const Point& p);  

  //! Return new extent mirrored according to planes
  Extent Mirror(int planes);

  //! Generates a continues, uniqe index per point contained within extent
  unsigned int Point2Offset(const Point& p);

  // operators
  bool operator==(const Extent& b) const {return equal(b);}
  bool operator!=(const Extent& b) const {return !equal(b);}
  
  void Shift(const Point& p);

 private:
  Point start_,end_;

  void set(const Point& p1, const Point& p2);
  bool equal(const Extent& b) const;

};

// 'global' functions using extent

DLLEXPORT_OST_IMG_BASE bool HasOverlap(const Extent& e1, const Extent& e2);
DLLEXPORT_OST_IMG_BASE Extent Overlap(const Extent& e1, const Extent& e2);

DLLEXPORT_OST_IMG_BASE std::ostream& operator<<(std::ostream& os, const img::Extent &b);

}} // namespace img


#endif
