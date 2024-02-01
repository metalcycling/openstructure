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
  iterate over a volume defined by an extent

  Author: Ansgar Philippsen
*/

#ifndef IMG_EXTENT_ITERATOR_H
#define IMG_EXTENT_ITERATOR_H

#include "point.hh"
#include "data.hh"

namespace ost { namespace img {

class Extent;

//! Extent iterator
/*!
  The extent iterator is implicitely convertible to Point, 
  and is therefore ideally suited for Get/Set methods which
  take a point as an argument.

  \code
  ImageHandle ih = get_some_image();
  for(ExtentIterator it(ih.GetExtent());!it.AtEnd();++it) {
    Real r = ih.GetReal(it);
    // do something with value
  }
  \endcode
  \sa \ref view_phase_diff.py "View Phase Diference Example"
*/
class DLLEXPORT_OST_IMG_BASE ExtentIterator {
 public:
  //! Initialize with Extent
  explicit ExtentIterator(const Extent& b, unsigned int domain=SPATIAL);

  //! copy ctor
  ExtentIterator(const ExtentIterator &i);

  //! assignement op
  ExtentIterator& operator=(const ExtentIterator& p);

  //! Prefix increase operator
  ExtentIterator& operator++();

  //! Postfix increase operator
  ExtentIterator operator++(int);

  //! Prefix decrease operator
  ExtentIterator& operator--(); 

  //! Postfix decrease operator
  ExtentIterator operator--(int);

  //! Skip ahead to next point
  void Inc();

  //! Go back to previous point
  void Dec();

  //! Reset to start
  ExtentIterator& ToStart();

  //! Reset to end
  ExtentIterator& ToEnd();

  //! Retrieve original extent
  Extent GetExtent() const;

  //! Retrieve start point
  const Point &GetStart() const;

  //! Retrieve end point
  const Point &GetEnd() const;

  //! Returns true if iterator has reached end
  bool AtEnd() const;

  //! Returns true if iterator is at the start
  bool AtStart() const;

  //! element access
  int operator[](int i) const;

  operator Point() const;

 private:
  Point start_, end_, current_,startpos_;
};

}} // namespace img

#endif
