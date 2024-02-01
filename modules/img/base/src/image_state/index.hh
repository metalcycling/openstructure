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
  value holder index

  Author: Ansgar Philippsen
*/

#ifndef IMG_INDEX_H
#define IMG_INDEX_H

#include <iosfwd>

#include <ost/base.hh>
#include <ost/img/size.hh>

namespace ost { namespace img { namespace image_state {

// encapsulated unsigned integer triplet
/*
  provides a compact mean to adress values in a ValueHolder
*/
struct DLLEXPORT_OST_IMG_BASE Index {
  Index(unsigned int uu,unsigned int vv, unsigned int ww):
    u(uu),v(vv),w(ww)
  {}

  bool equal(const Index& i) const;

  unsigned int u,v,w;
};

DLLEXPORT_OST_IMG_BASE bool operator==(const Index& i1, const Index& i2);
DLLEXPORT_OST_IMG_BASE bool operator!=(const Index& i1, const Index& i2);

DLLEXPORT_OST_IMG_BASE std::ostream& operator<<(std::ostream& out, const Index& i);

class DLLEXPORT_OST_IMG_BASE IndexIterator {
public:
  IndexIterator(int w, int h, int d);
  IndexIterator(const Size& s);

  bool AtStart() const;
  bool AtEnd() const;

  IndexIterator& ToStart();
  IndexIterator& ToEnd();

  // prefix
  IndexIterator& operator++();
  // postfix
  IndexIterator operator++(int);
  // prefix
  IndexIterator& operator--();
  // postfix
  IndexIterator operator--(int);

  operator Index() const;
private:
  unsigned int w_,h_,d_;
  Index indx_;
  bool at_end_;

  void increment();
  void decrement();
};

}}} // namespaces

#endif
