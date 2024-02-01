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
  Authors: Ansgar Philippsen, Giani Signorell
*/

#ifndef IMG_ALG_FILTER_LINE_ITERATOR_H
#define IMG_ALG_FILTER_LINE_ITERATOR_H

#include <ost/img/alg/module_config.hh>
#include <ost/img/extent.hh>

namespace ost { namespace img {

class ExtentIterator;

namespace alg { 

//! line iterator generator
/*!
  generates consecutive lines over given extent
*/
class DLLEXPORT_IMG_ALG LineIterator {
public:
  //! Initialization with extent and flag
  /*! Flag determines which axes are scanned by the iterator:
      (use the img::geometry::axis constants) 
  */


  LineIterator(const Extent& ext, unsigned int flag= (Axis::X | Axis::Y | Axis::Z) );
  //! return true if iterator has reached end
  bool AtEnd();
  //! increase method
  void Inc();
  //! retrieve current line (in form of an ExtentIterator)
  ExtentIterator GetLine();
  //! Retrieve Size of current line
  int GetSize() const;

  //! Prefix increase operator
  LineIterator& operator++();
  
  //! Postfix increase operator
  LineIterator operator++(int);

  //! convenience method, but preferred!
  ExtentIterator GetNextLine();
  
private:

  unsigned int flag_;
  Extent ext_;
  int current_dir_;
  Point start_,end_,indx_;
};


}}} //namespaces

#endif
