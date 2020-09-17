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
  data observer

  Author: Ansgar Philippsen
*/

#ifndef IMG_DATA_OBSERVER
#define IMG_DATA_OBSERVER


#include "image_handle.hh"

namespace ost { namespace img {


class DLLEXPORT InvalidObserver: public Error {
public:
  InvalidObserver(const String& s = String("unknown")):
    Error(String("InvalidObserver exception occurred: ") + s)
  {}
};

//! Abstract base class of data observer
/*
  Based on the pattern of same name (293). The expression Observer or
  Observed is part of the method names in purpose, in order to avoid
  conflict with other interfaces (such as in wxWidgets)

  ObserverUpdate() and ObserverRelease() are the abstract methods and
  must be implemented by a derived class. It is important to note that
  the reference to Data given to the ctor should not be saved in the
  derived class, but rather the GetObservedData method should be used.
  
  If the observed data goes out-of-scope, a call to ObserverRelease will
  be followed by a call to ObserverInvalidate, which will mark the data
  reference as invalidated. A further call to GetObserved Data will throw
  an InvalidObserver exception!
*/
class DLLEXPORT_OST_IMG_BASE DataObserver {
 public:
  //! constructor
  /*!
    requires reference to data, then attaches itself to data
  */
  DataObserver(const ImageHandle& d);

  DataObserver(const DataObserver& o);

  DataObserver& operator=(const DataObserver& o);

  virtual ~DataObserver();

  //! called upon data change, noop by default
  virtual void ObserverUpdate();

  //! special update, calls ObserverUpdate() by default
  virtual void ObserverUpdate(const Extent&);

  //! special update, calls ObserverUpdate() by default
  virtual void ObserverUpdate(const Point&);

  //! Called if data is released from memory
  virtual void ObserverRelease() = 0;

  //! Called after ObserverRelease
  /*
    This will automatically invalidate the observer so that
    a next call to GetData() will throw an InvalidObserver
    exception
  */
  void ObserverInvalidate();
  
  //! Returns data that is observed
  virtual const ImageHandle& GetObservedData() const;

  bool IsDataValid() const {return is_valid();}

 protected:

  //! Set new data to observe
  /*!
    Automatically unregister previously observed data
  */
  void SetObservedData(const ImageHandle& d);

  bool is_valid() const;

private:
  ImageHandle data_;

};

}} // namespace img

#endif

