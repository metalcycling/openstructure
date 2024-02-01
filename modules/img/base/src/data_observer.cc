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
  abstract observer implementation

  Author: Ansgar Philippsen
*/
#include <iostream>

#include "data_observer.hh"
#include "null_data.hh"

namespace ost { namespace img {

DataObserver::DataObserver(const ImageHandle& d):
  data_(d)
{
  data_.Attach(this);
}

DataObserver::DataObserver(const DataObserver& o):
  data_(o.data_) 
{
  data_.Attach(this);
}

DataObserver& DataObserver::operator=(const DataObserver& o)
{
  if(this!=&o) {
    data_.Detach(this);
    data_=o.data_;
    data_.Attach(this);
  }
  return *this;
}

DataObserver::~DataObserver() 
{
  data_.Detach(this);
}

void DataObserver::ObserverInvalidate()
{
  data_=ImageHandle();
}

const ImageHandle& DataObserver::GetObservedData() const
{
  if(!is_valid()) {
    std::cerr << "invalid observed data access" << std::endl;
    throw InvalidObserver("GetData called for invalidated observer");
  }
  return data_;
}

void DataObserver::SetObservedData(const ImageHandle& d)
{
  if(is_valid()) {
    data_.Detach(this);
  }
  data_ = d;
  data_.Attach(this);
}

void DataObserver::ObserverUpdate()
{
  // noop
}

void DataObserver::ObserverUpdate(const Extent&)
{
  this->ObserverUpdate();
}

void DataObserver::ObserverUpdate(const Point&)
{
  this->ObserverUpdate();
}

bool DataObserver::is_valid() const 
{
  return data_.IsValid();
}

}} // namespace
