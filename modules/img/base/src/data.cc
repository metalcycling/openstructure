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
  abstract data implementation

  Author: Ansgar Philippsen
*/

#include <ost/message.hh>

#include "data.hh"
#include "data_observer.hh"
#include "extent_iterator.hh"
#include "pixel_sampling.hh"

namespace ost { namespace img {;

// ConstData

ConstData::ConstData()
{
}

ConstData::ConstData(const ConstData& d)
{
}
  
ConstData::~ConstData()
{
}


ConstData& ConstData::operator=(const ConstData& d)
{
  return *this;
}




Vec3 ConstData::GetPixelSampling() const
{
  return Sampling().GetPixelSampling();
}


Vec3 ConstData::GetSpatialSampling() const
{
  return Sampling().GetSpatialSampling();
}

Vec3 ConstData::GetFrequencySampling() const
{
  return Sampling().GetFrequencySampling();
}

//Vec3 ConstData::IndexToCoord(const Point &i) const
//{
//  return Sampling().Point2Coord(i);
//}

//Vec3 ConstData::CoordToIndex(const Vec3& c) const
//{
//  return Sampling().Coord2Point(c);
//}

Real ConstData::OverallDifference(const ConstData& ref) const
{
  Real sum=0.0;
  for(ExtentIterator it(GetExtent()); !it.AtEnd(); ++it) {
    Real d=this->GetReal(it) - ref.GetReal(it);
    sum+=d*d;
  }
  return sqrt(sum)/(Real)GetExtent().GetSize().GetVol();
}

Real ConstData::NormDifference(const ConstData& ref) const
{
  // TODO: specialize for real and complex data...
  Real sum1=0.0;
  Real sum2=0.0;
  for(ExtentIterator it(GetExtent()); !it.AtEnd(); ++it) {
    sum1+=this->GetReal(it);
    sum2+=ref.GetReal(it);
  }
  sum1/=(Real)GetExtent().GetSize().GetVol();
  sum2/=(Real)GetExtent().GetSize().GetVol();

  Real sum=0.0;
  for(ExtentIterator it(GetExtent()); !it.AtEnd(); ++it) {
    Real d=(this->GetReal(it)-sum1) - (ref.GetReal(it)-sum2);
    sum+=d*d;
  }
  return sqrt(sum)/(Real)GetExtent().GetSize().GetVol();
}




// Data

Data::Data():
  ConstData()
{
}


Data::Data(const Data& d):
  ConstData(d)
{
}
  

Data::~Data()
{
}

Data& Data::operator=(const Data& d)
{
  ConstData::operator=(d);
  return *this;
}


void Data::SetPixelSampling(const Vec3& dim)
{
  Sampling().SetPixelSampling(dim);
  Notify();
}

void Data::SetSpatialSampling(const Vec3& dim)
{
  Sampling().SetSpatialSampling(dim);
  Notify();
}

}} //ns
