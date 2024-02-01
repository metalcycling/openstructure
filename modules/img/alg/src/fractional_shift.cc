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
  Authors: Giani Signorell, Andreas Schenk
*/

#include <cmath>
#include <ost/img/alg/fft.hh>
#include "fractional_shift.hh"

namespace ost { namespace img { namespace alg {

FractionalShift::FractionalShift(Real sx, Real sy, Real sz):
  ConstModIPAlgorithm("FractionalShift"),
  shift_(sx,sy,sz)
{
}
FractionalShift::FractionalShift(const Vec3& v):
  ConstModIPAlgorithm("FractionalShift"),
  shift_(v)
{
}

void FractionalShift::SetShift(Real sx, Real sy, Real sz)
{
  SetShift(Vec3(sx,sy,sz));
}

void FractionalShift::SetShift(const Vec3& v)
{
  shift_=v;
}

Vec3 FractionalShift::GetShift() const
{
  return shift_;
}

void FractionalShift::Visit(ImageHandle& ih) const
{  
  Vec3 ao=ih.GetAbsoluteOrigin();

  bool spatialflag=false;

  if(ih.GetDomain()==SPATIAL){
    spatialflag=true;
    ih.ApplyIP(FFT());
  }
  
  //calculate shift
  Size s=ih.GetSize();
  Vec3 phaseshift = 2.0 * M_PI * CompDivide(shift_,Vec3(static_cast<Real>(s[0]),
							static_cast<Real>(s[1]),
							static_cast<Real>(s[2])));

  //calculate shifted image
  for(ExtentIterator it=ih.GetIterator();!it.AtEnd();++it) {
    ih.SetComplex(it,ih.GetComplex(it)*std::polar<Real>(1.0,-Dot(phaseshift,Point(it).ToVec3())));
  }

  if(spatialflag){
    ih.ApplyIP(FFT());
  }
  // update spatial origin
  ih.SetAbsoluteOrigin(ao-shift_);
}

}}}//ns
