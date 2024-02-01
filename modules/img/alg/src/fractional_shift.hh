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

#ifndef FRACTIONAL_SHIFT_HH_
#define FRACTIONAL_SHIFT_HH_

#include <ost/img/algorithm.hh>

namespace ost { namespace img { namespace alg {

class DLLEXPORT_IMG_ALG FractionalShift: public ConstModIPAlgorithm
{
public:
  FractionalShift(Real sx=0.0, Real sy=0.0, Real sz=0.0);
  FractionalShift(const Vec3& v);

  void SetShift(Real sx=0.0, Real sy=0.0, Real sz=0.0);
  void SetShift(const Vec3& v);
  Vec3 GetShift() const;
  virtual void Visit(ImageHandle& ih) const;
private:
  Vec3 shift_;
  
};

}}}//ns

#endif /*FRACTIONAL_SHIFT_HH_*/
