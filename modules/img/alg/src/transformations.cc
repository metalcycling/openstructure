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
  specific transformation implementations

  Author: Ansgar Philippsen
*/

#include <ost/message.hh>

#include "transformations.hh"

namespace ost { namespace img { namespace alg {

Rotate2D::Rotate2D(Real angle)
{
  Real p = angle;
  Real c = cos(p);
  Real s = sin(p);
  Matrix() = Mat4(c,-s, 0.0, 0.0, 
		  s, c, 0.0, 0.0, 
		  0.0, 0.0, 1.0, 0.0,
		  0.0, 0.0, 0.0, 1.0);

  InverseMatrix() = Mat4(c, s, 0.0, 0.0, 
			 -s, c, 0.0, 0.0, 
			 0.0, 0.0, 1.0, 0.0,
			 0.0, 0.0, 0.0, 1.0);
}
  
Rotate2D90::Rotate2D90(int step)
{
  step=step%4;
  
  if(step==1) {
    // 90deg CW
    Matrix()  = Mat4(0.0,-1.0,0.0,0.0, 1.0,0.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0);
    InverseMatrix() = Mat4(0.0,1.0,0.0,0.0, -1.0,0.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0);
  } else if(step==2) {
    // 180 deg
    Matrix()  = Mat4(-1.0,0.0,0.0,0.0, 0.0,-1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0);
    InverseMatrix() = Matrix();
  } else if(step==3) {
    // 90deg CCW
    Matrix()  = Mat4(0.0,1.0,0.0,0.0, -1.0,0.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0);
    InverseMatrix() = Mat4(0.0,-1.0,0.0,0.0, 1.0,0.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0);
  } else {
    // identity
  }
}


Scale::Scale(Real scal)
{
  if(scal==0.0) throw(Error("scaling factor is zero"));
  Real iscal = 1.0/scal;
  Matrix() = Mat4(scal,0.0,0.0,0.0,
		  0.0,scal,0.0,0.0,
		  0.0,0.0,scal,0.0,
		  0.0,0.0,0.0,1.0);

  InverseMatrix() = Mat4(iscal,0.0,0.0,0.0,
			 0.0,iscal,0.0,0.0,
			 0.0,0.0,iscal,0.0,
			 0.0,0.0,0.0,1.0);
}


Translate::Translate(const Vec3& t)
{
  Matrix() = Mat4(1.0,0.0,0.0,t[0],
		  0.0,1.0,0.0,t[1],
		  0.0,0.0,1.0,t[2],
		  0.0,0.0,0.0,1.0);

  InverseMatrix() = Mat4(1.0,0.0,0.0,-t[0],
			 0.0,1.0,0.0,-t[1],
			 0.0,0.0,1.0,-t[2],
			 0.0,0.0,0.0,1.0);
}

}}} // namespaces
