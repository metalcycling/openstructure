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
  custom transformer

  Author: Ansgar Philippsen
*/

#ifndef IMG_CUSTOM_TRANSFORMER_IMPL_HH
#define IMG_CUSTOM_TRANSFORMER_IMPL_HH

#include <ost/message.hh>

#include "transformer_base.hh"

namespace ost { namespace img { namespace alg {

class DLLEXPORT CustomTransformerException: public Error {
public:
  virtual ~CustomTransformerException() throw() {}
  CustomTransformerException(const String& s):
    Error(s) {}
};


//! custom transformation
/*!
  Allows an arbitrary transformation to be specified. If only
  the forward matrix is given, the inverse will be calculated
  automatically.
*/
class DLLEXPORT_IMG_ALG CustomTransformer: public Transformer {
public:
  CustomTransformer(const Mat4& tmat);
  CustomTransformer(const Mat4& tmat,const Mat4& imat);

  virtual Mat4 GetMatrix() const {return tmat_;}
  virtual Mat4 GetInverseMatrix() const {return imat_;}

private:
  Mat4 tmat_;
  Mat4 imat_;
};

}}}

#endif
