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
  authors: Giani Signorell, Ansgar Philippsen
*/

#ifndef IMG_ALG_FILTER_ANISOTROPIC_H
#define IMG_ALG_FILTER_ANISOTROPIC_H

#include "gaussian.hh"

namespace ost { namespace img { namespace alg { 

class DLLEXPORT_IMG_ALG AnisotropicFilter: public GaussianFilter {
public:
  AnisotropicFilter(Real sigma=1.0);

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& s);

protected:
  template <typename T, class D>
  void lineTransformFW(ImageStateImpl<T, D>& s,
		       ExtentIterator& it, 
		       std::vector<T>& tmp_avg, 
		       std::vector<T>& tmp_min, 
		       std::vector<T>& tmp_max);

  template <typename T, class D>
  void lineTransformBW(ImageStateImpl<T,D>& s,
		       ExtentIterator& it,
		       std::vector<T>& tmp_avg,
		       std::vector<T>& tmp_min,
		       std::vector<T>& tmp_max );
  
private:
  int n_;
  Real b1_b0_, b2_b0_, b3_b0_;

  template <typename T>
  bool is_larger(const T& v1, const T& v2);
};

}}} // namespaces

#endif
