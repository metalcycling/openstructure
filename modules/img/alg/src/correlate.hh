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

#ifndef IMG_ALG_CORRELATE_HH
#define IMG_ALG_CORRELATE_HH

#include <ost/img/image_state.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

class DLLEXPORT_IMG_ALG ExplicitCorrelateFnc {
public:
  ExplicitCorrelateFnc();
  ExplicitCorrelateFnc(const ConstImageHandle& ref, bool wrap=true);
  
  template <typename V, class D>
  ImageStateBasePtr VisitState(const ImageStateImpl<V,D>& isi) const;

  String GetAlgorithmName() {return "Correlate";}

private:
  ConstImageHandle ref_;
  bool wrap_;
};

typedef ImageStateConstModOPAlgorithm<ExplicitCorrelateFnc> ExplicitCorrelate;

}

OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateConstModOPAlgorithm<alg::ExplicitCorrelateFnc>)

  }} // ns

#endif
