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
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_GUI_VIEWER_NORMALIZER_HH
#define IMG_GUI_VIEWER_NORMALIZER_HH


#include <boost/shared_ptr.hpp>
#include <ost/img/normalizer_impl.hh>

#include <ost/gui/module_config.hh>


namespace ost { namespace img { namespace gui {

class DLLEXPORT_OST_GUI ViewerNormalizer: public NormalizerImpl, 
                                          public RangeHandler
{
public:
  ViewerNormalizer();
  ViewerNormalizer(Real mininput, Real maxinput, Real minoutput, 
                   Real maxoutput, Real gamma, bool invert);
  ViewerNormalizer(const ViewerNormalizer& n);
  virtual ~ViewerNormalizer();

  
  // normalizer interface
  virtual Real BackConvert(Real v) const;
  virtual Complex BackConvert(Complex v) const;
  virtual Real Convert(Real v) const;
  virtual Complex Convert(Complex v) const;
  virtual NormalizerImpl* Clone() const {return new ViewerNormalizer((*this));}

  Real GetGamma();
  bool GetInvert();
  
private:
  Real p1_;
  Real p2_;
  Real p3_;
  
  Real gamma_;
  Real invert_sign_;
  Real clipmax_;
};

typedef boost::shared_ptr<ViewerNormalizer> ViewerNormalizerPtr;


}}}  //ns

#endif
