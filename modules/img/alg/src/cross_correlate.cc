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
  author: Andreas Schenk
*/

#include "cross_correlate.hh"
#include "dft.hh"
#include "stat.hh"
#include "conjugate.hh"


namespace ost { namespace img { namespace alg {
CrossCorrelate::CrossCorrelate(const ConstImageHandle& ref):
  ConstModIPAlgorithm("CrossCorrelate"),
  ref_(ref)
{}

void CrossCorrelate::Visit(ImageHandle& i) const
{
  DFT dftalg;
  Conj conjugate_alg;
  Stat statalg;
  i.Apply(statalg);
  i-=statalg.GetMean();
  i/=statalg.GetStandardDeviation();
  i.ApplyIP(dftalg);
  ref_.Apply(statalg);
  ImageHandle tmp=ref_-statalg.GetMean();
  tmp/=statalg.GetStandardDeviation();
  tmp.ApplyIP(dftalg);
  tmp.ApplyIP(conjugate_alg);
  i*=tmp;
  i.ApplyIP(dftalg);
  i/=i.GetSize().GetVolume();
}

}}} // ns
