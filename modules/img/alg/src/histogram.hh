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
#ifndef OST_IMG_ALG_HISTOGRAM_HH
#define OST_IMG_ALG_HISTOGRAM_HH

/*
  Author: Ansgar Philippsen
*/

#include <ost/base.hh>
#include <ost/message.hh>

#include <ost/img/algorithm.hh>
#include <ost/img/image_state.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

class DLLEXPORT HistogramError: public Error {
public:
    HistogramError(const String& s):
    Error(String("A histogram error occurred: ") + s)
  {}
};


class DLLEXPORT_IMG_ALG HistogramBase
{
public:
  typedef std::vector<int> Bins;

  // default ctor used for explicit instanciation. Don't use!
  HistogramBase();
  //! Initialize with number of bins to use
  HistogramBase(int bin_count, Real minimum, Real maximum);

  // image state algorithm interface
  template <typename T, class D>
  void VisitState(const ImageStateImpl<T,D>& isi);

  const Bins& GetBins() const;

  static String GetAlgorithmName() {return "Histogram";}

protected:

private:
  int bin_count_;
  Real min_,max_,cfac_;
  Bins bins_;
};

typedef ImageStateNonModAlgorithm<HistogramBase> Histogram;

}

OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateNonModAlgorithm<ost::img::alg::HistogramBase>)

}} // namespaces

#endif
