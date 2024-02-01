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
  Author: Ansgar Philippsen
*/

#ifndef IMG_TEST_UTILS
#define IMG_TEST_UTILS

#include <iostream>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

//#include <img/base.hh>
#include <ost/img/image_state.hh>
#include <ost/img/value_util.hh>

using namespace ost::img;

template<typename V,class DP>
void Randomize(image_state::ImageStateImpl<V,DP>& isi)
{
  for(image_state::IndexIterator it(isi.GetSize());!it.AtEnd();++it) {
    isi.Value(it) = Random<V>();
  }
}

template <typename V>
bool CheckClose(const V& v1, const V& v2, Real tol)
{
  return std::abs(v1-v2)<tol;
}

template <typename V, class DP>
bool CheckEqualStates(const image_state::ImageStateImpl<V,DP>& one,
		      const image_state::ImageStateImpl<V,DP>& two,
		      Real tol=1e-4)
{
  std::ostringstream msg;
  msg << one.GetSize() << " != " << two.GetSize();
  if(one.GetSize()!=two.GetSize()) {
    BOOST_TEST_MESSAGE(msg.str());
    return false;
  }

  for(image_state::IndexIterator it(one.GetSize()); !it.AtEnd(); ++it) {
    if(!CheckClose(one.Value(it),two.Value(it),tol)) {
      msg.str("");
      msg << " @" << it << ": " << one.Value(it) << " != " << two.Value(it);
      msg << " (delta=" << std::abs(one.Value(it)-two.Value(it)) << ")";
      BOOST_ERROR(msg.str());
      return false;
    }
  }

  return true;
}

#endif
