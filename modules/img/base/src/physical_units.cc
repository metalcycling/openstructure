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
#include <boost/static_assert.hpp>
#include <boost/preprocessor/seq.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/mpl/int.hpp>
#include "physical_units.hh"



namespace ost { namespace img { namespace physics{

namespace detail{

template <class Dimensions>
quantity<Dimensions>::quantity():value_(0.0){}

template <class Dimensions>
template<class D2>
quantity<Dimensions>::quantity(const quantity<D2> & rhs): value_(rhs.value())
{
  BOOST_STATIC_ASSERT((mpl::equal<Dimensions,D2>::type::value));
}

template <class Dimensions>
quantity<Dimensions>::quantity(Real x): value_(x)
{
}

template <class Dimensions>
Real quantity<Dimensions>::value() const 
{
  return value_;
}

template <class D>
quantity<D> operator+(quantity<D> x, quantity<D> y)
{
  return quantity<D>(x.value() + y.value());
}

template <class D>
quantity<D> operator-(quantity<D> x, quantity<D> y)
{
  return quantity<D>(x.value() - y.value());
}

template <class D1, class D2>
quantity<typename multiply_dimensions<D1,D2>::type> operator*(quantity<D1> x, quantity<D2> y)
{
   return quantity<typename multiply_dimensions<D1,D2>::type>(x.value() * y.value());
}

template <class D1, class D2>
quantity<typename divide_dimensions<D1,D2>::type> operator/(quantity<D1> x, quantity<D2> y)
{
   return quantity<typename divide_dimensions<D1,D2>::type>(x.value() / y.value());
}



}//detail ns
}}}//ns
