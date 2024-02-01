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
#ifndef IMG_PHYSICS_UNITS_HH
#define IMG_PHYSICS_UNITS_HH

#include <boost/mpl/vector_c.hpp>
#include<boost/mpl/placeholders.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/transform.hpp>
#include <ost/base.hh>

namespace ost { namespace img { namespace physics{

namespace detail{

namespace mpl = boost::mpl;

typedef mpl::vector_c<int,1,0,0,0,0,0,0> mass_dimension;
typedef mpl::vector_c<int,0,1,0,0,0,0,0> length_dimension; // or position 
typedef mpl::vector_c<int,0,0,1,0,0,0,0> time_dimension;
typedef mpl::vector_c<int,0,0,0,1,0,0,0> charge_dimension;
typedef mpl::vector_c<int,0,0,0,0,1,0,0> temperature_dimension;
typedef mpl::vector_c<int,0,0,0,0,0,1,0> intensity_dimension;
typedef mpl::vector_c<int,0,0,0,0,0,0,1> amount_dimension;
typedef mpl::vector_c<int,0,1,-1,0,0,0,0> velocity_dimension;     // l/t
typedef mpl::vector_c<int,0,1,-2,0,0,0,0> acceleration_dimension; // l/(t2)
typedef mpl::vector_c<int,1,1,-1,0,0,0,0> momentum_dimension;     // ml/t
typedef mpl::vector_c<int,1,1,-2,0,0,0,0> force_dimension;        // ml/(t2)
typedef mpl::vector_c<int,0,0,0,0,0,0,0> scalar_dimension;


template <class Dimensions>
class quantity
{
public:
  quantity();
  template <class D2>
  quantity(const quantity<D2>& rhs);
  explicit quantity(Real x);
  Real value() const;

protected:
  Real value_;
};

template <class D>
quantity<D> operator+(quantity<D> x, quantity<D> y);
template <class D>
quantity<D> operator-(quantity<D> x, quantity<D> y);


template <class D1, class D2>
struct multiply_dimensions: mpl::transform<D1,D2,mpl::plus<mpl::placeholders::_1,mpl::placeholders::_2> >
{
};
template <class D1, class D2>
quantity<typename multiply_dimensions<D1,D2>::type> operator*(quantity<D1> x, quantity<D2> y);

template <class D1, class D2>
struct divide_dimensions: mpl::transform<D1,D2,mpl::minus<mpl::placeholders::_1,mpl::placeholders::_2> >
{
};
template <class D1, class D2>
quantity<typename divide_dimensions<D1,D2>::type> operator/(quantity<D1> x, quantity<D2> y);


}// detail ns


typedef detail::quantity<detail::mass_dimension> Mass;
typedef detail::quantity<detail::length_dimension> Length;
typedef detail::quantity<detail::time_dimension> Time;
typedef detail::quantity<detail::charge_dimension> Charge;
typedef detail::quantity<detail::temperature_dimension> Temperature;
typedef detail::quantity<detail::intensity_dimension> Intensity;
typedef detail::quantity<detail::amount_dimension> Amount;
typedef detail::quantity<detail::velocity_dimension> Velocity;  
typedef detail::quantity<detail::momentum_dimension> Momentum;    
typedef detail::quantity<detail::force_dimension> Force;       


// provide implicit conversion from/to Real for Scalar
class Scalar : public detail::quantity<detail::scalar_dimension>
{
public:
  Scalar(Real x): detail::quantity<detail::scalar_dimension>(x)
  {
  }
  
  operator Real()
  {
    return value_;
  }
  
}; 
}}}//ns

#endif
