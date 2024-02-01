//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
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

#include <ost/geom/geom.hh>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "helper.hh"
using namespace geom;


BOOST_AUTO_TEST_SUITE( geom );

BOOST_AUTO_TEST_CASE(init_quat)
{
  // default
  Quat q1;
  BOOST_CHECK_CLOSE(Real(q1.w), Real(1.0), Real(1.0e-5));
  BOOST_CHECK_CLOSE(Real(q1.x), Real(0.0), Real(1.0e-5));
  BOOST_CHECK_CLOSE(Real(q1.y), Real(0.0), Real(1.0e-5));
  BOOST_CHECK_CLOSE(Real(q1.z), Real(0.0), Real(1.0e-5));

  Quat q2(2.0,3.0,4.0,5.0);
  BOOST_CHECK_CLOSE(Real(q2.w), Real(2.0), Real(1.0e-5));
  BOOST_CHECK_CLOSE(Real(q2.x), Real(3.0), Real(1.0e-5));
  BOOST_CHECK_CLOSE(Real(q2.y), Real(4.0), Real(1.0e-5));
  BOOST_CHECK_CLOSE(Real(q2.z), Real(5.0), Real(1.0e-5));
}

BOOST_AUTO_TEST_CASE(quat_rotate)
{
  Vec3 v(1,0,0);
  Quat q(30.0*M_PI/180.0,Vec3(0,1,0));
  Vec3 vrot=q.Rotate(v);
  BOOST_CHECK_CLOSE(Real(cos(30.0*M_PI/180.0)), Real(vrot[0]), Real(1e-5));
  BOOST_CHECK_SMALL(vrot[1],Real(1e-5));
  BOOST_CHECK_CLOSE(Real(-sin(30.0*M_PI/180.0)), Real(vrot[2]), Real(1e-5));
}


BOOST_AUTO_TEST_SUITE_END();

