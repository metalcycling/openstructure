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
/*
  Authors: Tobias Schmidt
 */
 #define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/geom/geom.hh>

using namespace geom;

BOOST_AUTO_TEST_SUITE( geom_base );

BOOST_AUTO_TEST_CASE(test_transform_essentials)
{
  Transform tf;

  BOOST_CHECK_EQUAL(tf.GetMatrix(), geom::Mat4());
  BOOST_CHECK_EQUAL(tf.GetRot(), geom::Mat3());
  BOOST_CHECK_EQUAL(tf.GetTrans(), geom::Vec3());
  BOOST_CHECK_EQUAL(tf.GetCenter(), geom::Vec3());

  Mat3 rot=AxisRotation(Vec3(0.4,1.0,-0.8),2.13253);
  Vec3 cen(0.5,1.0,2.2);
  Vec3 tra(10,20,30);
  tf.SetRot(rot);
  tf.SetTrans(tra);
  tf.SetCenter(cen);

  BOOST_CHECK_EQUAL(tf.GetRot(), rot);
  BOOST_CHECK_EQUAL(tf.GetCenter(), cen);
  BOOST_CHECK_EQUAL(tf.GetTrans(), tra);

  Transform tf2;
  tf2.SetMatrix(tf.GetMatrix());
  BOOST_CHECK_EQUAL(tf2.GetRot(), rot);
  BOOST_CHECK_EQUAL(tf2.GetCenter(), Vec3(0,0,0));
  BOOST_CHECK_EQUAL(tf2.GetTrans(), rot*-cen+tra);
}

BOOST_AUTO_TEST_CASE(test_transform_apply_transform)
{
  Transform t1,t2;
  Mat3 rm1=AxisRotation(Vec3(0.4,1.0,-0.8),2.13253);
  Mat3 rm2=AxisRotation(Vec3(2.4,-0.1,3.8),-1.123);
  t1.SetRot(rm1);
  t1.SetCenter(Vec3(0.5,1.0,2.2));
  t1.SetTrans(Vec3(10,20,30));
  t2.SetRot(rm2);
  t2.SetCenter(Vec3(1.3,2.7,-1.1));
  t2.SetTrans(Vec3(-40,-60,80));

  Mat4 mat1=t1.GetMatrix();
  Mat4 mat2=t2.GetMatrix();
  Mat4 mat12=mat2*mat1;
  
  Transform t3 = t1.Apply(t2);
  Mat4 mat3=t3.GetMatrix();
  for(int i=0;i<16;++i) BOOST_CHECK_CLOSE(mat3.Data()[i],mat12.Data()[i],1e-6);
            
}

BOOST_AUTO_TEST_SUITE_END();
