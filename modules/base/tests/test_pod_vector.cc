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
  Author: Marco Biasini
 */
#include <ost/pod_vector.hh>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using namespace ost;

struct Payload {
  Payload(int p): x(p) {}
  Payload(): x(0) {}
  int x;  
};

BOOST_AUTO_TEST_SUITE(base);

BOOST_AUTO_TEST_CASE(pod_vector_trivia)
{
  PodVector<Payload> vec;
  
  BOOST_CHECK_EQUAL(vec.begin(), vec.end());
  BOOST_CHECK_EQUAL(vec.size(), size_t(0));
  BOOST_CHECK_EQUAL(vec.capacity(), size_t(0));
}

BOOST_AUTO_TEST_CASE(pod_vector_push_pop)
{
  PodVector<Payload> vec;
  
  vec.push_back(Payload(0));
  vec.push_back(Payload(1));
  vec.push_back(Payload(2));
  BOOST_CHECK_EQUAL(vec.size(), size_t(3));
  BOOST_CHECK_EQUAL(vec.capacity(), size_t(4));
  BOOST_CHECK_EQUAL(vec[0].x, 0);
  BOOST_CHECK_EQUAL(vec[1].x, 1);
  BOOST_CHECK_EQUAL(vec[2].x, 2);
  vec.pop_back();
  BOOST_CHECK_EQUAL(vec.size(), size_t(2));
  BOOST_CHECK_EQUAL(vec.capacity(), size_t(4));  
}

BOOST_AUTO_TEST_SUITE_END();
