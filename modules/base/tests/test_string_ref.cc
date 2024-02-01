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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/string_ref.hh>
using namespace ost;

BOOST_AUTO_TEST_SUITE( base );


BOOST_AUTO_TEST_CASE(test_whitespace_split)
{
  std::string s=" 1 2  \r\n 3  4\t5 "; 
  StringRef sr(s.c_str(), s.size());
  std::vector<StringRef> srl=sr.split();
  BOOST_CHECK_EQUAL(srl.size(), static_cast<size_t>(5));
  BOOST_CHECK_EQUAL(srl[0], StringRef("1", 1));
  BOOST_CHECK_EQUAL(srl[1], StringRef("2", 1));
  BOOST_CHECK_EQUAL(srl[2], StringRef("3", 1));
  BOOST_CHECK_EQUAL(srl[3], StringRef("4", 1));
  BOOST_CHECK_EQUAL(srl[4], StringRef("5", 1));
}


BOOST_AUTO_TEST_CASE(test_string_ref)
{
  StringRef sr("12345", 5);
  BOOST_CHECK_EQUAL(sr.length(), size_t(5));
  BOOST_CHECK_EQUAL(sr.size(), size_t(5));  
  BOOST_CHECK(!sr.empty());
  BOOST_CHECK_EQUAL(sr[0], '1');  
  BOOST_CHECK_EQUAL(sr[1], '2');  
  BOOST_CHECK_EQUAL(sr[2], '3');  
  BOOST_CHECK_EQUAL(sr[3], '4');  
  BOOST_CHECK_EQUAL(sr[4], '5');
  
  // compare
  BOOST_CHECK_EQUAL(sr, sr);
  StringRef srx("1234", 4);
  BOOST_CHECK(sr!=srx);
  BOOST_CHECK(srx!=sr);
  
  // substr
  StringRef sr2=sr.substr(2);
  BOOST_CHECK_EQUAL(sr2.length(), size_t(3));  
  BOOST_CHECK_EQUAL(sr2[0], '3');  
  BOOST_CHECK_EQUAL(sr2[1], '4');  
  BOOST_CHECK_EQUAL(sr2[2], '5');
  
  StringRef sr3=sr.substr(2, 2);
  BOOST_CHECK_EQUAL(sr3.length(), size_t(2));  
  BOOST_CHECK_EQUAL(sr3[0], '3');  
  BOOST_CHECK_EQUAL(sr3[1], '4');  
  
  // rtrim
  BOOST_CHECK_EQUAL(StringRef("12", 2).rtrim(), StringRef("12", 2));
  BOOST_CHECK_EQUAL(StringRef("12 \t", 4).rtrim(), StringRef("12", 2));
  
  // to_int
  std::pair<bool, int> r=sr.to_int();
  BOOST_CHECK(r.first==true);
  BOOST_CHECK(r.second==12345);
  r=StringRef("x", 1).to_int();
  BOOST_CHECK(r.first==false);
  r=StringRef("12.34", 5).to_int();
  BOOST_CHECK(r.first==false);
  
  // to_float
  std::pair<bool, Real> r2=StringRef("1", 1).to_float();
  BOOST_CHECK(r2.first==true);
  BOOST_CHECK(r2.second==1.0);
  r2=StringRef("1.5", 3).to_float();
  BOOST_CHECK(r2.first==true);
  BOOST_CHECK(r2.second==1.5);  
  r2=StringRef("x", 1).to_float();
  BOOST_CHECK(r2.first==false);
  r2=StringRef("12.3.4", 6).to_float();
  BOOST_CHECK(r2.first==false);
  r2=StringRef("12.34.", 6).to_float();
  BOOST_CHECK(r2.first==false);

  // to_float_with_exp
  std::pair<bool, Real> r3=StringRef("1", 1).to_float();
  BOOST_CHECK(r3.first==true);
  BOOST_CHECK(r3.second==1.0);
  r3=StringRef("1.5", 3).to_float();
  BOOST_CHECK(r3.first==true);
  BOOST_CHECK(r3.second==1.5);
  r3=StringRef("x", 1).to_float();
  BOOST_CHECK(r3.first==false);
  r3=StringRef("12.3.4", 6).to_float();
  BOOST_CHECK(r3.first==false);
  r3=StringRef("12.34.", 6).to_float();
  BOOST_CHECK(r3.first==false);
  r3=StringRef("12.34e5", 7).to_float();
  BOOST_CHECK(r3.first==true);
  r3=StringRef("12e4", 4).to_float();
  BOOST_CHECK(r3.first==true);
  BOOST_CHECK_CLOSE(Real(120000), Real(r3.second), Real(1e-4));
  r3=StringRef("2e+4", 4).to_float();
  BOOST_CHECK(r3.first==true);
  BOOST_CHECK_CLOSE(Real(20000), Real(r3.second), Real(1e-4));
  r3=StringRef("2.3E+4", 6).to_float();
  BOOST_CHECK(r3.first==true);
  BOOST_CHECK_CLOSE(Real(23000), Real(r3.second), Real(1e-4));
  r3=StringRef("2.3E-4", 6).to_float();
  BOOST_CHECK(r3.first==true);
  BOOST_CHECK_CLOSE(Real(0.00023), Real(r3.second), Real(1e-4));
  r3=StringRef("2.010000e+00", 12).to_float();
  BOOST_CHECK(r3.first==true);
  BOOST_CHECK_CLOSE(Real(2.01), Real(r3.second), Real(1e-4));
  r3=StringRef("5e-34", 5).to_float();
  BOOST_CHECK(r3.first==true);
  r3=StringRef("5E-34", 5).to_float();
  BOOST_CHECK(r3.first==true);
  r3=StringRef("5.34e-34", 8).to_float();
  BOOST_CHECK(r3.first==true);
  r3=StringRef("5.34e-34e", 9).to_float();
  BOOST_CHECK(r3.first==false);
  r3=StringRef("5.34ee34", 8).to_float();
  BOOST_CHECK(r3.first==false);
  r3=StringRef("5.34e--34e", 10).to_float();
  BOOST_CHECK(r3.first==false);
  r3=StringRef("5.34e+3+4", 9).to_float();
  BOOST_CHECK(r3.first==false);
  r3=StringRef("5.34e-+34e", 10).to_float();
  BOOST_CHECK(r3.first==false);
  r3=StringRef("5.34e-3-4", 9).to_float();
  BOOST_CHECK(r3.first==false);
}

BOOST_AUTO_TEST_SUITE_END();
