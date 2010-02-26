//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
 *  Authors: Marco Biasini, Juergen Haas
 */
#include <ost/mol/mol.hh>
#include <ost/log.hh>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/message.hh>
#include <ost/geom/geom.hh>

using ost::Logger;



using namespace ost;
using namespace ost::mol;
struct Structure {
  Structure()
  {
    e=CreateEntity();
    ICSEditor editor=e.RequestICSEditor();
    c=editor.InsertChain("A");
    r=editor.AppendResidue(c, "ANGELIN");
    a1=editor.InsertAtom(r, "CC", geom::Vec3(-1.0, 0.0,  0.0));
    a2=editor.InsertAtom(r, "X1", geom::Vec3( 0.0, 0.0,  0.0));
    a3=editor.InsertAtom(r, "X2", geom::Vec3( 1.0, 0.0,  0.0));
    a4=editor.InsertAtom(r, "X3", geom::Vec3( 0.0, 0.0, -1.0));
    a5=editor.InsertAtom(r, "X4", geom::Vec3( 0.0, 0.0,  1.0));
    editor.Connect(a1, a2);
    editor.Connect(a2, a3);
    editor.Connect(a2, a4);
    editor.Connect(a2, a5);
  }
  EntityHandle  e;
  ChainHandle   c;
  ResidueHandle r;
  AtomHandle    a1;
  AtomHandle    a2;
  AtomHandle    a3;
  AtomHandle    a4;
  AtomHandle    a5;
};

void test_ics_angle_trivia() {
  Structure s;
  ICSEditor e=s.e.RequestICSEditor();
  BOOST_CHECK(!e.SetAngle(s.a1, s.a3, s.a2, 0));
  BOOST_CHECK(!e.SetAngle(s.a2, s.a3, s.a1, 0));
  BOOST_CHECK(e.SetAngle(s.a1, s.a2, s.a3, 0));

  BOOST_CHECK(!e.SetAngle(s.a1, s.a3, s.a2, 0));
  BOOST_CHECK(!e.SetAngle(s.a2, s.a3, s.a1, 0));
  BOOST_CHECK(e.SetAngle(s.a1, s.a2, s.a3, 0));

  BOOST_CHECK_THROW(s.e.GetAngle(s.a1, s.a3, s.a2), Error);
  BOOST_CHECK_THROW(s.e.GetAngle(s.a2, s.a3, s.a1), Error);
  BOOST_CHECK_NO_THROW(s.e.GetAngle(s.a1, s.a2, s.a3));
}

const static Real EPSILON=0.000001;

Real angle_xcs(AtomHandle a1, AtomHandle a2, AtomHandle a3) {
  return acos(Dot(geom::Normalize(a1.GetPos()-a2.GetPos()),
              geom::Normalize(a3.GetPos()-a2.GetPos())));
}

bool test_angle(Real a, Real e) {
  return std::abs(fmod(float(a-e), float(M_PI/2)))<EPSILON;
}

void test_ics_set_angle_sec() {
  Structure s;
  // test for set angle of two secondary connectors
  ICSEditor e=s.e.RequestICSEditor();  
  geom::Plane p(s.a3.GetPos(), s.a2.GetPos(), s.a4.GetPos());
  e.SetAngle(s.a3, s.a2, s.a4, M_PI/4);
  e.UpdateXCS();  
  BOOST_CHECK_MESSAGE(test_angle(M_PI/4, angle_xcs(s.a3, s.a2, s.a4)),
                      "expected " << M_PI/4 << " but " <<
                      angle_xcs(s.a3, s.a2, s.a4));

  BOOST_CHECK(geom::IsInPlane(p, s.a3.GetPos(), EPSILON));
  BOOST_CHECK(geom::IsInPlane(p, s.a2.GetPos(), EPSILON));
  BOOST_CHECK(geom::IsInPlane(p, s.a4.GetPos(), EPSILON));
}

void test_ics_set_angle_prim() {
  Structure s;
  // test for set angle of two secondary connectors
  ICSEditor e=s.e.RequestICSEditor();    
  e.SetAngle(s.a1, s.a2, s.a4, M_PI/4);
  e.UpdateXCS();
  geom::Plane p(s.a1.GetPos(), s.a2.GetPos(), s.a4.GetPos());
  BOOST_CHECK_MESSAGE(test_angle(M_PI/4, angle_xcs(s.a1, s.a2, s.a4)),
                      "expected " << M_PI/4 << " but " <<
                      angle_xcs(s.a1, s.a2, s.a4));

  BOOST_CHECK(geom::IsInPlane(p, s.a1.GetPos(), EPSILON));
  BOOST_CHECK(geom::IsInPlane(p, s.a2.GetPos(), EPSILON));
  BOOST_CHECK(geom::IsInPlane(p, s.a4.GetPos(), EPSILON));
}

void test_ics_set_angle() {
  //Logger::Instance().PushVerbosityLevel(Logger::DUMP);
  test_ics_set_angle_sec();
  test_ics_set_angle_prim();
}

void test_ics_get_angle() {
  Structure s;
  // test for get angle with secondary connectors only
  Real a=s.e.GetAngle(s.a3, s.a2, s.a4);
  Real e=angle_xcs(s.a3, s.a2, s.a4);
  BOOST_CHECK_MESSAGE(test_angle(a, e),
                      "expected " << e << " but " << a
                      << " found");

  a=s.e.GetAngle(s.a5, s.a2, s.a4);
  e=angle_xcs(s.a5, s.a2, s.a4);
  BOOST_CHECK_MESSAGE(test_angle(a, e),
                      "expected " << e << " but " << a
                      << " found");

  a=s.e.GetAngle(s.a5, s.a2, s.a5);
  e=angle_xcs(s.a5, s.a2, s.a5);
  BOOST_CHECK_MESSAGE(test_angle(a, e),
                      "expected " << e << " but " << a
                      << " found");
}

void test_ics() {
  // test trivial things first
  test_ics_angle_trivia();
  test_ics_get_angle();
  test_ics_set_angle();
}

BOOST_AUTO_TEST_SUITE( mol_base )

BOOST_AUTO_TEST_CASE(test_ics_angle_trivia) 
{
  test_ics_angle_trivia();
}

BOOST_AUTO_TEST_CASE(test_angle) 
{
  test_angle();
}

BOOST_AUTO_TEST_CASE(test_ics_set_angle_sec) 
{
  test_ics_set_angle_sec();
}

BOOST_AUTO_TEST_CASE(test_ics_set_angle_prim) 
{
  test_ics_set_angle_prim();
}

BOOST_AUTO_TEST_CASE(test_ics_set_angle) 
{
  test_ics_set_angle();
}

BOOST_AUTO_TEST_CASE(test_ics_get_angle) 
{
  test_ics_get_angle();
}

BOOST_AUTO_TEST_CASE(test_ics) 
{
  test_ics();
}

BOOST_AUTO_TEST_SUITE_END()