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
 *  Authors: Marco Biasini, Juergen Haas
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/mol/mol.hh>
#include <ost/message.hh>
#include "dummy_ent.hh"


using namespace ost;
using namespace ost::mol;

BOOST_AUTO_TEST_SUITE(mol_base);

BOOST_AUTO_TEST_CASE(switch_pos)
{
  test::DummyEnt ent;
  geom::Vec3 y=geom::Vec3(0.0, 1.0, 0.0);
  geom::Vec3 z=geom::Vec3(0.0, 0.0, 1.0);  
  // Should throw, since aa has no alt-location set.
  XCSEditor editor=ent.e.EditXCS();  
  BOOST_CHECK_THROW(editor.AddAltAtomPos("B", ent.aa, y),
                    Error);
  AtomHandle xx=editor.InsertAltAtom(ent.r, "ZZ", "A", y); 
  BOOST_CHECK_NO_THROW(editor.AddAltAtomPos("B", xx, z));
  BOOST_CHECK_EQUAL(ent.r.GetCurrentAltGroupName(), "A");
  BOOST_CHECK(ent.r.SwitchAtomPos("B"));
  BOOST_CHECK_EQUAL(xx.GetPos(), z);
  BOOST_CHECK(ent.r.SwitchAtomPos("A"));
  BOOST_CHECK_EQUAL(xx.GetPos(), y);
  editor.SetAtomPos(xx, geom::Vec3(0.0, 0.0, 0.0));
  BOOST_CHECK(ent.r.SwitchAtomPos("B"));  
  BOOST_CHECK(ent.r.SwitchAtomPos("A"));    
  BOOST_CHECK_EQUAL(xx.GetPos(), geom::Vec3(0.0, 0.0, 0.0));
}
BOOST_AUTO_TEST_CASE(delete_atom) 
{
  test::DummyEnt ent;
  XCSEditor editor=ent.e.EditXCS();       
  AtomHandle xx=editor.InsertAltAtom(ent.r, "ZZ", "A", geom::Vec3());
  AtomHandle zz=editor.InsertAltAtom(ent.r, "XX", "A", geom::Vec3());
  
  editor.AddAltAtomPos("B", xx, geom::Vec3());
  editor.AddAltAtomPos("B", zz, geom::Vec3());
  BOOST_CHECK_EQUAL(ent.r.GetCurrentAltGroupName(), "A");  
  editor.DeleteAtom(xx);  
  BOOST_CHECK(ent.r.SwitchAtomPos("B"));
  editor.DeleteAtom(zz);
  BOOST_CHECK_EQUAL(ent.r.GetCurrentAltGroupName(), "");
  BOOST_CHECK(ent.r.HasAltAtoms()==false);
}
BOOST_AUTO_TEST_CASE(atom_group_triv) 
{
  test::DummyEnt ent;
  BOOST_CHECK(ent.r.HasAltAtoms()==false);
  BOOST_CHECK(ent.r.HasAltAtomGroup("")==false);
  std::vector<String> names=ent.r.GetAltAtomGroupNames();
  BOOST_CHECK(names.empty());
  XCSEditor editor=ent.e.EditXCS();
  AtomHandle xx=editor.InsertAltAtom(ent.r, "A", "ZZ", geom::Vec3(0.0,1.0,0.0));   
  editor.AddAltAtomPos("B", xx, geom::Vec3(1.0, 0.0, 0.0));
  BOOST_CHECK(ent.r.HasAltAtoms()==true);
  BOOST_CHECK(ent.r.HasAltAtomGroup("B")==true);
  names=ent.r.GetAltAtomGroupNames();  
  BOOST_CHECK_EQUAL(names.size(), size_t(2));
  BOOST_CHECK_EQUAL(names[0], "B");
}

BOOST_AUTO_TEST_SUITE_END();
