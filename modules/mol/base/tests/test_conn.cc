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

using namespace ost;
using namespace ost::mol;

BOOST_AUTO_TEST_SUITE( mol_base );

BOOST_AUTO_TEST_CASE(connect) 
{
  EntityHandle e=CreateEntity();
  XCSEditor editor=e.EditXCS();  
  ChainHandle c=editor.InsertChain("A");
  ResidueHandle r=editor.AppendResidue(c, "GLY");

  AtomHandle n=editor.InsertAtom(r, "N", geom::Vec3(0,1,0));
  AtomHandle o=editor.InsertAtom(r, "O", geom::Vec3(0,0,0));
  AtomHandle x=editor.InsertAtom(r, "X", geom::Vec3(1,0,0));

  editor.Connect(n, o);
  editor.Connect(o, x);
  // test symmetry of operation.
  BOOST_CHECK(BondExists(n, o));
  BOOST_CHECK(BondExists(o, n));

  BOOST_CHECK(BondExists(o, x));
  BOOST_CHECK(BondExists(x, o));

  BOOST_CHECK(!BondExists(n, x));
  BOOST_CHECK(!BondExists(x, n));

  BOOST_CHECK(n.GetBondCount()==1);
  BOOST_CHECK(o.GetBondCount()==2);
  BOOST_CHECK(n.GetBondCount()==1);
}

BOOST_AUTO_TEST_SUITE_END();
