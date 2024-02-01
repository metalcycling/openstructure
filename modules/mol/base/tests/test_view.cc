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

struct FixtureView {
  FixtureView() {
    e=CreateEntity();
    XCSEditor editor=e.EditXCS();
    c=editor.InsertChain("A");
    r=editor.AppendResidue(c, "ANGELIN");
    aa=editor.InsertAtom(r, "X1", geom::Vec3( 0.0, 0.0,  0.0));
    ab=editor.InsertAtom(r, "X2", geom::Vec3( 1.0, 0.0,  0.0));
    ac=editor.InsertAtom(r, "X3", geom::Vec3( 0.0, 0.0, -1.0));
    ad=editor.InsertAtom(r, "X4", geom::Vec3( 0.0, 0.0,  1.0));
    editor.Connect(aa, ab);
    editor.Connect(ab, ac);
    editor.Connect(ac, ad);
  }
  EntityHandle e;
  ChainHandle  c;
  ResidueHandle r;
  AtomHandle aa;
  AtomHandle ab;
  AtomHandle ac;
  AtomHandle ad;
};

BOOST_AUTO_TEST_SUITE( mol_base );

BOOST_AUTO_TEST_CASE(gen_full_view) 
{
  FixtureView f;

  EntityView ev = f.e.CreateFullView();
  BOOST_CHECK_EQUAL(ev.GetChainCount(),f.e.GetChainCount());
  BOOST_CHECK_EQUAL(ev.GetResidueCount(),f.e.GetResidueCount());
  BOOST_CHECK_EQUAL(ev.GetAtomCount(),f.e.GetAtomCount());

  ChainView cv = (ev.GetChainList())[0];
  BOOST_CHECK_EQUAL(cv.GetHandle(),f.c);

  ResidueView rv = (cv.GetResidueList())[0];
  BOOST_CHECK_EQUAL(rv.GetHandle(),f.r);

  AtomViewList avl  = rv.GetAtomList();
  BOOST_CHECK_EQUAL(avl[0].GetHandle(),f.aa);
  BOOST_CHECK_EQUAL(avl[1].GetHandle(),f.ab);
  BOOST_CHECK_EQUAL(avl[2].GetHandle(),f.ac);
  BOOST_CHECK_EQUAL(avl[3].GetHandle(),f.ad);

  BOOST_CHECK_EQUAL(avl[0].GetBondCount(),1);
  BOOST_CHECK_EQUAL(avl[1].GetBondCount(),2);
  BOOST_CHECK_EQUAL(avl[2].GetBondCount(),2);
  BOOST_CHECK_EQUAL(avl[3].GetBondCount(),1);
  
  mol::AtomView av1=avl[0];
  mol::AtomView av2=avl[1];
  
  BOOST_CHECK(av1!=av2);
  BOOST_CHECK(av1.GetHashCode()!=av2.GetHashCode());
  BOOST_CHECK(av1.GetHandle().GetHashCode()!=av1.GetHashCode());
  
  mol::AtomView av3=av1;
  BOOST_CHECK(av1==av3);
  BOOST_CHECK_EQUAL(av1.GetHashCode(), av3.GetHashCode());  
}

BOOST_AUTO_TEST_SUITE_END();
