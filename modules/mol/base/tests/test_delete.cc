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
#include <ost/log.hh>
#include <ost/message.hh>
#include <ost/geom/geom.hh>


using namespace ost;
using namespace ost::mol;


class E : public EntityVisitor {
public:
  virtual bool VisitTorsion(const TorsionHandle& t) {
    throw Error("should never get here");
    return false;
  }
};

struct Fixture {
  Fixture() {
    e=CreateEntity();
    XCSEditor editor=e.EditXCS();    
    c=editor.InsertChain("A");
    r=editor.AppendResidue(c, "XXX");
    aa=editor.InsertAtom(r, "A", geom::Vec3());
    ab=editor.InsertAtom(r, "B", geom::Vec3());
    ac=editor.InsertAtom(r, "C", geom::Vec3());
    ad=editor.InsertAtom(r, "D", geom::Vec3());

    BondHandle b1=editor.Connect(aa, ab);
    BondHandle b2=editor.Connect(ab, ac);
    BondHandle b3=editor.Connect(ac, ad);

    t=editor.AddTorsion("TORSION", aa, ab, ac, ad);
  }
  EntityHandle e;
  ChainHandle  c;
  ResidueHandle r;
  AtomHandle aa;
  AtomHandle ab;
  AtomHandle ac;
  AtomHandle ad;
  TorsionHandle t;
};

BOOST_AUTO_TEST_SUITE( mol_base );

BOOST_AUTO_TEST_CASE(atom_delete) 
{
  Fixture f;
  XCSEditor e=f.e.EditXCS();
  e.DeleteAtom(f.ab);
  E ev;
  BOOST_CHECK_NO_THROW(f.e.Apply(ev));
}

BOOST_AUTO_TEST_SUITE_END();
