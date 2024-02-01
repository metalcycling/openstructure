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
 * Author Juergen Haas
 */
#include <ost/mol/alg/svd_superpose.hh>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/mol/mol.hh>
#include <ost/mol/builder.hh>
#include <ost/message.hh>

using namespace ost;
using namespace ost::mol;
using namespace ost::mol::alg;

struct Fixture {
  Fixture() {

    e = Builder()
           .Chain("A")
              .Residue("DUM")
                .Atom("A",geom::Vec3(-5,-5,-5))
                .Atom("B",geom::Vec3(-5, 5,-5))
                .Atom("C",geom::Vec3(-5, 5, 5))
                .Atom("D",geom::Vec3(-5,-5, 5))
                .Atom("E",geom::Vec3(5,-5, 5))
                .Atom("F",geom::Vec3(5,-5,-5))
                .Atom("G",geom::Vec3(5, 5,-5))
                .Atom("H",geom::Vec3(5, 5, 5));

  }
  EntityHandle e;

};

BOOST_AUTO_TEST_SUITE( mol_alg );

BOOST_AUTO_TEST_CASE(superposition_svd) 
{
  Fixture f1, f2;
  EntityView ev1 = f1.e.CreateFullView();
  XCSEditor ed=f1.e.EditXCS();
  EntityView ev2 = f2.e.CreateFullView();
  AtomHandleList atoms = f1.e.GetAtomList();
  for (AtomHandleList::const_iterator i = atoms.begin(), e = atoms.end(); i !=e ; ++i) {
      AtomHandle atom=*i;
      ed.SetAtomPos(atom, geom::EulerTransformation(0.1, 0, 0)*atom.GetPos()
                  +geom::Vec3(0, 5, 0));
  }
  SuperpositionResult res;
  res=SuperposeSVD(ev1, ev2, true);
  BOOST_CHECK(res.rmsd<0.0001);
}

BOOST_AUTO_TEST_CASE(iterative_superposition_svd) 
{
  Fixture f1, f2;
  EntityView ev1 = f1.e.CreateFullView();
  EntityView ev2 = f2.e.CreateFullView();
  XCSEditor ed=f1.e.EditXCS();
  AtomHandleList atoms = f1.e.GetAtomList();
  ed.SetAtomPos(atoms[0],geom::Vec3(-10,-10,-10));
  for (AtomHandleList::const_iterator i = atoms.begin(), e = atoms.end(); i !=e ; ++i) {
      AtomHandle atom=*i;
      ed.SetAtomPos(atom, geom::EulerTransformation(0.1, 0, 0)*atom.GetPos()
                  +geom::Vec3(0, 5, 0));
  }
  SuperpositionResult res_classic, res_1, res_2, res_3;
  res_1=IterativeSuperposeSVD(ev1, ev2, 1, 5.0, false);
  res_2=IterativeSuperposeSVD(ev1, ev2, 2, 5.0, false);
  res_3=IterativeSuperposeSVD(ev1, ev2, 10, 5.0, false);
  res_classic=SuperposeSVD(ev1,ev2, false); 
  BOOST_CHECK(res_2.rmsd_superposed_atoms<0.0001);
  BOOST_CHECK(res_2.rmsd>1.0);
  BOOST_CHECK(res_1.rmsd==res_classic.rmsd);
  BOOST_CHECK(res_1.ncycles==1);
  BOOST_CHECK(res_2.ncycles==2);
  BOOST_CHECK(res_1.fraction_superposed==1.0);
  BOOST_CHECK(res_2.fraction_superposed==7.0/8.0);
  BOOST_CHECK(res_3.ncycles==3);

  EntityHandle eh_2_atoms = Builder()
                              .Chain("A")
                                .Residue("D")
                                  .Atom("A",geom::Vec3(0,0,0))
                                  .Atom("B",geom::Vec3(0,0,0));

  //Error, that has to be thrown during initialization of superposition object
  BOOST_CHECK_THROW(IterativeSuperposeSVD(eh_2_atoms.CreateFullView(),
                                          eh_2_atoms.CreateFullView(),
                                          1,5.0,false), Error);

  //forces an error, since there won't be enough atoms within the distance 
  //threshold after the first superposition.
  ed.SetAtomPos(atoms[0],geom::Vec3(50,50,50));
  ed.SetAtomPos(atoms[1],geom::Vec3(50,50,50));
  ed.SetAtomPos(atoms[2],geom::Vec3(50,50,50));
  ed.SetAtomPos(atoms[3],geom::Vec3(50,50,50));

  BOOST_CHECK_THROW(IterativeSuperposeSVD(ev1, ev2, 5, 5.0, false), Error);
}




BOOST_AUTO_TEST_SUITE_END();
