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
#include <ost/mol/builder.hh>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using namespace ost;
using namespace ost::mol;

BOOST_AUTO_TEST_SUITE( mol_base_builder );

BOOST_AUTO_TEST_CASE(inserts_chains)
{
  EntityHandle ent = Builder()
                       .Chain("A")
                       .Chain("B");

  BOOST_CHECK_EQUAL(2, ent.GetChainCount());
  BOOST_CHECK(ent.FindChain("A").IsValid());
  BOOST_CHECK(ent.FindChain("B").IsValid());
}


BOOST_AUTO_TEST_CASE(appends_residues) {
  EntityHandle ent = Builder()
                        .Chain("A")
                           .Residue("GLY")
                           .Residue("GLN");

  BOOST_CHECK_EQUAL(1, ent.GetChainCount());
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 2);
  BOOST_CHECK(ent.FindResidue("A", 1).IsValid());
  BOOST_CHECK(ent.FindResidue("A", 2).IsValid());
}

BOOST_AUTO_TEST_CASE(inserts_atoms) {
  EntityHandle ent = Builder()
                        .Chain("A")
                           .Residue("GLY")
                               .Atom("N")
                               .Atom("CA")
                               .Atom("C")
                               .Atom("O")
                      ;

  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 4);
  BOOST_CHECK(ent.FindAtom("A", 1, "N").IsValid());
  BOOST_CHECK(ent.FindAtom("A", 1, "CA").IsValid());
  BOOST_CHECK(ent.FindAtom("A", 1, "C").IsValid());
  BOOST_CHECK(ent.FindAtom("A", 1, "O").IsValid());
}

BOOST_AUTO_TEST_SUITE_END();
