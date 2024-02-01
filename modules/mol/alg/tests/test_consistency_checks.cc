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
 * Author Valerio Mariani
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <ost/mol/mol.hh>
#include <ost/mol/alg/consistency_checks.hh>

using namespace ost::mol;
using namespace ost::mol::alg;


BOOST_AUTO_TEST_SUITE( mol_alg );

BOOST_AUTO_TEST_CASE(consistency_check) 
{

  EntityHandle a = CreateEntity();
  XCSEditor eda = a.EditXCS();
  ChainHandle ca = eda.InsertChain("A");
  ResidueHandle r1a = eda.AppendResidue(ca, "ALA", 1);
  ResidueHandle r2a = eda.AppendResidue(ca, "TYR", 2);
  ResidueHandle r3a = eda.AppendResidue(ca, "LEU", 3);
  ResidueHandle r4a = eda.AppendResidue(ca, "GLY", 4);
  ChainHandle ca2 = eda.InsertChain("C");
  ResidueHandle r1a2 = eda.AppendResidue(ca2, "PRO", 1);

  EntityHandle b = CreateEntity();
  XCSEditor edb = b.EditXCS();
  ChainHandle cb = edb.InsertChain("A");
  ResidueHandle r1b = edb.AppendResidue(cb, "ALA", 1);
  ResidueHandle r2b = edb.AppendResidue(cb, "TYR", 2);
  ResidueHandle r3b = edb.AppendResidue(cb, "LEU", 3);
  ResidueHandle r4b = edb.AppendResidue(cb, "GLY", 4);
  ResidueHandle r5b = edb.AppendResidue(cb, "PRO", 5);

  EntityHandle c = CreateEntity();
  XCSEditor edc = c.EditXCS();
  ChainHandle cc = edc.InsertChain("A");
  ResidueHandle r1c = edc.AppendResidue(cc, "ALA", 1);
  ResidueHandle r2c = edc.AppendResidue(cc, "PRO", 2);
  ResidueHandle r3c = edc.AppendResidue(cc, "LEU", 3);
  ResidueHandle r4c = edc.AppendResidue(cc, "GLY", 4);
   
  const EntityView eva = a.CreateFullView();
  const EntityView evb = b.CreateFullView();
  const EntityView evc = c.CreateFullView();
  BOOST_CHECK(ResidueNamesMatch(eva, evb) == true);
  BOOST_CHECK(ResidueNamesMatch(evb, eva) == true);
  BOOST_CHECK(ResidueNamesMatch(evc, evb) == false);
  BOOST_CHECK(ResidueNamesMatch(evb, evc) == false);
}

BOOST_AUTO_TEST_SUITE_END();
