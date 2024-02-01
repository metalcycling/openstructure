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
 *  Authors: Marco Biasini
 */

#include <ost/mol/builder.hh>
#include <ost/mol/transfer_connectivity.hh>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using namespace ost;
using namespace ost::mol;

BOOST_AUTO_TEST_SUITE( mol_transfer_conn );

BOOST_AUTO_TEST_CASE(test_transfer_conn)
{
  EntityHandle ent = Builder()
                       .Chain("A")
                         .Gly()
                         .Gly()
                         .ConnectToPrev()
                         .Gly()
                         .ConnectToPrev()
                   ;

  ent.EditXCS().Connect(ent.FindAtom("A", 1, "O"), ent.FindAtom("A", 3, "O"));

  EntityHandle ent2 = Builder()
                       .Chain("A")
                         .Gly(false)
                         .Gly(false)
                         .Gly(false)
                   ;
  std::map<ResidueHandle, ResidueHandle> to_from;
  ResidueHandleList r1 = ent.GetResidueList();
  ResidueHandleList r2 = ent2.GetResidueList();
  for (size_t i = 0;i<r1.size(); ++i) {
    to_from[r2[i]] = r1[i];
  }
  TransferConnectivity(ent2, to_from);
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 1, "N"), ent2.FindAtom("A", 1, "CA")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 1, "CA"), ent2.FindAtom("A", 1, "C")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 1, "C"), ent2.FindAtom("A", 1, "O")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 1, "C"), ent2.FindAtom("A", 2, "N")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 2, "N"), ent2.FindAtom("A", 2, "CA")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 2, "CA"), ent2.FindAtom("A", 2, "C")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 2, "C"), ent2.FindAtom("A", 2, "O")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 2, "C"), ent2.FindAtom("A", 3, "N")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 3, "N"), ent2.FindAtom("A", 3, "CA")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 3, "CA"), ent2.FindAtom("A", 3, "C")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 3, "C"), ent2.FindAtom("A", 3, "O")));
  BOOST_CHECK(BondExists(ent2.FindAtom("A", 1, "O"), ent2.FindAtom("A", 3, "O")));
  BOOST_CHECK_EQUAL(12, ent2.GetBondCount());
}

BOOST_AUTO_TEST_CASE(test_transfer_two_from_one_chain)
{
  // Test that if, like by PDBize, a single chain is used to connect two chains,
  // all the bonds stay in the right chain.
  EntityHandle src = CreateEntity();
  XCSEditor edi = src.EditXCS();
  ChainHandle fst_chain = edi.InsertChain(String("A"));
  ResidueHandle bma1 = edi.AppendResidue(fst_chain, "BMA");
  ResidueHandle man1 = edi.AppendResidue(fst_chain, "MAN");
  ResidueHandle man2 = edi.AppendResidue(fst_chain, "MAN");

  AtomHandle bma_o3;
  AtomHandle bma_o6;
  AtomHandle man_c1;

  bma_o3 = edi.InsertAtom(bma1, "O3", geom::Vec3(9.646, -29.415, 33.307));
  bma_o6 = edi.InsertAtom(bma1, "O6", geom::Vec3(8.642, -28.695, 27.405));

  man_c1 = edi.InsertAtom(man1, "C1", geom::Vec3(10.253, -30.563, 33.946));
  edi.Connect(bma_o3, man_c1);

  man_c1 = edi.InsertAtom(man2, "C1", geom::Vec3(8.091, -29.369, 26.233));
  edi.Connect(bma_o6, man_c1);

  EntityHandle dst = CreateEntity();
  edi = dst.EditXCS();
  fst_chain = edi.InsertChain(String("B"));
  bma1 = edi.AppendResidue(fst_chain, "BMA");
  man1 = edi.AppendResidue(fst_chain, "MAN");
  man2 = edi.AppendResidue(fst_chain, "MAN");

  edi.InsertAtom(bma1, "O3", geom::Vec3(9.646, -29.415, 33.307));
  edi.InsertAtom(bma1, "O6", geom::Vec3(8.642, -28.695, 27.405));

  edi.InsertAtom(man1, "C1", geom::Vec3(10.253, -30.563, 33.946));

  edi.InsertAtom(man2, "C1", geom::Vec3(8.091, -29.369, 26.233));

  ChainHandle snd_chain = edi.InsertChain(String("C"));
  bma1 = edi.AppendResidue(snd_chain, "BMA");
  man1 = edi.AppendResidue(snd_chain, "MAN");
  man2 = edi.AppendResidue(snd_chain, "MAN");

  edi.InsertAtom(bma1, "O3", geom::Vec3(-9.646, 29.415, 33.307));
  edi.InsertAtom(bma1, "O6", geom::Vec3(-8.642, 28.695, 27.405));

  edi.InsertAtom(man1, "C1", geom::Vec3(-10.253, 30.563, 33.946));

  edi.InsertAtom(man2, "C1", geom::Vec3(-8.091, 29.369, 26.233));

  std::map<ResidueHandle, ResidueHandle> to_from;
  ResidueHandleList r1 = src.GetResidueList();
  ResidueHandleList r2 = fst_chain.GetResidueList();
  for (size_t i = 0; i < r1.size(); ++i) {
    to_from[r2[i]] = r1[i];
  }
  r2 = snd_chain.GetResidueList();
  for (size_t i = 0; i < r1.size(); ++i) {
    to_from[r2[i]] = r1[i];
  }

  TransferConnectivity(dst, to_from);

  BOOST_CHECK(BondExists(dst.FindAtom("B", 1, "O3"),
                         dst.FindAtom("B", 2, "C1")));
  BOOST_CHECK(BondExists(dst.FindAtom("B", 1, "O6"),
                         dst.FindAtom("B", 3, "C1")));
  BOOST_CHECK(BondExists(dst.FindAtom("C", 1, "O3"),
                         dst.FindAtom("C", 2, "C1")));
  BOOST_CHECK(BondExists(dst.FindAtom("C", 1, "O6"),
                         dst.FindAtom("C", 3, "C1")));
  BOOST_CHECK_EQUAL(BondExists(dst.FindAtom("C", 1, "O6"),
                               dst.FindAtom("B", 3, "C1")), false);
  BOOST_CHECK_EQUAL(BondExists(dst.FindAtom("B", 1, "O6"),
                               dst.FindAtom("C", 3, "C1")), false);
}

BOOST_AUTO_TEST_CASE(test_transfer_two_combined_one_chain)
{
  // Test that if, like by PDBize, a single chain ends up two times in '_', the
  // connectivity does not get messy.
  EntityHandle src = CreateEntity();
  XCSEditor edi = src.EditXCS();
  ChainHandle fst_chain = edi.InsertChain(String("A"));
  ResidueHandle bma1 = edi.AppendResidue(fst_chain, "BMA");
  ResidueHandle man1 = edi.AppendResidue(fst_chain, "MAN");
  ResidueHandle man2 = edi.AppendResidue(fst_chain, "MAN");

  AtomHandle bma_o3;
  AtomHandle bma_o6;
  AtomHandle man_c1;

  bma_o3 = edi.InsertAtom(bma1, "O3", geom::Vec3(9.646, -29.415, 33.307));
  bma_o6 = edi.InsertAtom(bma1, "O6", geom::Vec3(8.642, -28.695, 27.405));

  man_c1 = edi.InsertAtom(man1, "C1", geom::Vec3(10.253, -30.563, 33.946));
  edi.Connect(bma_o3, man_c1);

  man_c1 = edi.InsertAtom(man2, "C1", geom::Vec3(8.091, -29.369, 26.233));
  edi.Connect(bma_o6, man_c1);

  EntityHandle dst = CreateEntity();
  edi = dst.EditXCS();
  fst_chain = edi.InsertChain(String("_"));
  bma1 = edi.AppendResidue(fst_chain, "BMA", ResNum(1, 'A'));
  man1 = edi.AppendResidue(fst_chain, "MAN", ResNum(1, 'B'));
  man2 = edi.AppendResidue(fst_chain, "MAN", ResNum(1, 'C'));

  edi.InsertAtom(bma1, "O3", geom::Vec3(9.646, -29.415, 33.307));
  edi.InsertAtom(bma1, "O6", geom::Vec3(8.642, -28.695, 27.405));
  edi.InsertAtom(man1, "C1", geom::Vec3(10.253, -30.563, 33.946));
  edi.InsertAtom(man2, "C1", geom::Vec3(8.091, -29.369, 26.233));

  bma1 = edi.AppendResidue(fst_chain, "BMA", ResNum(2, 'A'));
  man1 = edi.AppendResidue(fst_chain, "MAN", ResNum(2, 'B'));
  man2 = edi.AppendResidue(fst_chain, "MAN", ResNum(2, 'C'));

  edi.InsertAtom(bma1, "O3", geom::Vec3(-9.646, 29.415, 33.307));
  edi.InsertAtom(bma1, "O6", geom::Vec3(-8.642, 28.695, 27.405));
  edi.InsertAtom(man1, "C1", geom::Vec3(-10.253, 30.563, 33.946));
  edi.InsertAtom(man2, "C1", geom::Vec3(-8.091, 29.369, 26.233));

  std::map<ResidueHandle, ResidueHandle> to_from;
  ResidueHandleList r1 = src.GetResidueList();
  ResidueHandleList r2 = fst_chain.GetResidueList();
  for (size_t i = 0; i < r1.size(); ++i) {
    to_from[r2[i]] = r1[i];
  }
  for (size_t i = 0; i < r1.size(); ++i) {
    to_from[r2[i+3]] = r1[i];
  }

  TransferConnectivity(dst, to_from);

  BOOST_CHECK(BondExists(dst.FindAtom("_", ResNum(1, 'A'), "O3"),
                         dst.FindAtom("_", ResNum(1, 'B'), "C1")));
  BOOST_CHECK(BondExists(dst.FindAtom("_", ResNum(1, 'A'), "O6"),
                         dst.FindAtom("_", ResNum(1, 'C'), "C1")));
  BOOST_CHECK(BondExists(dst.FindAtom("_", ResNum(2, 'A'), "O3"),
                         dst.FindAtom("_", ResNum(2, 'B'), "C1")));
  BOOST_CHECK(BondExists(dst.FindAtom("_", ResNum(2, 'A'), "O6"),
                         dst.FindAtom("_", ResNum(2, 'C'), "C1")));
  BOOST_CHECK_EQUAL(BondExists(dst.FindAtom("B", ResNum(2, 'A'), "O6"),
                               dst.FindAtom("C", ResNum(1, 'C'), "C1")), false);
}

BOOST_AUTO_TEST_SUITE_END();
