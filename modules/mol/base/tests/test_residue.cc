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


using namespace ost;
using namespace ost::mol;

BOOST_AUTO_TEST_SUITE( mol_base );

BOOST_AUTO_TEST_CASE(test_in_sequence) 
{
  EntityHandle eh=CreateEntity();
  XCSEditor e=eh.EditXCS();
  ChainHandle ch1=e.InsertChain("A");
  ResidueHandle rA = e.AppendResidue(ch1, "A");
  ResidueHandle rB = e.AppendResidue(ch1, "B");
  ResidueHandle rC = e.AppendResidue(ch1, "C");  

  BOOST_CHECK(!InSequence(rA,rB));
  BOOST_CHECK(!InSequence(rB,rC));
  BOOST_CHECK(!InSequence(rA,rC));
}


BOOST_AUTO_TEST_CASE(test_res_index_bzdng227) 
{
  EntityHandle eh=CreateEntity();
  XCSEditor e=eh.EditXCS();
  ChainHandle ch1=e.InsertChain("A");
  ResidueHandle ra = e.AppendResidue(ch1, "A", 1);
  ResidueHandle rb = e.AppendResidue(ch1, "B", 2);
  ResidueHandle rc = e.AppendResidue(ch1, "C", 1);  
  
  BOOST_CHECK_EQUAL(ra.GetIndex(), 0);
  BOOST_CHECK_EQUAL(rb.GetIndex(), 1);
  BOOST_CHECK_EQUAL(rc.GetIndex(), 2);
}


BOOST_AUTO_TEST_CASE(throw_invalid_res_handle)
{
  ChainHandle chain;
  EntityHandle ent=CreateEntity();
  XCSEditor edi=ent.EditXCS();
  chain=edi.InsertChain("A");
  ResidueHandle res=chain.FindResidue(ResNum(1));
  BOOST_CHECK_THROW(CheckHandleValidity(res), InvalidHandle);
  edi.AppendResidue(chain, "GLY");
  res=chain.FindResidue(ResNum(1));
  BOOST_CHECK_NO_THROW(CheckHandleValidity(res));
}

BOOST_AUTO_TEST_CASE(throw_invalid_res_view)
{
  ChainHandle chain;
  EntityHandle ent=CreateEntity();
  XCSEditor edi=ent.EditXCS();
  chain=edi.InsertChain("A");
  ResidueView res;
  BOOST_CHECK_THROW(CheckHandleValidity(res), InvalidHandle);
  edi.AppendResidue(chain, "GLY");
  EntityView ent_view=ent.CreateFullView();
  res=ent_view.FindChain("A").FindResidue(1);
  BOOST_CHECK_NO_THROW(CheckHandleValidity(res));
}

// Test that .GetHashCode() throws on invalid residue handle
BOOST_AUTO_TEST_CASE(throw_hash_code_invalid_res_handle)
{
  ChainHandle chain;
  EntityHandle ent=CreateEntity();
  XCSEditor edi=ent.EditXCS();
  chain=edi.InsertChain("A");
  ResidueHandle res=chain.FindResidue(ResNum(1));
  BOOST_CHECK_THROW(res.GetHashCode(), InvalidHandle);
  edi.AppendResidue(chain, "GLY");
  res=chain.FindResidue(ResNum(1));
  BOOST_CHECK_NO_THROW(res.GetHashCode());
}

// Test that .GetHashCode() throws on invalid residue views
BOOST_AUTO_TEST_CASE(throw_hash_code_invalid_res_view)
{
  ChainHandle chain;
  EntityHandle ent=CreateEntity();
  XCSEditor edi=ent.EditXCS();
  chain=edi.InsertChain("A");
  ResidueView res;
  BOOST_CHECK_THROW(res.GetHashCode(), InvalidHandle);
  edi.AppendResidue(chain, "GLY");
  EntityView ent_view=ent.CreateFullView();
  res=ent_view.FindChain("A").FindResidue(1);
  BOOST_CHECK_NO_THROW(res.GetHashCode());
}

BOOST_AUTO_TEST_CASE(test_res_index)
{
  EntityHandle eh=CreateEntity();
  XCSEditor e=eh.EditXCS();
  ChainHandle ch1=e.InsertChain("A");
  ResidueHandle rA = e.AppendResidue(ch1, "A");
  ResidueHandle rB = e.AppendResidue(ch1, "B");
  ResidueHandle rC = e.AppendResidue(ch1, "C");
  EntityView v(eh);
  ChainView cv1(v, ch1);
  ResidueView rvA = cv1.AddResidue(rA);
  ResidueView rvB = cv1.AddResidue(rB);
  ResidueView rvC = cv1.AddResidue(rC);
  BOOST_CHECK_EQUAL(rvA.GetIndex(),0);
  BOOST_CHECK_EQUAL(rvA.GetIndex(),cv1.GetResidueIndex(rvA.GetNumber()));
  BOOST_CHECK_EQUAL(rvB.GetIndex(),1);
  BOOST_CHECK_EQUAL(rvB.GetIndex(),cv1.GetResidueIndex(rvB.GetNumber()));
  BOOST_CHECK_EQUAL(rvC.GetIndex(),2);
  BOOST_CHECK_EQUAL(rvC.GetIndex(),cv1.GetResidueIndex(rvC.GetNumber()));
  cv1.RemoveResidue(rvB);
  BOOST_CHECK_EQUAL(rvA.GetIndex(),0);
  BOOST_CHECK_EQUAL(rvA.GetIndex(),cv1.GetResidueIndex(rvA.GetNumber()));
  BOOST_CHECK_EQUAL(rvC.GetIndex(),1);
  BOOST_CHECK_EQUAL(rvC.GetIndex(),cv1.GetResidueIndex(rvC.GetNumber()));
}

BOOST_AUTO_TEST_CASE(rename_res)
{
   EntityHandle eh=CreateEntity();
   XCSEditor e=eh.EditXCS();
   ChainHandle ch1=e.InsertChain("A");
   ResidueHandle rA2B = e.AppendResidue(ch1, "A");
   e.RenameResidue(rA2B, "B");
   BOOST_CHECK_EQUAL(rA2B.GetName(), "B");
}

BOOST_AUTO_TEST_CASE(test_centralatom)
{
  // COOK UP ENTITY FOR TEST
  EntityHandle eh = CreateEntity();
  XCSEditor e = eh.EditXCS();
  ChainHandle ch = e.InsertChain("A");
  // decent peptide with all entries
  ResidueHandle rp1 = e.AppendResidue(ch, "A");
  e.InsertAtom(rp1, "CA", geom::Vec3(2, 0, 0));
  e.InsertAtom(rp1, "CB", geom::Vec3(1, 0, 0));
  e.InsertAtom(rp1, "C",  geom::Vec3(0, 0, 0));
  e.InsertAtom(rp1, "O",  geom::Vec3(0, 1, 0));
  rp1.SetChemClass(ChemClass(ChemClass::PEPTIDE_LINKING));
  // weird peptide with only CA and CB
  ResidueHandle rp2 = e.AppendResidue(ch, "B");
  e.InsertAtom(rp2, "CA", geom::Vec3(3, 0, 0));
  e.InsertAtom(rp2, "CB", geom::Vec3(3, 1, 0));
  rp2.SetChemClass(ChemClass(ChemClass::PEPTIDE_LINKING));
  // CA-only peptide
  ResidueHandle rp3 = e.AppendResidue(ch, "C");
  e.InsertAtom(rp3, "CA", geom::Vec3(4, 0, 0));
  rp3.SetChemClass(ChemClass(ChemClass::PEPTIDE_LINKING));
  // peptide with custom atoms
  ResidueHandle rp4 = e.AppendResidue(ch, "D");
  AtomHandle rp4_ax = e.InsertAtom(rp4, "XX", geom::Vec3(5, 0, 0));
  rp4.SetChemClass(ChemClass(ChemClass::PEPTIDE_LINKING));
  // nucleotide with all needed entries
  ResidueHandle rn1 = e.AppendResidue(ch, "E");
  e.InsertAtom(rn1, "P", geom::Vec3(6, 0, 0));
  e.InsertAtom(rn1, "OP1", geom::Vec3(6, 0.5, 0));
  e.InsertAtom(rn1, "OP2", geom::Vec3(6, 1.5, 0));
  rn1.SetChemClass(ChemClass(ChemClass::DNA_LINKING));
  // nucleotide with only P
  ResidueHandle rn2 = e.AppendResidue(ch, "F");
  e.InsertAtom(rn2, "P", geom::Vec3(7, 0, 0));
  rn2.SetChemClass(ChemClass(ChemClass::DNA_LINKING));
  // nucleotide with custom atoms
  ResidueHandle rn3 = e.AppendResidue(ch, "G");
  AtomHandle rn3_ax = e.InsertAtom(rn3, "XX", geom::Vec3(8, 0, 0));
  rn3.SetChemClass(ChemClass(ChemClass::DNA_LINKING));
  // unknown chem class
  ResidueHandle ru = e.AppendResidue(ch, "H");
  e.InsertAtom(ru, "P",  geom::Vec3(9, 0, 0));
  e.InsertAtom(ru, "CA", geom::Vec3(9, 1, 0));
  AtomHandle ru_ax = e.InsertAtom(ru, "XX", geom::Vec3(9, 2, 0));
  ru.SetChemClass(ChemClass(ChemClass::UNKNOWN));

  // CHECK CENTRAL ATOMS
  BOOST_CHECK(rp1.GetCentralAtom().IsValid());
  BOOST_CHECK_EQUAL(rp1.GetCentralAtom().GetQualifiedName(), "A.A1.CA");
  BOOST_CHECK(rp2.GetCentralAtom().IsValid());
  BOOST_CHECK_EQUAL(rp2.GetCentralAtom().GetQualifiedName(), "A.B2.CA");
  BOOST_CHECK(rp3.GetCentralAtom().IsValid());
  BOOST_CHECK_EQUAL(rp3.GetCentralAtom().GetQualifiedName(), "A.C3.CA");
  BOOST_CHECK(!rp4.GetCentralAtom().IsValid());
  BOOST_CHECK(rn1.GetCentralAtom().IsValid());
  BOOST_CHECK_EQUAL(rn1.GetCentralAtom().GetQualifiedName(), "A.E5.P");
  BOOST_CHECK(rn2.GetCentralAtom().IsValid());
  BOOST_CHECK_EQUAL(rn2.GetCentralAtom().GetQualifiedName(), "A.F6.P");
  BOOST_CHECK(!rn3.GetCentralAtom().IsValid());
  BOOST_CHECK(!ru.GetCentralAtom().IsValid());

  // CHECK NORMALS
  BOOST_CHECK_EQUAL(rp1.GetCentralNormal(), geom::Vec3(0, 1, 0));
  BOOST_CHECK_EQUAL(rp2.GetCentralNormal(), geom::Vec3(0, -1, 0));
  BOOST_CHECK_EQUAL(rp3.GetCentralNormal(), geom::Vec3(0, 0, 1));
  BOOST_CHECK_EQUAL(rp4.GetCentralNormal(), geom::Vec3(1, 0, 0));
  BOOST_CHECK_EQUAL(rn1.GetCentralNormal(), geom::Vec3(0, -1, 0));
  BOOST_CHECK_EQUAL(rn2.GetCentralNormal(), geom::Vec3(0, 0, 1));
  BOOST_CHECK_EQUAL(rn3.GetCentralNormal(), geom::Vec3(1, 0, 0));
  BOOST_CHECK_EQUAL(ru.GetCentralNormal(), geom::Vec3(1, 0, 0));

  // CHECK SETTING CENTRAL ATOMS
  rp4.SetCentralAtom(rp4_ax);
  BOOST_CHECK(rp4.GetCentralAtom().IsValid());
  BOOST_CHECK_EQUAL(rp4.GetCentralAtom().GetQualifiedName(), "A.D4.XX");
  BOOST_CHECK_EQUAL(rp4.GetCentralNormal(), geom::Vec3(0, 0, 1));
  rn3.SetCentralAtom(rn3_ax);
  BOOST_CHECK(rn3.GetCentralAtom().IsValid());
  BOOST_CHECK_EQUAL(rn3.GetCentralAtom().GetQualifiedName(), "A.G7.XX");
  BOOST_CHECK_EQUAL(rn3.GetCentralNormal(), geom::Vec3(0, 0, 1));
  ru.SetCentralAtom(ru_ax);
  BOOST_CHECK(ru.GetCentralAtom().IsValid());
  BOOST_CHECK_EQUAL(ru.GetCentralAtom().GetQualifiedName(), "A.H8.XX");
  // no normal for unknown residues
  BOOST_CHECK_EQUAL(ru.GetCentralNormal(), geom::Vec3(1, 0, 0));
}

BOOST_AUTO_TEST_SUITE_END();
