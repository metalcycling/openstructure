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
#include <ost/integrity_error.hh>

using namespace ost;
using namespace ost::mol;

void find_and_check_res(ChainHandle chain, const ResNum& n) {
  ResidueHandle fa=chain.FindResidue(n);
  BOOST_CHECK_MESSAGE(fa.IsValid(), "residue " << n << " not found");
  if (!fa.IsValid())
    return;
  BOOST_CHECK_MESSAGE(fa.GetNumber()==n,  
                      "residue " << n << " expected, but "
                      << fa.GetNumber() << " found");             
}

void no_find_res(ChainHandle chain, const ResNum& n) {
  ResidueHandle fa=chain.FindResidue(n);
  BOOST_CHECK_MESSAGE(fa.IsValid()==false, "residue "
                      << " found while looking for " << n
                      << " (should be invalid)");    
}

void test_res_pos() {
  EntityHandle eh=CreateEntity();
  XCSEditor e=eh.EditXCS();  
  ChainHandle ch1=e.InsertChain("A");
  e.AppendResidue(ch1, "A");
  e.AppendResidue(ch1, "B");
  e.AppendResidue(ch1, "C", ResNum(5));
  e.AppendResidue(ch1, "D");  
  e.AppendResidue(ch1, "E", ResNum(9));
  e.AppendResidue(ch1, "F");  
  e.AppendResidue(ch1, "E", ResNum(12));  
  e.AppendResidue(ch1, "G");    
  find_and_check_res(ch1, ResNum(1));
  find_and_check_res(ch1, ResNum(2));  
  no_find_res(ch1, ResNum(3));
  no_find_res(ch1, ResNum(4));  
  find_and_check_res(ch1, ResNum(5));    
  find_and_check_res(ch1, ResNum(6));    
  no_find_res(ch1, ResNum(7));    
  no_find_res(ch1, ResNum(8));    
  find_and_check_res(ch1, ResNum(9));
  find_and_check_res(ch1, ResNum(10));  
  no_find_res(ch1, ResNum(11));    
  find_and_check_res(ch1, ResNum(12));    
  find_and_check_res(ch1, ResNum(13));  
}

BOOST_AUTO_TEST_SUITE( mol_base );

BOOST_AUTO_TEST_CASE(test_comparison) 
{
  EntityHandle eh=CreateEntity();
  XCSEditor e=eh.EditXCS();  
  ChainHandle ch1=e.InsertChain("A");
  ChainHandle ch2=e.InsertChain("B");
  ChainHandle cc=ch1;
  BOOST_CHECK(ch1!=ch2);
  BOOST_CHECK(ch1==ch1);
  BOOST_CHECK(cc==ch1);
}

BOOST_AUTO_TEST_CASE(throw_invalid_chain_handle)
{
  EntityHandle ent=CreateEntity();
  ChainHandle ch=ent.FindChain("A");
  BOOST_CHECK_THROW(CheckHandleValidity(ch), InvalidHandle);
  XCSEditor edi=ent.EditXCS();
  edi.InsertChain("A");
  ch=ent.FindChain("A");
  BOOST_CHECK_NO_THROW(CheckHandleValidity(ch));


  EntityHandle eh1 = CreateEntity();
  XCSEditor e1=eh1.EditXCS();
  ChainHandle chain = e1.InsertChain("C");
  ResidueHandle res1 = e1.AppendResidue(chain, "FOO", ResNum(13));
  AtomHandle atom1 = e1.InsertAtom(res1, "X1",geom::Vec3());
  AtomHandle atom2 = e1.InsertAtom(res1, "X2",geom::Vec3());
  ResidueHandle res2 = e1.AppendResidue(chain, "FOO", ResNum(42));
  e1.InsertAtom(res2, "X1",geom::Vec3());
  e1.InsertAtom(res2, "X2",geom::Vec3());

  EntityHandle eh2=CreateEntity();
  XCSEditor e2=eh2.EditXCS();
  ChainHandle inserted_chain1 = e2.InsertChain("Q",chain);
  res1.SetIntProp("amazing",42);
  ResidueHandle inserted_residue1 = e2.AppendResidue(inserted_chain1,res1);
  BOOST_CHECK(inserted_residue1!=res1);
  BOOST_CHECK(eh2.GetResidueCount()==1);
  BOOST_CHECK(eh2.GetAtomCount()==0);
  BOOST_CHECK(inserted_residue1.HasProp("amazing"));
  ResidueHandle inserted_residue2 = e2.AppendResidue(inserted_chain1,res2,true);
  BOOST_CHECK(eh2.GetAtomCount()==2);
}

BOOST_AUTO_TEST_CASE(throw_invalid_chain_view)
{
  EntityHandle ent=CreateEntity();
  ChainHandle chain=ent.FindChain("A");
  ChainView chain_view;
  BOOST_CHECK_THROW(CheckHandleValidity(chain_view), InvalidHandle);
  XCSEditor edi=ent.EditXCS();
  edi.InsertChain("A");
  EntityView ent_view=ent.CreateFullView();  
  BOOST_CHECK_NO_THROW(CheckHandleValidity(ent_view.FindChain("A")));
}

BOOST_AUTO_TEST_CASE(res_pos) 
{
  EntityHandle eh=CreateEntity();
  XCSEditor e=eh.EditXCS();  
  ChainHandle ch1=e.InsertChain("A");
  e.AppendResidue(ch1, "A");
  e.AppendResidue(ch1, "B");
  e.AppendResidue(ch1, "C", ResNum(5));
  e.AppendResidue(ch1, "D");  
  e.AppendResidue(ch1, "E", ResNum(9));
  e.AppendResidue(ch1, "F");  
  e.AppendResidue(ch1, "E", ResNum(12));  
  e.AppendResidue(ch1, "G");    
  find_and_check_res(ch1, ResNum(1));
  find_and_check_res(ch1, ResNum(2));  
  no_find_res(ch1, ResNum(3));
  no_find_res(ch1, ResNum(4));  
  find_and_check_res(ch1, ResNum(5));    
  find_and_check_res(ch1, ResNum(6));    
  no_find_res(ch1, ResNum(7));    
  no_find_res(ch1, ResNum(8));    
  find_and_check_res(ch1, ResNum(9));
  find_and_check_res(ch1, ResNum(10));  
  no_find_res(ch1, ResNum(11));    
  find_and_check_res(ch1, ResNum(12));    
  find_and_check_res(ch1, ResNum(13));
}


BOOST_AUTO_TEST_CASE(res_pos_with_insertion_code) 
{
  EntityHandle eh=CreateEntity();
  XCSEditor e=eh.EditXCS();  
  ChainHandle ch1=e.InsertChain("A");
  e.AppendResidue(ch1, "A");
  e.AppendResidue(ch1, "B");
  // test for normal insertion
  e.AppendResidue(ch1, "C", ResNum(5, 'A'));
  e.AppendResidue(ch1, "D", ResNum(5, 'B'));  
  e.AppendResidue(ch1, "E", ResNum(6));  
  e.AppendResidue(ch1, "E", ResNum(7));
  // test for insertion plus gap right after
  e.AppendResidue(ch1, "E", ResNum(7, 'B'));  
  e.AppendResidue(ch1, "E", ResNum(9));      
  find_and_check_res(ch1, ResNum(1));
  find_and_check_res(ch1, ResNum(2));  
  no_find_res(ch1, ResNum(3));
  no_find_res(ch1, ResNum(4));  
  find_and_check_res(ch1, ResNum(5, 'A'));    
  find_and_check_res(ch1, ResNum(5, 'B'));    
  find_and_check_res(ch1, ResNum(6));        

  find_and_check_res(ch1, ResNum(7));      
  find_and_check_res(ch1, ResNum(7, 'B'));   
  no_find_res(ch1, ResNum(8));  
  find_and_check_res(ch1, ResNum(9));        
  no_find_res(ch1, ResNum(10));
}

BOOST_AUTO_TEST_CASE(prev_next) 
{
  EntityHandle eh=CreateEntity();
  XCSEditor e=eh.EditXCS();  
  ChainHandle ch1=e.InsertChain("A");
  e.AppendResidue(ch1, "A");
  e.AppendResidue(ch1, "B");
  e.AppendResidue(ch1, "C");  
  
  BOOST_REQUIRE(ch1.GetPrev(ch1.FindResidue(ResNum(2))).IsValid());
  BOOST_CHECK(ch1.GetPrev(ch1.FindResidue(ResNum(2))).GetKey()=="A");  
  BOOST_REQUIRE(ch1.GetNext(ch1.FindResidue(ResNum(2))).IsValid());
  BOOST_CHECK(ch1.GetNext(ch1.FindResidue(ResNum(2))).GetKey()=="C");
  BOOST_CHECK(!ch1.GetPrev(ch1.FindResidue(ResNum(1))).IsValid());  
  BOOST_CHECK(!ch1.GetNext(ch1.FindResidue(ResNum(3))).IsValid());    
  BOOST_CHECK(!ch1.GetPrev(ResidueHandle()).IsValid());
  BOOST_CHECK(!ch1.GetNext(ResidueHandle()).IsValid());
}

BOOST_AUTO_TEST_CASE(rename_chain)
{
   EntityHandle eh=CreateEntity();
   XCSEditor e=eh.EditXCS();
   ChainHandle ch1=e.InsertChain("A");
   ChainHandle ch2=e.InsertChain("B");
   e.RenameChain(ch1, "A"); // renaming chain with its current name should work
   BOOST_CHECK_EQUAL(ch1.GetName(), "A");
   BOOST_CHECK_THROW(e.RenameChain(ch1, "B"), IntegrityError);
   e.RenameChain(ch2, "C");
   BOOST_CHECK_EQUAL(ch2.GetName(), "C");
   BOOST_CHECK_EQUAL(eh.GetChainCount(), 2);
}

BOOST_AUTO_TEST_CASE(chain_type)
{
   EntityHandle eh = CreateEntity();
   XCSEditor e = eh.EditXCS();
   ChainHandle ch1 = e.InsertChain("A");

   // setting/ getting
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_UNKNOWN);
   BOOST_CHECK(!ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   e.SetChainType(ch1, CHAINTYPE_POLY);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_POLY);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   e.SetChainType(ch1, CHAINTYPE_NON_POLY);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_NON_POLY);
   BOOST_CHECK(!ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   e.SetChainType(ch1, CHAINTYPE_WATER);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_WATER);
   BOOST_CHECK(!ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   e.SetChainType(ch1, CHAINTYPE_POLY_PEPTIDE_D);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_POLY_PEPTIDE_D);
   e.SetChainType(ch1, CHAINTYPE_POLY_PEPTIDE_L);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_POLY_PEPTIDE_L);
   e.SetChainType(ch1, CHAINTYPE_POLY_DN);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(ch1.IsPolynucleotide());
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_POLY_DN);
   e.SetChainType(ch1, CHAINTYPE_POLY_RN);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(ch1.IsPolynucleotide());
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_POLY_RN);
   e.SetChainType(ch1, CHAINTYPE_POLY_SAC_D);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_POLY_SAC_D);
   e.SetChainType(ch1, CHAINTYPE_POLY_SAC_L);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_POLY_SAC_L);
   e.SetChainType(ch1, CHAINTYPE_POLY_DN_RN);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_POLY_DN_RN);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(ch1.IsPolynucleotide());
   e.SetChainType(ch1, CHAINTYPE_N_CHAINTYPES);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_N_CHAINTYPES);
   e.SetChainType(ch1, CHAINTYPE_UNKNOWN);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_UNKNOWN);
   e.SetChainType(ch1, CHAINTYPE_MACROLIDE);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_MACROLIDE);
   BOOST_CHECK(!ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   e.SetChainType(ch1, CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   e.SetChainType(ch1, CHAINTYPE_POLY_PEPTIDE_DN_RN);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_POLY_PEPTIDE_DN_RN);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   e.SetChainType(ch1, CHAINTYPE_BRANCHED);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_BRANCHED);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(!ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());
   e.SetChainType(ch1, CHAINTYPE_OLIGOSACCHARIDE);
   BOOST_CHECK(ch1.GetType() == CHAINTYPE_OLIGOSACCHARIDE);
   BOOST_CHECK(ch1.IsPolymer());
   BOOST_CHECK(!ch1.IsPolysaccharide());
   BOOST_CHECK(ch1.IsOligosaccharide());
   BOOST_CHECK(!ch1.IsPolypeptide());
   BOOST_CHECK(!ch1.IsPolynucleotide());

   // string -> chain type
   BOOST_CHECK(ChainTypeFromString("polymer") == CHAINTYPE_POLY);
   BOOST_CHECK(ChainTypeFromString("non-polymer") == CHAINTYPE_NON_POLY);
   BOOST_CHECK(ChainTypeFromString("water") == CHAINTYPE_WATER);
   BOOST_CHECK(ChainTypeFromString("polypeptide(D)") ==
               CHAINTYPE_POLY_PEPTIDE_D);
   BOOST_CHECK(ChainTypeFromString("polypeptide(L)") ==
               CHAINTYPE_POLY_PEPTIDE_L);
   BOOST_CHECK(ChainTypeFromString("polydeoxyribonucleotide") ==
               CHAINTYPE_POLY_DN);
   BOOST_CHECK(ChainTypeFromString("polyribonucleotide") ==
               CHAINTYPE_POLY_RN);
   BOOST_CHECK(ChainTypeFromString("polysaccharide(D)") ==
               CHAINTYPE_POLY_SAC_D);
   BOOST_CHECK(ChainTypeFromString("polysaccharide(L)") ==
               CHAINTYPE_POLY_SAC_L);
   BOOST_CHECK(ChainTypeFromString(
                      "polydeoxyribonucleotide/polyribonucleotide hybrid") ==
               CHAINTYPE_POLY_DN_RN);
   BOOST_CHECK(ChainTypeFromString("other") == CHAINTYPE_UNKNOWN);
   BOOST_CHECK(ChainTypeFromString("macrolide") == CHAINTYPE_MACROLIDE);
   BOOST_CHECK(ChainTypeFromString("cyclic-pseudo-peptide") ==
               CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE);
   BOOST_CHECK(ChainTypeFromString("peptide nucleic acid") ==
               CHAINTYPE_POLY_PEPTIDE_DN_RN);
   BOOST_CHECK(ChainTypeFromString("branched") == CHAINTYPE_BRANCHED);
   BOOST_CHECK(ChainTypeFromString("oligosaccharide") ==
               CHAINTYPE_OLIGOSACCHARIDE);
   BOOST_CHECK_THROW(ChainTypeFromString("supposed to fail"),
                     Error);

   // chain type -> string
   BOOST_CHECK(StringFromChainType(CHAINTYPE_POLY) == "polymer");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_NON_POLY) == "non-polymer");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_WATER) == "water");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_POLY_PEPTIDE_D) ==
               "polypeptide(D)");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_POLY_PEPTIDE_L) ==
               "polypeptide(L)");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_POLY_DN) ==
               "polydeoxyribonucleotide");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_POLY_RN) == "polyribonucleotide");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_POLY_SAC_D) ==
               "polysaccharide(D)");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_POLY_SAC_L) ==
               "polysaccharide(L)");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_POLY_DN_RN) ==
               "polydeoxyribonucleotide/polyribonucleotide hybrid");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_UNKNOWN) == "other");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_MACROLIDE) == "macrolide");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE) ==
               "cyclic-pseudo-peptide");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_POLY_PEPTIDE_DN_RN) ==
               "peptide nucleic acid");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_BRANCHED) == "branched");
   BOOST_CHECK(StringFromChainType(CHAINTYPE_OLIGOSACCHARIDE) ==
               "oligosaccharide");
   BOOST_CHECK_THROW(StringFromChainType(CHAINTYPE_N_CHAINTYPES),
                     Error);
}

BOOST_AUTO_TEST_CASE(chain_description)
{
  EntityHandle eh=CreateEntity();
  XCSEditor e = eh.EditXCS();
  ChainHandle ch1 = e.InsertChain("A");
  String description = "Very important information";
  e.SetChainDescription(ch1, description);
  BOOST_CHECK(ch1.GetDescription() == description);
}

BOOST_AUTO_TEST_SUITE_END();
