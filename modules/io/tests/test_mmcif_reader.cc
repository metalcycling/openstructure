//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2011 by the OpenStructure authors
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

#include <fstream>
#include <ost/platform.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/mol/mmcif_reader.hh>
#include <ost/conop/conop.hh>
#include <ost/conop/rule_based_builder.hh>

#define BOOST_AUTO_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


using namespace ost;
using namespace ost::io;

class TestMMCifParserProtected : MMCifParser {
public:
  TestMMCifParserProtected(std::istream& stream, mol::EntityHandle& ent_handle):
    MMCifParser(stream, ent_handle, IOProfile())
  { }

  TestMMCifParserProtected(std::istream& stream, mol::EntityHandle& ent_handle,
                           const IOProfile& profile):
    MMCifParser(stream, ent_handle, profile)
  { }

  TestMMCifParserProtected(const String& filename,
                           mol::EntityHandle& ent_handle):
    MMCifParser(filename, ent_handle, IOProfile())
  { }

  using MMCifParser::OnBeginLoop;
  using MMCifParser::OnEndData;
  using MMCifParser::IsValidPDBIdent;
  using MMCifParser::ParseAtomIdent;
  using MMCifParser::ParseAndAddAtom;
  using MMCifParser::ParseEntity;
  using MMCifParser::ParseEntityPoly;
  using MMCifParser::ParseCitation;
  using MMCifParser::ParseRefine;
  using MMCifParser::ParsePdbxStructAssemblyGen;
  using MMCifParser::ParsePdbxStructAssembly;
  using MMCifParser::ParsePdbxStructOperList;
  using MMCifParser::ParseStruct;
  using MMCifParser::TryStoreIdx;
  using MMCifParser::SetReadSeqRes;
  using MMCifParser::SetReadCanonicalSeqRes;
  using MMCifParser::ClearState;
  using MMCifParser::ConvertSEQRES;
  using MMCifParser::GetInfo;
  using MMCifParser::DetermineSecStructType;
  using MMCifParser::MMCifSecStructElement;
  using MMCifParser::MMCIF_HELIX;
  using MMCifParser::MMCIF_TURN;
  using MMCifParser::MMCIF_STRAND;
};

void SetAtomSiteHeader(StarLoopDesc* mmcif_h)
{
  mmcif_h->Clear();
  mmcif_h->SetCategory(StringRef("atom_site", 9));
  mmcif_h->Add(StringRef("auth_asym_id", 12));
  mmcif_h->Add(StringRef("id", 2));
  mmcif_h->Add(StringRef("label_alt_id", 12));
  mmcif_h->Add(StringRef("label_asym_id", 13));
  mmcif_h->Add(StringRef("label_atom_id", 13));
  mmcif_h->Add(StringRef("label_comp_id", 13));
  mmcif_h->Add(StringRef("label_entity_id", 15));
  mmcif_h->Add(StringRef("label_seq_id", 12));
  mmcif_h->Add(StringRef("type_symbol", 11));
  mmcif_h->Add(StringRef("Cartn_x", 7));
  mmcif_h->Add(StringRef("Cartn_y", 7));
  mmcif_h->Add(StringRef("Cartn_z", 7));
}

BOOST_AUTO_TEST_SUITE( io );

BOOST_AUTO_TEST_CASE(mmcif_isvalidpdbident)
{
  mol::EntityHandle eh=mol::CreateEntity();

  // on changing the tests for a PDB id in mmcif files, extend this unit test
  BOOST_MESSAGE("  Running mmcif_isvalidpdbident tests...");
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  TestMMCifParserProtected tmmcif_p(s, eh);
  StringRef id = StringRef("1FOO", 4);
  BOOST_MESSAGE("    Testing valid id ('"+ id.str() +"')...");
  BOOST_CHECK(tmmcif_p.IsValidPDBIdent(id));
  BOOST_MESSAGE("    done.");
  id = StringRef("this is to long for a PDB id", 28);
  BOOST_MESSAGE("    Testing oversized PDB id ('"+ id.str() +"')...");
  BOOST_CHECK(!tmmcif_p.IsValidPDBIdent(id));
  BOOST_MESSAGE("    done.");
  id = StringRef("nFOO", 4);
  BOOST_MESSAGE("    Testing PDB id with missing number ('"+ id.str() +"')...");
  BOOST_CHECK(!tmmcif_p.IsValidPDBIdent(id));
  BOOST_MESSAGE("    done.");
}

BOOST_AUTO_TEST_CASE(mmcif_trystoreidx)
{
  mol::EntityHandle eh = mol::CreateEntity();

  BOOST_MESSAGE("  Running mmcif_trystoreidx tests...");
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  TestMMCifParserProtected tmmcif_p(s, eh, IOProfile());
  StarLoopDesc mmcif_h;
  mmcif_h.SetCategory(StringRef("Foo", 3));
  // negative
  BOOST_CHECK_THROW(tmmcif_p.TryStoreIdx(0, "bar", mmcif_h), IOException);
  // positive
  mmcif_h.Add(StringRef("bar", 3));
  BOOST_CHECK_NO_THROW(tmmcif_p.TryStoreIdx(0, "bar", mmcif_h));
}

BOOST_AUTO_TEST_CASE(mmcif_convert_seqres)
{
  SetPrefixPath(getenv("OST_ROOT"));
  String lib_path=GetSharedDataPath()+"/compounds.chemlib";
  conop::CompoundLibPtr compound_lib=conop::CompoundLib::Load(lib_path);  
  if (!compound_lib) {
    std::cout << "WARNING: skipping SEQRES import unit test. " 
              << "Rule-based builder is required" << std::endl;
    return;    
  }
  conop::RuleBasedBuilderPtr rbb(new conop::RuleBasedBuilder(compound_lib));
  conop::Conopology::Instance().RegisterBuilder(rbb, "RBB");
  conop::Conopology::Instance().SetDefaultBuilder("RBB");
  mol::EntityHandle eh=mol::CreateEntity();
  
  TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
  std::vector<StringRef> columns;
  StarLoopDesc tmmcif_h;  
  BOOST_CHECK_EQUAL(tmmcif_p.ConvertSEQRES("A(MSE)Y", compound_lib), "AMY");
  BOOST_CHECK_THROW(tmmcif_p.ConvertSEQRES("A(MSEY", compound_lib), 
                    IOException);
  conop::Conopology::Instance().SetDefaultBuilder("HEURISTIC");
}

BOOST_AUTO_TEST_CASE(mmcif_onbeginloop)
{
  mol::EntityHandle eh=mol::CreateEntity();

  // add more tests on new mandatory items
  BOOST_MESSAGE("  Running mmcif_onbeginloop tests...");
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  MMCifParser mmcif_p(s, eh, IOProfile());
  StarLoopDesc mmcif_h;
  BOOST_MESSAGE("          testing atom_site items...");
  mmcif_h.SetCategory(StringRef("atom_site", 9));
  BOOST_MESSAGE("             auth_asym_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("auth_asym_id", 12));
  BOOST_MESSAGE("             id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("id", 2));
  BOOST_MESSAGE("             label_alt_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_alt_id", 12));
  BOOST_MESSAGE("             label_asym_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_asym_id", 13));
  BOOST_MESSAGE("             label_atom_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_atom_id", 13));
  BOOST_MESSAGE("             label_comp_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_comp_id", 13));
  BOOST_MESSAGE("             label_entity_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_entity_id", 15));
  BOOST_MESSAGE("             label_seq_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_seq_id", 12));
  BOOST_MESSAGE("             type_symbol");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("type_symbol", 11));
  BOOST_MESSAGE("             Cartn_x");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("Cartn_x", 7));
  BOOST_MESSAGE("             Cartn_y");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("Cartn_y", 7));
  BOOST_MESSAGE("             Cartn_z");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("Cartn_z", 7));
  BOOST_CHECK_NO_THROW(mmcif_p.OnBeginLoop(mmcif_h));
  BOOST_MESSAGE("          done.");
  mmcif_h.Clear();
  BOOST_MESSAGE("          testing entity items...");
  mmcif_h.SetCategory(StringRef("entity", 6));
  BOOST_MESSAGE("             id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("id", 2));
  BOOST_CHECK_NO_THROW(mmcif_p.OnBeginLoop(mmcif_h));
  BOOST_MESSAGE("          done.");
  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_parse_models)
{
  BOOST_MESSAGE("  Running mmcif_parse_models tests...");
  IOProfile profile;

  // positive w models
  BOOST_MESSAGE("          true positive test for models...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifParser mmcif_p("testfiles/mmcif/model_truepos.mmcif", eh, profile);
    BOOST_CHECK_NO_THROW(mmcif_p.Parse());
    BOOST_REQUIRE_EQUAL(eh.GetChainCount(),   2);
    BOOST_REQUIRE_EQUAL(eh.GetResidueCount(), 2);
    BOOST_REQUIRE_EQUAL(eh.GetAtomCount(),   26);
  }
  BOOST_MESSAGE("          done.");

  // positive wo models atom_site.mmcif
  BOOST_MESSAGE("          test absent atom_site.pdbx_PDB_model_num entry...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifParser mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);
    BOOST_CHECK_NO_THROW(mmcif_p.Parse());
  }
  BOOST_MESSAGE("          done.");
  // negative, more than 1 atom_site category
  BOOST_MESSAGE("          testing more than one atom_site block...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifParser mmcif_p("testfiles/mmcif/model_multi_atom_site.mmcif", eh,
                        profile);
    BOOST_CHECK_THROW(mmcif_p.Parse(), IOException);
  }
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifParser mmcif_p("testfiles/mmcif/model_multi_atom_site_inverted.mmcif",
                        eh, profile);
    BOOST_CHECK_THROW(mmcif_p.Parse(), IOException);
  }
  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("          testing single model with model no. entry...");
  {
    // build dummy header
    mol::EntityHandle eh = mol::CreateEntity();
    StarLoopDesc tmmcif_h;
    TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
    SetAtomSiteHeader(&tmmcif_h);
    tmmcif_p.OnBeginLoop(tmmcif_h);
    tmmcif_h.Add(StringRef("pdbx_PDB_model_num", 18));
    BOOST_CHECK_THROW(tmmcif_p.OnBeginLoop(tmmcif_h), IOException);
  }
  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_changing_label_entity_id)
{
  BOOST_MESSAGE("  Running mmcif_changing_label_entity_id tests...");
  IOProfile profile;

  // positive
  BOOST_MESSAGE("          true positive test...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifParser mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);
    BOOST_CHECK_NO_THROW(mmcif_p.Parse());
  }
  BOOST_MESSAGE("          done.");

  // negative
  BOOST_MESSAGE("          true negative test...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifParser mmcif_p("testfiles/mmcif/changing_label_entity_id.mmcif", eh,
                        profile);
    BOOST_CHECK_THROW(mmcif_p.Parse(), IOException);
  }
  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_unknown_entity_type)
{
  BOOST_MESSAGE("  Running mmcif_unknown_entity_type tests...");

  mol::EntityHandle eh = mol::CreateEntity();
  std::vector<StringRef> columns;
  TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
  StarLoopDesc tmmcif_h;

  // build dummy header
  tmmcif_h.SetCategory(StringRef("entity", 6));
  tmmcif_h.Add(StringRef("id", 2));
  tmmcif_h.Add(StringRef("type", 4));
  tmmcif_p.OnBeginLoop(tmmcif_h);

  // positive
  BOOST_MESSAGE("          known type...");
  // build datarow
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("polymer", 7));
  BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntity(columns));
  columns.pop_back();
  columns.push_back(StringRef("non-polymer", 11));
  BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntity(columns));
  columns.pop_back();
  columns.push_back(StringRef("water", 5));
  BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntity(columns));
  BOOST_MESSAGE("          done.");

  // negative
  BOOST_MESSAGE("          unknown type...");
  columns.pop_back();
  columns.push_back(StringRef("foo", 3));
  BOOST_CHECK_THROW(tmmcif_p.ParseEntity(columns), std::runtime_error);
  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_entity_tests)
{
  BOOST_MESSAGE("  Running mmcif_entity_tests...");
  mol::ChainHandle ch;
  IOProfile profile;

  // positive
  BOOST_MESSAGE("          fetching chain type...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifParser mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);
    mmcif_p.Parse();
    ch = eh.FindChain("A");
    BOOST_CHECK(ch.IsValid());
    BOOST_CHECK(ch.GetType() == mol::CHAINTYPE_POLY_PEPTIDE_L);
    ch = eh.FindChain("C");
    BOOST_CHECK(ch.IsValid());
    BOOST_CHECK(ch.GetType() == mol::CHAINTYPE_POLY_PEPTIDE_L);
    ch = eh.FindChain("O");
    BOOST_CHECK(ch.IsValid());
    BOOST_CHECK(ch.GetType() == mol::CHAINTYPE_WATER);
  }
  BOOST_MESSAGE("          done.");
  // negative: no entity description
  BOOST_MESSAGE("          check missing entity description...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifParser mmcif_p("testfiles/mmcif/model_truepos.mmcif",
                        eh,
                        profile);
    mmcif_p.Parse();
    ch = eh.FindChain("A");
    BOOST_CHECK(ch.IsValid());
    BOOST_CHECK(ch.GetType() == mol::CHAINTYPE_UNKNOWN);
    ch = eh.FindChain("B");
    BOOST_CHECK(ch.IsValid());
    BOOST_CHECK(ch.GetType() == mol::CHAINTYPE_UNKNOWN);
  }
  BOOST_MESSAGE("          done.");
  BOOST_MESSAGE("          fetch details...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifParser mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);
    mmcif_p.Parse();
    ch = eh.FindChain("A");
    BOOST_CHECK(ch.IsValid());
    BOOST_CHECK(ch.GetDescription() == "Very important information.");
  }
  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_entity_poly_tests)
{
  conop::Conopology::Instance().SetDefaultBuilder("RBB");
  BOOST_MESSAGE("  Running mmcif_entity_poly_tests...");
  mol::ChainHandle ch;
  IOProfile profile;
  StarLoopDesc tmmcif_h;

  mol::EntityHandle eh = mol::CreateEntity();
  MMCifParser mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);

  mmcif_p.SetReadSeqRes(true);
  mmcif_p.Parse();

  seq::SequenceList seqres = mmcif_p.GetSeqRes();
  seq::SequenceHandle curr = seqres.FindSequence("A");
  BOOST_CHECK(curr.GetString() == "VTI");

  BOOST_MESSAGE("          testing missing corresponding entity entry...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    std::vector<StringRef> columns;
    TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);

    tmmcif_h.SetCategory(StringRef("entity_poly", 11));
    tmmcif_h.Add(StringRef("entity_id", 9));
    tmmcif_p.OnBeginLoop(tmmcif_h);

    columns.push_back(StringRef("1", 1));
    BOOST_CHECK_THROW(tmmcif_p.ParseEntityPoly(columns), IOException);
  }
  BOOST_MESSAGE("          done.");
  BOOST_MESSAGE("          testing type recognition...");
  {
    TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
    std::vector<StringRef> columns;

    // create corresponding entity entry
    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("entity", 6));
    tmmcif_h.Add(StringRef("id", 2));
    tmmcif_h.Add(StringRef("type", 4));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("polymer", 7));
    tmmcif_p.ParseEntity(columns);
    columns.pop_back();
    columns.pop_back();

    // build dummy entity_poly header
    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("entity_poly", 11));
    tmmcif_h.Add(StringRef("entity_id", 9));
    tmmcif_h.Add(StringRef("type", 4));
    tmmcif_p.OnBeginLoop(tmmcif_h);

    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("polypeptide(D)", 14));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    columns.pop_back();
    columns.push_back(StringRef("polypeptide(L)", 14));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    columns.pop_back();
    columns.push_back(StringRef("polydeoxyribonucleotide", 23));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    columns.pop_back();
    columns.push_back(StringRef("polyribonucleotide", 18));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    columns.pop_back();
    columns.push_back(StringRef("polysaccharide(D)", 17));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    columns.pop_back();
    columns.push_back(StringRef("polysaccharide(L)", 17));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    columns.pop_back();
columns.push_back(StringRef("polydeoxyribonucleotide/polyribonucleotide hybrid",
                                49));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    columns.pop_back();
    columns.push_back(StringRef("other", 5));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    columns.pop_back();
    columns.pop_back();
    columns.push_back(StringRef("badbadprion", 11));
    BOOST_CHECK_THROW(tmmcif_p.ParseEntityPoly(columns), IOException);
    columns.pop_back();
  }
  BOOST_MESSAGE("          done.");
  BOOST_MESSAGE("          testing pdbx_seq_one_letter_code reading...");
  {
    TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
    std::vector<StringRef> columns;

    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("entity", 6));
    tmmcif_h.Add(StringRef("id", 2));
    tmmcif_h.Add(StringRef("type", 4));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("polymer", 7));
    tmmcif_p.ParseEntity(columns);
    columns.pop_back();
    columns.pop_back();

    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("entity_poly", 11));
    tmmcif_h.Add(StringRef("entity_id", 9));
    tmmcif_h.Add(StringRef("type", 4));
    tmmcif_h.Add(StringRef("pdbx_seq_one_letter_code", 24));
    tmmcif_p.OnBeginLoop(tmmcif_h);

    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("other", 5));
    columns.push_back(StringRef("ABRND", 5));
    tmmcif_p.SetReadSeqRes(true);
    tmmcif_p.SetReadCanonicalSeqRes(true);
    BOOST_CHECK_THROW(tmmcif_p.ParseEntityPoly(columns), IOException);
    tmmcif_p.SetReadCanonicalSeqRes(false);
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    BOOST_CHECK_THROW(tmmcif_p.ParseEntityPoly(columns), IOException);
  }
  BOOST_MESSAGE("          done.");
  BOOST_MESSAGE("          testing pdbx_seq_one_letter_code_can reading...");
  {
    TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
    std::vector<StringRef> columns;

    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("entity", 6));
    tmmcif_h.Add(StringRef("id", 2));
    tmmcif_h.Add(StringRef("type", 4));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("polymer", 7));
    tmmcif_p.ParseEntity(columns);
    columns.pop_back();
    columns.pop_back();

    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("entity_poly", 11));
    tmmcif_h.Add(StringRef("entity_id", 9));
    tmmcif_h.Add(StringRef("type", 4));
    tmmcif_h.Add(StringRef("pdbx_seq_one_letter_code_can", 28));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    tmmcif_p.SetReadCanonicalSeqRes(false);
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("other", 5));
    columns.push_back(StringRef("ABRND", 5));
    tmmcif_p.SetReadSeqRes(true);
    BOOST_CHECK_THROW(tmmcif_p.ParseEntityPoly(columns), IOException);
    tmmcif_p.SetReadCanonicalSeqRes(true);
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntityPoly(columns));
    BOOST_CHECK_THROW(tmmcif_p.ParseEntityPoly(columns), IOException);
  }
  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("  done.");
  conop::Conopology::Instance().SetDefaultBuilder("HEURISTIC");  
}

BOOST_AUTO_TEST_CASE(mmcif_citation_tests)
{
  BOOST_MESSAGE("  Running mmcif_citation_tests...");
  //build dummy citation
  mol::EntityHandle eh;
  TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
  StarLoopDesc tmmcif_h;
  std::vector<StringRef> columns;

  tmmcif_h.SetCategory(StringRef("citation", 8));
  tmmcif_h.Add(StringRef("id", 2));
  tmmcif_h.Add(StringRef("book_title", 10));
  tmmcif_h.Add(StringRef("journal_abbrev", 14));
  tmmcif_p.OnBeginLoop(tmmcif_h);

  columns.push_back(StringRef("Foo", 3));
  columns.push_back(StringRef("The Guide", 9));
  columns.push_back(StringRef(".", 1));

  BOOST_CHECK_NO_THROW(tmmcif_p.ParseCitation(columns));

  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef(".", 1));
  columns.push_back(StringRef("Hitch", 5));

  BOOST_CHECK_NO_THROW(tmmcif_p.ParseCitation(columns));

  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef("The Guide", 9));
  columns.push_back(StringRef("Hitch", 5));

  BOOST_CHECK_THROW(tmmcif_p.ParseCitation(columns), IOException);

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_citation_author_tests)
{
  BOOST_MESSAGE("  Running mmcif_citation_author_tests...");

  mol::EntityHandle eh = mol::CreateEntity();
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  IOProfile profile;
  MMCifParser mmcif_p(s, eh, profile);
  BOOST_CHECK_NO_THROW(mmcif_p.Parse());

  std::vector<String> authors =
    mmcif_p.GetInfo().GetCitations().back().GetAuthorList();

  BOOST_CHECK(authors.size() == 3);
  BOOST_CHECK(authors[0] == "Whiskers, P.D.");
  BOOST_CHECK(authors[1] == "McCheese, B.M.");
  BOOST_CHECK(authors[2] == "Van Hummel, J.F.");

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_refine_tests)
{
  BOOST_MESSAGE("  Running mmcif_refine_tests...");
  BOOST_MESSAGE("         positive test...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    std::ifstream s("testfiles/mmcif/atom_site.mmcif");
    IOProfile profile;
    MMCifParser mmcif_p(s, eh, profile);
    BOOST_CHECK_NO_THROW(mmcif_p.Parse());
    BOOST_CHECK_CLOSE(mmcif_p.GetInfo().GetResolution(), 2.0f, 0.001f);
  }
  BOOST_MESSAGE("         done.");
  BOOST_MESSAGE("         capturing fishy data lines...");
  {
    mol::EntityHandle eh;
    TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
    StarLoopDesc tmmcif_h;
    std::vector<StringRef> columns;
    
    tmmcif_h.SetCategory(StringRef("refine", 6));
    tmmcif_h.Add(StringRef("entry_id", 8));
    tmmcif_h.Add(StringRef("ls_d_res_high", 13));
    tmmcif_h.Add(StringRef("ls_d_res_low", 12));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    
    columns.push_back(StringRef("1Foo", 4));
    columns.push_back(StringRef("Foo", 3));
    columns.push_back(StringRef("1", 1));
    
    BOOST_CHECK_THROW(tmmcif_p.ParseRefine(columns), IOException);
  }
  BOOST_MESSAGE("         done.");
  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_biounit_tests)
{
  BOOST_MESSAGE("  Running mmcif_biounit_tests...");
  //build dummy biounit
  mol::EntityHandle eh = mol::CreateEntity();
  TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
  StarLoopDesc tmmcif_h;
  std::vector<StringRef> columns;

  tmmcif_h.SetCategory(StringRef("pdbx_struct_assembly_gen", 24));
  tmmcif_h.Add(StringRef("assembly_id", 11));
  tmmcif_h.Add(StringRef("asym_id_list", 12));
  tmmcif_h.Add(StringRef("oper_expression", 15));
  tmmcif_p.OnBeginLoop(tmmcif_h);

  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("A", 1));
  columns.push_back(StringRef("1", 1));

  BOOST_CHECK_NO_THROW(tmmcif_p.ParsePdbxStructAssemblyGen(columns));
  BOOST_CHECK_THROW(tmmcif_p.OnEndData(), IOException);

  tmmcif_h.Clear();
  tmmcif_h.SetCategory(StringRef("pdbx_struct_assembly", 20));
  tmmcif_h.Add(StringRef("id", 2));
  columns.pop_back();
  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef("1", 1));
  tmmcif_p.OnBeginLoop(tmmcif_h);
  tmmcif_p.ParsePdbxStructAssembly(columns);

  tmmcif_h.Clear();
  tmmcif_h.SetCategory(StringRef("pdbx_struct_assembly_gen", 24));
  tmmcif_h.Add(StringRef("assembly_id", 11));
  tmmcif_h.Add(StringRef("asym_id_list", 12));
  tmmcif_h.Add(StringRef("oper_expression", 15));
  tmmcif_p.OnBeginLoop(tmmcif_h);

  columns.pop_back();
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("A", 1));
  columns.push_back(StringRef("1", 1));

  BOOST_CHECK_NO_THROW(tmmcif_p.ParsePdbxStructAssemblyGen(columns));
  BOOST_CHECK_THROW(tmmcif_p.OnEndData(), IOException);

  tmmcif_h.Clear();
  tmmcif_h.SetCategory(StringRef("pdbx_struct_assembly_gen", 24));
  tmmcif_h.Add(StringRef("assembly_id", 11));
  tmmcif_h.Add(StringRef("asym_id_list", 12));
  tmmcif_h.Add(StringRef("oper_expression", 15));
  tmmcif_p.OnBeginLoop(tmmcif_h);

  columns.pop_back();
  columns.push_back(StringRef("1-Z", 3));
  BOOST_CHECK_THROW(tmmcif_p.ParsePdbxStructAssemblyGen(columns), IOException);

  columns.pop_back();
  columns.push_back(StringRef("--", 3));
  BOOST_CHECK_THROW(tmmcif_p.ParsePdbxStructAssemblyGen(columns), IOException);

  columns.pop_back();
  columns.push_back(StringRef("A-3", 3));
  BOOST_CHECK_THROW(tmmcif_p.ParsePdbxStructAssemblyGen(columns), IOException);

  tmmcif_h.Clear();
  tmmcif_h.SetCategory(StringRef("pdbx_struct_oper_list", 21));
  tmmcif_h.Add(StringRef("id", 2));
  tmmcif_h.Add(StringRef("type", 4));
  tmmcif_h.Add(StringRef("vector[1]", 9));
  tmmcif_h.Add(StringRef("vector[2]", 9));
  tmmcif_h.Add(StringRef("vector[3]", 9));
  tmmcif_h.Add(StringRef("matrix[1][1]", 12));
  tmmcif_h.Add(StringRef("matrix[1][2]", 12));
  tmmcif_h.Add(StringRef("matrix[1][3]", 12));
  tmmcif_h.Add(StringRef("matrix[2][1]", 12));
  tmmcif_h.Add(StringRef("matrix[2][2]", 12));
  tmmcif_h.Add(StringRef("matrix[2][3]", 12));
  tmmcif_h.Add(StringRef("matrix[3][1]", 12));
  tmmcif_h.Add(StringRef("matrix[3][2]", 12));
  tmmcif_h.Add(StringRef("matrix[3][3]", 12));

  tmmcif_p.OnBeginLoop(tmmcif_h);

  columns.pop_back();
  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("Nan", 3));
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("A", 1));
  columns.push_back(StringRef("3", 1));
  BOOST_CHECK_THROW(tmmcif_p.ParsePdbxStructOperList(columns), IOException);

  columns.pop_back();
  columns.pop_back();
  columns.pop_back();
  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("Nan", 3));
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("2", 1));
  columns.push_back(StringRef("3", 1));
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("2", 1));
  columns.push_back(StringRef("3", 1));
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("2", 1));
  columns.push_back(StringRef("3", 1));
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("A", 1));
  columns.push_back(StringRef("3", 1));
  BOOST_CHECK_THROW(tmmcif_p.ParsePdbxStructOperList(columns), IOException);

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_struct_tests)
{
  BOOST_MESSAGE("  Running mmcif_struct_tests...");

  mol::EntityHandle eh = mol::CreateEntity();
  TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
  StarLoopDesc tmmcif_h;
  std::vector<StringRef> columns;

  tmmcif_h.SetCategory(StringRef("struct", 6));
  tmmcif_h.Add(StringRef("entry_id", 8));
  tmmcif_h.Add(StringRef("pdbx_CASP_flag", 14));
  tmmcif_h.Add(StringRef("pdbx_model_details", 18));
  tmmcif_h.Add(StringRef("pdbx_model_type_details", 23));
  tmmcif_h.Add(StringRef("pdbx_formula_weight", 19));
  tmmcif_p.OnBeginLoop(tmmcif_h);

  columns.push_back(StringRef("1BAR", 4));
  columns.push_back(StringRef("?", 1));
  columns.push_back(StringRef("?", 1));
  columns.push_back(StringRef("?", 1));
  columns.push_back(StringRef("1.0", 3));
  BOOST_CHECK_NO_THROW(tmmcif_p.ParseStruct(columns));

  MMCifInfoStructDetails sd = MMCifInfoStructDetails();
  sd = tmmcif_p.GetInfo().GetStructDetails();

  BOOST_CHECK(sd.GetCASPFlag() == '\0');
  BOOST_CHECK(sd.GetModelDetails() == "");
  BOOST_CHECK(sd.GetModelTypeDetails() == "");

  columns.pop_back();
  columns.push_back(StringRef("A", 1));
  BOOST_CHECK_THROW(tmmcif_p.ParseStruct(columns), IOException);

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_struct_conf_tests)
{
  BOOST_MESSAGE("  Running mmcif_struct_conf_tests...");
  mol::EntityHandle eh = mol::CreateEntity();
  TestMMCifParserProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);

  BOOST_MESSAGE("          testing type validation");
  StringRef type = StringRef("HELX_P", 6);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_OT_P", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_P", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_OT_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_AL_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_GA_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_OM_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_PI_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_27_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_3T_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_PP_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_P",     9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_OT_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_AL_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_GA_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_OM_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_PI_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_27_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
 type = StringRef("HELX_LH_3T_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_PP_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_N", 6);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_OT_N", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_N", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_OT_N", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_A_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_B_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_Z_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_N", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_OT_N", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_A_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_B_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_Z_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_HELIX);
  type = StringRef("TURN_P", 6);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_TURN);
  type = StringRef("TURN_OT_P", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_TURN);
  type = StringRef("TURN_TY1_P", 10);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_TURN);
  type = StringRef("TURN_TY1P_P", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_TURN);
  type = StringRef("TURN_TY2_P", 10);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_TURN);
  type = StringRef("TURN_TY2P_P", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_TURN);
  type = StringRef("TURN_TY3_P", 10);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_TURN);
  type = StringRef("TURN_TY3P_P", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_TURN);
  type = StringRef("STRN", 4);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifParserProtected::MMCIF_STRAND);
  type = StringRef("Foo", 3);
  BOOST_CHECK_THROW(tmmcif_p.DetermineSecStructType(type), IOException);

  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_parseatomident)
{
  BOOST_MESSAGE("  Running mmcif_parseatomident tests...");

  mol::EntityHandle eh = mol::CreateEntity();

  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  IOProfile profile;
  StarLoopDesc tmmcif_h;
  TestMMCifParserProtected tmmcif_p(s, eh, profile);
  std::vector<StringRef> columns;
  String chain_name;
  StringRef res_name;
  mol::ResNum resnum(0);
  bool valid_res_num = false;
  char alt_loc;
  StringRef atom_name;

  BOOST_MESSAGE("          testing valid line");
  //tmmcif_p.ParseAtomIdent();
  BOOST_MESSAGE("          done.");
  // negative
  //cols.push_back(StringRef("ATOM", 4));
  //BOOST_CHECK_THROW(tmmcif_p.ParseAtomIdent(cols,
  //                                          chain_name,
  //                                          res_name,
  //                                          resnum,
  //                                          atom_name,
  //                                          alt_loc), IOException);
  // positive
  //StarLoopDesc tmmcif_h;
  //tmmcif_h.SetCategory(StringRef("atom_site", 9));
  // build header
  //mmcif_h.Add(StringRef("AUTH_ASYM_ID", 12));
  /*
    this->TryStoreIdx(AUTH_ASYM_ID,    "auth_asym_id",    header);
    this->TryStoreIdx(ID,              "id",              header);
    this->TryStoreIdx(LABEL_ALT_ID,    "label_alt_id",    header);
    this->TryStoreIdx(LABEL_ASYM_ID,   "label_asym_id",   header);
    this->TryStoreIdx(LABEL_ATOM_ID,   "label_atom_id",   header);
    this->TryStoreIdx(LABEL_COMP_ID,   "label_comp_id",   header);
    this->TryStoreIdx(LABEL_ENTITY_ID, "label_entity_id", header);
    this->TryStoreIdx(LABEL_SEQ_ID,    "label_seq_id",    header);
    this->TryStoreIdx(TYPE_SYMBOL,     "type_symbol",     header);
    this->TryStoreIdx(CARTN_X, "Cartn_x", header);
    this->TryStoreIdx(CARTN_Y, "Cartn_y", header);
    this->TryStoreIdx(CARTN_Z, "Cartn_z", header);
*/

  BOOST_MESSAGE("          testing profile to read calpha only");
  {
    profile.calpha_only = true;
    // set up dummy header to pre-set indices
    SetAtomSiteHeader(&tmmcif_h);
    tmmcif_p.OnBeginLoop(tmmcif_h);
    // create CA dummy line
    columns.push_back(StringRef("A", 1));
    columns.push_back(StringRef("2", 1));
    columns.push_back(StringRef(".", 1));
    columns.push_back(StringRef("A", 1));
    columns.push_back(StringRef("CA", 2));
    columns.push_back(StringRef("VAL", 3));
    columns.push_back(StringRef("1", 1));      // label_entity_id
    columns.push_back(StringRef("11", 2));     // label_seq_id
    columns.push_back(StringRef("C", 1));      // type_symbol
    columns.push_back(StringRef("25.369", 6)); // Cartn_x
    columns.push_back(StringRef("30.691", 6)); // Cartn_y
    columns.push_back(StringRef("11.795", 6)); // Cartn_z
    BOOST_CHECK_EQUAL(tmmcif_p.ParseAtomIdent(columns, chain_name, res_name,
                                              resnum, valid_res_num, atom_name,
                                              alt_loc), true);
    columns.pop_back();
    columns.pop_back();
    columns.pop_back();
    columns.pop_back();
    columns.pop_back();
    columns.pop_back();
    columns.pop_back();
    columns.pop_back();
    columns.push_back(StringRef("CB", 2));
    columns.push_back(StringRef("VAL", 3));
    columns.push_back(StringRef("1", 1));      // label_entity_id
    columns.push_back(StringRef("11", 2));     // label_seq_id
    columns.push_back(StringRef("C", 1));      // type_symbol
    columns.push_back(StringRef("25.369", 6)); // Cartn_x
    columns.push_back(StringRef("30.691", 6)); // Cartn_y
    columns.push_back(StringRef("11.795", 6)); // Cartn_z
    BOOST_CHECK_EQUAL(tmmcif_p.ParseAtomIdent(columns, chain_name, res_name,
                                              resnum, valid_res_num, atom_name,
                                              alt_loc), false);
  }
  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_parseandaddatom)
{
  mol::EntityHandle eh = mol::CreateEntity();

  BOOST_MESSAGE("  Running mmcif_parseandaddatom tests...");
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  TestMMCifParserProtected tmmcif_p(s, eh, IOProfile());
  std::vector<StringRef> cols;

  //BOOST_MESSAGE("    testing short atom_site entry");
  //cols.push_back(StringRef("ATOM", 4));
  //BOOST_CHECK_THROW(tmmcif_p.ParseAndAddAtom(cols), IOException);
  //BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_testreader)
{
  BOOST_MESSAGE("  Running mmcif_testreader tests...");
  mol::EntityHandle eh = mol::CreateEntity();
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  IOProfile profile;
  MMCifParser mmcif_p(s, eh, profile);

  BOOST_MESSAGE("          testing Parse()...");
  BOOST_CHECK_NO_THROW(mmcif_p.Parse());

  BOOST_REQUIRE_EQUAL(eh.GetChainCount(),    3);
  BOOST_REQUIRE_EQUAL(eh.GetResidueCount(), 14);
  BOOST_REQUIRE_EQUAL(eh.GetAtomCount(),    35);

  mol::ChainHandle ch = eh.FindChain("A");
  BOOST_CHECK(ch.IsValid());
  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("          testing numbering water...");
  ch = eh.FindChain("O");
  BOOST_CHECK(ch.IsValid());
  mol::ResidueHandleList rl = ch.GetResidueList();
  mol::ResidueHandleList::const_iterator rs;
  int i = 1;
  for (rs = rl.begin(); rs != rl.end(); ++rs, ++i) {
    BOOST_CHECK_EQUAL(rs->GetNumber().GetNum(), i);
  }

  // add checking of struct_conf info, here

  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("          reading data fields which should not fail...");
  BOOST_CHECK(mmcif_p.GetInfo().GetMethod().str() == "Deep-fry");
  BOOST_CHECK(mmcif_p.GetInfo().GetBioUnits().back().GetDetails() ==
              "author_defined_assembly");
  BOOST_CHECK(mmcif_p.GetInfo().GetBioUnits().back().GetChainList().back() ==
              "F");
  MMCifInfoBioUnit bu = mmcif_p.GetInfo().GetBioUnits().back();
  BOOST_CHECK(bu.GetOperations().back().back()->GetType() ==
              "identity operation");
  MMCifInfoStructDetails sd = mmcif_p.GetInfo().GetStructDetails();
  BOOST_CHECK(sd.GetEntryID() == "1BAR");
  BOOST_CHECK(sd.GetTitle() == "A Title");
  BOOST_CHECK(sd.GetCASPFlag() == 'Y');
  BOOST_CHECK(sd.GetDescriptor() == "ADENYLATE KINASE");
  BOOST_CHECK_CLOSE(sd.GetMass(), 1.0f, 0.001f);
  BOOST_CHECK(sd.GetMassMethod() == "Good Guess");
  BOOST_CHECK(sd.GetModelDetails() == "Even better guessing");
  BOOST_CHECK(sd.GetModelTypeDetails() == "Guess");
  BOOST_MESSAGE("          done.");

  BOOST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_SUITE_END();
