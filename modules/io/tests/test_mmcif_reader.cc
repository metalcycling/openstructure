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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <ost/platform.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/mol/mmcif_reader.hh>
#include <ost/conop/conop.hh>



using namespace ost;
using namespace ost::io;

class TestMMCifReaderProtected : MMCifReader {
public:
  TestMMCifReaderProtected(std::istream& stream, mol::EntityHandle& ent_handle):
    MMCifReader(stream, ent_handle, IOProfile())
  { }

  TestMMCifReaderProtected(std::istream& stream, mol::EntityHandle& ent_handle,
                           const IOProfile& profile):
    MMCifReader(stream, ent_handle, profile)
  { }

  TestMMCifReaderProtected(const String& filename,
                           mol::EntityHandle& ent_handle):
    MMCifReader(filename, ent_handle, IOProfile())
  { }

  using MMCifReader::OnBeginLoop;
  using MMCifReader::OnEndData;
  using MMCifReader::ParseAtomIdent;
  using MMCifReader::ParseAndAddAtom;
  using MMCifReader::ParseEntity;
  using MMCifReader::ParseEntityPoly;
  using MMCifReader::ParseCitation;
  using MMCifReader::ParseRefine;
  using MMCifReader::ParsePdbxStructAssemblyGen;
  using MMCifReader::ParsePdbxStructAssembly;
  using MMCifReader::ParsePdbxStructOperList;
  using MMCifReader::ParseStruct;
  using MMCifReader::ParseStructConf;
  using MMCifReader::ParseStructSheetRange;
  using MMCifReader::ParsePdbxEntityBranch;
  using MMCifReader::ParsePdbxEntityBranchLink;
  using MMCifReader::TryStoreIdx;
  using MMCifReader::SetRestrictChains;
  using MMCifReader::SetReadSeqRes;
  using MMCifReader::SetReadCanonicalSeqRes;
  using MMCifReader::ClearState;
  using MMCifReader::ConvertSEQRES;
  using MMCifReader::GetInfo;
  using MMCifReader::DetermineSecStructType;
  using MMCifReader::MMCifSecStructElement;
  using MMCifReader::MMCIF_HELIX;
  using MMCifReader::MMCIF_TURN;
  using MMCifReader::MMCIF_STRAND;
  using MMCifReader::SetAuthChainID;
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

conop::CompoundLibPtr SetDefaultCompoundLib() {
  // return NULL if not successful, else return newly set default lib
  // REQ: OST_ROOT to be set
  char * ost_root = getenv("OST_ROOT");
  if (!ost_root) return conop::CompoundLibPtr();
  SetPrefixPath(ost_root);
  String lib_path = GetSharedDataPath() + "/compounds.chemlib";
  conop::CompoundLibPtr compound_lib = conop::CompoundLib::Load(lib_path);
  if (compound_lib) {
    conop::Conopology::Instance().SetDefaultLib(compound_lib);
  }
  return compound_lib;
}

BOOST_AUTO_TEST_SUITE( io );

BOOST_AUTO_TEST_CASE(mmcif_trystoreidx)
{
  mol::EntityHandle eh = mol::CreateEntity();

  BOOST_TEST_MESSAGE("  Running mmcif_trystoreidx tests...");
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  TestMMCifReaderProtected tmmcif_p(s, eh, IOProfile());
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
  conop::CompoundLibPtr compound_lib = SetDefaultCompoundLib();
  if (!compound_lib) {
    std::cout << "WARNING: skipping mmcif_convert_seqres unit test. " 
              << "Compound library is required" << std::endl;
    return;
  }

  mol::EntityHandle eh=mol::CreateEntity();
  
  TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
  BOOST_CHECK_EQUAL(tmmcif_p.ConvertSEQRES("A(MSE)Y", compound_lib), "AMY");
  BOOST_CHECK_THROW(tmmcif_p.ConvertSEQRES("A(MSEY", compound_lib), 
                    IOException);
}

BOOST_AUTO_TEST_CASE(mmcif_onbeginloop)
{
  mol::EntityHandle eh=mol::CreateEntity();

  // add more tests on new mandatory items
  BOOST_TEST_MESSAGE("  Running mmcif_onbeginloop tests...");
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  MMCifReader mmcif_p(s, eh, IOProfile());
  StarLoopDesc mmcif_h;
  BOOST_TEST_MESSAGE("          testing atom_site items...");
  mmcif_h.SetCategory(StringRef("atom_site", 9));
  BOOST_TEST_MESSAGE("             auth_asym_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("auth_asym_id", 12));
  BOOST_TEST_MESSAGE("             id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("id", 2));
  BOOST_TEST_MESSAGE("             label_alt_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_alt_id", 12));
  BOOST_TEST_MESSAGE("             label_asym_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_asym_id", 13));
  BOOST_TEST_MESSAGE("             label_atom_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_atom_id", 13));
  BOOST_TEST_MESSAGE("             label_comp_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_comp_id", 13));
  BOOST_TEST_MESSAGE("             label_entity_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_entity_id", 15));
  BOOST_TEST_MESSAGE("             label_seq_id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("label_seq_id", 12));
  BOOST_TEST_MESSAGE("             type_symbol");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("type_symbol", 11));
  BOOST_TEST_MESSAGE("             Cartn_x");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("Cartn_x", 7));
  BOOST_TEST_MESSAGE("             Cartn_y");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("Cartn_y", 7));
  BOOST_TEST_MESSAGE("             Cartn_z");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("Cartn_z", 7));
  BOOST_CHECK_NO_THROW(mmcif_p.OnBeginLoop(mmcif_h));
  BOOST_TEST_MESSAGE("          done.");
  mmcif_h.Clear();
  BOOST_TEST_MESSAGE("          testing entity items...");
  mmcif_h.SetCategory(StringRef("entity", 6));
  BOOST_TEST_MESSAGE("             id");
  BOOST_CHECK_THROW(mmcif_p.OnBeginLoop(mmcif_h), IOException);
  mmcif_h.Add(StringRef("id", 2));
  BOOST_CHECK_NO_THROW(mmcif_p.OnBeginLoop(mmcif_h));
  BOOST_TEST_MESSAGE("          done.");
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_parse_models)
{
  BOOST_TEST_MESSAGE("  Running mmcif_parse_models tests...");
  IOProfile profile;

  // positive w models
  BOOST_TEST_MESSAGE("          true positive test for models...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifReader mmcif_p("testfiles/mmcif/model_truepos.mmcif", eh, profile);
    BOOST_CHECK_NO_THROW(mmcif_p.Parse());
    BOOST_REQUIRE_EQUAL(eh.GetChainCount(),   2);
    BOOST_REQUIRE_EQUAL(eh.GetResidueCount(), 2);
    BOOST_REQUIRE_EQUAL(eh.GetAtomCount(),   26);
  }
  BOOST_TEST_MESSAGE("          done.");

  // positive wo models atom_site.mmcif
  BOOST_TEST_MESSAGE("          test absent atom_site.pdbx_PDB_model_num "
                     "entry...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifReader mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);
    BOOST_CHECK_NO_THROW(mmcif_p.Parse());
  }
  BOOST_TEST_MESSAGE("          done.");
  // negative, more than 1 atom_site category
  BOOST_TEST_MESSAGE("          testing more than one atom_site block...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifReader mmcif_p("testfiles/mmcif/model_multi_atom_site.mmcif", eh,
                        profile);
    BOOST_CHECK_THROW(mmcif_p.Parse(), IOException);
  }
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifReader mmcif_p("testfiles/mmcif/model_multi_atom_site_inverted.mmcif",
                        eh, profile);
    BOOST_CHECK_THROW(mmcif_p.Parse(), IOException);
  }
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("          testing single model with model no. entry...");
  {
    // build dummy header
    mol::EntityHandle eh = mol::CreateEntity();
    StarLoopDesc tmmcif_h;
    TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
    SetAtomSiteHeader(&tmmcif_h);
    tmmcif_p.OnBeginLoop(tmmcif_h);
    tmmcif_h.Add(StringRef("pdbx_PDB_model_num", 18));
    BOOST_CHECK_THROW(tmmcif_p.OnBeginLoop(tmmcif_h), IOException);
  }
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_changing_label_entity_id)
{
  BOOST_TEST_MESSAGE("  Running mmcif_changing_label_entity_id tests...");
  IOProfile profile;

  // positive
  BOOST_TEST_MESSAGE("          true positive test...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifReader mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);
    BOOST_CHECK_NO_THROW(mmcif_p.Parse());
  }
  BOOST_TEST_MESSAGE("          done.");

  // negative
  BOOST_TEST_MESSAGE("          true negative test...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifReader mmcif_p("testfiles/mmcif/changing_label_entity_id.mmcif", eh,
                        profile);
    BOOST_CHECK_THROW(mmcif_p.Parse(), IOException);
  }
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_unknown_entity_type)
{
  BOOST_TEST_MESSAGE("  Running mmcif_unknown_entity_type tests...");

  mol::EntityHandle eh = mol::CreateEntity();
  std::vector<StringRef> columns;
  TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
  StarLoopDesc tmmcif_h;

  // build dummy header
  tmmcif_h.SetCategory(StringRef("entity", 6));
  tmmcif_h.Add(StringRef("id", 2));
  tmmcif_h.Add(StringRef("type", 4));
  tmmcif_p.OnBeginLoop(tmmcif_h);

  // positive
  BOOST_TEST_MESSAGE("          known type...");
  // build datarow
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("polymer", 7));
  BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntity(columns));
  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef("2", 1));
  columns.push_back(StringRef("non-polymer", 11));
  BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntity(columns));
  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef("3", 1));  
  columns.push_back(StringRef("water", 5));
  BOOST_CHECK_NO_THROW(tmmcif_p.ParseEntity(columns));
  BOOST_TEST_MESSAGE("          done.");

  // negative
  BOOST_TEST_MESSAGE("          unknown type...");
  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef("4", 1));
  columns.push_back(StringRef("foo", 3));
  BOOST_CHECK_THROW(tmmcif_p.ParseEntity(columns), Error);
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_entity_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_entity_tests...");
  mol::ChainHandle ch;
  IOProfile profile;

  // positive
  BOOST_TEST_MESSAGE("          fetching chain type...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifReader mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);
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
  BOOST_TEST_MESSAGE("          done.");
  // negative: no entity description
  BOOST_TEST_MESSAGE("          check missing entity description...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifReader mmcif_p("testfiles/mmcif/model_truepos.mmcif",
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
  BOOST_TEST_MESSAGE("          done.");
  BOOST_TEST_MESSAGE("          fetch details...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    MMCifReader mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);
    mmcif_p.Parse();
    ch = eh.FindChain("A");
    BOOST_CHECK(ch.IsValid());
    BOOST_CHECK(ch.GetDescription() == "Very important information.");
  }
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_entity_poly_tests)
{
  if (!SetDefaultCompoundLib()) {
    std::cout << "WARNING: skipping mmcif_entity_poly_tests unit test. " 
              << "Compound library is required" << std::endl;
    return;
  }
  
  BOOST_TEST_MESSAGE("  Running mmcif_entity_poly_tests...");
  mol::ChainHandle ch;
  IOProfile profile;
  StarLoopDesc tmmcif_h;

  mol::EntityHandle eh = mol::CreateEntity();
  MMCifReader mmcif_p("testfiles/mmcif/atom_site.mmcif", eh, profile);

  mmcif_p.SetReadSeqRes(true);
  mmcif_p.Parse();

  seq::SequenceList seqres = mmcif_p.GetSeqRes();
  seq::SequenceHandle curr = seqres.FindSequence("A");
  BOOST_CHECK(curr.GetString() == "VTI");

  BOOST_TEST_MESSAGE("          testing type recognition...");
  {
    TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
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
    columns.push_back(StringRef("badbadprion", 11));
    BOOST_CHECK_THROW(tmmcif_p.ParseEntityPoly(columns), ost::Error);
    columns.pop_back();
  }
  BOOST_TEST_MESSAGE("          done.");
  BOOST_TEST_MESSAGE("          testing pdbx_seq_one_letter_code reading...");
  {
    TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
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
  BOOST_TEST_MESSAGE("          done.");
  BOOST_TEST_MESSAGE("          testing pdbx_seq_one_letter_code_can "
                     "reading...");
  {
    TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
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
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_citation_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_citation_tests...");
  //build dummy citation
  mol::EntityHandle eh;
  TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
  StarLoopDesc tmmcif_h;
  std::vector<StringRef> columns;

  tmmcif_h.SetCategory(StringRef("citation", 8));
  tmmcif_h.Add(StringRef("id", 2));
  tmmcif_h.Add(StringRef("year", 4));
  tmmcif_h.Add(StringRef("book_publisher_city", 19));
  tmmcif_h.Add(StringRef("book_title", 10));
  tmmcif_h.Add(StringRef("book_publisher", 14));
  tmmcif_h.Add(StringRef("journal_abbrev", 14));
  tmmcif_p.OnBeginLoop(tmmcif_h);

  // ensure that we use book_title if no journal given (no RCSB use of this)
  columns.push_back(StringRef("Foo", 3));
  columns.push_back(StringRef("1979", 4));
  columns.push_back(StringRef("The restaurant", 14));
  columns.push_back(StringRef("The Guide", 9));
  columns.push_back(StringRef("Doug", 4));
  columns.push_back(StringRef(".", 1));

  BOOST_CHECK_NO_THROW(tmmcif_p.ParseCitation(columns));
  const MMCifInfoCitation& cit = tmmcif_p.GetInfo().GetCitations().back();
  BOOST_CHECK_EQUAL(cit.GetID(), String("Foo"));
  BOOST_CHECK_EQUAL(cit.GetYear(), 1979);
  BOOST_CHECK_EQUAL(cit.GetPublishedIn(), String("The Guide"));
  BOOST_CHECK_EQUAL(cit.GetBookPublisher(), String("Doug"));
  BOOST_CHECK_EQUAL(cit.GetBookPublisherCity(), String("The restaurant"));
  BOOST_CHECK_EQUAL(cit.IsCitationTypeBook(), true);

  // ensure that we override book_title if not properly given
  columns.pop_back();
  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef(".", 1));
  columns.push_back(StringRef("Doug", 4));
  columns.push_back(StringRef("Hitch", 5));

  BOOST_CHECK_NO_THROW(tmmcif_p.ParseCitation(columns));
  BOOST_CHECK_EQUAL(tmmcif_p.GetInfo().GetCitations().back().GetPublishedIn(),
                    String("Hitch"));

  // ensure that we override journal if book_title given
  columns.pop_back();
  columns.pop_back();
  columns.pop_back();
  columns.push_back(StringRef("The Guide", 9));
  columns.push_back(StringRef("Doug", 4));
  columns.push_back(StringRef("Hitch", 5));

  BOOST_CHECK_NO_THROW(tmmcif_p.ParseCitation(columns));
  BOOST_CHECK_EQUAL(tmmcif_p.GetInfo().GetCitations().back().GetPublishedIn(),
                    String("The Guide"));

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_citation_author_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_citation_author_tests...");

  mol::EntityHandle eh = mol::CreateEntity();
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  IOProfile profile;
  MMCifReader mmcif_p(s, eh, profile);
  BOOST_REQUIRE_NO_THROW(mmcif_p.Parse());

  std::vector<String> authors =
    mmcif_p.GetInfo().GetCitations().back().GetAuthorList();

  BOOST_CHECK(authors.size() == 3);
  BOOST_CHECK(authors[0] == "Whiskers, P.D.");
  BOOST_CHECK(authors[1] == "McCheese, B.M.");
  BOOST_CHECK(authors[2] == "Van Hummel, J.F.");

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_struct_ref)
{
  mol::EntityHandle eh = mol::CreateEntity();
  std::ifstream s("testfiles/mmcif/struct_ref.cif");
  IOProfile profile;
  MMCifReader mmcif_p(s, eh, profile);
  mmcif_p.Parse();
  MMCifInfoStructRefs refs=mmcif_p.GetInfo().GetStructRefs();
  BOOST_CHECK_EQUAL(refs.size(), size_t(1));
  MMCifInfoStructRefPtr sr1=refs[0];
  BOOST_CHECK_EQUAL(sr1->GetDBName(), "UNP");
  BOOST_CHECK_EQUAL(sr1->GetDBID(), "BLA2_BACCE");
  BOOST_CHECK_EQUAL(sr1->GetDBAccess(), "P04190");
  BOOST_CHECK_EQUAL(sr1->GetID(), "1");
  MMCifInfoStructRefSeqs seqs=sr1->GetAlignedSeqs();
  BOOST_CHECK_EQUAL(seqs.size(), size_t(2));
  BOOST_CHECK_EQUAL(seqs[0]->GetID(), "1");
  BOOST_CHECK_EQUAL(seqs[0]->GetChainName(), "A");
  BOOST_CHECK_EQUAL(seqs[0]->GetSeqBegin(), 1);
  BOOST_CHECK_EQUAL(seqs[0]->GetSeqEnd(), 19);
  BOOST_CHECK_EQUAL(seqs[0]->GetDBBegin(), 31);
  BOOST_CHECK_EQUAL(seqs[0]->GetDBEnd(), 49);
  BOOST_CHECK_EQUAL(seqs[1]->GetID(), "13");
  BOOST_CHECK_EQUAL(seqs[1]->GetChainName(), "B");
  BOOST_CHECK_EQUAL(seqs[1]->GetSeqBegin(), 1);
  BOOST_CHECK_EQUAL(seqs[1]->GetSeqEnd(), 19);
  BOOST_CHECK_EQUAL(seqs[1]->GetDBBegin(), 31);
  BOOST_CHECK_EQUAL(seqs[1]->GetDBEnd(), 49);
  MMCifInfoStructRefSeqDifs diffs=seqs[0]->GetDifs();
  BOOST_CHECK_EQUAL(diffs.size(), size_t(1));
  BOOST_CHECK_EQUAL(diffs[0]->GetSeqRNum(), 91);
  BOOST_CHECK_EQUAL(diffs[0]->GetDetails(), "ENGINEERED MUTATION");
}

BOOST_AUTO_TEST_CASE(mmcif_refine_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_refine_tests...");
  BOOST_TEST_MESSAGE("         positive test...");
  {
    mol::EntityHandle eh = mol::CreateEntity();
    std::ifstream s("testfiles/mmcif/atom_site.mmcif");
    IOProfile profile;
    MMCifReader mmcif_p(s, eh, profile);
    BOOST_REQUIRE_NO_THROW(mmcif_p.Parse());
    BOOST_CHECK_CLOSE(mmcif_p.GetInfo().GetResolution(), Real(2), Real(0.001));
    BOOST_CHECK_CLOSE(mmcif_p.GetInfo().GetRFree(), Real(0.229), Real(0.01));
    BOOST_CHECK_CLOSE(mmcif_p.GetInfo().GetRWork(), Real(0.200), Real(0.01));
  }
  BOOST_TEST_MESSAGE("         done.");
  BOOST_TEST_MESSAGE("         capturing fishy data lines...");
  {
    mol::EntityHandle eh;
    TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
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
  BOOST_TEST_MESSAGE("         done.");
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_biounit_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_biounit_tests...");
  //build dummy biounit
  mol::EntityHandle eh = mol::CreateEntity();
  TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
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

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_struct_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_struct_tests...");

  mol::EntityHandle eh = mol::CreateEntity();
  TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
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

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_struct_conf_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_struct_conf_tests...");
  mol::EntityHandle eh = mol::CreateEntity();
  TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);

  BOOST_TEST_MESSAGE("          testing type validation");
  StringRef type = StringRef("HELX_P", 6);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_OT_P", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_P", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_OT_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_AL_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_GA_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_OM_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_PI_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_27_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_3T_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_PP_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_P",     9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_OT_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_AL_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_GA_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_OM_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_PI_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_27_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
 type = StringRef("HELX_LH_3T_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_PP_P", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_N", 6);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_OT_N", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_N", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_OT_N", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_A_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_B_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_RH_Z_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_N", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_OT_N", 12);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_A_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_B_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("HELX_LH_Z_N", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_HELIX);
  type = StringRef("TURN_P", 6);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_TURN);
  type = StringRef("TURN_OT_P", 9);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_TURN);
  type = StringRef("TURN_TY1_P", 10);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_TURN);
  type = StringRef("TURN_TY1P_P", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_TURN);
  type = StringRef("TURN_TY2_P", 10);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_TURN);
  type = StringRef("TURN_TY2P_P", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_TURN);
  type = StringRef("TURN_TY3_P", 10);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_TURN);
  type = StringRef("TURN_TY3P_P", 11);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_TURN);
  type = StringRef("STRN", 4);
  BOOST_CHECK(tmmcif_p.DetermineSecStructType(type) ==
              TestMMCifReaderProtected::MMCIF_STRAND);
  type = StringRef("Foo", 3);
  BOOST_CHECK_THROW(tmmcif_p.DetermineSecStructType(type), IOException);

  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("          testing auth_chain_id switch...");
  StarLoopDesc tmmcif_h;
  {
    std::vector<StringRef> columns;
    tmmcif_h.SetCategory(StringRef("struct_conf", 11));
    tmmcif_h.Add(StringRef("beg_label_asym_id", 17));
    tmmcif_h.Add(StringRef("beg_label_comp_id", 17));
    tmmcif_h.Add(StringRef("beg_label_seq_id", 16));
    tmmcif_h.Add(StringRef("conf_type_id", 12));
    tmmcif_h.Add(StringRef("end_label_asym_id", 17));
    tmmcif_h.Add(StringRef("end_label_comp_id", 17));
    tmmcif_h.Add(StringRef("end_label_seq_id", 16));
    tmmcif_h.Add(StringRef("id", 2));
    tmmcif_h.Add(StringRef("beg_auth_asym_id", 16));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    columns.push_back(StringRef("A", 1));
    columns.push_back(StringRef("ARG", 3));
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("HELX_RH_AL_P", 12));
    columns.push_back(StringRef("A", 1));
    columns.push_back(StringRef("ARG", 3));
    columns.push_back(StringRef("2", 1));
    columns.push_back(StringRef("DHLX1", 5));
    columns.push_back(StringRef("A", 1));
    tmmcif_p.SetAuthChainID(true);
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseStructConf(columns));
    tmmcif_p.SetAuthChainID(false);
    BOOST_CHECK_NO_THROW(tmmcif_p.ParseStructConf(columns));
    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("struct_conf", 11));
    tmmcif_h.Add(StringRef("beg_label_asym_id", 17));
    tmmcif_h.Add(StringRef("beg_label_comp_id", 17));
    tmmcif_h.Add(StringRef("beg_label_seq_id", 16));
    tmmcif_h.Add(StringRef("conf_type_id", 12));
    tmmcif_h.Add(StringRef("end_label_asym_id", 17));
    tmmcif_h.Add(StringRef("end_label_comp_id", 17));
    tmmcif_h.Add(StringRef("end_label_seq_id", 16));
    tmmcif_h.Add(StringRef("id", 2));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    columns.pop_back();
    tmmcif_p.SetAuthChainID(true);
    BOOST_CHECK_THROW(tmmcif_p.ParseStructConf(columns), IOException);
  }
  tmmcif_p.SetAuthChainID(false);
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_struct_sheet_range_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_struct_sheet_range_tests...");
  mol::EntityHandle eh = mol::CreateEntity();
  TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);

  BOOST_TEST_MESSAGE("          testing auth_chain_id switch...");
  StarLoopDesc tmmcif_h;
  std::vector<StringRef> columns;
  tmmcif_h.SetCategory(StringRef("struct_sheet_range", 18));
  tmmcif_h.Add(StringRef("sheet_id", 8));
  tmmcif_h.Add(StringRef("beg_label_asym_id", 17));
  tmmcif_h.Add(StringRef("beg_label_comp_id", 17));
  tmmcif_h.Add(StringRef("beg_label_seq_id", 16));
  tmmcif_h.Add(StringRef("end_label_asym_id", 17));
  tmmcif_h.Add(StringRef("end_label_comp_id", 17));
  tmmcif_h.Add(StringRef("end_label_seq_id", 16));
  tmmcif_h.Add(StringRef("id", 2));
  tmmcif_h.Add(StringRef("beg_auth_asym_id", 16));
  tmmcif_p.OnBeginLoop(tmmcif_h);
  columns.push_back(StringRef("Sheet1", 6));
  columns.push_back(StringRef("A", 1));
  columns.push_back(StringRef("ARG", 3));
  columns.push_back(StringRef("1", 1));
  columns.push_back(StringRef("A", 1));
  columns.push_back(StringRef("ARG", 3));
  columns.push_back(StringRef("2", 1));
  columns.push_back(StringRef("DSTRAND", 7));
  columns.push_back(StringRef("A", 1));
  tmmcif_p.SetAuthChainID(true);
  BOOST_CHECK_NO_THROW(tmmcif_p.ParseStructSheetRange(columns));
  tmmcif_p.SetAuthChainID(false);
  BOOST_CHECK_NO_THROW(tmmcif_p.ParseStructSheetRange(columns));
  tmmcif_h.Clear();
  tmmcif_h.SetCategory(StringRef("struct_sheet_range", 11));
  tmmcif_h.Add(StringRef("sheet_id", 8));
  tmmcif_h.Add(StringRef("beg_label_asym_id", 17));
  tmmcif_h.Add(StringRef("beg_label_comp_id", 17));
  tmmcif_h.Add(StringRef("beg_label_seq_id", 16));
  tmmcif_h.Add(StringRef("end_label_asym_id", 17));
  tmmcif_h.Add(StringRef("end_label_comp_id", 17));
  tmmcif_h.Add(StringRef("end_label_seq_id", 16));
  tmmcif_h.Add(StringRef("id", 2));
  tmmcif_p.OnBeginLoop(tmmcif_h);
  columns.pop_back();
  tmmcif_p.SetAuthChainID(true);
  BOOST_CHECK_THROW(tmmcif_p.ParseStructSheetRange(columns), IOException);
  tmmcif_p.SetAuthChainID(false);
  BOOST_TEST_MESSAGE("          done.");



  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_parseatomident)
{
  BOOST_TEST_MESSAGE("  Running mmcif_parseatomident tests...");

  mol::EntityHandle eh = mol::CreateEntity();

  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  IOProfile profile;
  StarLoopDesc tmmcif_h;
  TestMMCifReaderProtected tmmcif_p(s, eh, profile);
  std::vector<StringRef> columns;
  String cif_chain_name;
  String auth_chain_name;
  StringRef res_name;
  mol::ResNum resnum(0);
  bool valid_res_num = false;
  char alt_loc;
  StringRef atom_name;

  BOOST_TEST_MESSAGE("          testing valid line");
  BOOST_TEST_MESSAGE("          done.");
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

  BOOST_TEST_MESSAGE("          testing profile to read calpha only");
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
    BOOST_CHECK_EQUAL(tmmcif_p.ParseAtomIdent(columns, auth_chain_name,
    			                                    cif_chain_name, res_name,
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
    BOOST_CHECK_EQUAL(tmmcif_p.ParseAtomIdent(columns, auth_chain_name, 
    			                                    cif_chain_name, res_name,
                                              resnum, valid_res_num, atom_name,
                                              alt_loc), false);
  }
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("  done.");
}

/*
BOOST_AUTO_TEST_CASE(mmcif_parseandaddatom)
{
  mol::EntityHandle eh = mol::CreateEntity();

  BOOST_TEST_MESSAGE("  Running mmcif_parseandaddatom tests...");
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  TestMMCifReaderProtected tmmcif_p(s, eh, IOProfile());
  std::vector<StringRef> cols;

  //BOOST_TEST_MESSAGE("    testing short atom_site entry");
  //cols.push_back(StringRef("ATOM", 4));
  //BOOST_CHECK_THROW(tmmcif_p.ParseAndAddAtom(cols), IOException);
  //BOOST_TEST_MESSAGE("  done.");
}
*/

BOOST_AUTO_TEST_CASE(mmcif_testreader)
{
  BOOST_TEST_MESSAGE("  Running mmcif_testreader tests...");
  mol::EntityHandle eh = mol::CreateEntity();
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  IOProfile profile;
  MMCifReader mmcif_p(s, eh, profile);

  mmcif_p.SetRestrictChains("A O C");

  BOOST_TEST_MESSAGE("          testing Parse()...");
  BOOST_REQUIRE_NO_THROW(mmcif_p.Parse());

  BOOST_REQUIRE_EQUAL(eh.GetChainCount(),    3);
  BOOST_REQUIRE_EQUAL(eh.GetResidueCount(), 14);
  BOOST_REQUIRE_EQUAL(eh.GetAtomCount(),    35);

  mol::ChainHandle ch = eh.FindChain("A");
  BOOST_CHECK(ch.IsValid());
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("          testing numbering water...");
  ch = eh.FindChain("O");
  BOOST_CHECK(ch.IsValid());
  mol::ResidueHandleList rl = ch.GetResidueList();
  mol::ResidueHandleList::const_iterator rs;
  int i = 1;
  for (rs = rl.begin(); rs != rl.end(); ++rs, ++i) {
    BOOST_CHECK_EQUAL(rs->GetNumber().GetNum(), i);
  }
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("          testing secondary structure...");
  // pick chains, iterate residues, check for correct sec.struct.
  ch = eh.FindChain("A");
  rl = ch.GetResidueList();
  BOOST_CHECK_EQUAL(rl[0].GetSecStructure().IsHelical(), true);
  BOOST_CHECK_EQUAL(rl[1].GetSecStructure().IsHelical(), true);
  BOOST_CHECK_EQUAL(rl[2].GetSecStructure().IsExtended(), true);
  ch = eh.FindChain("C");
  rl = ch.GetResidueList();
  BOOST_CHECK_EQUAL(rl[0].GetSecStructure().IsExtended(), true);
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("          reading data fields which should not fail...");
  BOOST_TEST_MESSAGE(mmcif_p.GetInfo().GetBioUnits().back().GetID());
  BOOST_CHECK(mmcif_p.GetInfo().GetMethod().str() == "Deep-fry");
  BOOST_CHECK(mmcif_p.GetInfo().GetBioUnits().back().GetDetails() ==
              "author_defined_assembly");
  BOOST_CHECK(mmcif_p.GetInfo().GetBioUnits().back().GetMethodDetails() == "?");
  BOOST_CHECK(mmcif_p.GetInfo().GetBioUnits().back().GetID() == "2");
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
  #if OST_DOUBLE_PRECISION
  BOOST_CHECK_CLOSE(sd.GetMass(), 1.0, 0.001);
  #else
  BOOST_CHECK_CLOSE(sd.GetMass(), 1.0f, 0.001f);
  #endif
  BOOST_CHECK(sd.GetMassMethod() == "Good Guess");
  BOOST_CHECK(sd.GetModelDetails() == "Even better guessing");
  BOOST_CHECK(sd.GetModelTypeDetails() == "Guess");
  MMCifInfoObsolete obs = mmcif_p.GetInfo().GetObsoleteInfo();
  BOOST_CHECK(obs.GetDate() == "2011-08-31");
  BOOST_CHECK(obs.GetID() == "Obsolete");
  BOOST_CHECK(obs.GetPDBID() == "1FOO");
  BOOST_CHECK(obs.GetReplacedPDBID() == "2BAR");
  BOOST_TEST_MESSAGE("          done.");

  BOOST_TEST_MESSAGE("  done.");
}

// helper for mmcif_test_chain_mappings
inline void CheckChainMap(mol::EntityHandle eh, const MMCifInfo& info,
                          const String& cif_name, const String& pdb_name,
                          bool check_info_map) {
  // check chain
  mol::ChainHandle ch = eh.FindChain(cif_name);
  BOOST_CHECK(ch.IsValid());
  BOOST_CHECK(ch.HasProp("pdb_auth_chain_name"));
  BOOST_CHECK(ch.GetStringProp("pdb_auth_chain_name") == pdb_name);
  // info mapping
  if (check_info_map) {
    BOOST_CHECK(info.GetMMCifPDBChainTr(cif_name) == pdb_name);
    BOOST_CHECK(info.GetPDBMMCifChainTr(pdb_name) == cif_name);
  }
}

BOOST_AUTO_TEST_CASE(mmcif_test_chain_mappings)
{
  BOOST_TEST_MESSAGE("  Running mmcif_test_chain_mappings tests...");
  
  // check compound lib
  bool compound_lib_available = static_cast<bool>(SetDefaultCompoundLib());

  // load data
  mol::EntityHandle eh = mol::CreateEntity();
  std::ifstream s("testfiles/mmcif/atom_site.mmcif");
  IOProfile profile;
  MMCifReader mmcif_p(s, eh, profile);
  if (compound_lib_available) {
    mmcif_p.SetReadSeqRes(true);
  }
  BOOST_REQUIRE_NO_THROW(mmcif_p.Parse());
  const MMCifInfo& info = mmcif_p.GetInfo();
  
  // check 1-to-1 mappings
  CheckChainMap(eh, info, "A", "A", compound_lib_available);
  CheckChainMap(eh, info, "C", "C", compound_lib_available);
  CheckChainMap(eh, info, "O", "B", false); // water
  CheckChainMap(eh, info, "Z", "Z", compound_lib_available);

  // check entity ID mapping
  BOOST_CHECK(info.GetMMCifEntityIdTr("A") == "1");
  BOOST_CHECK(info.GetMMCifEntityIdTr("C") == "1");
  BOOST_CHECK(info.GetMMCifEntityIdTr("O") == "5");
  BOOST_CHECK(info.GetMMCifEntityIdTr("Z") == "1");
  
  // check non-existent mappings
  BOOST_CHECK(info.GetMMCifPDBChainTr("B") == "");
  BOOST_CHECK(info.GetPDBMMCifChainTr("O") == "");
  BOOST_CHECK(info.GetMMCifEntityIdTr("B") == "");

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_test_revisions_old)
{
  BOOST_TEST_MESSAGE("  Running mmcif_test_revisions_old tests...");
  
  // load data
  mol::EntityHandle eh = mol::CreateEntity();
  std::ifstream s("testfiles/mmcif/3IMJ_rev_old.cif");
  IOProfile profile;
  MMCifReader mmcif_p(s, eh, profile);
  BOOST_REQUIRE_NO_THROW(mmcif_p.Parse());
  const MMCifInfo& info = mmcif_p.GetInfo();
  const MMCifInfoRevisions& revs = info.GetRevisions();
  
  // check items (note: version numbers only for new mmCIF versions)
  BOOST_CHECK_EQUAL(revs.GetSize(), size_t(3));
  BOOST_CHECK_EQUAL(revs.GetNum(0), 1);
  BOOST_CHECK_EQUAL(revs.GetDate(0), String("2009-11-17"));
  BOOST_CHECK_EQUAL(revs.GetStatus(0), String("full release"));
  BOOST_CHECK_EQUAL(revs.GetMajor(0), -1);
  BOOST_CHECK_EQUAL(revs.GetMinor(0), -1);
  BOOST_CHECK_EQUAL(revs.GetNum(1), 2);
  BOOST_CHECK_EQUAL(revs.GetDate(1), String("2011-07-13"));
  BOOST_CHECK_EQUAL(revs.GetStatus(1), String("?"));
  BOOST_CHECK_EQUAL(revs.GetMajor(1), -1);
  BOOST_CHECK_EQUAL(revs.GetMinor(1), -1);
  BOOST_CHECK_EQUAL(revs.GetNum(2), 3);
  BOOST_CHECK_EQUAL(revs.GetDate(2), String("2012-12-12"));
  BOOST_CHECK_EQUAL(revs.GetStatus(2), String("?"));
  BOOST_CHECK_EQUAL(revs.GetMajor(2), -1);
  BOOST_CHECK_EQUAL(revs.GetMinor(2), -1);
  // check rest
  BOOST_CHECK_EQUAL(revs.GetDateOriginal(), String("2009-08-10"));
  BOOST_CHECK_EQUAL(revs.GetLastDate(), String("2012-12-12"));
  BOOST_CHECK_EQUAL(revs.GetLastMajor(), -1);
  BOOST_CHECK_EQUAL(revs.GetLastMinor(), -1);
  BOOST_CHECK_EQUAL(revs.GetFirstRelease(), size_t(1));

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_test_revisions_new)
{
  BOOST_TEST_MESSAGE("  Running mmcif_test_revisions_new tests...");
  
  // load data
  mol::EntityHandle eh = mol::CreateEntity();
  std::ifstream s("testfiles/mmcif/3IMJ_rev_new.cif");
  IOProfile profile;
  MMCifReader mmcif_p(s, eh, profile);
  BOOST_REQUIRE_NO_THROW(mmcif_p.Parse());
  const MMCifInfo& info = mmcif_p.GetInfo();
  const MMCifInfoRevisions& revs = info.GetRevisions();
  
  // check items
  BOOST_CHECK_EQUAL(revs.GetSize(), size_t(3));
  BOOST_CHECK_EQUAL(revs.GetNum(0), 1);
  BOOST_CHECK_EQUAL(revs.GetDate(0), String("2009-11-17"));
  BOOST_CHECK_EQUAL(revs.GetStatus(0), String("Initial release"));
  BOOST_CHECK_EQUAL(revs.GetMajor(0), 1);
  BOOST_CHECK_EQUAL(revs.GetMinor(0), 0);
  BOOST_CHECK_EQUAL(revs.GetNum(1), 2);
  BOOST_CHECK_EQUAL(revs.GetDate(1), String("2011-07-13"));
  BOOST_CHECK_EQUAL(revs.GetStatus(1), String("?"));
  BOOST_CHECK_EQUAL(revs.GetMajor(1), 1);
  BOOST_CHECK_EQUAL(revs.GetMinor(1), 1);
  BOOST_CHECK_EQUAL(revs.GetNum(2), 3);
  BOOST_CHECK_EQUAL(revs.GetDate(2), String("2012-12-12"));
  BOOST_CHECK_EQUAL(revs.GetStatus(2), String("?"));
  BOOST_CHECK_EQUAL(revs.GetMajor(2), 2);
  BOOST_CHECK_EQUAL(revs.GetMinor(2), 0);
  // check rest
  BOOST_CHECK_EQUAL(revs.GetDateOriginal(), String("2009-08-10"));
  BOOST_CHECK_EQUAL(revs.GetLastDate(), String("2012-12-12"));
  BOOST_CHECK_EQUAL(revs.GetLastMajor(), 2);
  BOOST_CHECK_EQUAL(revs.GetLastMinor(), 0);
  BOOST_CHECK_EQUAL(revs.GetFirstRelease(), size_t(1));

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_pdbx_entity_branch_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_pdbx_entity_branch_tests...");
  StarLoopDesc tmmcif_h;
  mol::EntityHandle eh = mol::CreateEntity();

  BOOST_TEST_MESSAGE("          testing chain type recognition...");
  {
    TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
    std::vector<StringRef> columns;

    // create corresponding entity entry
    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("entity", 6));
    tmmcif_h.Add(StringRef("id", 2));
    tmmcif_h.Add(StringRef("type", 4));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("branched", 8));
    tmmcif_p.ParseEntity(columns);
    columns.pop_back();
    columns.pop_back();

    // build dummy pdbx_entity_branch header
    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("pdbx_entity_branch", 18));
    tmmcif_h.Add(StringRef("entity_id", 9));
    tmmcif_h.Add(StringRef("type", 4));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("oligosaccharide", 15));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParsePdbxEntityBranch(columns));
    columns.pop_back();

    columns.push_back(StringRef("ordinarysugar", 13));
    BOOST_CHECK_THROW(tmmcif_p.ParsePdbxEntityBranch(columns), ost::Error);
  }
  BOOST_TEST_MESSAGE("          done.");
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_pdbx_entity_branch_link_tests)
{
  BOOST_TEST_MESSAGE("  Running mmcif_pdbx_entity_branch_link_tests...");
  mol::EntityHandle eh = mol::CreateEntity();
  StarLoopDesc tmmcif_h;

  BOOST_TEST_MESSAGE("          testing link retrieval...");
  {
    TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
    std::vector<StringRef> columns;

    // build dummy pdbx_entity_branch_link header
    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("pdbx_entity_branch_link", 23));
    tmmcif_h.Add(StringRef("entity_id", 9));
    tmmcif_h.Add(StringRef("atom_id_1", 9));
    tmmcif_h.Add(StringRef("atom_id_2", 9));
    tmmcif_h.Add(StringRef("comp_id_1", 9));
    tmmcif_h.Add(StringRef("comp_id_2", 9));
    tmmcif_h.Add(StringRef("entity_branch_list_num_1", 24));
    tmmcif_h.Add(StringRef("entity_branch_list_num_2", 24));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("C1", 2));
    columns.push_back(StringRef("O3", 2));
    columns.push_back(StringRef("MAN", 3));
    columns.push_back(StringRef("BMA", 3));
    columns.push_back(StringRef("2", 1));
    columns.push_back(StringRef("1", 1));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParsePdbxEntityBranchLink(columns));
  }
  BOOST_TEST_MESSAGE("          done.");
  BOOST_TEST_MESSAGE("          testing link with bond order...");
  {
    TestMMCifReaderProtected tmmcif_p("testfiles/mmcif/atom_site.mmcif", eh);
    std::vector<StringRef> columns;

    // build dummy pdbx_entity_branch_link header
    tmmcif_h.Clear();
    tmmcif_h.SetCategory(StringRef("pdbx_entity_branch_link", 23));
    tmmcif_h.Add(StringRef("entity_id", 9));
    tmmcif_h.Add(StringRef("atom_id_1", 9));
    tmmcif_h.Add(StringRef("atom_id_2", 9));
    tmmcif_h.Add(StringRef("comp_id_1", 9));
    tmmcif_h.Add(StringRef("comp_id_2", 9));
    tmmcif_h.Add(StringRef("entity_branch_list_num_1", 24));
    tmmcif_h.Add(StringRef("entity_branch_list_num_2", 24));
    tmmcif_h.Add(StringRef("value_order", 11));
    tmmcif_p.OnBeginLoop(tmmcif_h);
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("C1", 2));
    columns.push_back(StringRef("O3", 2));
    columns.push_back(StringRef("MAN", 3));
    columns.push_back(StringRef("BMA", 3));
    columns.push_back(StringRef("2", 1));
    columns.push_back(StringRef("1", 1));
    columns.push_back(StringRef("2", 1));
    BOOST_CHECK_NO_THROW(tmmcif_p.ParsePdbxEntityBranchLink(columns));
  }
  BOOST_TEST_MESSAGE("          done.");
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_SUITE_END();
