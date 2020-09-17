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

#include <ost/io/io_exception.hh>
#include <ost/io/mol/mmcif_info.hh>
#include <ost/mol/mol.hh>

using namespace ost;
using namespace ost::io;

BOOST_AUTO_TEST_SUITE( io );

BOOST_AUTO_TEST_CASE(mmcif_info_obsolete)
{
  BOOST_TEST_MESSAGE("  Running mmcif_info_obsolete tests...");

  MMCifInfoObsolete obs = MMCifInfoObsolete();

  obs.SetDate("2011-08-31");
  obs.SetID(StringRef("OBSLTE", 6));
  obs.SetPDBID("1FOO");
  obs.SetReplacedPDBID("1BAR");

  BOOST_CHECK(obs.GetDate() == "2011-08-31");
  BOOST_CHECK(obs.GetID() == "Obsolete");
  BOOST_CHECK(obs.GetPDBID() == "1FOO");
  BOOST_CHECK(obs.GetReplacedPDBID() == "1BAR");

  obs.SetID(StringRef("SPRSDE", 6));
  BOOST_CHECK(obs.GetID() == "Supersede");

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_info_citation)
{
  BOOST_TEST_MESSAGE("  Running mmcif_info_citation tests...");

  MMCifInfoCitation cit = MMCifInfoCitation();
  std::vector<String> author_list;

  author_list.push_back("Kabel, H.");

  cit.SetID("ID");
  cit.SetCAS("FOO");
  cit.SetISBN("0-0-0-0-0");
  cit.SetPublishedIn("Journal of Uncanny Science");
  cit.SetVolume("3");
  cit.SetPageFirst("1");
  cit.SetPageLast("10");
  cit.SetDOI("HERE");
  cit.SetPubMed(815);
  cit.SetYear(815);
  cit.SetTitle("Foo");
  cit.SetBookPublisher("Brackelmann and Sons");
  cit.SetBookPublisherCity("Stenkelfeld");
  cit.SetAuthorList(author_list);
  cit.SetCitationType(MMCifInfoCitation::JOURNAL);
  author_list.clear();

  BOOST_CHECK(cit.GetID() == "ID");
  BOOST_CHECK(cit.GetCAS() == "FOO");
  BOOST_CHECK(cit.GetISBN() == "0-0-0-0-0");
  BOOST_CHECK(cit.GetPublishedIn() == "Journal of Uncanny Science");
  BOOST_CHECK(cit.GetVolume() == "3");
  BOOST_CHECK(cit.GetPageFirst() == "1");
  BOOST_CHECK(cit.GetPageLast() == "10");
  BOOST_CHECK(cit.GetDOI() == "HERE");
  BOOST_CHECK(cit.GetPubMed() == 815);
  BOOST_CHECK(cit.GetYear() == 815);
  BOOST_CHECK(cit.GetTitle() == "Foo");
  BOOST_CHECK(cit.GetBookPublisher() == "Brackelmann and Sons");
  BOOST_CHECK(cit.GetBookPublisherCity() == "Stenkelfeld");
  BOOST_CHECK(cit.GetCitationType() == MMCifInfoCitation::JOURNAL);
  BOOST_CHECK(cit.IsCitationTypeJournal() == true);
  BOOST_CHECK(cit.IsCitationTypeBook() == false);
  BOOST_CHECK(cit.IsCitationTypeUnknown() == false);
  author_list = cit.GetAuthorList();
  BOOST_CHECK(author_list.back() == "Kabel, H.");

  // checking all possible variants of citation type
  cit.SetCitationTypeJournal();
  BOOST_CHECK(cit.IsCitationTypeJournal() == true);
  cit.SetCitationTypeBook();
  BOOST_CHECK(cit.IsCitationTypeBook() == true);
  cit.SetCitationTypeUnknown();
  BOOST_CHECK(cit.IsCitationTypeUnknown() == true);

  BOOST_TEST_MESSAGE("  done.");
  BOOST_TEST_MESSAGE("  trying to add everything to an info object");
  MMCifInfo info = MMCifInfo();
  info.AddCitation(cit);
  std::vector<MMCifInfoCitation> citations = info.GetCitations();
  BOOST_CHECK(citations.size() == 1);
  BOOST_CHECK(citations.back() == cit);
  BOOST_CHECK_THROW(info.AddAuthorsToCitation(StringRef("Foo", 3),
                                              author_list),
                    IOException);
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_info_biounit)
{
  BOOST_TEST_MESSAGE("  Running mmcif_info_biounit tests...");

  MMCifInfoBioUnit bu = MMCifInfoBioUnit();

  bu.SetDetails("author_defined_assembly");
  bu.SetMethodDetails("PISA");
  bu.SetID("1");
  bu.AddChain("B");
  bu.AddChain("A");

  BOOST_CHECK(bu.GetDetails() == "author_defined_assembly");
  BOOST_CHECK(bu.GetMethodDetails() == "PISA");
  BOOST_CHECK(bu.GetChainList().back() == "A");
  std::vector<std::pair<int, int> > tr = bu.GetChainIntervalList();
  BOOST_CHECK(tr[0].first == 0);
  BOOST_CHECK(tr[0].second == 2);

  std::vector<String> chains;
  chains.push_back("B");
  bu.SetChainList(chains);
  BOOST_CHECK(bu.GetChainList().back() == "B");
  tr = bu.GetChainIntervalList();
  BOOST_CHECK(tr[0].first == 0);
  BOOST_CHECK(tr[0].second == 1);

  MMCifInfo info = MMCifInfo();
  info.AddBioUnit(bu);
  std::vector<MMCifInfoBioUnit> biounits = info.GetBioUnits();
  BOOST_CHECK(biounits.size() == 1);
  BOOST_CHECK(biounits.back() == bu);

  info.AddBioUnit(bu);
  bu.SetID("2");
  info.AddBioUnit(bu);
  biounits = info.GetBioUnits();
  BOOST_CHECK(biounits.size() == 2);
  tr = biounits.front().GetChainIntervalList();
  BOOST_CHECK(tr.size() == 2);
  BOOST_CHECK(tr[0].first == 0);
  BOOST_CHECK(tr[0].second == 1);
  BOOST_CHECK(tr[1].first == 1);
  BOOST_CHECK(tr[1].second == 2);

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_info_transoperation)
{
  BOOST_TEST_MESSAGE("  Running mmcif_info_transoperation tests...");

  MMCifInfoTransOpPtr op(new MMCifInfoTransOp);

  op->SetID("1");
  op->SetType("identity operation");
  op->SetVector(1, 2, 3);
  op->SetMatrix(1, 2, 3, 4, 5, 6, 7 ,8, 9);

  BOOST_CHECK(op->GetID() == "1");
  BOOST_CHECK(op->GetVector() == geom::Vec3(1, 2, 3));
  BOOST_CHECK(op->GetMatrix() == geom::Mat3(1, 2, 3, 4, 5, 6, 7, 8, 9));
  BOOST_CHECK(op->GetType() == "identity operation");

  MMCifInfo info = MMCifInfo();
  info.AddOperation(op);
  BOOST_CHECK(info.GetOperations().back() == op);

  std::vector<MMCifInfoTransOpPtr> ops;
  ops.push_back(*(info.GetOperations().begin()));
  MMCifInfoBioUnit bu = MMCifInfoBioUnit();
  bu.AddOperations(ops);
  BOOST_CHECK((*(bu.GetOperations().begin()->begin())) == op);
  std::vector<std::pair<int, int> > tr = bu.GetOperationsIntervalList();
  BOOST_CHECK(tr.size() == 1);
  BOOST_CHECK(tr.back().first == 0);
  BOOST_CHECK(tr.back().second == 1);

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_info_structdetails)
{
  BOOST_TEST_MESSAGE("  Running mmcif_info_structdetails tests...");

  MMCifInfoStructDetails sd = MMCifInfoStructDetails();

  sd.SetEntryID("1BAR");
  sd.SetTitle("More than a structure");
  sd.SetCASPFlag('Y');
  sd.SetDescriptor("ADENYLATE KINASE");
  sd.SetMass(1.0);
  sd.SetMassMethod("Good Guess");
  sd.SetModelDetails("Even more guessing");
  sd.SetModelTypeDetails("MINIMIZED AVERAGE");

  BOOST_CHECK(sd.GetEntryID() == "1BAR");
  BOOST_CHECK(sd.GetTitle() == "More than a structure");
  BOOST_CHECK(sd.GetCASPFlag() == 'Y');
  BOOST_CHECK(sd.GetDescriptor() == "ADENYLATE KINASE");
  #if OST_DOUBLE_PRECISION
  BOOST_CHECK_CLOSE(sd.GetMass(), 1.0, 0.001);
  #else
  BOOST_CHECK_CLOSE(sd.GetMass(), 1.0f, 0.001f);
  #endif
  BOOST_CHECK(sd.GetMassMethod() == "Good Guess");
  BOOST_CHECK(sd.GetModelDetails() == "Even more guessing");
  BOOST_CHECK(sd.GetModelTypeDetails() == "MINIMIZED AVERAGE");

  MMCifInfo info = MMCifInfo();
  info.SetStructDetails(sd);
  BOOST_CHECK(info.GetStructDetails() == sd);

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_info_revisions)
{
  BOOST_TEST_MESSAGE("  Running mmcif_info_revisions tests...");

  MMCifInfoRevisions rev = MMCifInfoRevisions();

  BOOST_CHECK(rev.GetDateOriginal() == "?");
  BOOST_CHECK(rev.GetFirstRelease() == 0);
  BOOST_CHECK(rev.GetSize() == 0);
  BOOST_CHECK(rev.GetLastDate() == "?");
  BOOST_CHECK(rev.GetLastMajor() == -1);
  BOOST_CHECK(rev.GetLastMinor() == -1);
  BOOST_CHECK_THROW(rev.GetDate(0), std::out_of_range);
  BOOST_CHECK_THROW(rev.GetNum(0), std::out_of_range);
  BOOST_CHECK_THROW(rev.GetStatus(0), std::out_of_range);

  rev.SetDateOriginal("2012-05-03");
  rev.AddRevision(1, "2012-05-04", "in preparation");
  rev.AddRevision(2, "2012-05-05", "full release", 1, 2);

  BOOST_CHECK(rev.GetSize() == 2);
  BOOST_CHECK(rev.GetDateOriginal() == "2012-05-03");
  BOOST_CHECK(rev.GetDate(0) == "2012-05-04");
  BOOST_CHECK(rev.GetNum(0) == 1);
  BOOST_CHECK(rev.GetStatus(0) == "in preparation");
  BOOST_CHECK(rev.GetMajor(0) == -1);
  BOOST_CHECK(rev.GetMinor(0) == -1);
  BOOST_CHECK(rev.GetDate(1) == "2012-05-05");
  BOOST_CHECK(rev.GetDate(1) == rev.GetLastDate());
  BOOST_CHECK(rev.GetFirstRelease() == 2);
  BOOST_CHECK(rev.GetNum(1) == 2);
  BOOST_CHECK(rev.GetStatus(1) == "full release");
  BOOST_CHECK(rev.GetMajor(1) == 1);
  BOOST_CHECK(rev.GetMinor(1) == 2);
  BOOST_CHECK(rev.GetMajor(1) == rev.GetLastMajor());
  BOOST_CHECK(rev.GetMinor(1) == rev.GetLastMinor());
  BOOST_CHECK_THROW(rev.GetDate(2), std::out_of_range);
  BOOST_CHECK_THROW(rev.GetNum(2), std::out_of_range);
  BOOST_CHECK_THROW(rev.GetStatus(2), std::out_of_range);
  BOOST_CHECK_THROW(rev.GetMajor(2), std::out_of_range);
  BOOST_CHECK_THROW(rev.GetMinor(2), std::out_of_range);

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_info_branch)
{
  BOOST_TEST_MESSAGE("  Running mmcif_info_branch tests...");

  // create a dummy entity to start an editor...
  mol::EntityHandle eh = mol::CreateEntity();
  mol::XCSEditor editor = eh.EditXCS();
  mol::ChainHandle ch = editor.InsertChain("A");
  mol::ResidueHandle res1 = editor.AppendResidue(ch, "NAG");
  mol::ResidueHandle res2 = editor.AppendResidue(ch, "NAG");
  // create AtomHandles for testing
  mol::AtomHandle atom1 = editor.InsertAtom(res2, "C1",geom::Vec3());
  mol::AtomHandle atom2 = editor.InsertAtom(res1, "O4",geom::Vec3());

  MMCifInfoEntityBranchLink branch1(atom1, atom2, 1);
  BOOST_CHECK(branch1.GetAtom1().GetQualifiedName() == "A.NAG2.C1");
  BOOST_CHECK(branch1.GetAtom2().GetQualifiedName() == "A.NAG1.O4");
  BOOST_CHECK(branch1.GetBondOrder() == 1);

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_info)
{
  BOOST_TEST_MESSAGE("  Running mmcif_info tests...");

  MMCifInfo info = MMCifInfo();

  info.SetMethod("Cooking.");
  info.SetResolution(1.9f);

  BOOST_CHECK(info.GetMethod() == StringRef("Cooking.", 8));
  #if OST_DOUBLE_PRECISION
  BOOST_CHECK_CLOSE(info.GetResolution(), 1.9, 0.001);
  #else
  BOOST_CHECK_CLOSE(info.GetResolution(), 1.9f, 0.001f);
  #endif

  info.AddMMCifPDBChainTr("A", "B");
  BOOST_CHECK_THROW(info.AddMMCifPDBChainTr("A", "B"), IOException);
  BOOST_CHECK("B" == info.GetMMCifPDBChainTr("A"));
  BOOST_CHECK("" == info.GetMMCifPDBChainTr("C"));

  info.AddPDBMMCifChainTr("A", "B");
  BOOST_CHECK_THROW(info.AddPDBMMCifChainTr("A", "B"), IOException);
  BOOST_CHECK("B" == info.GetPDBMMCifChainTr("A"));
  BOOST_CHECK("" == info.GetPDBMMCifChainTr("C"));

  info.AddMMCifEntityIdTr("A", "1");
  BOOST_CHECK_THROW(info.AddMMCifEntityIdTr("A", "1"), IOException);
  BOOST_CHECK("1" == info.GetMMCifEntityIdTr("A"));
  BOOST_CHECK("" == info.GetMMCifEntityIdTr("C"));

  BOOST_CHECK(info.GetRevisions().GetSize() == 0);

  // simple check that we can add branch links
  mol::EntityHandle eh = mol::CreateEntity();
  mol::XCSEditor editor = eh.EditXCS();
  mol::ChainHandle ch1 = editor.InsertChain("A");
  mol::ResidueHandle res11 = editor.AppendResidue(ch1, "NAG");
  mol::ResidueHandle res12 = editor.AppendResidue(ch1, "NAG");
  // create AtomHandles for testing
  mol::AtomHandle atom11 = editor.InsertAtom(res12, "C1",geom::Vec3());
  mol::AtomHandle atom12 = editor.InsertAtom(res11, "O4",geom::Vec3());
  mol::ChainHandle ch2 = editor.InsertChain("B");
  mol::ResidueHandle res21 = editor.AppendResidue(ch2, "BMA");
  mol::ResidueHandle res22 = editor.AppendResidue(ch2, "MAN");
  // create AtomHandles for testing
  mol::AtomHandle atom21 = editor.InsertAtom(res22, "C1",geom::Vec3());
  mol::AtomHandle atom22 = editor.InsertAtom(res21, "O3",geom::Vec3());
  info.AddEntityBranchLink(ch1.GetName(), atom11, atom12, 1);
  info.AddEntityBranchLink(ch2.GetName(), atom21, atom22, 1);
  std::vector<MMCifInfoEntityBranchLink> blinks = info.GetEntityBranchLinks();

  BOOST_CHECK(blinks.size() == 2);
  BOOST_CHECK(blinks[0].GetAtom1().GetQualifiedName() == "A.NAG2.C1");
  BOOST_CHECK(blinks[0].GetAtom2().GetQualifiedName() == "A.NAG1.O4");
  BOOST_CHECK(blinks[0].GetBondOrder() == 1);
  BOOST_CHECK(blinks[1].GetAtom1().GetQualifiedName() == "B.MAN2.C1");
  BOOST_CHECK(blinks[1].GetAtom2().GetQualifiedName() == "B.BMA1.O3");
  BOOST_CHECK(blinks[1].GetBondOrder() == 1);

  // check that branch links get bonds
  info.ConnectBranchLinks();

  BOOST_CHECK(atom11.GetBondPartners()[0] == atom12);
  BOOST_CHECK(atom22.GetBondPartners()[0] == atom21);

  // check chain(name) retrieval works
  std::vector<String> chain_names = info.GetEntityBranchChainNames();
  BOOST_CHECK(chain_names[0] == "A");
  BOOST_CHECK(chain_names[1] == "B");

  // check chain(handle) retrieval works
  mol::ChainHandleList chains = info.GetEntityBranchChains();
  BOOST_CHECK(chains[0].GetName() == "A");
  BOOST_CHECK(chains[1].GetName() == "B");

  // check retrieval of links by chain name
  std::vector<MMCifInfoEntityBranchLink> cblinks =
    info.GetEntityBranchByChain("A");
  BOOST_CHECK(cblinks.size() == 1);
  cblinks = info.GetEntityBranchByChain("C");
  BOOST_CHECK(cblinks.size() == 0);

  BOOST_TEST_MESSAGE("  done.");
}


BOOST_AUTO_TEST_SUITE_END();
