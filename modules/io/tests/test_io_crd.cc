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

#include <ost/test_utils/compare_files.hh>
#include <ost/mol/mol.hh>
#include <ost/io/mol/entity_io_crd_handler.hh>
#include <ost/io/mol/save_entity.hh>

using namespace ost;
using namespace ost::io;


BOOST_AUTO_TEST_SUITE( io );

BOOST_AUTO_TEST_CASE(test_io_crd_import_handler)
{
  String fname("testfiles/crd/test_in.crd");

  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOCRDHandler crdh;

  BOOST_CHECK(EntityIOCRDHandler::ProvidesImport("","crd"));
  BOOST_CHECK(EntityIOCRDHandler::ProvidesImport(fname));
  BOOST_CHECK(EntityIOCRDHandler::ProvidesImport("test_in.CRD"));

  BOOST_CHECK(EntityIOCRDHandler::ProvidesExport("","crd"));
  BOOST_CHECK(EntityIOCRDHandler::ProvidesExport(fname));
  BOOST_CHECK(EntityIOCRDHandler::ProvidesExport("test_in.CRD"));

}

/// check standard format input
BOOST_AUTO_TEST_CASE(test_io_crd) 
{
  const String fname("testfiles/crd/test_in.crd");

  mol::EntityHandle eh=mol::CreateEntity();
  CRDReader crdr(fname);

  crdr.Import(eh);
  std::vector<mol::AtomHandle> alist = crdr.GetSequentialAtoms();

  BOOST_CHECK_EQUAL(eh.GetAtomCount(), 9);
  BOOST_CHECK_EQUAL(eh.GetResidueCount(), 2);
  BOOST_CHECK_EQUAL(eh.GetChainCount(), 2);
  BOOST_CHECK_EQUAL(eh.GetResidueList()[0].GetNumber().GetNum(), 3);
  BOOST_CHECK_EQUAL(eh.GetResidueList()[1].GetNumber().GetNum(), 4);
  BOOST_CHECK_EQUAL(eh.GetResidueList()[0].GetName(), "ALA");
  BOOST_CHECK_EQUAL(eh.GetResidueList()[1].GetName(), "P2AC");
  BOOST_CHECK_EQUAL(eh.GetAtomList()[0].GetName(), "N");
  BOOST_CHECK_EQUAL(eh.GetAtomList()[7].GetName(), "HA51");
  BOOST_CHECK_EQUAL(eh.GetChainList()[0].GetName(), "PEPT");
  BOOST_CHECK_EQUAL(eh.GetChainList()[1].GetName(), "RNA");
}

/// check extended format input
BOOST_AUTO_TEST_CASE(test_io_crd_ext)
{
  const String fname("testfiles/crd/multisegment.crd");

  mol::EntityHandle eh=mol::CreateEntity();
  CRDReader crdr(fname);

  crdr.Import(eh);
  std::vector<mol::AtomHandle> alist = crdr.GetSequentialAtoms();

  BOOST_CHECK_EQUAL(eh.GetAtomCount(), 43);
  BOOST_CHECK_EQUAL(eh.GetResidueCount(), 6);
  BOOST_CHECK_EQUAL(eh.GetChainCount(), 3);
  BOOST_CHECK_EQUAL(eh.GetResidueList()[0].GetNumber().GetNum(), 1);
  BOOST_CHECK_EQUAL(eh.GetResidueList()[1].GetNumber().GetNum(), 2);
  BOOST_CHECK_EQUAL(eh.GetResidueList()[2].GetNumber().GetNum(), 1);
  BOOST_CHECK_EQUAL(eh.GetResidueList()[0].GetName(), "GLY");
  BOOST_CHECK_EQUAL(eh.GetResidueList()[1].GetName(), "GLU");
  BOOST_CHECK_EQUAL(eh.GetResidueList()[2].GetName(), "G3AC");
  BOOST_CHECK_EQUAL(eh.GetAtomList()[0].GetName(), "N");
  BOOST_CHECK_EQUAL(eh.GetAtomList()[20].GetName(), "OE1");
  BOOST_CHECK_EQUAL(eh.GetChainList()[0].GetName(), "PROTEIN");
  BOOST_CHECK_EQUAL(eh.GetChainList()[1].GetName(), "RNA");
  BOOST_CHECK_EQUAL(eh.GetChainList()[2].GetName(), "SAM");
}

/// check that empty lines at end of file are ignored
BOOST_AUTO_TEST_CASE(test_io_crd_empty_line)
{
  const String fname("testfiles/crd/empty_line_at_end.crd");
  mol::EntityHandle eh=mol::CreateEntity();
  CRDReader reader(fname);
  reader.Import(eh);
}

/// check extended format with multisegment file output
BOOST_AUTO_TEST_CASE(test_io_crd_writer_ext)
{
  const String fname("testfiles/crd/multisegment.crd");
  {
    mol::EntityHandle eh=mol::CreateEntity();
    CRDReader reader(fname);
    reader.Import(eh);

    CRDWriter writer(String("testfiles/crd/multisegment-out.crd"), true);
    writer.Write(eh);
  }
  BOOST_CHECK(compare_files("testfiles/crd/multisegment-out.crd",
                            "testfiles/crd/multisegment.crd"));
}

/// check standard format output
BOOST_AUTO_TEST_CASE(test_io_crd_writer_std)
{
  const String fname("testfiles/crd/test_in.crd");
  {
    mol::EntityHandle eh=mol::CreateEntity();
    CRDReader reader(fname);
    reader.Import(eh);

    CRDWriter writer(String("testfiles/crd/test-out.crd"));
    writer.Write(eh);
  }
  BOOST_CHECK(compare_files("testfiles/crd/test-out.crd",
                            "testfiles/crd/test_in.crd"));
}

/// check standard format is written by default
BOOST_AUTO_TEST_CASE(test_io_crd_default_format)
{
  const String fname("testfiles/crd/test_in.crd");
  {
    mol::EntityHandle eh=mol::CreateEntity();
    CRDReader reader(fname);
    reader.Import(eh);

    ost::io::SaveEntity(eh, "testfiles/crd/test-out.crd");
  }
  BOOST_CHECK(compare_files("testfiles/crd/test-out.crd",
                            "testfiles/crd/test_in.crd"));
}


BOOST_AUTO_TEST_SUITE_END();
