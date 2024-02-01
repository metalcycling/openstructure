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
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
using boost::unit_test_framework::test_suite;

#include <ost/test_utils/compare_files.hh>
#include <ost/mol/mol.hh>
#include <ost/io/mol/entity_io_sdf_handler.hh>
#include <ost/io/mol/save_entity.hh>
#include <ost/io/io_exception.hh>

using namespace ost;
using namespace ost::io;

BOOST_AUTO_TEST_SUITE( io );

BOOST_AUTO_TEST_CASE(simple_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;
  sdfh.Import(eh,"testfiles/sdf_v3000/simple.sdf");

  // check compounds/atoms/bonds count
  BOOST_CHECK_EQUAL(eh.GetChainCount(), 1);
  BOOST_CHECK_EQUAL(eh.GetAtomCount(),  6);
  BOOST_CHECK_EQUAL(eh.GetBondCount(),  6);
  BOOST_CHECK_CLOSE(eh.GetMass(), Real(121.546997), Real(1e-4));

  // check atom/bond types
  mol::AtomHandle ah=eh.GetAtomList()[0];
  mol::AtomHandle ah2=eh.GetAtomList()[5];

  BOOST_CHECK_EQUAL(ah.GetElement(),  "N");
  BOOST_CHECK_EQUAL(ah2.GetElement(), "CL");
  BOOST_CHECK_CLOSE(ah.GetRadius(),  Real(1.55), Real(1e-2));
  BOOST_CHECK_CLOSE(ah2.GetRadius(), Real(1.75), Real(1e-2));
  BOOST_CHECK_CLOSE(ah.GetMass(),  Real(14.007), Real(1e-4));
  BOOST_CHECK_CLOSE(ah2.GetMass(), Real(35.453), Real(1e-3));
  BOOST_CHECK_EQUAL(ah.GetBondCount(),  3);
  BOOST_CHECK_EQUAL(ah2.GetBondCount(), 1);
  BOOST_CHECK_EQUAL(ah.GetCharge(),  1);
  BOOST_CHECK_EQUAL(ah2.GetCharge(), 0);

  mol::BondHandle bh=ah.GetBondList()[0];
  BOOST_CHECK_EQUAL(bh.GetBondOrder(), 2);
}

BOOST_AUTO_TEST_CASE(multiple_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;
  sdfh.Import(eh,"testfiles/sdf_v3000/multiple.sdf");

  // check number of compounds
  BOOST_CHECK_EQUAL(eh.GetChainCount(), 4);
}

BOOST_AUTO_TEST_CASE(properties_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;
  sdfh.Import(eh,"testfiles/sdf_v3000/properties.sdf");

  // check number of compounds
  mol::ChainHandleList chl=eh.GetChainList();
  int count=1;
  for (mol::ChainHandleList::iterator i=chl.begin();i!=chl.end();++i,count++)
  {
    BOOST_REQUIRE(i->HasProp("prop_one"));
    BOOST_REQUIRE(i->HasProp("prop_two"));
    BOOST_CHECK_CLOSE(boost::lexical_cast<Real>(i->GetStringProp("prop_one")),
                      Real(count),Real(1e-4));
    BOOST_CHECK_CLOSE(boost::lexical_cast<Real>(i->GetStringProp("prop_two")),
                      Real(count*(-2.2)),Real(1e-4));
  }
}

BOOST_AUTO_TEST_CASE(read_sdf_v3000)
{
  const String fname("testfiles/sdf_v3000/compound.sdf");

  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;

  // check import
  sdfh.Import(eh,"testfiles/sdf_v3000/compound.sdf");

  // check atoms/bonds
  BOOST_CHECK_EQUAL(eh.GetChainCount(), 4);
  BOOST_CHECK_EQUAL(eh.GetAtomCount(), 180);
  BOOST_CHECK_EQUAL(eh.GetBondCount(), 188);

  // check molecule name
  mol::ChainHandle ch=eh.FindChain("00003_Test Ligand");
  BOOST_CHECK(ch.IsValid());

  // check properties
  BOOST_CHECK(ch.HasProp("r_i_glide_rmsd"));
  BOOST_CHECK_EQUAL(boost::lexical_cast<float>(boost::trim_copy
                     (ch.GetStringProp("r_i_glide_rmsd"))),
                     0.543804f);
}

BOOST_AUTO_TEST_CASE(wrong_atomcount_error_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;
  BOOST_CHECK_THROW(sdfh.Import(eh,"testfiles/sdf_v3000/wrong_atomcount.sdf"), IOException);
}

BOOST_AUTO_TEST_CASE(wrong_bondcount_error_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;
  BOOST_CHECK_THROW(sdfh.Import(eh,"testfiles/sdf_v3000/wrong_bondcount.sdf"), IOException);
}

BOOST_AUTO_TEST_CASE(wrong_atomlinelength_error_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;

  BOOST_CHECK_THROW(sdfh.Import(eh,"testfiles/sdf_v3000/wrong_atomlinelength.sdf"), IOException);
}

BOOST_AUTO_TEST_CASE(wrong_atompos_error_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;

  BOOST_CHECK_THROW(sdfh.Import(eh,"testfiles/sdf_v3000/wrong_atompos.sdf"), IOException);
}

BOOST_AUTO_TEST_CASE(wrong_charge_error_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;

  BOOST_CHECK_THROW(sdfh.Import(eh,"testfiles/sdf_v3000/wrong_charge.sdf"), IOException);
}

BOOST_AUTO_TEST_CASE(wrong_bondlinelength_error_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;

  BOOST_CHECK_THROW(sdfh.Import(eh,"testfiles/sdf_v3000/wrong_bondlinelength.sdf"), IOException);
}

BOOST_AUTO_TEST_CASE(wrong_bondtype_error_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;

  BOOST_CHECK_THROW(sdfh.Import(eh,"testfiles/sdf_v3000/wrong_bondtype.sdf"), IOException);
}

BOOST_AUTO_TEST_CASE(wrong_bondatomnumber_error_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;

  BOOST_CHECK_THROW(sdfh.Import(eh,"testfiles/sdf_v3000/wrong_bondatomnumber.sdf"), IOException);
}

BOOST_AUTO_TEST_CASE(wrong_keyword_error_sdf_v3000)
{
  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOSDFHandler sdfh;

  BOOST_CHECK_THROW(sdfh.Import(eh,"testfiles/sdf_v3000/wrong_keyword.sdf"), IOException);
}


BOOST_AUTO_TEST_SUITE_END();
