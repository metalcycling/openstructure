
// -----------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>

// Copyright (C) 2008-2020 by the OpenStructure authors

// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.

// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// -----------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <ost/mol/mol.hh>
#include <ost/platform.hh>
#include <ost/log.hh>
#include <ost/conop/rule_based.hh>
#include <ost/conop/conop.hh>
#include <ost/mol/builder.hh>
#include "helper.hh"

using boost::unit_test_framework::test_suite;
using namespace ost;
using namespace ost::conop;
//using namespace ost::io;
using namespace ost::mol;

CompoundLibPtr load_lib()
{
  if (!getenv("OST_ROOT")) {
    LOG_ERROR("OST_ROOT environment variable not set. Can't load "
              "compound library without a proper OST_ROOT");
    return CompoundLibPtr();
  }
  SetPrefixPath(getenv("OST_ROOT"));
  String lib_path=GetSharedDataPath()+"/compounds.chemlib";
  CompoundLibPtr compound_lib=CompoundLib::Load(lib_path);  
  return compound_lib;
}

BOOST_AUTO_TEST_SUITE(conop);

BOOST_AUTO_TEST_CASE(rule_based_init_check)
{
  CompoundLibPtr lib; // null ptr is return value when library loading failed
  BOOST_CHECK_THROW(RuleBasedProcessor rbc1(lib), ost::Error);
  BOOST_CHECK_THROW(RuleBasedProcessor rbc2(lib, true, false, CONOP_WARN,
                                            CONOP_WARN, false, true, true, true,
                                            true, CONOP_WARN), ost::Error);
  lib = load_lib();
  if (!lib) { return; }
  BOOST_CHECK_NO_THROW(RuleBasedProcessor rbc3(lib));
  BOOST_CHECK_NO_THROW(RuleBasedProcessor rbc4(lib, true, false, CONOP_WARN,
                                               CONOP_WARN, false, true, true,
                                               true, true, CONOP_WARN));
}

BOOST_AUTO_TEST_CASE(rule_based_set_get_flags)
{
   CompoundLibPtr lib = load_lib();
   if (!lib) { return; }
   RuleBasedProcessor rbc(lib);
   // check the defaults
   BOOST_CHECK_EQUAL(rbc.GetConnect(), true);
   BOOST_CHECK_EQUAL(rbc.GetCheckBondFeasibility(), false);
   BOOST_CHECK_EQUAL(rbc.GetUnkAtomTreatment(), CONOP_WARN);
   BOOST_CHECK_EQUAL(rbc.GetUnkResidueTreatment(), CONOP_WARN);
   BOOST_CHECK_EQUAL(rbc.GetStrictHydrogens(), false);
   rbc.SetConnect(false);
   rbc.SetStrictHydrogens(true);
   rbc.SetCheckBondFeasibility(true);
   rbc.SetUnkResidueTreatment(CONOP_FATAL);
   rbc.SetUnkAtomTreatment(CONOP_REMOVE);
   BOOST_CHECK_EQUAL(rbc.GetConnect(), false);
   BOOST_CHECK_EQUAL(rbc.GetStrictHydrogens(), true);
   BOOST_CHECK_EQUAL(rbc.GetCheckBondFeasibility(),true);
   BOOST_CHECK_EQUAL(rbc.GetUnkResidueTreatment(),CONOP_FATAL);
   BOOST_CHECK_EQUAL(rbc.GetUnkAtomTreatment(),CONOP_REMOVE);
}


BOOST_AUTO_TEST_CASE(rule_based_connect)
{
  CompoundLibPtr lib = load_lib();
  if (!lib) { return; }
  RuleBasedProcessor rbc(lib);
  EntityHandle ent = CreateEntity();
  ChainHandle ch=ent.EditXCS().InsertChain("A");
  ResidueHandle r = make_cytosine(ch);
  rbc.SetConnect(false); 
  rbc.Process(ent);
  BOOST_CHECK_EQUAL(static_cast<size_t>(0),
                    ent.GetBondList().size());
  rbc.SetConnect(true); 
  rbc.Process(ent);
  verify_nucleotide_connectivity(r);
}


BOOST_AUTO_TEST_CASE(rule_based_unk_atoms)
{
  CompoundLibPtr lib = load_lib();
  if (!lib) { return; }
  RuleBasedProcessor rbc(lib);
  EntityHandle ent = CreateEntity();
  XCSEditor edi=ent.EditXCS();
  ChainHandle ch = edi.InsertChain("A");
  ResidueHandle gly = edi.AppendResidue(ch, "GLY");
  edi.InsertAtom(gly, "N", geom::Vec3());
  edi.InsertAtom(gly, "CA", geom::Vec3());
  edi.InsertAtom(gly, "C", geom::Vec3());
  edi.InsertAtom(gly, "O", geom::Vec3());
  edi.InsertAtom(gly, "CB", geom::Vec3());
  EntityHandle ent2 = ent.Copy();
  rbc.SetUnkAtomTreatment(CONOP_SILENT);
  DiagnosticsPtr diags=rbc.Process(ent2);
  BOOST_CHECK_EQUAL(diags->diag_count(), static_cast<size_t>(0));
  rbc.SetUnkAtomTreatment(CONOP_WARN);
  diags = rbc.Process(ent2);
  BOOST_CHECK_EQUAL(diags->diag_count(), static_cast<size_t>(1));
  ent2 = ent.Copy();
  rbc.SetUnkAtomTreatment(CONOP_REMOVE);
  rbc.Process(ent2);
  BOOST_CHECK(!ent2.FindAtom("A", 1, "CB"));
  ent2 = ent.Copy();
  rbc.SetUnkAtomTreatment(CONOP_REMOVE_RESIDUE);
  rbc.Process(ent2);
  BOOST_CHECK_EQUAL(ent2.GetResidueCount(), 0);
}

BOOST_AUTO_TEST_CASE(guesses_elements_of_unknown_atoms)
{
  CompoundLibPtr lib = load_lib();
  if (!lib) { return; }
  RuleBasedProcessor rbc(lib);
  EntityHandle ent = CreateEntity();
  XCSEditor edi=ent.EditXCS();
  ChainHandle ch = edi.InsertChain("A");
  ResidueHandle gly = edi.AppendResidue(ch, "???");
  ResidueHandle gly2 = edi.AppendResidue(ch, "GLY");
  AtomHandle ca1 = edi.InsertAtom(gly, "CA", geom::Vec3(0,0,0));
  AtomHandle ca2 = edi.InsertAtom(gly2, "CX", geom::Vec3(0,0,0));
  DiagnosticsPtr diags=rbc.Process(ent);
  BOOST_CHECK_EQUAL(ca1.GetElement(), "C");
  BOOST_CHECK_EQUAL(ca2.GetElement(), "C");
}

BOOST_AUTO_TEST_CASE(fills_properties_of_unknown_residues)
{
  CompoundLibPtr lib = load_lib();
  if (!lib) { return; }
  RuleBasedProcessor rbc(lib);
  EntityHandle ent = CreateEntity();
  XCSEditor edi=ent.EditXCS();
  ChainHandle ch = edi.InsertChain("A");
  ResidueHandle gly = edi.AppendResidue(ch, "???");
  gly.SetOneLetterCode('G');
  gly.SetChemClass(ChemClass(ChemClass::PEPTIDE_LINKING));
  DiagnosticsPtr diags=rbc.Process(ent);
  BOOST_CHECK_EQUAL(gly.GetOneLetterCode(), '?');
  BOOST_CHECK_EQUAL(gly.GetName(), "???");
  BOOST_CHECK_EQUAL(gly.GetChemClass(), ChemClass(ChemClass::UNKNOWN));
}

BOOST_AUTO_TEST_CASE(connects_atoms_of_unknown_residues_based_on_distance) 
{

  CompoundLibPtr lib = load_lib();
  if (!lib) { return; }
  RuleBasedProcessor rbc(lib);
  EntityHandle ent = CreateEntity();
  XCSEditor edi=ent.EditXCS();
  ChainHandle ch = edi.InsertChain("A");
  ResidueHandle gly = edi.AppendResidue(ch, "???");
  AtomHandle ca = edi.InsertAtom(gly, "CA", geom::Vec3(0,0,0));
  AtomHandle c = edi.InsertAtom(gly, "C", geom::Vec3(1.54,0,0));
  AtomHandle o = edi.InsertAtom(gly, "O", geom::Vec3(1.54,3,0));
  DiagnosticsPtr diags=rbc.Process(ent);
  BOOST_CHECK(BondExists(ca, c));
  BOOST_CHECK(!BondExists(c, o));
}

BOOST_AUTO_TEST_CASE(rule_based_unk_res)
{
  CompoundLibPtr lib = load_lib();
  if (!lib) { return; }
  RuleBasedProcessor rbc(lib);
  EntityHandle ent = CreateEntity();
  XCSEditor edi=ent.EditXCS();
  ChainHandle ch = edi.InsertChain("A");
  ResidueHandle gly = edi.AppendResidue(ch, "???");
  EntityHandle ent2 = ent.Copy();
  rbc.SetUnkResidueTreatment(CONOP_SILENT); 
  DiagnosticsPtr diags=rbc.Process(ent2);
  BOOST_CHECK_EQUAL(diags->diag_count(), static_cast<size_t>(0));
  rbc.SetUnkResidueTreatment(CONOP_WARN);
  diags = rbc.Process(ent2);
  BOOST_CHECK_EQUAL(diags->diag_count(), static_cast<size_t>(1));
  ent2 = ent.Copy();
  rbc.SetUnkResidueTreatment(CONOP_REMOVE);
  rbc.Process(ent2);
  BOOST_CHECK(!ent2.FindResidue("A", 1));
}

BOOST_AUTO_TEST_CASE(hetatom_connect_rule_based) {

  CompoundLibPtr lib = load_lib();
  if (!lib) { return; }

  // STEP 1: Specify two atoms as hetatoms - they should be connected
  // by the processor
  mol::EntityHandle e = Builder()
    .Chain("A")
      .Residue("GLY")
        .Atom("N", geom::Vec3(-8.22, 35.20, 22.39))
        .Atom("CA", geom::Vec3(-8.28, 36.36, 21.49))
        .Atom("C", geom::Vec3(-8.59, 35.93, 20.06))
        .Atom("O", geom::Vec3(-7.88, 36.30, 19.12))
        .Atom("CB", geom::Vec3(-6.96, 37.11, 21.53))
  ;

  ost::mol::AtomHandleList atoms = e.GetAtomList();
  atoms[0].SetHetAtom(true); // N
  atoms[1].SetHetAtom(true); // CA
  RuleBasedProcessor proc(lib);
  proc.Process(e);
  BOOST_CHECK(mol::BondExists(e.FindAtom("A", 1, "N"), e.FindAtom("A", 1, "CA")));
  BOOST_CHECK(mol::BondExists(e.FindAtom("A", 1, "CA"), e.FindAtom("A", 1, "C")));
  
  // STEP 2: Same thing again but we tell the processor NOT to connect
  // hetatoms
  e = Builder()
    .Chain("A")
      .Residue("GLY")
        .Atom("N", geom::Vec3(-8.22, 35.20, 22.39))
        .Atom("CA", geom::Vec3(-8.28, 36.36, 21.49))
        .Atom("C", geom::Vec3(-8.59, 35.93, 20.06))
        .Atom("O", geom::Vec3(-7.88, 36.30, 19.12))
        .Atom("CB", geom::Vec3(-6.96, 37.11, 21.53))
  ;

  atoms = e.GetAtomList();
  atoms[0].SetHetAtom(true); // N
  atoms[1].SetHetAtom(true); // CA
  proc.SetConnectHetatm(false);
  proc.Process(e);
  BOOST_CHECK(!mol::BondExists(e.FindAtom("A", 1, "N"), e.FindAtom("A", 1, "CA")));
  BOOST_CHECK(mol::BondExists(e.FindAtom("A", 1, "CA"), e.FindAtom("A", 1, "C")));
}

BOOST_AUTO_TEST_SUITE_END();
