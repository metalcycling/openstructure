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
using boost::unit_test_framework::test_suite;

#include <ost/test_utils/compare_files.hh>
#include <ost/geom/vec_mat_predicates.hh>

#include <ost/platform.hh>
#include <ost/dyn_cast.hh>
#include <ost/mol/mol.hh>
#include <ost/conop/conop.hh>
#include <ost/conop/heuristic.hh>

#include <ost/io/mol/entity_io_pdb_handler.hh>
#include <ost/io/mol/pdb_reader.hh>
#include <ost/io/mol/pdb_writer.hh>
#include <ost/log.hh>

#include <ost/io/io_exception.hh>

using namespace ost;
using namespace ost::io;

BOOST_AUTO_TEST_SUITE( io );

BOOST_AUTO_TEST_CASE(test_pdb_import_handler) 
{
  String fname("testfiles/test_in.pdb");

  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOPDBHandler pdbh;

  BOOST_CHECK(EntityIOPDBHandler::ProvidesImport("","pdb"));
  BOOST_CHECK(EntityIOPDBHandler::ProvidesImport(fname));
  BOOST_CHECK(EntityIOPDBHandler::ProvidesImport("test_in.PDB"));
  
  BOOST_CHECK(EntityIOPDBHandler::ProvidesExport("","pdb"));
  BOOST_CHECK(EntityIOPDBHandler::ProvidesExport(fname));
  BOOST_CHECK(EntityIOPDBHandler::ProvidesExport("test_in.PDB"));
    
  pdbh.Import(eh,"testfiles/pdb/simple.pdb");
}

BOOST_AUTO_TEST_CASE(test_parse_compnd_record) 
{
  String fname("testfiles/pdb/compnd.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_REQUIRE_EQUAL(ent.GetChainCount(), 28);
  BOOST_REQUIRE_EQUAL(ent.GetResidueCount(), 28);

  mol::ChainHandle ch = ent.FindChain("A");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
  
  ch = ent.FindChain("B");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==2);

  ch = ent.FindChain("W");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==2);

  
}

BOOST_AUTO_TEST_CASE(test_parse_compnd_record2) 
{
  String fname("testfiles/pdb/1AKE.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);

  mol::ChainHandle ch = ent.FindChain("A");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
}

BOOST_AUTO_TEST_CASE(test_parse_compnd_record3) 
{
  String fname("testfiles/pdb/1AKE.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent, "B");

  mol::ChainHandle ch = ent.FindChain("B");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
}

BOOST_AUTO_TEST_CASE(test_parse_compnd_record4) 
{
  String fname("testfiles/pdb/3mk3.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);

  mol::ChainHandle ch = ent.FindChain("B");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
  
  ch = ent.FindChain("y");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
  
  ch = ent.FindChain("3");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
  
  ch = ent.FindChain("U");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
  
  ch = ent.FindChain("i");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
}

//COMPND CHAIN record misses B chain
BOOST_AUTO_TEST_CASE(test_parse_compnd_record5) 
{
  Logger::Instance().PushVerbosityLevel(0);

  String fname("testfiles/pdb/1AKE_noBchain.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  
  BOOST_CHECK_NO_THROW(reader.Import(ent));
  
  mol::ChainHandle ch = ent.FindChain("A");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
  
  ch = ent.FindChain("B");
  BOOST_CHECK(ch.HasProp("mol_id")==false);
}

//chain I in MOL_ID record but no chain I
BOOST_AUTO_TEST_CASE(test_parse_compnd_record6) 
{
  Logger::Instance().PushVerbosityLevel(0);
  String fname("testfiles/pdb/1oax.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
  
}

// has an empy MOLECULE record (unsupported ATM anyway, but crashed ost)
BOOST_AUTO_TEST_CASE(test_parse_compnd_record7) 
{
  Logger::Instance().PushVerbosityLevel(0);
  String fname("testfiles/pdb/2p6a.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  
  reader.Import(ent);
  mol::ChainHandle ch = ent.FindChain("A");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==1);
  ch = ent.FindChain("D");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==2);

  ch = ent.FindChain("E");
  BOOST_CHECK(ch.HasProp("mol_id")==true);
  BOOST_CHECK(ch.GetIntProp("mol_id")==3);
}


//COMPND   3 CHAIN: A V; 
BOOST_AUTO_TEST_CASE(test_parse_compnd_missing_comma_chain_record) 
{
  String fname("testfiles/pdb/compnd_missing_comma_chain_record.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

//COMPND   3 CHAIN A, V; 
BOOST_AUTO_TEST_CASE(test_parse_compnd_missing_colon_chain_record) 
{
  String fname("testfiles/pdb/compnd_missing_colon_chain_record.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

//COMPND    MOL_ID: ;
BOOST_AUTO_TEST_CASE(test_parse_compnd_missing_mol_id_complete_mol_id_record) 
{
  String fname("testfiles/pdb/compnd_missing_mol_id_complete_mol_id_record.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

//COMPND    MOL_ID 1;
BOOST_AUTO_TEST_CASE(test_parse_compnd_missing_colon_mol_id_record) 
{
  String fname("testfiles/pdb/compnd_missing_colon_mol_id_record.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

//COMPND    MOL_ID: 1 
BOOST_AUTO_TEST_CASE(test_parse_compnd_missing_semicolon_mol_id_record) 
{
  String fname("testfiles/pdb/compnd_missing_semicolon_mol_id_record.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

//COMPND   3 CHAIN: A, V  
BOOST_AUTO_TEST_CASE(test_parse_compnd_missing_semicolon_chain_record) 
{
  String fname("testfiles/pdb/compnd_missing_semicolon_chain_record.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

//COMPND   3 CHAIN: 
BOOST_AUTO_TEST_CASE(test_parse_compnd_missing_chains_chain_record) 
{
  String fname("testfiles/pdb/compnd_missing_chains_chain_record.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

//no chain record
BOOST_AUTO_TEST_CASE(test_parse_compnd_missing_complete_chain_record) 
{
  String fname("testfiles/pdb/compnd_missing_complete_chain_record.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}
//no MOL_ID record
BOOST_AUTO_TEST_CASE(test_parse_compnd_missing_complete_mol_id_record) 
{
  String fname("testfiles/pdb/compnd_missing_complete_mol_id_record.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

BOOST_AUTO_TEST_CASE(atom_record)
{
  String fname("testfiles/pdb/atom.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_REQUIRE_EQUAL(ent.GetChainCount(), 2);
  BOOST_REQUIRE_EQUAL(ent.GetResidueCount(), 2);
  BOOST_REQUIRE_EQUAL(ent.GetAtomCount(), 2);
  mol::AtomHandle a1=ent.FindAtom("A", mol::ResNum(1), "N");
  BOOST_REQUIRE(a1.IsValid());
  BOOST_CHECK_EQUAL(a1.GetName(), "N");
  BOOST_CHECK_EQUAL(a1.GetResidue().GetName(), "MET");
  BOOST_CHECK_EQUAL(a1.GetResidue().GetChain().GetName(), "A");  
  
  BOOST_CHECK_EQUAL(a1.GetPos(), geom::Vec3(16.0, 64.0, 8.0));
  BOOST_CHECK_EQUAL(a1.GetBFactor(), 1.0);
  BOOST_CHECK_EQUAL(a1.GetOccupancy(), 0.5);  
  BOOST_CHECK_EQUAL(a1.GetElement(), "N");
  BOOST_CHECK_EQUAL(a1.IsHetAtom(), false);
  mol::AtomHandle a2=ent.FindAtom(" ", mol::ResNum(1), "CA");
  BOOST_REQUIRE(a2.IsValid());  
  BOOST_CHECK_EQUAL(a2.GetName(), "CA");
  BOOST_CHECK_EQUAL(a2.GetResidue().GetName(), "MET");
  BOOST_CHECK_EQUAL(a2.GetResidue().GetChain().GetName(), " ");  
  
  BOOST_CHECK_EQUAL(a2.GetPos(), geom::Vec3(32.0, -128.0, -2.5));
  BOOST_CHECK_EQUAL(a2.GetBFactor(), 128.0);
  BOOST_CHECK_EQUAL(a2.GetOccupancy(), 1.0);  
  BOOST_CHECK_EQUAL(a2.GetElement(), "C");
  BOOST_CHECK_EQUAL(a2.IsHetAtom(), true);
}

BOOST_AUTO_TEST_CASE(end_record)
{
  String fname("testfiles/pdb/end.pdb");
  PDBReader reader(fname, IOProfile()); 
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 1);
}

BOOST_AUTO_TEST_CASE(join_spread_records_on)
{
  String fname("testfiles/pdb/join-spread-records.pdb");
  IOProfile profile;
  profile.join_spread_atom_records=true;
  PDBReader reader(fname, profile);
  
  mol::EntityHandle ent=mol::CreateEntity();

  reader.Import(ent);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 2);
  mol::ResidueHandle res1=ent.FindResidue("A", mol::ResNum(1));
  BOOST_CHECK(res1.IsValid());
  BOOST_CHECK_EQUAL(res1.GetAtomCount(), 2);
  BOOST_CHECK(res1.FindAtom("N"));
  BOOST_CHECK(res1.FindAtom("CE"));
  mol::ResidueHandle res2=ent.FindResidue("A", mol::ResNum(2));  
  BOOST_CHECK(res2.IsValid());
  BOOST_CHECK_EQUAL(res2.GetAtomCount(), 1);
  BOOST_CHECK(res2.FindAtom("N"));
}

BOOST_AUTO_TEST_CASE(join_spread_records_off)
{
  String fname("testfiles/pdb/join-spread-records.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 3);
  mol::ResidueHandle res1=ent.FindResidue("A", mol::ResNum(1));
  BOOST_CHECK(res1.IsValid());
  BOOST_CHECK_EQUAL(res1.GetAtomCount(), 1);
  BOOST_CHECK(res1.FindAtom("N"));
  mol::ResidueHandle res2=ent.FindResidue("A", mol::ResNum(2));  
  BOOST_REQUIRE(res2.IsValid());
  BOOST_CHECK_EQUAL(res2.GetAtomCount(), 1);
  BOOST_CHECK(res2.FindAtom("N"));
}

BOOST_AUTO_TEST_CASE(calpha_only_import_on)
{
  String fname("testfiles/pdb/calpha.pdb");
  IOProfile profile;
  profile.calpha_only=true;
  PDBReader reader(fname, profile);
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 2);
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 2);
}

BOOST_AUTO_TEST_CASE(het_import)
{
  String fname("testfiles/pdb/het.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_CHECK_EQUAL(ent.Select("ligand=true").GetResidueCount(), 1);
}

BOOST_AUTO_TEST_CASE(calpha_only_import_off)
{
  String fname("testfiles/pdb/calpha.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 2);
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 4);
}

BOOST_AUTO_TEST_CASE(anisou_record)
{
  String fname("testfiles/pdb/anisou.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_REQUIRE(ent.GetAtomCount()==1);
  mol::AtomHandle a1=ent.FindAtom("A", mol::ResNum(7), "N");
  BOOST_REQUIRE(a1.IsValid());
  const geom::Mat3& anisou=a1.GetAnisou();
  BOOST_CHECK_CLOSE(Real( 0.0100), anisou(0, 0), Real(1e-4));
  BOOST_CHECK_CLOSE(Real(-0.0016), anisou(1, 0), Real(1e-4));
  BOOST_CHECK_CLOSE(Real(-0.0026), anisou(2, 0), Real(1e-4));
  BOOST_CHECK_CLOSE(Real(-0.0016), anisou(0, 1), Real(1e-4));
  BOOST_CHECK_CLOSE(Real( 0.0110), anisou(1, 1), Real(1e-4));
  BOOST_CHECK_CLOSE(Real(-0.0054), anisou(2, 1), Real(1e-4));
  BOOST_CHECK_CLOSE(Real(-0.0026), anisou(0, 2), Real(1e-4));
  BOOST_CHECK_CLOSE(Real(-0.0054), anisou(1, 2), Real(1e-4));
  BOOST_CHECK_CLOSE(Real( 0.0120), anisou(2, 2), Real(1e-4));    
}

BOOST_AUTO_TEST_CASE(only_66_cols)
{
  String fname("testfiles/pdb/short.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
}

BOOST_AUTO_TEST_CASE(no_endmdl_record)
{
  String fname("testfiles/pdb/model.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

BOOST_AUTO_TEST_CASE(deuterium_import)
{
  String fname("testfiles/pdb/val-with-deuterium.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  // we use conopology to mark amino acids as peptide-linking. 
  conop::HeuristicProcessor heu_proc;
  heu_proc.Process(ent);
  // this check makes sure that we correctly detect deal with the deuterium
  // atoms in the residue.
  BOOST_CHECK(ent.FindResidue("A", 297).IsPeptideLinking());
}

BOOST_AUTO_TEST_CASE(bzdng_318)
{
  String fname("testfiles/pdb/bzdng-318.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  // we use conopology to mark amino acids as peptide-linking.
  conop::HeuristicProcessor heu_proc;
  heu_proc.Process(ent);
  {
    PDBWriter writer(std::string("testfiles/pdb/bzdng-318-out.pdb"),
                     IOProfile());
    writer.Write(ent);
  }

  BOOST_CHECK(compare_files("testfiles/pdb/bzdng-318.pdb",
                            "testfiles/pdb/bzdng-318-out.pdb"));
}

BOOST_AUTO_TEST_CASE(faulty_lines)
{
  String fname("testfiles/pdb/faulty.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
  IOProfile profile;
  profile.fault_tolerant=true;
  PDBReader reader2(fname, profile);
  reader2.Import(ent);
}

BOOST_AUTO_TEST_CASE(write_atom)
{
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  mol::ResidueHandle r=edi.AppendResidue(ch, "GLY");

  mol::AtomHandle a=edi.InsertAtom(r, "CA", geom::Vec3(32.0, -128.0, -2.5), "C");
  a.SetOccupancy(1.0);
  a.SetBFactor(128.0);
  writer.Write(ent);
  String s=out.str();
  BOOST_CHECK_EQUAL(s.substr(0, 54), 
                    "ATOM      1  CA  GLY A   1      32.000-128.000  -2.500");
  BOOST_CHECK_EQUAL(s.substr(54, 26), 
                    "  1.00128.00           C  ");
}

BOOST_AUTO_TEST_CASE(write_atom_100000)
{
  char c_names[] = "ABCDEFGHIJ";
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch;
  mol::ResidueHandle r;
  mol::AtomHandle a;

  for (unsigned long i = 0; i < 10; ++i) {
    ch=edi.InsertChain(String(1, c_names[i]));
    for (unsigned long j = 1; j < 1000; ++j) {
      r=edi.AppendResidue(ch, "ARG");
      a=edi.InsertAtom(r,"N",   geom::Vec3(26.861, 50.841, 38.803), "N");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"CA",  geom::Vec3(27.437, 49.969, 37.786), "C");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"C",   geom::Vec3(26.336, 48.959, 37.429), "C");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"O",   geom::Vec3(25.745, 48.313, 38.312), "O");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"CB",  geom::Vec3(28.653, 49.266, 38.349), "C");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"CG",  geom::Vec3(29.870, 50.188, 38.416), "C");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"CD",  geom::Vec3(31.033, 49.532, 39.173), "C");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"NE",  geom::Vec3(32.318, 50.244, 39.125), "N");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"CZ",  geom::Vec3(33.462, 49.750, 39.679), "C");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"NH1", geom::Vec3(33.522, 48.572, 40.308), "N");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
      a=edi.InsertAtom(r,"NH2", geom::Vec3(34.610, 50.427, 39.597), "N");
      a.SetOccupancy(1.0);
      a.SetBFactor(128.0);
    }
  }

  writer.Write(ent);
  String s=out.str();
  BOOST_CHECK_EQUAL(s.substr(8099844, 5), "99999");
  BOOST_CHECK_EQUAL(s.substr(8099925, 5), "*****");
}

BOOST_AUTO_TEST_CASE(write_hetatom)
{
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  mol::ResidueHandle r=edi.AppendResidue(ch, "CA");
  mol::AtomHandle a=edi.InsertAtom(r, "CA", geom::Vec3(32.0, -128.0, -2.5), 
                                   "CA");
  a.SetHetAtom(true);
  a.SetMass(40.01);
  a.SetOccupancy(1.0);
  a.SetBFactor(40.75);
  writer.Write(ent);
  String s=out.str();
  BOOST_CHECK_EQUAL(s.substr(0, 54), 
                    "HETATM    1 CA    CA A   1      32.000-128.000  -2.500");
  BOOST_CHECK_EQUAL(s.substr(54, 26), 
                    "  1.00 40.75          CA  ");
}

BOOST_AUTO_TEST_CASE(no_endmdl_record_fault_tolerant)
{
  String fname("testfiles/pdb/model.pdb");
  IOProfile profile;
  profile.fault_tolerant=true;
  PDBReader reader(fname, profile);
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 1);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 1);  
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 3);
  ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 1);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 1);  
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 2);  
}

BOOST_AUTO_TEST_CASE(alt_loc_import_export)
{
  String fname("testfiles/pdb/alt-loc.pdb");  
  // this scope is required to force the writer stream to be closed before 
  // opening the file again in compare_files. Avoids a race condition.
  {
    PDBReader reader(fname, IOProfile());
    PDBWriter writer(String("testfiles/pdb/alt-loc-out.pdb"), IOProfile());
    
    mol::EntityHandle ent=mol::CreateEntity();
    reader.Import(ent);
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/alt-loc.pdb", 
                            "testfiles/pdb/alt-loc-out.pdb"));
}

BOOST_AUTO_TEST_CASE(write_ter)
{
  String fname("testfiles/pdb/ter.pdb");  
  // this scope is required to force the writer stream to be closed before 
  // opening the file again in compare_files. Avoids a race condition.
  {
    PDBReader reader(fname, IOProfile());
    PDBWriter writer(String("testfiles/pdb/ter-out.pdb"), IOProfile());
    
    mol::EntityHandle ent=mol::CreateEntity();
    reader.Import(ent);
    // we use conopology to mark amino acids as peptide-linking. this is require 
    // for proper TER output
    conop::HeuristicProcessor heu_proc;

    heu_proc.Process(ent);
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/ter.pdb", 
                            "testfiles/pdb/ter-out.pdb"));
}

BOOST_AUTO_TEST_CASE(write_ter2)
{
  String fname("testfiles/pdb/ter2.pdb");  
  // this scope is required to force the writer stream to be closed before 
  // opening the file again in compare_files. Avoids a race condition.
  {
    PDBReader reader(fname, IOProfile());
    PDBWriter writer(String("testfiles/pdb/ter2-out.pdb"), IOProfile());
    
    mol::EntityHandle ent=mol::CreateEntity();
    reader.Import(ent);
    // we use conopology to mark amino acids as peptide-linking. this is 
    // require for proper TER output
    conop::HeuristicProcessor heu_proc;
    heu_proc.Process(ent);
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/ter2.pdb", 
                            "testfiles/pdb/ter2-out.pdb"));
}

BOOST_AUTO_TEST_CASE(write_ter3)
{
  String fname("testfiles/pdb/ter3.pdb");  
  // this scope is required to force the writer stream to be closed before 
  // opening the file again in compare_files. Avoids a race condition.
  {
    PDBReader reader(fname, IOProfile());
    PDBWriter writer(String("testfiles/pdb/ter3-out.pdb"), IOProfile());
    
    mol::EntityHandle ent=mol::CreateEntity();
    reader.Import(ent);
    // we use conopology to mark amino acids as peptide-linking. this is 
    // require for proper TER output
    conop::HeuristicProcessor heu_proc;
    heu_proc.Process(ent);
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/ter3.pdb", 
                            "testfiles/pdb/ter3-out.pdb"));
}

BOOST_AUTO_TEST_CASE(write_ter4)
{
  {
    PDBWriter writer(String("testfiles/pdb/ter_emptychain-out.pdb"), IOProfile());

    mol::EntityHandle ent=mol::CreateEntity();
    mol::XCSEditor edi=ent.EditXCS();
    mol::ChainHandle ch=edi.InsertChain("A");
    mol::ResidueHandle r=edi.AppendResidue(ch, "ALA");
    mol::AtomHandle a=edi.InsertAtom(r, "N", geom::Vec3(32.0, -128.0, -2.5));
    mol::ResidueHandle r2=edi.AppendResidue(ch, "GLY");
    mol::AtomHandle a2=edi.InsertAtom(r2, "N", geom::Vec3(35.0, -99.0, -10.5));
    mol::ChainHandle ch2=edi.InsertChain("B");
    conop::HeuristicProcessor heu_proc;
    heu_proc.Process(ent);
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/ter_emptychain.pdb",
                            "testfiles/pdb/ter_emptychain-out.pdb"));
}

BOOST_AUTO_TEST_CASE(write_ter5)
{
  {
    PDBWriter writer(String("testfiles/pdb/ter_view-emptychain-out.pdb"), IOProfile());

    mol::EntityHandle ent=mol::CreateEntity();
    mol::XCSEditor edi=ent.EditXCS();
    mol::ChainHandle ch=edi.InsertChain("A");
    mol::ResidueHandle r=edi.AppendResidue(ch, "ALA");
    mol::AtomHandle a=edi.InsertAtom(r, "N", geom::Vec3(32.0, -128.0, -2.5));
    mol::ResidueHandle r2=edi.AppendResidue(ch, "GLY");
    mol::AtomHandle a2=edi.InsertAtom(r2, "N", geom::Vec3(35.0, -99.0, -10.5));
    mol::ChainHandle ch2=edi.InsertChain("B");
    conop::HeuristicProcessor heu_proc;
    heu_proc.Process(ent);
    writer.Write(ent.Select(""));
  }
  BOOST_CHECK(compare_files("testfiles/pdb/ter_view-emptychain.pdb",
                            "testfiles/pdb/ter_view-emptychain-out.pdb"));
}

BOOST_AUTO_TEST_CASE(write_ter6)
{
  String fname("testfiles/pdb/ter4.pdb");  
  // this scope is required to force the writer stream to be closed before 
  // opening the file again in compare_files. Avoids a race condition.
  {
    PDBReader reader(fname, IOProfile());
    PDBWriter writer(String("testfiles/pdb/ter4-out.pdb"), IOProfile());
    
    mol::EntityHandle ent=mol::CreateEntity();
    reader.Import(ent);
    conop::HeuristicProcessor heu_proc;
    heu_proc.Process(ent);
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/ter4.pdb", 
                            "testfiles/pdb/ter4-out.pdb"));
}

BOOST_AUTO_TEST_CASE(read_write_conect)
{
  // this scope is required to force the writer stream to be closed before
  // opening the file again in compare_files. Avoids a race condition.
  {
    IOProfile profile;
    profile.read_conect=true;
    PDBReader reader(String("testfiles/pdb/conect.pdb"), profile);
    PDBWriter writer(String("testfiles/pdb/conect-out.pdb"), profile);
    mol::EntityHandle ent=mol::CreateEntity();
    reader.Import(ent);
    conop::HeuristicProcessor heu_proc;
    heu_proc.Process(ent);
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/conect.pdb",
                            "testfiles/pdb/conect-out.pdb"));
}

BOOST_AUTO_TEST_CASE(alt_loc_tf)
{
  String fname("testfiles/pdb/alt-loc.pdb");  
  // this scope is required to force the writer stream to be closed before 
  // opening the file again in compare_files. Avoids a race condition.
  mol::EntityHandle ent=mol::CreateEntity();  
  PDBReader reader(fname, IOProfile());
  reader.Import(ent);
  String out_name("testfiles/pdb/alt-loc-tf-out.pdb");
  {
    PDBWriter writer(out_name, IOProfile());
    geom::Mat4 shift;
    shift.PasteTranslation(geom::Vec3(10,20,30));
    ent.EditXCS().ApplyTransform(shift);
    writer.Write(ent);
  }
  PDBReader r2(out_name, IOProfile());
  mol::EntityHandle ent2=mol::CreateEntity();
  r2.Import(ent2);
  mol::ResidueHandle res1=ent2.FindResidue("A", mol::ResNum(1));
  mol::AtomHandle a1=res1.FindAtom("N");
  BOOST_CHECK_EQUAL(res1.GetAltAtomPos(a1, "A"), geom::Vec3(26,84,30));
  BOOST_CHECK_EQUAL(res1.GetAltAtomPos(a1, "B"), geom::Vec3(18,-108,30));  
}

BOOST_AUTO_TEST_CASE(res_name_too_long)
{
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  mol::ResidueHandle r=edi.AppendResidue(ch, "CALCIUM");
  mol::AtomHandle a=edi.InsertAtom(r, "CA", geom::Vec3(32.0, -128.0, -2.5));
  BOOST_CHECK_THROW(writer.Write(ent), IOException);
}

BOOST_AUTO_TEST_CASE(res_num_too_low)
{
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch = edi.InsertChain("A");
  mol::ResidueHandle r = edi.AppendResidue(ch, "ARG");
  edi.InsertAtom(r,"N",   geom::Vec3(26.861, 50.841, 38.803), "N");
  edi.InsertAtom(r,"CA",  geom::Vec3(27.437, 49.969, 37.786), "C");
  edi.InsertAtom(r,"C",   geom::Vec3(26.336, 48.959, 37.429), "C");
  edi.InsertAtom(r,"O",   geom::Vec3(25.745, 48.313, 38.312), "O");
  edi.InsertAtom(r,"CB",  geom::Vec3(28.653, 49.266, 38.349), "C");
  edi.InsertAtom(r,"CG",  geom::Vec3(29.870, 50.188, 38.416), "C");
  edi.InsertAtom(r,"CD",  geom::Vec3(31.033, 49.532, 39.173), "C");
  edi.InsertAtom(r,"NE",  geom::Vec3(32.318, 50.244, 39.125), "N");
  edi.InsertAtom(r,"CZ",  geom::Vec3(33.462, 49.750, 39.679), "C");
  edi.InsertAtom(r,"NH1", geom::Vec3(33.522, 48.572, 40.308), "N");
  edi.InsertAtom(r,"NH2", geom::Vec3(34.610, 50.427, 39.597), "N");

  edi.SetResidueNumber(r, mol::ResNum(-1000));
  BOOST_CHECK_THROW(writer.Write(ent), IOException);
}

BOOST_AUTO_TEST_CASE(res_num_too_high)
{
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch = edi.InsertChain("A");
  mol::ResidueHandle r = edi.AppendResidue(ch, "ARG");
  edi.InsertAtom(r,"N",   geom::Vec3(26.861, 50.841, 38.803), "N");
  edi.InsertAtom(r,"CA",  geom::Vec3(27.437, 49.969, 37.786), "C");
  edi.InsertAtom(r,"C",   geom::Vec3(26.336, 48.959, 37.429), "C");
  edi.InsertAtom(r,"O",   geom::Vec3(25.745, 48.313, 38.312), "O");
  edi.InsertAtom(r,"CB",  geom::Vec3(28.653, 49.266, 38.349), "C");
  edi.InsertAtom(r,"CG",  geom::Vec3(29.870, 50.188, 38.416), "C");
  edi.InsertAtom(r,"CD",  geom::Vec3(31.033, 49.532, 39.173), "C");
  edi.InsertAtom(r,"NE",  geom::Vec3(32.318, 50.244, 39.125), "N");
  edi.InsertAtom(r,"CZ",  geom::Vec3(33.462, 49.750, 39.679), "C");
  edi.InsertAtom(r,"NH1", geom::Vec3(33.522, 48.572, 40.308), "N");
  edi.InsertAtom(r,"NH2", geom::Vec3(34.610, 50.427, 39.597), "N");

  edi.SetResidueNumber(r, mol::ResNum(10000));
  BOOST_CHECK_THROW(writer.Write(ent), IOException);
}


BOOST_AUTO_TEST_CASE(res_name_mismatch_alt_loc)
{
  String fname("testfiles/pdb/arg-glu-gln.pdb");
  IOProfile profile;
  profile.fault_tolerant=true;
  PDBReader reader(fname, profile);
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 1);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 1);  
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 11);
}

BOOST_AUTO_TEST_CASE(res_name_mismatch_alt_loc_pedantic)
{
  String fname("testfiles/pdb/arg-glu-gln.pdb");
  IOProfile profile;
  PDBReader reader(fname, profile);
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_NO_THROW(reader.Import(ent));
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 1);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 1);  
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 11);
}

BOOST_AUTO_TEST_CASE(res_name_mismatch_pedantic)
{
  String fname("testfiles/pdb/more-than-one-name.pdb");
  IOProfile profile;
  PDBReader reader(fname, profile);
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

BOOST_AUTO_TEST_CASE(duplicate_atom_strict)
{
  String fname("testfiles/pdb/more-than-one-name.pdb");
  IOProfile profile;
  PDBReader reader(fname, profile);
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_THROW(reader.Import(ent), IOException);
}

BOOST_AUTO_TEST_CASE(duplicate_atom_tolerant)
{
  String fname("testfiles/pdb/duplicate-atom.pdb");
  IOProfile profile;
  profile.fault_tolerant=true;
  PDBReader reader(fname, profile);
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_NO_THROW(reader.Import(ent));
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 1);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 1);  
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 2);
}

BOOST_AUTO_TEST_CASE(res_name_mismatch_tolerant)
{
  String fname("testfiles/pdb/more-than-one-name.pdb");
  IOProfile profile;
  profile.fault_tolerant=true;
  PDBReader reader(fname, profile);
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_NO_THROW(reader.Import(ent));
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 1);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 1);  
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 4);
}

BOOST_AUTO_TEST_CASE(seqres_import)
{
  char * ost_root=getenv("OST_ROOT");
  if(!ost_root){
    std::cout << "WARNING: skipping SEQRES import unit test. " 
              << "Rule-based processor is required" << std::endl;
    return;
  }
  SetPrefixPath(ost_root);
  String lib_path=GetSharedDataPath()+"/compounds.chemlib";
  conop::CompoundLibPtr compound_lib=conop::CompoundLib::Load(lib_path);  
  if (!compound_lib) {
    std::cout << "WARNING: skipping SEQRES import unit test. Compound lib is " 
              << "required" << std::endl;
    return;    
  }
  conop::Conopology::Instance().SetDefaultLib(compound_lib);
  String fname("testfiles/pdb/seqres.pdb");
  IOProfile profile;
  PDBReader reader(fname, profile);
  reader.SetReadSeqRes(true);
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_NO_THROW(reader.Import(ent));
  seq::SequenceList seqres=reader.GetSeqRes();
  BOOST_CHECK_EQUAL(seqres.GetCount(), 1);
  BOOST_CHECK_EQUAL(seqres[0].GetName(), String("A"));
  BOOST_CHECK_EQUAL(seqres[0].GetString(), 
                    String("GSNIGETLGEKWKSRLNALGKSEFQIYKKSGIQEVDRTLAKEGIKRGETDHH"
                           "AVSRGSAKLRWFVERNLVTPEGKVVDLGCGRGGWSYYCGGLKNVREVKGLTK"
                           "GGPGHEEPIPMSTYGWNLVRLQSGVDVFFIPPERCDTLLCDIGESSPNPTVE"
                           "AGRTLRVLNLVENWLSNNTQFCVKVLNPYMPSVIEKMEALQRKHGGALVRNP"
                           "LSRNSTHEMYWVSNASGNIVSSVNMISRMLINRFTMRHKKATYEPDVDLGSG"
                           "TRNIGIESETPNLDIIGKRIEKIKQEHETSWHYDQ"));
}


BOOST_AUTO_TEST_CASE(chain_name_too_long)
{
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("AB");
  mol::ResidueHandle r=edi.AppendResidue(ch, "CA");
  mol::AtomHandle a=edi.InsertAtom(r, "CA", geom::Vec3(32.0, -128.0, -2.5));
  BOOST_CHECK_THROW(writer.Write(ent), IOException);
}

BOOST_AUTO_TEST_CASE(atom_name_too_long)
{
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  mol::ResidueHandle r=edi.AppendResidue(ch, "CA");
  mol::AtomHandle a=edi.InsertAtom(r, "CALCIUM", geom::Vec3(32.0, -128.0, -2.5));
  BOOST_CHECK_THROW(writer.Write(ent), IOException);
}

BOOST_AUTO_TEST_CASE(anisou_BZDNG_274)
{
  String fname("testfiles/pdb/bzdng274.pdb");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  BOOST_CHECK_NO_THROW(reader.Import(ent));
}

BOOST_AUTO_TEST_CASE(charmm_rname)
{
  {
    PDBWriter writer(String("testfiles/pdb/charmm_rname-out.pdb"),
                     IOProfile("CHARMM", false, false, false, false, false, false));

    mol::EntityHandle ent=mol::CreateEntity();
    mol::XCSEditor edi=ent.EditXCS();
    mol::ChainHandle ch=edi.InsertChain("A");
    mol::ResidueHandle r=edi.AppendResidue(ch, "ALA");
    mol::AtomHandle a=edi.InsertAtom(r, "N", geom::Vec3(32.0, -128.0, -2.5));
    mol::ResidueHandle r2=edi.AppendResidue(ch, "G3AG");
    mol::AtomHandle a2=edi.InsertAtom(r2, "N", geom::Vec3(32.0, -128.0, -2.5));
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/charmm_rname.pdb",
                            "testfiles/pdb/charmm_rname-out.pdb"));
}

BOOST_AUTO_TEST_CASE(charmm_longcname)
{
  {
    PDBWriter writer(String("testfiles/pdb/charmm_longcname-out.pdb"),
                     IOProfile("CHARMM", false, false, false, false, false, false));

    mol::EntityHandle ent=mol::CreateEntity();
    mol::XCSEditor edi=ent.EditXCS();
    mol::ChainHandle ch=edi.InsertChain("PROT");
    mol::ResidueHandle r=edi.AppendResidue(ch, "ALA");
    mol::AtomHandle a=edi.InsertAtom(r, "N", geom::Vec3(32.0, -128.0, -2.5));
    mol::ResidueHandle r2=edi.AppendResidue(ch, "GLY");
    mol::AtomHandle a2=edi.InsertAtom(r2, "N", geom::Vec3(35.0, -99.0, -10.5));
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/charmm_longcname.pdb",
                            "testfiles/pdb/charmm_longcname-out.pdb"));
}

BOOST_AUTO_TEST_CASE(write_charmm_ter)
{
  {
    PDBWriter writer(String("testfiles/pdb/charmm_ter-out.pdb"),
                     IOProfile("CHARMM", false, false, false, false, false, false));

    mol::EntityHandle ent=mol::CreateEntity();
    mol::XCSEditor edi=ent.EditXCS();
    mol::ChainHandle ch=edi.InsertChain("A");
    mol::ResidueHandle r=edi.AppendResidue(ch, "ALA");
    mol::AtomHandle a=edi.InsertAtom(r, "N", geom::Vec3(32.0, -128.0, -2.5));
    mol::ResidueHandle r2=edi.AppendResidue(ch, "GLY");
    mol::AtomHandle a2=edi.InsertAtom(r2, "N", geom::Vec3(35.0, -99.0, -10.5));
    conop::HeuristicProcessor heu_proc;
    heu_proc.Process(ent);
    writer.Write(ent);
  }
  BOOST_CHECK(compare_files("testfiles/pdb/charmm_ter.pdb",
                            "testfiles/pdb/charmm_ter-out.pdb"));
}

BOOST_AUTO_TEST_CASE(test_pqr_import_handler)
{
  String fname("testfiles/test_in.pqr");

  mol::EntityHandle eh=mol::CreateEntity();
  EntityIOPDBHandler pdbh;

  BOOST_CHECK(EntityIOPDBHandler::ProvidesImport(fname));
  BOOST_CHECK(EntityIOPDBHandler::ProvidesImport("test_in.PQR"));

  BOOST_CHECK(EntityIOPDBHandler::ProvidesExport(fname));
  BOOST_CHECK(EntityIOPDBHandler::ProvidesExport("test_in.PQR"));

  pdbh.Import(eh,"testfiles/pdb/simple.pqr");
}


BOOST_AUTO_TEST_CASE(test_pqr_read_atom)
{
  String fname("testfiles/pdb/simple.pqr");
  PDBReader reader(fname, IOProfile());
  mol::EntityHandle ent=mol::CreateEntity();
  reader.Import(ent);
  BOOST_REQUIRE_EQUAL(ent.GetChainCount(), 1);
  BOOST_REQUIRE_EQUAL(ent.GetResidueCount(), 5);
  BOOST_REQUIRE_EQUAL(ent.GetAtomCount(), 91);
  mol::AtomHandle a1=ent.FindAtom(" ", mol::ResNum(1), "N");
  BOOST_REQUIRE(a1.IsValid());
  BOOST_CHECK_EQUAL(a1.GetName(), "N");
  BOOST_CHECK_EQUAL(a1.GetResidue().GetName(), "MET");
  BOOST_CHECK_EQUAL(a1.GetResidue().GetChain().GetName(), " ");

  BOOST_CHECK(vec3_is_close(a1.GetPos(), geom::Vec3(21.6, 35.3, 56.7),Real(0.1)));
  BOOST_CHECK_EQUAL(a1.GetElement(), "");
  BOOST_CHECK_EQUAL(a1.IsHetAtom(), false);
  BOOST_CHECK_CLOSE(a1.GetCharge(), Real(-0.3755), Real(1e-4));
  BOOST_CHECK_CLOSE(a1.GetRadius(), Real(2.0005), Real(1e-4));

  mol::AtomHandle a2=ent.FindAtom(" ", mol::ResNum(2), "CZ");
  BOOST_REQUIRE(a2.IsValid());
  BOOST_CHECK_EQUAL(a2.GetName(), "CZ");
  BOOST_CHECK_EQUAL(a2.GetResidue().GetName(), "ARG");
  BOOST_CHECK_EQUAL(a2.GetResidue().GetChain().GetName(), " ");

  BOOST_CHECK(vec3_is_close(a2.GetPos(), geom::Vec3(23.9, 28.7, 56.5),Real(0.1)));
  BOOST_CHECK_EQUAL(a2.GetElement(), "");
  BOOST_CHECK_EQUAL(a2.IsHetAtom(), false);
  BOOST_CHECK_CLOSE(a2.GetCharge(), Real(0.2507), Real(1e-4));
  BOOST_CHECK_CLOSE(a2.GetRadius(), Real(1.7503), Real(1e-4));
}

BOOST_AUTO_TEST_CASE(checks_for_atom_pos_overflow)
{
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  writer.SetIsPQR(true);

  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  mol::ResidueHandle r=edi.AppendResidue(ch, "GLY");

  mol::AtomHandle a=edi.InsertAtom(r, "CA", geom::Vec3(0, -1000,0), "C");
  
  BOOST_CHECK_THROW(writer.Write(ent), IOException);
  edi.SetAtomPos(a, geom::Vec3(10000,0,0));
  BOOST_CHECK_THROW(writer.Write(ent), IOException);
}

BOOST_AUTO_TEST_CASE(test_pqr_write_atom)
{
  std::stringstream out;
  PDBWriter writer(out, IOProfile());
  writer.SetIsPQR(true);

  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  mol::ResidueHandle r=edi.AppendResidue(ch, "GLY");

  mol::AtomHandle a=edi.InsertAtom(r, "CA", geom::Vec3(32.0, -128.0, -2.5), "C");
  a.SetOccupancy(1.0);
  a.SetBFactor(128.0);
  a.SetCharge(-0.6543);
  a.SetRadius(1.2345);
  mol::AtomHandle a2=edi.InsertAtom(r, "CA", geom::Vec3(32.0, -128.0, -2.5), "C");
  a2.SetOccupancy(1.0);
  a2.SetBFactor(128.0);
  a2.SetCharge(0.1234);
  a2.SetRadius(2.5432);
  writer.Write(ent);
  String s=out.str();
  BOOST_CHECK_EQUAL(s.substr(0, 54),
                    "ATOM      1  CA  GLY A   1      32.000-128.000  -2.500");
  BOOST_CHECK_EQUAL(s.substr(54, 26),
                    " -0.6543 1.2345        C  ");

  PDBWriter fwriter(String("testfiles/pdb/pqr_atom-out.pqr"), IOProfile());
  BOOST_CHECK_EQUAL(fwriter.IsPQR(), true);
}

BOOST_AUTO_TEST_SUITE_END();
