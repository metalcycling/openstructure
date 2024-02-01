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
#include <ost/geom/vec_mat_predicates.hh>
#include <ost/mol/chem_class.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/property_id.hh>
#include <cmath>

#define CHECK_TRANSFORMED_ATOM_POSITION(ATOM,TARGET) \
   BOOST_CHECK(vec3_is_close(ATOM.GetPos(), TARGET,Real(0.1)))

#define CHECK_ORIGINAL_ATOM_POSITION(ATOM,TARGET) \
   BOOST_CHECK(vec3_is_close(ATOM.GetOriginalPos(), TARGET,Real(0.1)))

#define CHECK_ALTERNATE_ATOM_POSITION(ATOM,TARGET,GROUP) \
   BOOST_CHECK(vec3_is_close(ATOM.GetAltPos(GROUP), TARGET,Real(0.1)))

using namespace ost;
using namespace ost::mol;


EntityHandle make_test_entity()
{
  EntityHandle eh = CreateEntity();
  XCSEditor e=eh.EditXCS();
  ChainHandle chain = e.InsertChain("A");
  ResidueHandle res1 = e.AppendResidue(chain, "MET");

  eh.SetName("TestEntity");

  e.InsertAtom(res1, "N",geom::Vec3(21.609,35.384,56.705));
  e.InsertAtom(res1, "CA",geom::Vec3(20.601,35.494,57.793));
  e.InsertAtom(res1, "C",geom::Vec3(19.654,34.300,57.789));
  e.InsertAtom(res1, "O",geom::Vec3(18.447,34.456,57.595));
  e.Connect(res1.FindAtom("N"), res1.FindAtom("CA"));
  e.Connect(res1.FindAtom("CA"), res1.FindAtom("C"));
  e.Connect(res1.FindAtom("C"), res1.FindAtom("O"));    
  ResidueHandle res2 = e.AppendResidue(chain, "ARG");
  e.InsertAtom(res2, "N",geom::Vec3(20.202,33.112,58.011));
  e.InsertAtom(res2, "CA",geom::Vec3(19.396,31.903,58.033));
  e.InsertAtom(res2, "C",geom::Vec3(18.608,31.739,59.328));
  e.InsertAtom(res2, "O",geom::Vec3(17.651,30.965,59.381));
  e.Connect(res1.FindAtom("C"), res2.FindAtom("N"));
  e.Connect(res2.FindAtom("N"), res2.FindAtom("CA"));
  e.Connect(res2.FindAtom("CA"), res2.FindAtom("C"));
  e.Connect(res2.FindAtom("C"), res2.FindAtom("O"));
  res1.SetChemClass(ChemClass(ChemClass::L_PEPTIDE_LINKING));
  res2.SetChemClass(ChemClass(ChemClass::L_PEPTIDE_LINKING));  
  e.AddTorsion("PHI", res1.FindAtom("C"), res2.FindAtom("N"), 
               res2.FindAtom("CA"), res2.FindAtom("C"));
  return eh;
}

BOOST_AUTO_TEST_SUITE( mol_base );


BOOST_AUTO_TEST_CASE(throw_invalid_ent_handle)
{
  EntityHandle ent;
  BOOST_CHECK_THROW(CheckHandleValidity(ent), InvalidHandle);  
  ent=CreateEntity();
  BOOST_CHECK_NO_THROW(CheckHandleValidity(ent));
}

BOOST_AUTO_TEST_CASE(throw_invalid_ent_view)
{
  EntityView ent_view;
  BOOST_CHECK_THROW(CheckHandleValidity(ent_view), InvalidHandle);  
  EntityHandle ent=CreateEntity();
  BOOST_CHECK_NO_THROW(CheckHandleValidity(ent.CreateFullView()));
  BOOST_CHECK_NO_THROW(CheckHandleValidity(ent.CreateEmptyView()));
}

BOOST_AUTO_TEST_CASE(bzdng250)
{
  EntityHandle eh=make_test_entity();
  eh.EditXCS().ApplyTransform(geom::Mat4(-1, 0, 0, 0,
                                          0, 1, 0, 0,
                                          0, 0,-1, 0,
                                          0, 0, 0, 1));
  BOOST_CHECK_NO_THROW(eh.FindResidue("A", 1).FindTorsion("PHI").GetAngle());
}

BOOST_AUTO_TEST_CASE(entity_creator) 
{
  EntityHandle eh = CreateEntity();
  eh.SetName("TestEntity");
  BOOST_CHECK(eh.GetName()=="TestEntity");

  XCSEditor e=eh.EditXCS();
  ChainHandle chain = e.InsertChain("C");
  BOOST_CHECK(chain.GetName()=="C");

  ResidueHandle res = e.AppendResidue(chain, "FOO", ResNum(13));
  BOOST_CHECK(chain.GetName()==res.GetChain().GetName());
  BOOST_CHECK(res.GetNumber()==ResNum(13));
  BOOST_CHECK(res.GetKey()=="FOO");

  AtomHandle atom1 = e.InsertAtom(res, "X1",geom::Vec3());
  BOOST_CHECK(res==atom1.GetResidue());
  BOOST_CHECK(atom1.GetName()=="X1");
  BOOST_CHECK(atom1.GetIndex()==0);

  AtomHandle atom2 = e.InsertAtom(res, "X2",geom::Vec3());
  BOOST_CHECK(res==atom2.GetResidue());
  BOOST_CHECK(atom2.GetName()=="X2");
  BOOST_CHECK(atom2.GetIndex()==1);

  BondHandle bond = e.Connect(atom1, atom2);

  BOOST_CHECK(eh.GetChainCount()==1);
  BOOST_CHECK(eh.GetResidueCount()==1);
  BOOST_CHECK(eh.GetAtomCount()==2);
  
  BOOST_CHECK(eh.GetAtomCount()==chain.GetAtomCount());
  BOOST_CHECK(eh.GetResidueCount()==chain.GetResidueCount());

  chain.SetIntProp("amazing",42);
  EntityHandle eh2 = CreateEntity();
  XCSEditor e2 = eh2.EditXCS();
  ChainHandle inserted_chain=e2.InsertChain("Q",chain);
  BOOST_CHECK(eh2.GetChainCount()==1);
  BOOST_CHECK(eh2.GetResidueCount()==0);
  BOOST_CHECK(eh2.GetAtomCount()==0);
  BOOST_CHECK(inserted_chain.HasProp("amazing"));
  BOOST_CHECK(eh2.FindChain("Q").IsValid());

  EntityHandle eh3 = CreateEntity();
  XCSEditor e3 = eh3.EditXCS();
  e3.InsertChain("Q",chain,true);
  BOOST_CHECK(eh3.GetResidueCount()==1);
  BOOST_CHECK(eh3.GetAtomCount()==2);

  EntityVisitor v;
  eh.Apply(v);
}

BOOST_AUTO_TEST_CASE(spatial_organizer) 
{
  EntityHandle eh=CreateEntity();
  XCSEditor e=eh.EditXCS();  
  ChainHandle chain = e.InsertChain("A");
  ResidueHandle res = e.AppendResidue(chain, "X");

  AtomHandle a1=e.InsertAtom(res, "X1",geom::Vec3(0.0,0.0,0.0));
  AtomHandle a2=e.InsertAtom(res, "X2",geom::Vec3(1.0,0.0,0.0));
  AtomHandle a3=e.InsertAtom(res, "X3",geom::Vec3(0.0,-1.0,0.0));
  AtomHandle a4=e.InsertAtom(res, "X4",geom::Vec3(17.0,3.0,-33.0));
  AtomHandle a5=e.InsertAtom(res, "X5",geom::Vec3(-2.0,3.0,1.0));

  AtomHandleList ahv = eh.FindWithin(geom::Vec3(0.0,0.0,0.0),5.0);
  BOOST_CHECK_EQUAL(ahv.size(), size_t(4));
  BOOST_CHECK_EQUAL(std::count(ahv.begin(), ahv.end(), a1), 1);
  BOOST_CHECK_EQUAL(std::count(ahv.begin(), ahv.end(), a2), 1);
  BOOST_CHECK_EQUAL(std::count(ahv.begin(), ahv.end(), a3), 1);
  BOOST_CHECK_EQUAL(std::count(ahv.begin(), ahv.end(), a5), 1);
}

BOOST_AUTO_TEST_CASE(transformation) 
{
  EntityHandle eh = CreateEntity();
  XCSEditor e=eh.EditXCS();
  ChainHandle chain = e.InsertChain("A");
  ResidueHandle res = e.AppendResidue(chain, "X");

  AtomHandle atom1 = e.InsertAtom(res, "X1",geom::Vec3(1.0,0.0,0.0));
  AtomHandle atom2 = e.InsertAtom(res, "X2",geom::Vec3(0.0,2.0,0.0));
  AtomHandle atom3 = e.InsertAtom(res, "X3",geom::Vec3(0.0,0.0,-3.0));

  AtomHandleList within_list1 = eh.FindWithin(atom1.GetPos(),2.5);
  BOOST_CHECK(within_list1.size()==2);
  BOOST_CHECK(within_list1[0]==atom1);
  BOOST_CHECK(within_list1[1]==atom2);

  BOOST_CHECK(eh.HasTransform()==false);

  geom::Transform trans;
  trans.ApplyZAxisRotation(90.0);
  e.ApplyTransform(trans);
  BOOST_CHECK(eh.HasTransform()==true);

  geom::Vec3 orig_atom1=geom::Vec3(1.0,0.0,0.0);
  geom::Vec3 orig_atom2=geom::Vec3(0.0,2.0,0.0);
  geom::Vec3 orig_atom3=geom::Vec3(0.0,0.0,-3.0);

  geom::Vec3 tr_atom1=geom::Vec3(0.0,-1.0,0.0);
  geom::Vec3 tr_atom2=geom::Vec3(2.0,0.0,0.0);
  geom::Vec3 tr_atom3=geom::Vec3(0.0,0.0,-3.0);

  CHECK_TRANSFORMED_ATOM_POSITION(atom1,tr_atom1);
  CHECK_ORIGINAL_ATOM_POSITION(atom1,orig_atom1);
  CHECK_TRANSFORMED_ATOM_POSITION(atom2,tr_atom2);
  CHECK_ORIGINAL_ATOM_POSITION(atom2,orig_atom2);
  CHECK_TRANSFORMED_ATOM_POSITION(atom3,tr_atom3);
  CHECK_ORIGINAL_ATOM_POSITION(atom3,orig_atom3);

  geom::Transform trans2;
  trans2.ApplyXAxisTranslation(3.5);
  e.ApplyTransform(trans2);

  BOOST_CHECK(eh.HasTransform()==true);

  tr_atom1=geom::Vec3(3.5,-1.0,0.0);
  tr_atom2=geom::Vec3(5.5,0.0,0.0);
  tr_atom3=geom::Vec3(3.5,0.0,-3.0);

  CHECK_TRANSFORMED_ATOM_POSITION(atom1,tr_atom1);
  CHECK_ORIGINAL_ATOM_POSITION(atom1,orig_atom1);
  CHECK_TRANSFORMED_ATOM_POSITION(atom2,tr_atom2);
  CHECK_ORIGINAL_ATOM_POSITION(atom2,orig_atom2);
  CHECK_TRANSFORMED_ATOM_POSITION(atom3,tr_atom3);
  CHECK_ORIGINAL_ATOM_POSITION(atom3,orig_atom3);

  AtomHandleList within_list2 = eh.FindWithin(atom1.GetPos(),2.5);
  BOOST_CHECK_EQUAL(within_list2.size(), size_t(2));
  BOOST_CHECK_EQUAL(std::count(within_list2.begin(), within_list2.end(), atom1), 1);
  BOOST_CHECK_EQUAL(std::count(within_list2.begin(), within_list2.end(), atom2), 1);  
  e.SetAtomPos(atom3,geom::Vec3(3.5,0.0,0.0));
  e.SetAtomPos(atom1,geom::Vec3(3.5,-0.9,0.0));
  tr_atom3=geom::Vec3(3.5,0.0,0.0);
  orig_atom3=geom::Vec3(0.0,0.0,0.0);
  tr_atom1=geom::Vec3(3.5,-0.9,0.0);
  orig_atom1=geom::Vec3(0.9,0.0,0.0);
  CHECK_TRANSFORMED_ATOM_POSITION(atom3,tr_atom3);
  CHECK_ORIGINAL_ATOM_POSITION(atom3,orig_atom3);
  CHECK_TRANSFORMED_ATOM_POSITION(atom1,tr_atom1);
  CHECK_ORIGINAL_ATOM_POSITION(atom1,orig_atom1);

  AtomHandleList within_list3 = eh.FindWithin(atom1.GetPos(),2.5);
  BOOST_CHECK(within_list3.size()==3);
  BOOST_CHECK_EQUAL(std::count(within_list3.begin(), within_list3.end(), atom1), 1);
  BOOST_CHECK_EQUAL(std::count(within_list3.begin(), within_list3.end(), atom2), 1);
  BOOST_CHECK_EQUAL(std::count(within_list3.begin(), within_list3.end(), atom3), 1);    

  geom::Vec3 new_orig1 = geom::Vec3(0.0,0.0,0.0);
  geom::Vec3 new_orig2 = geom::Vec3(1.5,0.0,0.0);
  geom::Vec3 new_orig3 = geom::Vec3(0.0,0.0,-2.0);

  e.SetAtomOriginalPos(atom1,new_orig1);
  e.SetAtomOriginalPos(atom2,new_orig2);
  e.SetAtomOriginalPos(atom3,new_orig3);

  geom::Mat4 identity;
  e.SetTransform(identity);
  //BOOST_CHECK(eh.HasTransform()==false);

  BondHandle bond1 = e.Connect(atom1,atom2);
  BondHandle bond2 = e.Connect(atom1,atom3);

  BOOST_CHECK(bond1.GetLength()==1.5);
  BOOST_CHECK(bond2.GetLength()==2.0);

  e.SetTransform(trans);
  
  //BOOST_CHECK(eh.HasTransform()==true);

  BOOST_CHECK(bond1.GetLength()==1.5);
  BOOST_CHECK(bond2.GetLength()==2.0);
  ICSEditor i = eh.EditICS();
  i.SetBondLength(bond1,1.0);
  i.SetBondLength(bond2,3.0);

  BOOST_CHECK(bond1.GetLength()==1.0);
  BOOST_CHECK(bond2.GetLength()==3.0);

  geom::Vec3 atom1_af = geom::Vec3(0.0,0.0,0.0);
  geom::Vec3 atom2_af = geom::Vec3(1.0,0.0,0.0);
  geom::Vec3 atom3_af = geom::Vec3(0.0,0.0,-3.0);

  CHECK_ORIGINAL_ATOM_POSITION(atom1,atom1_af);
  CHECK_ORIGINAL_ATOM_POSITION(atom2,atom2_af);
  CHECK_ORIGINAL_ATOM_POSITION(atom3,atom3_af);

  geom::Vec3 atom1_tr = geom::Vec3(0.0,0.0,0.0);
  geom::Vec3 atom2_tr = geom::Vec3(0.0,-1.0,0.0);
  geom::Vec3 atom3_tr = geom::Vec3(0.0,0.0,-3.0);

  CHECK_TRANSFORMED_ATOM_POSITION(atom1,atom1_tr);
  CHECK_TRANSFORMED_ATOM_POSITION(atom2,atom2_tr);
  CHECK_TRANSFORMED_ATOM_POSITION(atom3,atom3_tr);
  i.SetAngle(atom2,atom1,atom3,M_PI);

  atom1_af = geom::Vec3(0.0,0.0,0.0);
  atom2_af = geom::Vec3(1.0,0.0,0.0);
  atom3_af = geom::Vec3(-3.0,0.0,0.0);

  CHECK_ORIGINAL_ATOM_POSITION(atom1,atom1_af);
  CHECK_ORIGINAL_ATOM_POSITION(atom2,atom2_af);
  CHECK_ORIGINAL_ATOM_POSITION(atom3,atom3_af);

  atom1_tr = geom::Vec3(0.0,0.0,0.0);
  atom2_tr = geom::Vec3(0.0,-1.0,0.0);
  atom3_tr = geom::Vec3(0.0,3.0,0.0);

  CHECK_TRANSFORMED_ATOM_POSITION(atom1,atom1_tr);
  CHECK_TRANSFORMED_ATOM_POSITION(atom2,atom2_tr);
  CHECK_TRANSFORMED_ATOM_POSITION(atom3,atom3_tr);

  geom::Vec3 orig_atom4 = geom::Vec3(0.0,0.0,2.0);
  geom::Vec3 tr_atom4 = geom::Vec3(0.0,0.0,2.0);
  geom::Vec3 alt_atom4 = geom::Vec3(1.0,0.0,2.0);
  geom::Vec3 tr_alt_atom4 = geom::Vec3(0.0,-1.0,2.0);


  AtomHandle atom4 = e.InsertAltAtom(res, "X4","Set1",geom::Vec3(0.0,0.0,2.0));
  e.AddAltAtomPos("Set2",atom4, alt_atom4);

  CHECK_ORIGINAL_ATOM_POSITION(atom4,orig_atom4);
  CHECK_TRANSFORMED_ATOM_POSITION(atom4,tr_atom4);
  CHECK_ALTERNATE_ATOM_POSITION(atom4,alt_atom4,"Set2");

  res.SwitchAtomPos("Set2");

  CHECK_ORIGINAL_ATOM_POSITION(atom4,alt_atom4);
  CHECK_TRANSFORMED_ATOM_POSITION(atom4,tr_alt_atom4);
  CHECK_ALTERNATE_ATOM_POSITION(atom4,orig_atom4,"Set1");
}

BOOST_AUTO_TEST_CASE(copy) 
{
  EntityHandle eh=make_test_entity();
  EntityHandle cp=eh.Copy();
  BOOST_CHECK_EQUAL(cp.GetName(),eh.GetName());

  ResidueHandle x2=eh.FindResidue("A", mol::ResNum(2));
  BOOST_CHECK_EQUAL(cp.GetAtomCount(),eh.GetAtomCount());
  BOOST_CHECK_EQUAL(cp.GetResidueCount(),eh.GetResidueCount());
  BOOST_CHECK_EQUAL(cp.GetChainCount(),eh.GetChainCount());
  BOOST_CHECK_EQUAL(cp.GetBondCount(),eh.GetBondCount());
  ResidueHandle r1=cp.FindResidue("A", mol::ResNum(1));
  ResidueHandle r2=cp.FindResidue("A", mol::ResNum(2));  
  BOOST_CHECK(mol::BondExists(r1.FindAtom("N"), r1.FindAtom("CA")));
  BOOST_CHECK(mol::BondExists(r1.FindAtom("CA"), r1.FindAtom("C")));
  BOOST_CHECK(mol::BondExists(r1.FindAtom("C"), r1.FindAtom("O")));
  
  BOOST_CHECK(mol::BondExists(r1.FindAtom("C"), r2.FindAtom("N")));   
  
  BOOST_CHECK(mol::BondExists(r2.FindAtom("N"), r2.FindAtom("CA")));
  BOOST_CHECK(mol::BondExists(r2.FindAtom("CA"), r2.FindAtom("C")));
  BOOST_CHECK(mol::BondExists(r2.FindAtom("C"), r2.FindAtom("O")));
  BOOST_CHECK(r2.GetPhiTorsion());
}

BOOST_AUTO_TEST_CASE(copy_residue_props)
{
  EntityHandle ent=mol::CreateEntity();
  XCSEditor edi=ent.EditXCS();
  ChainHandle ch=edi.InsertChain("A");
  ResidueHandle res=edi.AppendResidue(ch, "DUMMY", mol::ResNum(666, '6'));
  res.SetOneLetterCode('X');
  res.SetIsProtein(true);
  res.SetIsLigand(true);
  ChemClass cl(ChemClass::L_PEPTIDE_LINKING);  
  res.SetSecStructure(SecStructure(SecStructure::ALPHA_HELIX));
  res.SetChemClass(cl);
  EntityHandle copy=ent.Copy();
  ResidueHandle res2=copy.GetResidueList()[0];
  BOOST_CHECK_EQUAL(res2.GetName(), String("DUMMY"));
  BOOST_CHECK_EQUAL(res2.GetOneLetterCode(), 'X');

  BOOST_CHECK_EQUAL(res2.GetChemClass(), cl);
  BOOST_CHECK_EQUAL(res.GetSecStructure(), SecStructure::ALPHA_HELIX);
  BOOST_CHECK(res2.IsProtein()==true);
  BOOST_CHECK(res2.IsLigand()==true);  
  BOOST_CHECK_EQUAL(res2.GetNumber(), ResNum(666, '6'));
}

BOOST_AUTO_TEST_CASE(copy_atom_props)
{
  EntityHandle ent=mol::CreateEntity();
  XCSEditor edi=ent.EditXCS();
  ChainHandle ch=edi.InsertChain("A");
  ResidueHandle res=edi.AppendResidue(ch, "DUMMY");
  AtomHandle   atom=edi.InsertAtom(res, "X", geom::Vec3(1,2,3), "C");
  atom.SetMass(100.0);
  atom.SetBFactor(200.0);
  atom.SetOccupancy(300.0);
  atom.SetHetAtom(true);
  atom.SetRadius(500);
  atom.SetCharge(800);
  atom.SetAnisou(geom::Mat3(100,200,300));
  EntityHandle copy=ent.Copy();
  AtomHandle atom2=copy.GetAtomList()[0];
  BOOST_CHECK_EQUAL(atom2.GetPos(), geom::Vec3(1,2,3));
  BOOST_CHECK_EQUAL(atom2.GetName(), String("X"));
  BOOST_CHECK_EQUAL(atom2.GetElement(), String("C"));
  BOOST_CHECK_EQUAL(atom2.GetMass(), Real(100.0));
  BOOST_CHECK_EQUAL(atom2.GetCharge(), Real(800.0));  
  BOOST_CHECK_EQUAL(atom2.GetBFactor(), Real(200.0)); 
  BOOST_CHECK_EQUAL(atom2.IsHetAtom(), true);
  BOOST_CHECK_EQUAL(atom2.GetAnisou(), geom::Mat3(100,200,300));
  BOOST_CHECK_EQUAL(atom2.GetRadius(), Real(500.0));
}

BOOST_AUTO_TEST_CASE(rename_atom)
{
   EntityHandle ent=CreateEntity();
   XCSEditor edi=ent.EditXCS();
   ChainHandle ch1=edi.InsertChain("A");
   ResidueHandle res = edi.AppendResidue(ch1, "A");
   AtomHandle   atom=edi.InsertAtom(res, "A", geom::Vec3(1,2,3), "C");
   edi.RenameAtom(atom, "B");
   BOOST_CHECK_EQUAL(atom.GetName(), "B");
}

BOOST_AUTO_TEST_CASE(minmax)
{
  EntityHandle eh=make_test_entity();
  EntityView ev = eh.CreateFullView();
  mol::AtomViewList avl = ev.GetAtomList();
  mol::AtomViewList::iterator i;
  int n=0.0;
  for (i=avl.begin(); i!=avl.end(); ++i, ++n) {
    i->SetFloatProp("test", n);
  }
  std::pair<float,float> minmax = ev.GetMinMax("test", Prop::ATOM);
  BOOST_CHECK_EQUAL(minmax.first, 0.0);
  BOOST_CHECK_EQUAL(minmax.second, 7.0);
}

BOOST_AUTO_TEST_CASE(prune) {
  EntityHandle ent=mol::CreateEntity();
  XCSEditor edi = ent.EditXCS();

  ChainHandle ch1 = edi.InsertChain("A");
  ChainHandle ch2 = edi.InsertChain("B");

  ResidueHandle res1_1 = edi.AppendResidue(ch1, "DUMMY");
  ResidueHandle res1_2 = edi.AppendResidue(ch1, "DUMMY");
  AtomHandle   atom1_1_1 = edi.InsertAtom(res1_1, "X", geom::Vec3(1,2,3), "C");
  AtomHandle   atom1_2_1 = edi.InsertAtom(res1_2, "X", geom::Vec3(1,2,3), "C");

  ResidueHandle res2_1 = edi.AppendResidue(ch2, "DUMMY");
  ResidueHandle res2_2 = edi.AppendResidue(ch2, "DUMMY");
  AtomHandle   atom2_1_1 = edi.InsertAtom(res2_1, "X", geom::Vec3(1,2,3), "C");
  AtomHandle   atom2_2_1 = edi.InsertAtom(res2_2, "X", geom::Vec3(1,2,3), "C");

  BOOST_CHECK_EQUAL(ent.GetChainCount(), 2);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 4);
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 4);

  // shouldn't do anything
  edi.Prune();
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 2);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 4);
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 4);

  // remove one atom, the subsequent prune call should then remove the
  // according residue but the number of chains should remain the same
  edi.DeleteAtom(atom1_1_1);
  edi.Prune();
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 2);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 3);
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 3);

  // remove the second atom from the first chain, the whole chain should then
  // be removed by Prune
  edi.DeleteAtom(atom1_2_1);
  edi.Prune();
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 1);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 2);
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 2);

  // remove both residues from the second chain. The entity should be empty
  // afterwards
  edi.DeleteResidue(res2_1);
  edi.DeleteResidue(res2_2);
  edi.Prune();
  BOOST_CHECK_EQUAL(ent.GetChainCount(), 0);
  BOOST_CHECK_EQUAL(ent.GetResidueCount(), 0);
  BOOST_CHECK_EQUAL(ent.GetAtomCount(), 0);
}

BOOST_AUTO_TEST_SUITE_END();
