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
#include <iostream>

using namespace ost;
using namespace ost::mol;

template <typename HANDLE, typename ITER>
void test_res_handle_iterator_b(HANDLE ch) {
  ITER r1=ch.ResiduesBegin();
  ITER r2=ch.ResiduesEnd();  
  BOOST_CHECK(r1!=r2);    
  BOOST_CHECK_EQUAL((*r1).GetName(), "D");
  
  ++r1;  
  BOOST_CHECK(r1!=r2);  
  BOOST_CHECK_EQUAL((*r1).GetName(), "E");
  
  ++r1;  
  BOOST_CHECK(r1!=r2);  
  BOOST_CHECK_EQUAL((*r1).GetName(), "F");     
  ++r1; 
  BOOST_CHECK(r1==r2);
}

template <typename HANDLE, typename ITER>
void test_res_handle_iterator_a(HANDLE eh) {
  ITER r1=eh.ResiduesBegin();
  ITER r2=eh.ResiduesEnd();  
  BOOST_CHECK(r1!=r2);
  BOOST_CHECK_EQUAL((*r1).GetName(), "A");
  
  ++r1;
  BOOST_CHECK(r1!=r2);  
  BOOST_CHECK_EQUAL((*r1).GetName(), "B");
  
  ++r1;
  BOOST_CHECK(r1!=r2);  
  BOOST_CHECK_EQUAL((*r1).GetName(), "C");
  
  ++r1;  
  BOOST_CHECK(r1!=r2);    
  BOOST_CHECK_EQUAL((*r1).GetName(), "D");
  
  ++r1;  
  BOOST_CHECK(r1!=r2);  
  BOOST_CHECK_EQUAL((*r1).GetName(), "E");
  
  ++r1;  
  BOOST_CHECK(r1!=r2);  
  BOOST_CHECK_EQUAL((*r1).GetName(), "F");      
  
  ++r1;
  BOOST_CHECK(r1==r2);
}

template <typename HANDLE, typename ITER>
void test_atom_handle_iterator_c(HANDLE en) {
  ITER a1=en.AtomsBegin();
  ITER a2=en.AtomsEnd();  
  BOOST_CHECK(a1!=a2);
  BOOST_CHECK_EQUAL((*a1).GetName(), "X");
  ++a1;
  BOOST_CHECK(a1!=a2);  
  BOOST_CHECK_EQUAL((*a1).GetName(), "Y");
  ++a1;
  BOOST_CHECK(a1!=a2);  
  BOOST_CHECK_EQUAL((*a1).GetName(), "U");
  ++a1;  
  BOOST_CHECK(a1!=a2);    
  BOOST_CHECK_EQUAL((*a1).GetName(), "V");
  ++a1;
  BOOST_CHECK(a1!=a2);  
  BOOST_CHECK_EQUAL((*a1).GetName(), "M");
  ++a1;  
  BOOST_CHECK(a1!=a2);    
  BOOST_CHECK_EQUAL((*a1).GetName(), "N");  
  ++a1;
  BOOST_CHECK(a1==a2);  
}

template <typename HANDLE, typename ITER>
void test_atom_handle_iterator_a(HANDLE ch) {
  ITER a1=ch.AtomsBegin();
  ITER a2=ch.AtomsEnd();  
  BOOST_CHECK(a1!=a2);
  BOOST_CHECK_EQUAL((*a1).GetName(), "X");
  ++a1;
  BOOST_CHECK(a1!=a2);  
  BOOST_CHECK_EQUAL((*a1).GetName(), "Y");
  ++a1;
  BOOST_CHECK(a1!=a2);  
  BOOST_CHECK_EQUAL((*a1).GetName(), "U");
  ++a1;  
  BOOST_CHECK(a1!=a2);    
  BOOST_CHECK_EQUAL((*a1).GetName(), "V");
  ++a1;
  BOOST_CHECK(a1==a2);
}

template <typename HANDLE, typename ITER>
void test_atom_handle_iterator_b(HANDLE rr) {
  ITER a1=rr.AtomsBegin();
  ITER a2=rr.AtomsEnd();  
  BOOST_CHECK(a1!=a2);
  BOOST_CHECK_EQUAL((*a1).GetName(), "X");
  ++a1;
  BOOST_CHECK(a1!=a2);  
  BOOST_CHECK_EQUAL((*a1).GetName(), "Y");
  
  ++a1;
  BOOST_CHECK(a1==a2); 
}

BOOST_AUTO_TEST_SUITE( mol_base );

BOOST_AUTO_TEST_CASE(chain_handle_iterator) 
{
  EntityHandle ent=CreateEntity();
  XCSEditor editor=ent.EditXCS();  
  ChainHandle ch1=editor.InsertChain("A");
  ChainHandle ch2=editor.InsertChain("B");
  ChainHandle ch3=editor.InsertChain("C");
  ChainHandleIter ch_it=ent.ChainsBegin();
  ChainHandleIter ch_end=ent.ChainsEnd();
  BOOST_CHECK(ch_it!=ch_end);
  BOOST_CHECK_EQUAL((*ch_it).GetName(), "A");
  ++ch_it;
  BOOST_CHECK(ch_it!=ch_end);  
  BOOST_CHECK_EQUAL((*ch_it).GetName(), "B");
  ++ch_it;
  BOOST_CHECK(ch_it!=ch_end);  
  BOOST_CHECK_EQUAL((*ch_it).GetName(), "C");
  ++ch_it;
  BOOST_CHECK(ch_it==ch_end);
}

BOOST_AUTO_TEST_CASE(test_res_handle_iterator) 
{
  EntityHandle eh=CreateEntity();
  XCSEditor editor=eh.EditXCS();  
  ChainHandle ch1=editor.InsertChain("A");
  editor.AppendResidue(ch1, "A");
  editor.AppendResidue(ch1, "B");
  editor.AppendResidue(ch1, "C");  
  
  // make sure the iterator knows how to deal with empty chains
  ChainHandle ch3=editor.InsertChain("B");
  
  ChainHandle ch2=editor.InsertChain("C");
  editor.AppendResidue(ch2, "D");
  editor.AppendResidue(ch2, "E");
  editor.AppendResidue(ch2, "F");  
  EntityView v1=eh.CreateFullView();
  test_res_handle_iterator_b<ChainHandle, ResidueHandleIter>(ch2);
  test_res_handle_iterator_a<EntityHandle, ResidueHandleIter>(eh);
  test_res_handle_iterator_a<EntityView, ResidueViewIter>(v1);  
  // make the sure the iterator knows how to deal with an empty chain at 
  // the end
  ChainHandle ch4=editor.InsertChain("D");
  test_res_handle_iterator_a<EntityHandle, ResidueHandleIter>(eh);  
  EntityView v2=eh.CreateFullView();  
  test_res_handle_iterator_a<EntityView, ResidueViewIter>(v2);    
  test_res_handle_iterator_b<ChainHandle, ResidueHandleIter>(ch2);    
  editor.AppendResidue(ch4, "D");
  editor.AppendResidue(ch4, "E");
  editor.AppendResidue(ch4, "F");  
  test_res_handle_iterator_b<ChainHandle, ResidueHandleIter>(ch2);
}

BOOST_AUTO_TEST_CASE(test_empty_iterators) 
{
  EntityHandle eh=CreateEntity();
  BOOST_CHECK(eh.ChainsBegin()==eh.ChainsEnd());
  BOOST_CHECK(eh.ResiduesBegin()==eh.ResiduesEnd());  
  BOOST_CHECK(eh.AtomsBegin()==eh.AtomsEnd());    
  XCSEditor editor=eh.EditXCS();
  ChainHandle ch1=editor.InsertChain("A");  
  BOOST_CHECK(ch1.ResiduesBegin()==ch1.ResiduesEnd());
  BOOST_CHECK(ch1.AtomsBegin()==ch1.AtomsEnd());  
  ResidueHandle r1=editor.AppendResidue(ch1, "A");
  BOOST_CHECK(r1.AtomsBegin()==r1.AtomsEnd());    
  
  eh=CreateEntity();
  editor=eh.EditXCS();  
  ch1=editor.InsertChain("A");
  r1=editor.AppendResidue(ch1, "A");
  ResidueHandle r2=editor.AppendResidue(ch1, "B");
  editor.InsertAtom(r1, "X", geom::Vec3());
  editor.InsertAtom(r1, "Y", geom::Vec3()); 
  editor.InsertAtom(r2, "U", geom::Vec3());
  editor.InsertAtom(r2, "V", geom::Vec3());
  
  // testcase for ticket #120
  EntityView ev=eh.CreateFullView();
  AtomViewList avl=ev.GetAtomList();
  for (AtomViewList::const_iterator i=avl.begin(), e=avl.end(); i!=e; ++i) {
    ev.RemoveAtom(*i);
  }
  BOOST_CHECK(ev.AtomsBegin()==ev.AtomsEnd());
  BOOST_CHECK(ev.GetAtomList().empty());
  //BOOST_CHECK(ev.FindChain("A").AtomsBegin()==ev.FindChain("A").AtomsEnd());
}


BOOST_AUTO_TEST_CASE(test_atom_handle_iterator) 
{
  EntityHandle eh=CreateEntity();
  XCSEditor editor=eh.EditXCS();  
  ChainHandle ch1=editor.InsertChain("A");
  ResidueHandle r1=editor.AppendResidue(ch1, "A");
  ResidueHandle r2=editor.AppendResidue(ch1, "B");

  editor.InsertAtom(r1, "X", geom::Vec3());
  editor.InsertAtom(r1, "Y", geom::Vec3()); 
  editor.InsertAtom(r2, "U", geom::Vec3());
  editor.InsertAtom(r2, "V", geom::Vec3());
  
  ChainHandle ch2=editor.InsertChain("B");  
  ResidueHandle r4=editor.AppendResidue(ch2, "D");  
  editor.InsertAtom(r4, "M", geom::Vec3());
  editor.InsertAtom(r4, "N", geom::Vec3());  
  test_atom_handle_iterator_a<ChainHandle, AtomHandleIter>(ch1);
  test_atom_handle_iterator_b<ResidueHandle, AtomHandleIter>(r1);

  test_atom_handle_iterator_c<EntityHandle, AtomHandleIter>(eh);
  editor.InsertChain("C");
  // make sure the iterator knows how to deal with an empty chain at the end  
  EntityView v=eh.CreateFullView();
  test_atom_handle_iterator_c<EntityHandle, AtomHandleIter>(eh);  

  test_atom_handle_iterator_c<EntityView, AtomViewIter>(v);
}
BOOST_AUTO_TEST_SUITE_END();
