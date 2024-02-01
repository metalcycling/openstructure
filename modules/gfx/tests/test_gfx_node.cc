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
  Author: Marco Biasini
 */
#include <ost/gfx/gfx_node.hh>
#include <ost/gfx/gfx_object.hh>
#include <ost/gfx/scene.hh>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


using boost::unit_test_framework::test_suite;
using namespace ost;
using namespace ost::gfx;


struct Observer : public SceneObserver {
  Observer(): added_count(0),removed_count(0) {}
  
  virtual void NodeAdded(const GfxNodeP& node)
  {
    added_count+=1;
  }
  
  virtual void NodeRemoved(const GfxNodeP& node)
  {
    removed_count+=1;
  }
  
  int added_count;
  int removed_count;
};

BOOST_AUTO_TEST_SUITE(gfx_node)

BOOST_AUTO_TEST_CASE(gfx_node_add) 
{
  GfxNodeP n1(new GfxNode("1"));
  BOOST_CHECK_EQUAL(n1->GetParent(), GfxNodeP());
  BOOST_CHECK_EQUAL(n1->GetChildCount(), size_t(0));
  GfxNodeP n2(new GfxNode("2"));
  n1->Add(n2);
  BOOST_CHECK_EQUAL(n1->GetParent(), GfxNodeP());
  BOOST_CHECK_EQUAL(n1->GetChildCount(), size_t(1));   
  BOOST_CHECK_EQUAL(n2->GetParent(), n1);
  BOOST_CHECK_EQUAL(n2->GetChildCount(), size_t(0));
  
  // "move" node 2 from 1 to 3
  GfxNodeP n3(new GfxNode("3"));
  n3->Add(n2);
  BOOST_CHECK_EQUAL(n1->GetChildCount(), size_t(0));
  BOOST_CHECK_EQUAL(n2->GetParent(), n3);
  BOOST_CHECK_EQUAL(n3->GetChildCount(), size_t(1));
}

BOOST_AUTO_TEST_CASE(gfx_node_add_duplicate) 
{
  GfxNodeP n1(new GfxNode("1"));
  GfxNodeP n2(new GfxNode("2"));
  GfxNodeP n3(new GfxNode("2"));
  BOOST_CHECK_NO_THROW(n1->Add(n2));
  BOOST_CHECK_THROW(n1->Add(n3), Error);
}

BOOST_AUTO_TEST_CASE(gfx_node_remove) 
{
  GfxNodeP n1(new GfxNode("1"));
  BOOST_CHECK_EQUAL(n1->GetParent(), GfxNodeP());
  BOOST_CHECK_EQUAL(n1->GetChildCount(), size_t(0));
  GfxNodeP n2(new GfxNode("2"));
  n1->Add(n2);
  n1->Remove(n2);
  BOOST_CHECK_EQUAL(n1->GetChildCount(), size_t(0));
  BOOST_CHECK_EQUAL(n2->GetParent(), GfxNodeP());
}

BOOST_AUTO_TEST_CASE(gfx_node_remove_all) 
{
  GfxNodeP n1(new GfxNode("1"));
  BOOST_CHECK_EQUAL(n1->GetParent(), GfxNodeP());
  BOOST_CHECK_EQUAL(n1->GetChildCount(), size_t(0));
  GfxNodeP n2(new GfxNode("2"));
  GfxNodeP n3(new GfxNode("3"));
  GfxNodeP n4(new GfxNode("4"));
  n1->Add(n2);
  n1->Add(n3);
  n1->Add(n4);
  n1->RemoveAll();
  BOOST_CHECK_EQUAL(n1->GetChildCount(), size_t(0));
  BOOST_CHECK_EQUAL(n2->GetParent(), GfxNodeP());
}

BOOST_AUTO_TEST_CASE(is_attached_to_scene) 
{
  Scene::Instance().RemoveAll();  
  GfxNodeP n1(new GfxNode("1"));
  GfxNodeP n2(new GfxNode("2"));
  n1->Add(n2);
  BOOST_CHECK(!n1->IsAttachedToScene());
  BOOST_CHECK(!n2->IsAttachedToScene());
  
  Scene::Instance().Add(n1);
  BOOST_CHECK(n1->IsAttachedToScene());
  BOOST_CHECK(n2->IsAttachedToScene());
  
  // root node is always connected to scene!
  BOOST_CHECK(Scene::Instance().GetRootNode()->IsAttachedToScene());
}

BOOST_AUTO_TEST_CASE(observe_added_removed)
{
  Observer o1;
  Scene::Instance().RemoveAll();
  Scene::Instance().AttachObserver(&o1);
  GfxNodeP n1(new GfxNode("1"));
  GfxNodeP n2(new GfxNode("2"));
  GfxNodeP n3(new GfxNode("3"));
  n1->Add(n2);
  BOOST_CHECK_EQUAL(o1.added_count, 0);
  BOOST_CHECK_EQUAL(o1.removed_count, 0);
  Scene::Instance().Add(n1);
  BOOST_CHECK_EQUAL(o1.added_count, 1);
  BOOST_CHECK_EQUAL(o1.removed_count, 0);
  n2->Add(n3);
  BOOST_CHECK_EQUAL(o1.added_count, 2);
  BOOST_CHECK_EQUAL(o1.removed_count, 0);
  n2->Remove(n3);
  BOOST_CHECK_EQUAL(o1.added_count, 2);
  BOOST_CHECK_EQUAL(o1.removed_count, 1);
}



BOOST_AUTO_TEST_SUITE_END()
