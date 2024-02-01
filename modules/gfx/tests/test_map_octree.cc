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

#include <ost/gfx/impl/map_octree.hh>

#include <ost/img/image_handle.hh>
#include <ost/img/image_factory.hh>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


using boost::unit_test_framework::test_suite;
using namespace ost;
using namespace ost::gfx::impl;

namespace {

struct Pow2Vis {
  Pow2Vis():
    leaf_count(0), node_count(0)
  { }
  bool VisitNode(const OctreeNode& node, uint8_t level, 
                 const img::Extent& ext)
  {
    BOOST_CHECK_EQUAL(level, 0);
    BOOST_CHECK_EQUAL(node.GetMax(), 64.0f);
    BOOST_CHECK_EQUAL(node.GetMin(), 0.5f);    
    node_count++;
    return true;
  }
  void VisitLeaf(img::RealSpatialImageState* map, 
                 const img::Point& point) 
  {
    leaf_count++;
  }
  int leaf_count;
  int node_count;
};

struct NonPow2Vis {
  NonPow2Vis():
    leaf_count(0), node_count(0)
  { }
  bool VisitNode(const OctreeNode& node, uint8_t level, 
                 const img::Extent& ext)
  {
    node_count++;    
    if (level==0) {
      BOOST_CHECK_EQUAL(node.GetMax(), 64.0f);
      BOOST_CHECK_EQUAL(node.GetMin(), 0.0f);      
      return true;
    } 
    if (level==1) {
      switch(node_count) {
        case 2:
          BOOST_CHECK_EQUAL(ext, img::Extent(img::Point(0,0,0),
                                              img::Point(3,1,1)));
          break;
        case 5:
          BOOST_CHECK_EQUAL(ext, img::Extent(img::Point(4,0,0),
                                              img::Point(4,1,1)));
          break;
        default:
          BOOST_FAIL("unexpected node");
      }
    }
    return true;
  }
  void VisitLeaf(img::RealSpatialImageState* map, 
                 const img::Point& point) 
  {
    leaf_count++;
    // z-coordinate of pixels change fastest, followed by y and then x.
    switch (leaf_count) {
      case 1:
        BOOST_CHECK_EQUAL(point, img::Point(0,0,0));
        break;
      case 2:
        BOOST_CHECK_EQUAL(point, img::Point(0,0,1));
        break;
      case 3:
        BOOST_CHECK_EQUAL(point, img::Point(0,1,0));
        break;
      case 4:
        BOOST_CHECK_EQUAL(point, img::Point(0,1,1));
        break;
      case 5:
        BOOST_CHECK_EQUAL(point, img::Point(1,0,0));
        break;
      case 6:
        BOOST_CHECK_EQUAL(point, img::Point(1,0,1));
        break;
      case 7:
        BOOST_CHECK_EQUAL(point, img::Point(1,1,0));
        break;
      case 8:
        BOOST_CHECK_EQUAL(point, img::Point(1,1,1));
        break;
      case 9:
        BOOST_CHECK_EQUAL(point, img::Point(2,0,0));
        break;
      case 10:
        BOOST_CHECK_EQUAL(point, img::Point(2,0,1));
        break;
      case 11:
        BOOST_CHECK_EQUAL(point, img::Point(2,1,0));
        break;
      case 12:
        BOOST_CHECK_EQUAL(point, img::Point(2,1,1));
        break;
      case 13:
        BOOST_CHECK_EQUAL(point, img::Point(3,0,0));
        break;
      case 14:
        BOOST_CHECK_EQUAL(point, img::Point(3,0,1));
        break;
      case 15:
        BOOST_CHECK_EQUAL(point, img::Point(3,1,0));
        break;
      case 16:
        BOOST_CHECK_EQUAL(point, img::Point(3,1,1));
        break;
      case 17:
        BOOST_CHECK_EQUAL(point, img::Point(4,0,0));
        break;
      case 18:
        BOOST_CHECK_EQUAL(point, img::Point(4,0,1));
        break;
      case 19:
        BOOST_CHECK_EQUAL(point, img::Point(4,1,0));
        break;
      case 20:
        BOOST_CHECK_EQUAL(point, img::Point(4,1,1));
        break;
      default:
        BOOST_FAIL("should never get here");
    }
  }
  int leaf_count;
  int node_count;
};

}

BOOST_AUTO_TEST_SUITE(gfx)

BOOST_AUTO_TEST_CASE(octree_power_of_two) 
{
  img::ImageHandle img=img::CreateImage(img::Size(2, 2, 2));
  img.SetReal(img::Point(0, 0, 0), 1.0f);
  img.SetReal(img::Point(1, 0, 0), 2.0f);
  img.SetReal(img::Point(1, 1, 0), 4.0f);  
  img.SetReal(img::Point(1, 0, 1), 8.0f); 
  img.SetReal(img::Point(1, 1, 1), 16.0f);
  img.SetReal(img::Point(0, 1, 1), 32.0f);   
  img.SetReal(img::Point(0, 1, 0), 64.0f);
  img.SetReal(img::Point(0, 0, 1), 0.5f);
  Pow2Vis v;
  MapOctree octree(img);
  octree.Initialize();
  octree.VisitDF(v);  
  BOOST_CHECK_EQUAL(v.leaf_count, 8);
  BOOST_CHECK_EQUAL(v.node_count, 1);
}

BOOST_AUTO_TEST_CASE(octree_non_power_of_two)
{
  img::ImageHandle img=img::CreateImage(img::Size(5, 2, 2));
  img.SetReal(img::Point(0, 0, 0), 1.0f);
  img.SetReal(img::Point(1, 0, 0), 2.0f);
  img.SetReal(img::Point(1, 1, 0), 4.0f);  
  img.SetReal(img::Point(1, 0, 1), 8.0f); 
  img.SetReal(img::Point(1, 1, 1), 16.0f);
  img.SetReal(img::Point(0, 1, 1), 32.0f);   
  img.SetReal(img::Point(0, 1, 0), 64.0f);
  img.SetReal(img::Point(0, 0, 1), 0.5f);
  NonPow2Vis v;
  MapOctree octree(img);
  octree.Initialize();
  octree.VisitDF(v);  
}

BOOST_AUTO_TEST_SUITE_END()
