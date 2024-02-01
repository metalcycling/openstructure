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
#include <ost/mol/coord_group.hh>

using namespace ost;
using namespace ost::mol;

BOOST_AUTO_TEST_SUITE( io );

BOOST_AUTO_TEST_CASE(coord_group) 
{
  EntityHandle e=CreateEntity();
  XCSEditor editor=e.EditXCS();  
  ChainHandle c=editor.InsertChain("A");
  ResidueHandle r=editor.AppendResidue(c, "XXX");
  AtomHandle aa=editor.InsertAtom(r, "A", geom::Vec3(0,1,2));
  AtomHandle ab=editor.InsertAtom(r, "B", geom::Vec3(3,4,5));
  AtomHandle ac=editor.InsertAtom(r, "C", geom::Vec3(6,7,8));
  AtomHandle ad=editor.InsertAtom(r, "D", geom::Vec3(9,10,11));
  editor.UpdateICS();
  std::vector<AtomHandle> alist;
  alist.push_back(aa);
  alist.push_back(ab);
  alist.push_back(ac);

  CoordGroupHandle cg = CreateCoordGroup(alist);

  BOOST_CHECK(cg.GetAtomCount()==3);
  BOOST_CHECK(cg.GetFrameCount()==0);

  std::vector<geom::Vec3> clist(3);
  clist[0]=geom::Vec3(-1,-2,-3);
  clist[1]=geom::Vec3(-4,-5,-6);
  clist[2]=geom::Vec3(-7,-8,-9);
  cg.AddFrame(clist);
  BOOST_CHECK(cg.GetFrameCount()==1);  
  clist[0]=geom::Vec3(-11,-12,-13);
  clist[1]=geom::Vec3(-14,-15,-16);
  clist[2]=geom::Vec3(-17,-18,-19);
  cg.AddFrame(clist);
  BOOST_CHECK(cg.GetFrameCount()==2);

  BOOST_CHECK(aa.GetPos()==geom::Vec3(0,1,2));
  BOOST_CHECK(ab.GetPos()==geom::Vec3(3,4,5));
  BOOST_CHECK(ac.GetPos()==geom::Vec3(6,7,8));
  BOOST_CHECK(ad.GetPos()==geom::Vec3(9,10,11));

  cg.CopyFrame(0);
  BOOST_CHECK(aa.GetPos()==geom::Vec3(-1,-2,-3));
  BOOST_CHECK(ab.GetPos()==geom::Vec3(-4,-5,-6));
  BOOST_CHECK(ac.GetPos()==geom::Vec3(-7,-8,-9));
  BOOST_CHECK(ad.GetPos()==geom::Vec3(9,10,11));

  cg.CopyFrame(1);
  BOOST_CHECK(aa.GetPos()==geom::Vec3(-11,-12,-13));
  BOOST_CHECK(ab.GetPos()==geom::Vec3(-14,-15,-16));
  BOOST_CHECK(ac.GetPos()==geom::Vec3(-17,-18,-19));
  BOOST_CHECK(ad.GetPos()==geom::Vec3(9,10,11));

  geom::Transform tf;
  tf.ApplyXAxisRotation(17.0);
  tf.ApplyYAxisRotation(-135.0);
  tf.ApplyZAxisRotation(234.0);
  tf.SetCenter(geom::Vec3(-7.3,1.2,5.5));
  tf.SetTrans(geom::Vec3(14.5,-87.1,22.2));

  cg.ApplyTransform(tf);
  for(size_t i=0;i<3;++i) {
    BOOST_CHECK_CLOSE(geom::Distance(cg.GetAtomPos(1,alist[i]),tf.Apply(clist[i])),Real(0),1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END();
