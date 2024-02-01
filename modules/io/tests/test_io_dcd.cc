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

#include <ost/io/mol/dcd_io.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/mol/coord_group.hh>

using namespace ost;
using namespace ost::io;

#include <boost/random.hpp>
namespace {
  boost::mt19937 RandomGenerator(time(NULL));
  boost::uniform_01<boost::mt19937> UniformRandom(RandomGenerator);
}

BOOST_AUTO_TEST_SUITE( io );

BOOST_AUTO_TEST_CASE(test_io_dcd_charmm_frames)
{
  mol::EntityHandle eh=mol::CreateEntity();
  mol::XCSEditor ed=eh.EditXCS();
  mol::ChainHandle chain=ed.InsertChain("A");
  mol::ResidueHandle res=ed.AppendResidue(chain,mol::ResidueKey("UNK"));

  static unsigned int natoms=13;

  // create atoms and frame coords at the same time
  mol::AtomHandleList atoms(natoms);
  geom::Vec3List atom_pos(natoms);
  std::ostringstream aname;
  geom::Vec3 apos1,apos2;
  for(size_t i=0;i<natoms;++i) {
    for(size_t k=0;k<3;++k) {
      apos1[k]=UniformRandom();
      apos2[k]=UniformRandom();
    }
    aname.str("");
    aname << "X" << i;
    mol::AtomHandle atom=ed.InsertAtom(res,aname.str(),apos1);
    atom_pos[i]=apos2;
    atoms[i]=atom;
  }
  geom::Vec3 cell_size(UniformRandom(),UniformRandom(),UniformRandom());
  geom::Vec3 cell_angles(M_PI*UniformRandom(),M_PI*UniformRandom(),M_PI*UniformRandom());

  mol::CoordGroupHandle cg=mol::CreateCoordGroup(atoms);
  cg.AddFrame(atom_pos,cell_size,cell_angles);

  SaveCHARMMTraj(cg,"test_io_dcd_out.pdb","test_io_dcd_out.dcd");
  mol::CoordGroupHandle cg2=LoadCHARMMTraj(eh,"test_io_dcd_out.dcd");
  BOOST_CHECK_EQUAL(cg2.GetAtomCount(),natoms);
  BOOST_CHECK_EQUAL(cg2.GetFrameCount(),uint(1));

  mol::CoordFramePtr cf2 = cg2.GetFrame(0);
  BOOST_CHECK(geom::Distance(cf2->GetCellSize(),cell_size)<1e-5);
  BOOST_CHECK(geom::Distance(cf2->GetCellAngles(),cell_angles)<1e-5);

  geom::Vec3List atom_pos2=cg2.GetFramePositions(0);

  for(size_t i=0;i<natoms;++i) {
    BOOST_CHECK(geom::Distance(atom_pos[i],atom_pos2[i])<1e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END();
