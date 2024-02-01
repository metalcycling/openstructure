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

#include <ost/mol/mm/observer.hh>
#include <OpenMM.h>

namespace ost { namespace mol{ namespace mm{

void TrajObserver::Init(boost::shared_ptr<OpenMM::Context> c, 
                        TopologyPtr top,
                        ost::mol::EntityHandle& ent){

  if(registered_){
    throw ost::Error("Can register observer to only one simulation!");
  }

  registered_ = true;

  context_ = c;
  c_group_ = ost::mol::CreateCoordGroup(ent.GetAtomList());
}

void TrajObserver::Notify(){

  geom::Vec3List positions;
  StateExtractor::ExtractPositions(context_,positions,true);
  c_group_.AddFrame(positions);
}


void TrajWriter::Init(boost::shared_ptr<OpenMM::Context> c, 
                      TopologyPtr top,
                      ost::mol::EntityHandle& ent){

  if(registered_){
    throw ost::Error("Can register observer to only one simulation!");
  }

  registered_ = true;

  context_ = c;
  ost::io::IOProfile profile("CHARMM",false,false,false,false,false,false);
  ost::io::PDBWriter writer(pdb_filename_, profile);
  writer.Write(ent.GetAtomList());

  stream_.open(dcd_filename_.c_str(), std::ios::binary);

  x.resize(top->GetNumParticles());
  y.resize(top->GetNumParticles());
  z.resize(top->GetNumParticles());


  // size of first header block in bytes
  int32_t magic_number=84;
  stream_.write(reinterpret_cast<char*>(&magic_number), 4);

  // magic string
  char crd[]={'C', 'O', 'R', 'D'};
  stream_.write(crd, 4);

  // icntrl[0], NSET, number of frames
  // we don't know the number of frames yet, will be written
  // when the finalize function is called.
  int32_t zero=0;
  stream_.write(reinterpret_cast<char*>(&zero), 4);

  // icntrl[1], ISTART, starting timestep
  stream_.write(reinterpret_cast<char*>(&zero), 4);

  // icntrl[2], NSAVC, timesteps between DCD saves
  int32_t one=1;
  stream_.write(reinterpret_cast<char*>(&one), 4);

  // icntrl[3] to icntrl[7], unused
  for (int i=3; i<=7; ++i) {
    stream_.write(reinterpret_cast<char*>(&zero), 4);
  }

  // icntrl[8], NAMNF, number of fixed atoms
  stream_.write(reinterpret_cast<char*>(&zero), 4);

  // icntrl[9], DELTA, timestep as float for CHARMM format
  float delta=1.0;
  stream_.write(reinterpret_cast<char*>(&delta), 4);

  // icntrl[10], CHARMM format: ucell per frame
  stream_.write(reinterpret_cast<char*>(&one), 4);

  // icntrl[11] to icntrl[18], unused
  for (int i=11; i<=18; ++i) {
    stream_.write(reinterpret_cast<char*>(&zero), 4);
  }

  // icntrl[19], charmm version
  int32_t charmm_version=24;
  stream_.write(reinterpret_cast<char*>(&charmm_version), 4);
  // bracket first header block
  stream_.write(reinterpret_cast<char*>(&magic_number), 4);

  //  no titles in title block
  int32_t four=4;
  stream_.write(reinterpret_cast<char*>(&four), 4);
  stream_.write(reinterpret_cast<char*>(&zero), 4);
  stream_.write(reinterpret_cast<char*>(&four), 4);

  // atom count block
  stream_.write(reinterpret_cast<char*>(&four), 4);
  int32_t atom_count=top->GetNumParticles();  
  stream_.write(reinterpret_cast<char*>(&atom_count), 4);
  stream_.write(reinterpret_cast<char*>(&four), 4);
}

void TrajWriter::Notify(){

  if(!registered_){
    throw ost::Error("Trajectory writer is not registered to any simulation!");
  }

  if(!stream_){
    throw ost::Error("Trajectory writer does not contain open filestream. Did you already call the finalize funciton before?");
  }

  geom::Vec3List positions;
  StateExtractor::ExtractPositions(context_,positions,true,true);
  OpenMM::State openmm_state = context_->getState(0); //create minimal state
                                                      //to only extract periodic
                                                      //box information
  OpenMM::Vec3 ucell_a,ucell_b,ucell_c;

  openmm_state.getPeriodicBoxVectors(ucell_a,ucell_b,ucell_c);

  int32_t out_n=positions.size()*4;

  // ucell
  int32_t bsize=48; // ucell block size, 6 doubles

  double a_length = 10 * sqrt(ucell_a.dot(ucell_a)); //directly transform
  double b_length = 10 * sqrt(ucell_b.dot(ucell_b)); //stuff to Angstrom
  double c_length = 10 * sqrt(ucell_c.dot(ucell_c));

  //note, that following angles have a compatibility issue
  //CHARMM seems to write out the cosines of the actual 
  //unit cell angles, while NAMD seems to write the angles itself...
  //we stick with the CHARMM way of doing things

  double alpha = ucell_b.dot(ucell_c)/(b_length*c_length);
  double beta = ucell_c.dot(ucell_a)/(c_length*a_length);
  double gamma = ucell_a.dot(ucell_b)/(a_length*b_length);

  double ucell[]={a_length, alpha, b_length,
                  beta, gamma, c_length};

  stream_.write(reinterpret_cast<char*>(&bsize),4);
  stream_.write(reinterpret_cast<char*>(ucell),bsize);
  stream_.write(reinterpret_cast<char*>(&bsize),4);

  int k=0;
  for (geom::Vec3List::iterator i=positions.begin(); 
       i!=positions.end(); ++i, ++k) {
    x[k]=float((*i)[0]);
    y[k]=float((*i)[1]);
    z[k]=float((*i)[2]);
  }
  stream_.write(reinterpret_cast<char*>(&out_n), 4);
  stream_.write(reinterpret_cast<char*>(&x[0]), out_n);
  stream_.write(reinterpret_cast<char*>(&out_n), 4);
  stream_.write(reinterpret_cast<char*>(&out_n), 4);
  stream_.write(reinterpret_cast<char*>(&y[0]), out_n);
  stream_.write(reinterpret_cast<char*>(&out_n), 4);
  stream_.write(reinterpret_cast<char*>(&out_n), 4);
  stream_.write(reinterpret_cast<char*>(&z[0]), out_n);
  stream_.write(reinterpret_cast<char*>(&out_n), 4);

  ++frames_;

  stream_.seekp(8,std::ios_base::beg);
  stream_.write(reinterpret_cast<char*>(&frames_),4);
  stream_.seekp(0,std::ios_base::end);
}

void TrajWriter::Finalize(){
  if(!registered_){
    throw ost::Error("Trajectory writer is not registered to any simulation!");
  }
  if(!stream_){
    throw ost::Error("Trajectory writer does not contain open filestream. Did you already call the finalize function before?");
  }
  stream_.close();
}

}}}// ns
