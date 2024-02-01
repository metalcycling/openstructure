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

#ifndef OST_MM_OBSERVER_HH
#define OST_MM_OBSERVER_HH


#include <boost/shared_ptr.hpp>

#include <ost/io/mol/io_profile.hh>
#include <ost/io/mol/pdb_writer.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/coord_frame.hh>
#include <ost/mol/coord_group.hh>
#include <ost/mol/mm/state_extractor.hh>
#include <ost/mol/mm/topology.hh>
#include <ost/mol/mm/modeller.hh>
#include <ost/stdint.hh>

namespace OpenMM{
  class Context; //hacky way of telling the Context is around.
                 //will be included in source file to avoid
                 //dependencies on external libraries
}

namespace ost { namespace mol{ namespace mm{

class Observer;
class TrajObserver;
class TrajWriter;
typedef boost::shared_ptr<Observer> ObserverPtr;
typedef boost::shared_ptr<TrajObserver> TrajObserverPtr;
typedef boost::shared_ptr<TrajWriter> TrajWriterPtr;


class Observer{
public:

  virtual void Notify() = 0;

  virtual void Init(boost::shared_ptr<OpenMM::Context> c,
                    TopologyPtr top,
                    ost::mol::EntityHandle& ent) = 0;

  virtual int Rhythm() = 0;

};


class TrajObserver : public Observer{

public:  

  TrajObserver(int rhythm): rhythm_(rhythm), registered_(false) { }

  void Init(boost::shared_ptr<OpenMM::Context> c, 
            TopologyPtr top,
            ost::mol::EntityHandle& ent);

  void Notify();

  int Rhythm() { return rhythm_; }

  CoordGroupHandle GetTraj() { return c_group_; }

private:

  ost::mol::CoordGroupHandle c_group_;
  boost::shared_ptr<OpenMM::Context> context_;
  int rhythm_;
  bool registered_;
};

class TrajWriter : public Observer{

//Note, that this code is highly redundant to the code in io/mol/dcd_io.hh!!!

public:

  TrajWriter(int rhythm, const String& pdb_filename, const String& dcd_filename)
            : rhythm_(rhythm), pdb_filename_(pdb_filename),
              dcd_filename_(dcd_filename), stream_(), registered_(false),
              frames_(0) { }

  void Init(boost::shared_ptr<OpenMM::Context> c, 
            TopologyPtr top,
            ost::mol::EntityHandle& ent);

  void Notify();

  int Rhythm() { return rhythm_; }

  void Finalize();

private:

  TrajWriter(const TrajWriter& writer) { } //output filestream is unique!


  boost::shared_ptr<OpenMM::Context> context_;
  int rhythm_;
  String pdb_filename_;
  String dcd_filename_;
  std::ofstream stream_;
  bool registered_;
  uint32_t frames_;
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
};


}}}//ns

#endif
