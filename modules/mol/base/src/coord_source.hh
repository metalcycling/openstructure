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
#ifndef OST_MOL_COORD_SOURCE_HH
#define OST_MOL_COORD_SOURCE_HH

/*
  Author: Marco Biasini
 */
#include <boost/shared_ptr.hpp>
#include <ost/mol/module_config.hh>
#include <ost/mol/coord_frame.hh>
#include <ost/mol/entity_handle.hh>

namespace geom {
  // we should really have fw headers...
  class Transform;
}

namespace ost { namespace mol {

class CoordSource;

typedef boost::shared_ptr<CoordSource> CoordSourcePtr;

/// \brief coordinate source
/// 
/// The coord source implements a strategy to provide coordinates, i.e from a 
/// trajectory. In it's simplest incarnaction, InMemCoordSource, the frame
/// coordinates are already present in memory. 
class DLLEXPORT_OST_MOL CoordSource {
public:
  CoordSource(const AtomHandleList& atoms);
  
  virtual ~CoordSource() {}
  
  virtual uint GetFrameCount() const =0;

  virtual CoordFramePtr GetFrame(uint frame_id) const = 0;
 
  CoordSourcePtr Extract(int start=0, int stop=-1, int step=1);
  
  // time in ps between frames
  float GetFrameDelta() const {return delta_;}
  void SetFrameDelta(float d) {delta_=d;}

  // start time in ps
  float GetStartTime() const {return start_time_;}
  void SetStartTime(float t) {start_time_=t;}
  
  int GetAtomCount() const;
  
  EntityHandle GetEntity() const;
  
  const AtomHandleList& GetAtomList() const;

  void SetAtomPos(uint frame, AtomHandle atom, const geom::Vec3& pos);

  geom::Vec3 GetAtomPos(uint frame, AtomHandle atom) const;

  /// \brief assign the coordinates in the given frame to the atoms
  void CopyFrame(uint frame);
  bool IsMutable() const;
  
  void Capture();
  void CaptureInto(int pos);
  void Capture(uint f);
  void SetFramePositions(uint frame, const std::vector<geom::Vec3>& clist);
  void SetFrameCellSize(uint frame, const geom::Vec3& size);

  virtual void AddFrame(const std::vector<geom::Vec3>& coords) = 0;
  virtual void AddFrame(const std::vector<geom::Vec3>& coords,const geom::Vec3& cell_size,const geom::Vec3& cell_angles) = 0;
  virtual void InsertFrame(int pos, const std::vector<geom::Vec3>& coords) = 0;

  void ApplyTransform(const geom::Transform& tf);

protected:
  void SetMutable(bool flag);
private:
  AtomHandleList atoms_;
  EntityHandle   entity_;
  bool           mutable_;
  std::map<long,uint> atom_dict_;
  float delta_,start_time_;
};

}}

#endif
