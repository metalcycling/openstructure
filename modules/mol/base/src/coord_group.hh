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
#ifndef OST_COORD_GROUP_HH
#define OST_COORD_GROUP_HH

/*
  Author: Ansgar Philippsen, Marco Biasini
*/

#include <vector>
#include <boost/shared_array.hpp>
#include "atom_handle.hh"
#include "coord_source.hh"

namespace geom {
  // we should really have fw headers...
  class Transform;
}

namespace ost { namespace mol {

/// \brief coordinate group, for trajectories and such
class DLLEXPORT_OST_MOL CoordGroupHandle {
public:
  /// \brief create empty, invalid handle
  CoordGroupHandle();
  /// \brief used internally
  CoordGroupHandle(CoordSourcePtr source);

  /// \brief return trajectories entity handle
  EntityHandle GetEntity() const;

  //EntityGroupHandle ExtractFromTo(int start, int end, int step=1);
  //EntityGroupHandle ExtractEvery(int step, int start=0);

  /// \brief number of atoms per frame
  uint GetAtomCount() const;

  /// \brief number of frames
  uint GetFrameCount() const;

  /// \brief time in ps between frames, default 0.0
  float GetDelta() const;

  /// \brief set time between each frame in ps
  void SetDelta(float d);

  /// \brief start time in ps, default 0.0
  float GetStartTime() const;

  /// \brief set start time in ps
  void SetStartTime(float t);

  /// \brief assign positions to the given frame - order and count must match 
  ///      initial atomlist
  void SetFramePositions(uint frame, const geom::Vec3List& clist);

  /// \brief set cell size for the given frame
  void SetFrameCellSize(uint frame, const geom::Vec3& size);

  /// \brief get the positions of all the atoms in the given frame
  geom::Vec3List GetFramePositions(uint frame);  
  
  /// \brief copy atom positions of given frame to stored atoms in entity
  void CopyFrame(uint frame);
  
  /// \brief store current atom positions in new frame
  void Capture();

  void CaptureInto(int pos);

  /// \brief store current atom positions in given frame (must exist)
  void Capture(uint frame);
  
  /// \brief add frame 
  //void AddFrame(const std::vector<geom::Vec3>& clist);
  void AddFrame(const geom::Vec3List& clist);
  void AddFrame(const geom::Vec3List& clist,const geom::Vec3& cell_size,const geom::Vec3& cell_angles);
  
  void AddFrames(const CoordGroupHandle& cg);
  /// \brief set an indidivial atom position in the given frame
  void SetAtomPos(uint frame, AtomHandle atom, const geom::Vec3& pos);

  /// \brief retrieve an indidivial atom position in the given frame
  geom::Vec3 GetAtomPos(uint frame, AtomHandle atom) const;
  
  /// \brief check for handle validity
  bool IsValid() const;

  /// \brief whether the handle is valid
  operator bool() const;
  
  AtomHandleList GetAtomList() const;
  CoordFramePtr GetFrame(uint frame) const;
  CoordFrame GetFrame2(uint frame);
  
  /// \brief return a filtered coord group, containing only the atoms in the 
  ///     view
  CoordGroupHandle Filter(const EntityView& selected,int first=0,int last=-1,int stride=1) const;

  /// \brief apply in-place transform to each coordinate in each frame
  void ApplyTransform(const geom::Transform& tf);


private:
  void CheckValidity() const;
  CoordSourcePtr source_;
};

// factory method
// create with a reference set of atoms and zero frames
DLLEXPORT_OST_MOL CoordGroupHandle CreateCoordGroup(const AtomHandleList& atoms);

}} // ns

#endif
