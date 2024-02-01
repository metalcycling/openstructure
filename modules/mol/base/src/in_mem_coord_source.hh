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
#ifndef OST_MOL_IN_MEM_COORD_SOURCE_HH
#define OST_MOL_IN_MEM_COORD_SOURCE_HH

/*
  Author: Marco Biasini
 */

#include <ost/mol/coord_source.hh>
#include <ost/mol/handle_type_fw.hh>

namespace ost { namespace mol {


class InMemCoordSource;

typedef boost::shared_ptr<InMemCoordSource> InMemCoordSourcePtr;

/// \brief in-memory coordinate source
class DLLEXPORT_OST_MOL InMemCoordSource : public CoordSource {
public:
  InMemCoordSource(const AtomHandleList& atoms);
  
  virtual uint GetFrameCount() const;
  
  virtual CoordFramePtr GetFrame(uint frame_id) const;
  
  void AddFrame(const CoordFramePtr& frame);
  virtual void AddFrame(const std::vector<geom::Vec3>& coords);
  virtual void AddFrame(const std::vector<geom::Vec3>& coords,const geom::Vec3& cell_size,const geom::Vec3& cell_angles);
  virtual void InsertFrame(int pos, const std::vector<geom::Vec3>& coords);
  
private:
  CoordFrameList  frames_;
};

}}

#endif
