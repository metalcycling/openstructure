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

#ifndef OST_MM_STATE_EXTRACTOR_HH
#define OST_MM_STATE_EXTRACTOR_HH

#include <vector>

#include <boost/shared_ptr.hpp>
#include <ost/geom/vec3.hh>

namespace OpenMM{
  class Context; //hacky way of telling the Context is around.
                 //will be included in source file to avoid
                 //dependencies on external libraries

}

namespace ost { namespace mol{ namespace mm{

class StateExtractor{

public:

  static void ExtractPotentialEnergy(boost::shared_ptr<OpenMM::Context> context,
  	                                 Real& energy);

  static void ExtractKineticEnergy(boost::shared_ptr<OpenMM::Context> context,
  	                               Real& energy);

  static void ExtractEnergy(boost::shared_ptr<OpenMM::Context> context,
  	                        Real& energy);

  static void ExtractPositions(boost::shared_ptr<OpenMM::Context> context,
  	                           geom::Vec3List& positions,
                               bool enforce_periodic_box = false,
                               bool in_angstrom = true);

  static void ExtractVelocities(boost::shared_ptr<OpenMM::Context> context,
  	                            geom::Vec3List& velocities);

  static void ExtractForces(boost::shared_ptr<OpenMM::Context> context,
  	                        geom::Vec3List& forces);

};


}}} //ns

#endif
