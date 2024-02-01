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
#ifndef OST_IO_ENTITY_DCD_IO_HH
#define OST_IO_ENTITY_DCD_IO_HH

/*
  CHARMM trajectory IO
  Authors: Ansgar Philippsen, Marco Biasini
 */

#include <ost/io/module_config.hh>
#include <ost/mol/coord_group.hh>
#include <ost/io/mol/io_profile.hh>
namespace ost { namespace io {



/*! \brief import a CHARMM trajectory in dcd format with an existing entity
    requires the existing entity and the trajectory file - obviously the
    atom layout of the entity must match the trajectory file
*/
mol::CoordGroupHandle DLLEXPORT_OST_IO LoadCHARMMTraj(const mol::EntityHandle& ent,
                                                      const String& trj_filename,
                                                      unsigned int stride=1,
                                                      bool lazy_load=false,
                                                      bool detect_swap=true,
                                                      bool byte_swap=false);


/*! \brief export coord group as PDB file and DCD trajectory
    if the pdb filename is an empty string, it won't be exported
    the optional stride parameter will cause every nth frame to be exported
 */
void DLLEXPORT_OST_IO SaveCHARMMTraj(const mol::CoordGroupHandle& coord_group, 
                                     const String& pdb_filename, 
                                     const String& dcd_filename,
                                     unsigned int stride=1,
                                     const IOProfile& profile=IOProfile());


}} // ns

#endif
