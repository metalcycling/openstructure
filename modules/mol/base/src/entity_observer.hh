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
#ifndef OST_ENTITY_OBSERVER_HA
#define OST_ENTITY_OBSERVER_HA

#include "entity_observer_fw.hh"

#include <ost/mol/module_config.hh>

namespace ost { namespace mol {

/*
  the messages may grow more sophisticated to indicate the nature
  of the change, such as coordinate modifications, deletions/additions
  of structural elements, updates in the coonectivity, etc
*/

class DLLEXPORT_OST_MOL EntityObserver {
public:
  virtual ~EntityObserver();
  // main signal routine
  virtual void OnUpdate();
  // called when observed entity is about to go out of scope
  virtual void OnDestroy();
  
};

}} // ns

#endif
