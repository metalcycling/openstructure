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
#ifndef OST_HANDLE_TYPE_FW_HH
#define OST_HANDLE_TYPE_FW_HH

#include <vector>
#include <list>
namespace ost { namespace mol {

class ChainHandle;
typedef std::vector<ChainHandle> ChainHandleList;
class ResidueHandle;
typedef std::vector<ResidueHandle> ResidueHandleList;
class AtomHandle;
typedef std::vector<AtomHandle> AtomHandleList;
class EntityHandle;

class TorsionHandle;
typedef std::vector<TorsionHandle> TorsionHandleList;

class BondHandle;
typedef std::vector<BondHandle> BondHandleList;

class EntityHandle;
typedef std::vector<EntityHandle> EntityHandleList;
}}
#endif

