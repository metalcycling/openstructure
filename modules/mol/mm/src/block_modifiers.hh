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

#ifndef OST_MM_BLOCK_MODIFIERS
#define OST_MM_BLOCK_MODIFIERS

#include <ost/mol/residue_handle.hh>
#include <ost/mol/mm/buildingblock.hh>
#include <boost/shared_ptr.hpp>
#include <ost/mol/xcs_editor.hh>
#include <ost/io/binary_data_source.hh>
#include <ost/io/binary_data_sink.hh>

namespace ost{ namespace mol{ namespace mm{

class HydrogenConstructor;
class TerminiConstructor;
class BlockModifier;

typedef boost::shared_ptr<HydrogenConstructor> HydrogenConstructorPtr;
typedef boost::shared_ptr<TerminiConstructor> TerminiConstructorPtr;
typedef boost::shared_ptr<BlockModifier> BlockModifierPtr;

typedef enum {
  GromacsBlockModifiers,
  HeuristicBlockModifiers
} BlockModifierType;

class HydrogenConstructor{
public:
  HydrogenConstructor() { }
  virtual ~HydrogenConstructor() { }
  virtual void ApplyOnBuildingBlock(BuildingBlockPtr p) = 0;
  virtual void ApplyOnResidue(ost::mol::ResidueHandle& res, ost::mol::XCSEditor& ed) = 0;
  virtual void OnSave(ost::io::BinaryDataSink& ds) = 0;
  virtual BlockModifierType GetBlockModifierType() = 0;
};

class TerminiConstructor{
public:
  TerminiConstructor() { }
  virtual ~TerminiConstructor() { }
  virtual void ApplyOnBuildingBlock(BuildingBlockPtr p) = 0;
  virtual void ApplyOnResidue(ost::mol::ResidueHandle& res, ost::mol::XCSEditor& ed) = 0;
  virtual void OnSave(ost::io::BinaryDataSink& ds) = 0;
  virtual BlockModifierType GetBlockModifierType() = 0;
};

class BlockModifier{
public:
  BlockModifier() { }
  virtual ~BlockModifier() { }
  virtual void ApplyOnBuildingBlock(BuildingBlockPtr p) = 0;
  virtual void ApplyOnResidue(ost::mol::ResidueHandle& res, ost::mol::XCSEditor& ed) = 0;
  virtual void OnSave(ost::io::BinaryDataSink& ds) = 0;
  virtual BlockModifierType GetBlockModifierType() = 0;
};


}}} //ns

#endif
