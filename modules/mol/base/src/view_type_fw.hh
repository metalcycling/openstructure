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
#ifndef OST_VIEW_TYPE_FW_HH
#define OST_VIEW_TYPE_FW_HH

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>

namespace ost { namespace mol {
class AtomView;
class AtomViewData;
typedef boost::shared_ptr<AtomViewData> AtomViewDataPtr;
typedef boost::weak_ptr<AtomViewData> AtomViewDataWeakPtr;
typedef std::vector<AtomView> AtomViewList;

class ResidueView;
class ResidueViewData;
typedef boost::shared_ptr<ResidueViewData> ResidueViewDataPtr;
typedef boost::weak_ptr<ResidueViewData> ResidueViewDataWeakPtr;
typedef std::vector<ResidueView> ResidueViewList;

class ChainView;
class ChainViewData;
typedef boost::shared_ptr<ChainViewData> ChainViewDataPtr;
typedef boost::weak_ptr<ChainViewData> ChainViewDataWeakPtr;
typedef std::vector<ChainView> ChainViewList;

class EntityView;
class EntityViewData;
typedef boost::shared_ptr<EntityViewData> EntityViewDataPtr;
typedef boost::weak_ptr<EntityViewData> EntityViewDataWeakPtr;
typedef std::vector<EntityView> EntityViewList;

class QueryViewWrapper;
typedef std::vector<QueryViewWrapper> QueryViewWrapperList;

/// \brief Flags to control the behaviour when adding handles to views.
struct ViewAddFlag {
  typedef enum {
    /// Include all atoms when adding a residue handle to a view    
    INCLUDE_ATOMS     =1,
    /// Include all residues when adding a chain to a view
    INCLUDE_RESIDUES  =2,
    /// Include all chains when creating a new entity view
    INCLUDE_CHAINS    =4,
    /// Convenience flags to include all substructures.
    INCLUDE_ALL       =INCLUDE_ATOMS|INCLUDE_RESIDUES|INCLUDE_CHAINS,
  
    /// If set, it will be checked that no duplicates are created when adding 
    /// a new handle.
    CHECK_DUPLICATES  =8
 } Type;
};
typedef int ViewAddFlags;

}} // ns

#endif

