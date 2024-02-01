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
#ifndef OST_VIEW_OP_HH
#define OST_VIEW_OP_HH

/*
  Author: Marco Biasini
 */
#include <ost/mol/module_config.hh>
#include <ost/mol/view_type_fw.hh>
#include <ost/mol/handle_type_fw.hh>
#include <ost/integrity_error.hh>
#include <ost/mol/entity_handle.hh>
namespace ost { namespace mol {
  


/// \brief return union of two entity views
/// 
/// a view containing elements that are either in ev1 or ev2 are included.
/// no duplicates are generated.
///
/// \relates EntityView
EntityView DLLEXPORT_OST_MOL Union(const EntityView& ev1, 
                                    const EntityView& ev2);
                                    
/// \brief return difference between two entity views
/// 
/// a view containing elements in ev1 but not in ev2 is returned.
///
/// \relates EntityView
EntityView DLLEXPORT_OST_MOL Difference(const EntityView& ev1, 
                                         const EntityView& ev2);  
                                         
/// \brief return intersection of two entity views
///
/// a view containing elements in ev1 and ev2 is returned
/// \relates EntityView
EntityView DLLEXPORT_OST_MOL Intersection(const EntityView& ev1, 
                                          const EntityView& ev2);
                                 
/// \brief create view from atom list
/// \throw IntegrityError if atoms belonging to different entities are 
///      encountered
/// \relates EntityView
EntityView DLLEXPORT_OST_MOL CreateViewFromAtomList(const AtomHandleList& atoms);

/// \brief create view from atom list
/// 
/// \throw IntegrityError if atoms belonging to different entities are 
///      encountered
/// \relates EntityView
EntityView DLLEXPORT_OST_MOL CreateViewFromAtomList(const AtomViewList& atoms);

/// \brief compare two entity views
///
/// \throw IntegrityError if trying to compare two views that do not point to
///     the same entity.
/// \return a pair of views. \c first contains the elements that are in \p view1,
///     but not in \p view2, \c second contains the elements that are in \p view2
///     but not in \p view1
std::pair<EntityView, EntityView> DLLEXPORT_OST_MOL 
CompareViews(const EntityView& view1, const EntityView& view2);

/// \brief create new entity handle from entity view
/// 
/// This function behaves exactly like EntityHandle::Copy(), except that only
/// atoms, residues, chains and bonds that are present in the view will be 
/// copied.
/// 
/// \param view is the view to be converted to a handle
/// \param include_exlusive_atoms if true, atoms that are part of an exclusive
///     bond (only one bond partner is included in the view) will also be included
///     in the new entity handle.
/// \param handle If invalid a new entity will be created. If valid, the atoms, 
///     residues, chains, bonds and torsions will be added to handle. This is 
///     useful to combine several entities into one.
EntityHandle DLLEXPORT_OST_MOL 
CreateEntityFromView(const EntityView& view, bool include_exlusive_atoms,
                     EntityHandle handle=EntityHandle());

}} // ns
#endif
