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

#include <ost/gui/gosty_app.hh>
#include <ost/mol/query_view_wrapper.hh>
#include <ost/gui/scene_win/scene_win_model.hh>
#include <ost/gfx/scene.hh>
#include <ost/gfx/gfx_node.hh>
#include <ost/gfx/entity.hh>

#include <ost/gui/scene_win/scene_win.hh>

#include "current_selection_node.hh"

#include <QFont>
namespace ost { namespace gui {

CurrentSelectionNode::CurrentSelectionNode(gfx::EntityP entity, 
                                           SceneNode* parent):
   EntityPartNode("Current Selection", entity, 
                  mol::QueryViewWrapper(entity->GetSelection()),parent),
   wrapper_(mol::QueryViewWrapper(entity->GetSelection())){
}

void CurrentSelectionNode::SetQueryView(mol::QueryViewWrapper part)
{
  //Do Nothing
}

mol::QueryViewWrapper CurrentSelectionNode::GetQueryView() const
{
  wrapper_ = mol::QueryViewWrapper(this->GetEntity()->GetSelection());
  return wrapper_;
}

}}

