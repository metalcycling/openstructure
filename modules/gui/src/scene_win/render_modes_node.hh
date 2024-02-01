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
#ifndef OST_GUI_SCENE_WIN_RENDER_MODES_NODE_HH
#define OST_GUI_SCENE_WIN_RENDER_MODES_NODE_HH

#include <ost/gfx/entity_fw.hh>
#include <ost/gfx/gfx_node_fw.hh>

#include <ost/gui/module_config.hh>
#include <ost/gui/scene_win/label_node.hh>

#include "render_mode_node.hh"

#include <QObject>
#include <QVariant>
#include <QSet>
/*
  Author: Stefan Scheuber
 */

namespace ost { namespace gui {

class DLLEXPORT_OST_GUI RenderModesNode : public LabelNode {
  Q_OBJECT
public:
  RenderModesNode(gfx::EntityP entity, SceneNode* node_parent );

  virtual void RenderModeChanged();
  gfx::GfxNodeP GetGfxNode() const;

private:
  void Update();

  gfx::GfxNodeP node_;
  QMap<gfx::RenderMode::Type,RenderModeNode*> render_types_;
};

}}

#endif
