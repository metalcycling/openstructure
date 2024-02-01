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
#ifndef OST_GUI_SCENE_WIN_ENTITY_PART_NODE_HH
#define OST_GUI_SCENE_WIN_ENTITY_PART_NODE_HH


#include <ost/mol/query_view_wrapper.hh>

#include <ost/gfx/entity_fw.hh>

#include <ost/gui/module_config.hh>
#include <ost/gui/scene_win/scene_node.hh>

#include <QObject>
#include <QVariant>
#include <QModelIndex>
/*
  Author: Stefan Scheuber
 */

namespace ost { namespace gui {

class DLLEXPORT_OST_GUI EntityPartNode : public SceneNode {
  Q_OBJECT
public:
  EntityPartNode(QString name, gfx::EntityP entity, mol::QueryViewWrapper part, SceneNode* node_parent );

  virtual QVariant GetData(int column, int role);
  virtual Qt::ItemFlags Flags(int column) const;
  virtual int GetColumnCount() const;

  virtual void SetQueryView(mol::QueryViewWrapper part);
  virtual mol::QueryViewWrapper GetQueryView() const;

  virtual gfx::EntityP GetEntity() const;

  virtual void SetName(QString name);
  virtual const QString& GetName() const;

protected:
  virtual void SetVisible(bool visible);
  virtual bool GetVisible() const;

private:
  QString name_;
  gfx::EntityP entity_;
  mol::QueryViewWrapper query_view_;
  bool visible_;
};

}}

#endif
