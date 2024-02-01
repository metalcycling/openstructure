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
#include <ost/gui/tools/tool_manager.hh>

#include <ost/gui/scene_win/scene_win.hh>

#include "tool_bar.hh"

#include <QAction>
#include <QDebug>
#include <QActionGroup>

namespace ost { namespace gui {

ToolBar::ToolBar(QWidget* parent):
  QToolBar(parent), group_(new QActionGroup(this))
{ 
  group_->setExclusive(true);
  this->setIconSize(QSize(16,16));
  QWidget* scene_win=GostyApp::Instance()->GetSceneWin();       
  connect(group_, SIGNAL(triggered(QAction*)), this,
          SLOT(ActivateTool(QAction*)));  
  connect(this, SIGNAL(actionTriggered(QAction*)),
          this,SLOT(OnAction(QAction*)));
  connect(scene_win, SIGNAL(ActiveNodesChanged(gfx::NodePtrList,gfx::EntityP,mol::QueryViewWrapperList)),
          this,SLOT(OnActiveGfxNodeChanged(gfx::NodePtrList)));
  connect(&ToolManager::Instance(), SIGNAL(ToolAdded(const Tool*&)),
          this, SLOT(OnToolAdded(const Tool*&)));          
}


void ToolBar::OnAction(QAction* action)
{

}

void ToolBar::OnActiveGfxNodeChanged(gfx::NodePtrList nodes)
{
  // disable/enable QActions in toolbar based on capabilities of the Tool.
  ToolManager& man=ToolManager::Instance();
  QList<QAction*> actions=this->GetActions();
  for (QList<QAction*>::iterator i=actions.begin(), 
       e=actions.end(); i!=e; ++i) {
    QAction* a=*i;
    Tool* tool=man.FindToolByName(a->toolTip().toStdString());
    if (tool) {
      bool can_operate=tool->CanOperateOn(nodes);
      a->setEnabled(can_operate);
      if (!can_operate && a->isChecked()) {
        a->setChecked(false);
      }
    }
  }
}

QAction* ToolBar::AddAction(const QString& text, QIcon& icon)
{
  QString str=" "+QString(1, text[0])+" ";  
  QAction* action=group_->addAction(str);
  action->setIcon(icon);
  action->setToolTip(text);
  action->setCheckable(true);  
  action->setEnabled(false);
  this->addAction(action);
  return action;
}

QList<QAction*> ToolBar::GetActions()
{
  return group_->actions();
}

void ToolBar::ActivateTool(QAction* action)
{
  ToolManager& tm=ToolManager::Instance();
  Tool* tool=tm.FindToolByName(action->toolTip().toStdString());
  if (tool) {
    tm.SetActiveTool(tool);
  }
}

void ToolBar::OnToolAdded(const Tool*& tool)
{
  QIcon icon  = tool->GetIcon();
  this->AddAction(QString(tool->GetName().c_str()), icon);
}

}}
