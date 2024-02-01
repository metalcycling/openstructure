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
/*
  Author: Marco Biasini
 */
#include "selection_tool.hh"


#include <ost/platform.hh>
#include <ost/gfx/scene.hh>

#include <ost/gui/gl_canvas.hh>

#include <QDir>
namespace ost { namespace gui {

SelectionTool::SelectionTool()
  : Tool("Selection Tool")
{
  ToolOptions* opts=this->GetToolOptions();
  ToolOptionEnum* sel_mode(new ToolOptionEnum("sel_mode", "Selection Mode"));
  sel_mode->Add("Atom", 1);
  sel_mode->Add("Residue", 2);  
  sel_mode->Add("Chain", 3); 
  sel_mode->Add("Bond", 4);  
  sel_mode->SetIndex(0);
  opts->AddOption(sel_mode);
}

int SelectionTool::GetSelectionMode() const
{
  ToolOption* opt=this->GetToolOptions()->GetOption("", "sel_mode");
  ToolOptionEnum* sel_mode=dynamic_cast<ToolOptionEnum*>(opt);
  return sel_mode->GetValue();
}

void SelectionTool::Click(const MouseEvent& event)
{
  int sel_mode=this->GetSelectionMode();
  gfx::Scene& scene=gfx::Scene::Instance();
  scene.SetSelectionMode(sel_mode);
  scene.Pick(event.GetPos().x(), 
             scene.GetViewport().height-event.GetPos().y(), 8);
}


bool SelectionTool::CanOperateOn(gfx::NodePtrList nodes)
{
  return true;
}

QIcon SelectionTool::GetIcon() const
{
  QDir icon_path(GetSharedDataPath().c_str());
  icon_path.cd("gui");
  icon_path.cd("icons");
  QIcon icon = QIcon(icon_path.absolutePath()+QDir::separator()+QString("selection_icon.png"));
  return icon;
}

}}
