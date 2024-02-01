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
#include <ost/log.hh>
#include <ost/config.hh>
#include <ost/gfx/scene.hh>

#include <ost/gui/gl_win.hh>
#include <ost/gui/gl_canvas.hh>
#include <ost/gui/tools/tool_manager.hh>
#include <ost/gui/tools/selection_tool.hh>
#include <ost/gui/tools/manipulator.hh>
#include <ost/gui/tools/measure_tool.hh>
#include <ost/gui/perspective.hh>
#include <ost/gui/tools/tool.hh>
#include <ost/gui/tools/tool_bar.hh>
#include <ost/gui/scene_menu.hh>
#include <ost/gui/tools/map_tool.hh>

/*
  Authors: Ansgar Philippsen, Marco Biasini
 */

namespace ost { namespace gui {

GLWin::GLWin(QWidget* p, bool try_stereo): Widget(NULL, p) {


  gl_canvas_ = new GLCanvas;

  if(try_stereo) {
    gl_canvas_->SetStereoFormat();
  } else{
    gl_canvas_->SetDefaultFormat();
  }
  
  gl_canvas_widget_ = QWidget::createWindowContainer(gl_canvas_);
  gl_canvas_widget_->setFocusPolicy(Qt::StrongFocus);
  gl_canvas_widget_->setMouseTracking(true);

  toolbar_ = new ToolBar(this); 
  ToolManager::Instance().AddTool(new SelectionTool);
  ToolManager::Instance().AddTool(new Manipulator);
  ToolManager::Instance().AddTool(new MeasureTool);
  ToolManager::Instance().AddTool(new MapTool);

  scene_menu_ = new SceneMenu();

  QMainWindow* main = new QMainWindow;
  main->setCentralWidget(gl_canvas_widget_);
  main->addToolBar(Qt::LeftToolBarArea, toolbar_);
  this->SetInternalWidget(main);
  gfx::Scene::Instance().AttachObserver(this);

  connect(gl_canvas_, SIGNAL(CustomContextMenuRequested(const QPoint&)), 
          this, SLOT(RequestContextMenu(const QPoint&)));
  connect(&ToolManager::Instance(), SIGNAL(ActiveToolChanged(Tool*)), 
          this, SLOT(ActiveToolChanged(Tool*)));
}

void GLWin::ActiveToolChanged(Tool* t) {
  gfx::Scene::Instance().RequestRedraw();
}

GLWin::~GLWin() {
  gfx::Scene::Instance().DetachObserver(this);
}

void GLWin::SetTestMode(bool f) {
  gl_canvas_->SetTestMode(f);
}

void GLWin::StatusMessage(const String& m) {
  GostyApp::Instance()->GetPerspective()->StatusMessage(m);
}

bool GLWin::Restore(const QString& prefix) {
  return true;
}

bool GLWin::Save(const QString& prefix) {
  return true;
}

void GLWin::RequestContextMenu(const QPoint& pos) {
  QPoint pick_pos(pos.x(), gl_canvas_->height()-pos.y());
  scene_menu_->SetPickPoint(pick_pos);
  scene_menu_->ShowMenu(gl_canvas_->mapToGlobal(pos));
}

}} // ns
