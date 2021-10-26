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
#include "gosty_app.hh"

#include <ost/gui/python_shell/python_shell.hh>
#include <ost/gui/gl_win.hh>
#include <ost/gui/tools/tool_options_win.hh>
#include <ost/gui/perspective.hh>
#include <ost/gui/main_area.hh>
#include <ost/gui/python_shell/python_interpreter.hh>

#include <ost/gui/scene_win/scene_win.hh>
#include <ost/gui/sequence_viewer/sequence_viewer.hh>
#include <ost/gui/messages/message_widget.hh>
#include <ost/gui/data_viewer/data_viewer.hh>

#include <QApplication>
#include <QMainWindow>
#include <QMdiSubWindow>
#include <QMenuBar>
#include <QDesktopWidget>
/*
  Author: Marco Biasini, Andreas Schenk, Stefan Scheuber
 */

namespace ost { namespace gui {


GostyApp* GostyApp::app_=NULL;


GostyApp::GostyApp():
  py_shell_(NULL), gl_win_(NULL), scene_win_(NULL), message_widget_(NULL), seq_viewer_(NULL),
  tool_options_win_(NULL), main_(new GostyMainWindow),
  perspective_(NULL), external_widgets_(QMap<QString,WidgetGeomHandler *>()),
  try_stereo_(false)
{
  assert(GostyApp::app_==NULL);
  GostyApp::app_=this;
  main_->show();
}

GostyApp* GostyApp::Instance() {
  if (!GostyApp::app_) {
    GostyApp::app_=new GostyApp;
  }
  return GostyApp::app_;
}

bool GostyApp::ValidInstance() {
  return GostyApp::app_ != NULL;
}

void GostyApp::SetAppTitle(const QString& app_title)
{
  main_->setWindowTitle(app_title);
}


void GostyApp::OnQuit()
{
  QMapIterator<QString, WidgetGeomHandler *> i(external_widgets_);
  while (i.hasNext()) {
    i.next();
    i.value()->SaveGeom("ui/external_widgets/");
  }
}

ToolOptionsWin* GostyApp::GetToolOptionsWin()
{
  if (tool_options_win_==NULL) {
    tool_options_win_=new ToolOptionsWin;
    tool_options_win_->SetDestroyOnClose(false);
  }
  return tool_options_win_;  
}

SceneWin* GostyApp::GetSceneWin()
{
  if (scene_win_==NULL) {
    scene_win_=new SceneWin;
    scene_win_->SetDestroyOnClose(false);    
  }
  return scene_win_;  
}

SequenceViewer* GostyApp::GetSequenceViewer()
{
  if (seq_viewer_==NULL) {
    seq_viewer_=new SequenceViewer(false, true);
    seq_viewer_->SetDestroyOnClose(false);
  }
  return seq_viewer_;
}

ost::img::gui::DataViewer* GostyApp::CreateDataViewer(const ost::img::ImageHandle& d, const QString& name, const int noparent)
{
  ost::img::gui::DataViewer* viewer=NULL;
  if (noparent)
  {
    viewer=new ost::img::gui::DataViewer(NULL,d,name);
  }
  else
  {
    viewer=new ost::img::gui::DataViewer(NULL,d,name);
    this->GetPerspective()->GetMainArea()->AddWidget(name, viewer);
  }

  return viewer;

}
  
PythonShell* GostyApp::GetPyShell()
{
  if (py_shell_==NULL) {
    py_shell_=new PythonShell;
    py_shell_->SetDestroyOnClose(false);
  }
  return py_shell_;
}


void GostyApp::StopScript() 
{
  PythonInterpreter::Instance().StopScript();
}


GLWin* GostyApp::GetGLWin()
{
  if (gl_win_==NULL) {
    gl_win_=new GLWin(main_,try_stereo_);
    gl_win_->SetDestroyOnClose(false);    
  }
  return gl_win_;  
}

void GostyApp::ProcessEvents()
{
  QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents,100);
}

Perspective* GostyApp::GetPerspective()
{
  if (perspective_==NULL) {
    perspective_=new Perspective(main_);
  }
  return perspective_;
}

MessageWidget* GostyApp::GetMessageWidget()
{
  if (message_widget_==NULL) {
    message_widget_=new MessageWidget;
  }
  return message_widget_;
}

void GostyApp::AddWidgetToApp(const QString& ident, QWidget* widget)
{
  external_widgets_[ident] = new WidgetGeomHandler(ident,widget);
  external_widgets_[ident]->LoadGeom("ui/external_widgets/");
}

QWidget* GostyApp::GetWidget(const QString& ident)
{
  if (external_widgets_.contains(ident)) {
    return external_widgets_[ident]->GetWidget();
  }
  return NULL;
}

void GostyApp::RemoveWidgetFromApp(const QString& ident){
  if(external_widgets_.contains(ident)){
    external_widgets_[ident]->SaveGeom("ui/external_widgets/");
    external_widgets_.remove(ident);
  }
}

}}
