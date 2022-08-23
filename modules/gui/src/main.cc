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
#include "main.hh"

#include <sstream>


#include <ost/log.hh>
#include <ost/gui/python_shell/python_shell.hh>
#include <ost/gui/file_loader.hh>

#include "gl_win.hh"
//#include "entity_explorer/entity_explorer.hh"
#include "tools/tool_options_win.hh"


#include <QApplication>
#include <QDesktopWidget>
#include <QScreen>
#include <QAction>
#include <QMenu>
#include <QMenuBar>
#include <QFileDialog>
#include <QDockWidget>
#include <QSettings>
#include <QMdiSubWindow>
#include <QMessageBox>
#include <QPixmap>
#include <QDragEnterEvent>
#include <QList>
#include <QUrl>
#include <QMdiArea>
#include <QAction>
#include <QFile>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>

namespace ost { namespace gui {

QSize GostyMainWindow::GetDefaultSize()
{
  QRect desktop_geom=QGuiApplication::primaryScreen()->availableGeometry();
  if(desktop_geom.height() > desktop_geom.width()){
	  return QSize(desktop_geom.width(),desktop_geom.width());
  }
  else{
	  return QSize(desktop_geom.width()*0.8,desktop_geom.height());
  }
}
GostyMainWindow::GostyMainWindow():
  WidgetStateSaver<QMainWindow>("windows/main", QPoint(20,20),
                                GostyMainWindow::GetDefaultSize(), 0)
{
  setAcceptDrops(true);
}

void GostyMainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if(event->mimeData()->hasUrls()){
    event->acceptProposedAction();
  }
}

void GostyMainWindow::dropEvent(QDropEvent *event)
{
  if(event->mimeData()->hasUrls()){
    QList<QUrl> urls =  event->mimeData()->urls();
    for(int i=0;i<urls.size(); i++){
      FileLoader::LoadObject(urls[i].toLocalFile());
    }
    event->acceptProposedAction();
  }
}

void GostyMainWindow::OnQuit()
{
  //SetFullScreen(false);
  GostyApp::Instance()->OnQuit();
  deleteLater();
  // exit has to be called on the instance not the class,
  // otherwise the aboutToQuit singal doesn't get emitted
  QApplication::instance()->exit(0);
}

void GostyMainWindow::closeEvent(QCloseEvent * event)
{
  this->OnQuit();
}

}} // ns
