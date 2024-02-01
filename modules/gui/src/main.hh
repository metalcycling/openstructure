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
#ifndef OST_GUI_OST_MAIN_HH
#define OST_GUI_OST_MAIN_HH


#include <ost/gui/python_shell/text_logger.hh>
#include <boost/shared_ptr.hpp>
#include <ost/gui/module_config.hh>
#include <map>

#include "widget_state_saver.hh"


// Qt headers must come last
#include <QMainWindow>

class QDropEvent;
class QDragEnterEvent;
class QCloseEvent;

namespace ost { namespace gui {

class DLLEXPORT_OST_GUI GostyMainWindow : public WidgetStateSaver<QMainWindow>
{
  Q_OBJECT
public:
  GostyMainWindow();
public slots:
  void OnQuit();
protected:
  static QSize GetDefaultSize();
  virtual void dragEnterEvent (QDragEnterEvent * event);
  virtual void dropEvent(QDropEvent *event);
  virtual void closeEvent(QCloseEvent* event);
};

}} // ns

#endif

