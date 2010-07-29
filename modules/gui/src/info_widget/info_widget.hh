//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
#ifndef OST_GUI_INFO_WIDGET_INFO_WIDGET_HH
#define OST_GUI_INFO_WIDGET_INFO_WIDGET_HH

#include <QListView>
#include <QMessageBox>
#include <QStandardItemModel>

#include <ost/gui/widget.hh>
#include <ost/gui/module_config.hh>

/*
  Author: Stefan Scheuber
*/

namespace ost { namespace gui {

// the display window for all graphical objects
class DLLEXPORT_OST_GUI InfoWidget: public Widget
{
  Q_OBJECT;
public:
  InfoWidget(QWidget* parent=NULL);
  ~InfoWidget();

public:
  virtual void LogMessage(const QString& message, QMessageBox::Icon icon=QMessageBox::Information);
  virtual void LogMessage(QStandardItem* item);
  virtual void LogMessage(const QString& message, QIcon icon);

  virtual void SetMaxMessages(int max){}
  virtual int GetMaxMessages(){return 0;}

  virtual int GetMessagesCount(){return 0;}

  virtual bool Save(const QString& prefix) { return true; }
  virtual bool Restore(const QString& prefix) { return true; }

  ActionList GetActions();
public slots:
  void Clear();
  void RemoveSelected();
  void Update();

private:
  QPixmap GetIcon(QMessageBox::Icon icon, QWidget* widget);


  QStandardItemModel* model_;
  QListView* view_;

  ActionList actions_;
};

}} // ns

#endif
