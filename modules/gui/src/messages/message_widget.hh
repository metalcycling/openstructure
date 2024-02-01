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
#ifndef OST_GUI_MESSAGES_MESSAGE_WIDGET_HH
#define OST_GUI_MESSAGES_MESSAGE_WIDGET_HH


#include <ost/gui/widget.hh>
#include <ost/gui/module_config.hh>

#include <QListView>
#include <QMessageBox>
#include <QStandardItemModel>

/*
  Author: Stefan Scheuber
*/

namespace ost { namespace gui {

// the display window for Log Messages
class DLLEXPORT_OST_GUI MessageWidget: public Widget
{
  Q_OBJECT;
public:
  MessageWidget(QWidget* parent=NULL);
  ~MessageWidget();

public:
  virtual void LogMessage(const QString& message, QMessageBox::Icon icon=QMessageBox::Information);
  virtual void LogMessage(QStandardItem* item);
  virtual void LogMessage(const QString& message, QIcon icon);

  virtual int GetMessagesCount(QMessageBox::Icon icon=QMessageBox::NoIcon);
  virtual int GetTotalMessagesCount();
  virtual bool Save(const QString& prefix) { return true; }
  virtual bool Restore(const QString& prefix) { return true; }

  QPixmap GetIcon(QMessageBox::Icon icon, QWidget* widget);

  ActionList GetActions();

signals:
  void MessageCountChanged(QMessageBox::Icon);
  void AllCleared();

public slots:
  void ErrorOutput(unsigned int id,const QString& output);
  void Clear();
  void RemoveSelected();
  void Update();

private slots:
  void ContextMenuRequested(const QPoint& pos);

private:
  void Increase(QMessageBox::Icon icon);
  void Decrease(QMessageBox::Icon icon);

  QStandardItemModel* model_;
  QListView* view_;

  ActionList actions_;
  QMap<int,int> count_map_;
};

}} // ns

#endif
