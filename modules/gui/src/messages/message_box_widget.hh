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
#ifndef OST_GUI_MESSAGES_BOX_MESSAGE_WIDGET_HH
#define OST_GUI_MESSAGES_BOX_MESSAGE_WIDGET_HH


#include <ost/gui/module_config.hh>

#include <QWidget>
#include <QMap>
#include <QMessageBox>
/*
  Author: Stefan Scheuber
*/

namespace ost { namespace gui {

class DLLEXPORT_OST_GUI MessageLevel: public QWidget
{
  Q_OBJECT;
public:
  MessageLevel(QPixmap icon, QWidget* parent);

  void SetMessageCount(int message_count);
  int GetMessageCount();
private:
  int message_count_;
  QLabel* count_label_;
};

// a summarized version for the the Log Messages
class DLLEXPORT_OST_GUI MessageBoxWidget: public QWidget
{
  Q_OBJECT;
public:
  MessageBoxWidget(QWidget* parent=NULL);
  ~MessageBoxWidget();

public slots:
  void Update(QMessageBox::Icon icon);
  void UpdateAll();

protected:
  void mouseDoubleClickEvent(QMouseEvent* event);

private:
  QMap<QMessageBox::Icon,MessageLevel*> level_map_;

};

}} // ns

#endif
