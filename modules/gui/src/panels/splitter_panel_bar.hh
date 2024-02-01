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
#ifndef OST_GUI_SPLITTER_PANEL_BAR
#define OST_GUI_SPLITTER_PANEL_BAR


#include <ost/gui/module_config.hh>
#include <ost/gui/widget_pool.hh>
#include <ost/gui/widget.hh>
#include <ost/gui/thin_splitter.hh>

#include "panel_widget_container.hh"
#include "panel_bar_widget_holder.hh"

#include <QHBoxLayout>
#include <QWidget>
#include <QString>
#include <QAction>
namespace ost { namespace gui {

/// \brief bottom Bar
class DLLEXPORT_OST_GUI SplitterPanelBar : public PanelWidgetContainer {
  Q_OBJECT
public:

  SplitterPanelBar(PanelBar* parent, Qt::Orientation orientation);
  ~SplitterPanelBar();
  
  virtual bool Save(const QString& prefix);
  virtual bool Restore(const QString& prefix);

  void WidgetMoved(Widget* widget, int position);
  QString GetName();


private slots:
  void SplitterMoved(int pos, int index);

private:
  virtual int IndexOf(Widget* widget);
  virtual void ShowWidget(Widget* widget, int pos, bool show);

  QHBoxLayout* layout_;
  ThinSplitter* splitter_;
  QList<int> sizes_;
  QList<PanelBarWidgetHolder*> widget_holders_;
};

}}

#endif
