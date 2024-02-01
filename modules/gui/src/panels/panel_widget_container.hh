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
#ifndef OST_GUI_PANEL_BAR_PANEL_WIDGET_CONTAINER
#define OST_GUI_PANEL_BAR_PANEL_WIDGET_CONTAINER


#include <ost/gui/module_config.hh>
#include <ost/gui/widget.hh>

#include "panel_bar.hh"

#include <QHBoxLayout>
#include <QWidget>
#include <QString>
#include <QAction>
#include <QMap>
namespace ost { namespace gui {

/// \brief A PanelWidgetContainer can display various widgets
/// A PanelWidgetContainer is a abstract class which is used by the PanelBar class as view mode.
class DLLEXPORT_OST_GUI PanelWidgetContainer : public Widget {
  Q_OBJECT
public:
  PanelWidgetContainer(PanelBar* parent);
  /// \brief Get readable name
  /// This method must return a QString which identifies this WidgetContainer. This String will be visible in the gui to select the view mode of a \ref PanelBar .
  virtual QString GetName()= 0;
  /// \brief move widget to given position
  /// Move the given widget to the index position.
  /// \param widget widget to be moved to the position.
  /// \param index new index of the given widget. If the index is out of range add the widget will be added to the front or at the end.

  virtual void WidgetMoved(Widget* widget, int index) = 0;
  /// \brief display a widget at the given position
  /// Display the given widget at the index position.
  /// \param widget widget to be moved to the position.
  /// \param pos index of the given widget. If the index is out of range add the widget will be added to the front or at the end.
  /// \param show if show if false hide the widget even the position is not correct.
  virtual void ShowWidget(Widget* widget,int pos,bool show = true) = 0;
  /// \brief Begin with dragging widgets
  /// This method is initially called in the \ref PanelManager class in the \ref PanelManager::StartDrag() method.

  virtual void StartDrag();
  /// \brief Stop dragging widgets
  /// This method is initially called in the \ref PanelManager class in the \ref PanelManager::EndDrag() method.
  virtual void EndDrag();

  /// \brief Get the PanelBar
  /// \return Returns the PanelBar in which this PanelWidgetContainer is being displayed.
  PanelBar* GetPanelBar();

private:
  PanelBar* parent_;
};

}}

Q_DECLARE_METATYPE(ost::gui::PanelWidgetContainer*)

#endif
