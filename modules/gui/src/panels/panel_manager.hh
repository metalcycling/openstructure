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
#ifndef OST_GUI_PANEL_BAR_PANEL_MANAGER
#define OST_GUI_PANEL_BAR_PANEL_MANAGER


#include <ost/gui/widget.hh>
#include <ost/gui/widget_pool.hh>
#include <ost/gui/thin_splitter.hh>

#include "panel_widget_container.hh"
#include "button_box.hh"

#include <QWidget>
#include <QString>
#include <QList>
#include <QMenu>
#include <QDragEnterEvent>
namespace ost { namespace gui {

enum PanelPosition
{
  LEFT_PANEL = 0,
  RIGHT_PANEL,
  BOTTOM_PANEL,
  NONE
};

/// \brief Class which organizes all widgets which are in the side panels
/// This class handles all side bar widgets. It can be used to display, hide or move a widget to a PaneBar. There are three Bars (left, bottom, right) which are organized by this class.
/// Whenever a widget is being removed or added it checks first if the widget type is known and if there are available instances.
class DLLEXPORT_OST_GUI PanelManager : public Widget {
  Q_OBJECT
public:
  PanelManager(QWidget* widget);
  ~PanelManager(){ }
  virtual bool Save(const QString& prefix);
  virtual bool Restore(const QString& prefix);

  /// \brief Add a widget to the widget pool
  ///
  /// The widget must already be in the WidgetRegistry.
  /// If you are not sure if the Widget is in the WidgetRegistry, use the \ref PanelManager::AddWidgetToPool(const QString& name, Widget* widget) Method instead.
  /// \param class_name class name of class which should be added to WidgetRegistry.
  /// \param limit amount of parallel instances allowed (-1 if infinite)
  virtual void AddWidgetToPool(const QString& class_name, int limit=-1);
  /// \brief Add a widget to the widget pool
  ///
  /// Same as \ref AddWidgetToPool(const QString& class_name, int limit)
  /// \param name Name which is displayed in the gui.
  /// \param widget Widget which will be added to the WidgetPool of this class and the WidgetRegistry.
  virtual void AddWidgetToPool(const QString& name, Widget* widget);

  /// \brief Display a Widget in a PanelBar
  /// With Method you can add a widget to the given PanelBar. The widget which finally will be added to the gui will be created from the WidgetRegistry.
  /// If the WidgetPool does not know the class name of the given widget or if there are no instances left, nothing will happen.
  /// \param pos Indicates which PanelBar is affected
  /// \param widget the widget will not directly added to the PanelBar. The class_name will be used to identify the widget in the WidgetRegistry which will return a fresh instance of this class.
  /// \param is_hidden marks if the class should be displayed in the gui. Default the widget will be shown.
  void AddWidget(PanelPosition pos, Widget* widget, bool is_hidden=false);
  /// \brief Display a Widget in a PanelBar
  /// Same as \ref PanelManager::AddWidget(PanelPosition pos, Widget* widget, bool is_hidden)
  /// \param pos Indicates which PanelBar is affected
  /// \param class_name the class_name of the widget you would like to add.
  /// \param is_hidden marks if the class should be displayed in the gui. Default the widget will be shown.
  void AddWidgetByName(PanelPosition pos, const QString& class_name, bool is_hidden=false);

  /// \brief Replace a Widget with another Widget
  /// With this Method you can replace a Widget which is in a PanelBar with another Widget.
  /// If the Widget is not found in any of the PanelBar or there is no instance left in the pool of the second Widget nothing happens.
  /// \param w1 The widget which will be replaced
  /// \param w2 The new Widget which should replace the old one.
  void ReplaceWidget(Widget* w1, Widget* w2);
  /// Same as \ref PanelManager::ReplaceWidget(Widget* w1, Widget* w2)
  /// \param w1 The widget which will be replaced
  /// \param class_name The class_name which from which a instance will be created and then replaces the old one.
  void ReplaceWidget(Widget* w1, QString& class_name);

  /// \brief Get Menu
  /// The GetMenu method returns a QMenu reference, which contains various actions. The action states will be updated automatically.
  /// \return Reference to a QMenu which can be used for example in a QMenuBar
  QMenu* GetMenu();
  /// \brief Get a List with available Widgets
  /// Get all the Widgets which have available resources in the WidgetPool.
  /// \return QList With all available Widgets
  QList<QString> GetAvailableWidgets();

  /// \brief Begin with dragging widgets
  /// This method should be called, when someone tries to drop a widget on a PanelBar. It makes all PanelBar visible (even they are empty).
  /// Beware, you must call \ref EndDrag() after drag event is over! Otherwise, the gui might not look like you expect.
  void StartDrag();
  /// \brief End with dragging widgets
  /// This method should be called always when a drag event is over. Normally \ref StartDrag() is called before calling this method.
  void EndDrag();

public slots:
  /// \brief get a QMenu with available Widget
  /// Same as \ref GetMenu()
  /// \return Returns a QMenu with all available Widgets as QActions.
  QMenu* GetAvailableWidgetsMenu();

  /// \brief Toggle between the ViewMode 's of the Given Panel
  /// \param pos Panel which ViewMode should be toggled
  void ToggleViewMode(PanelPosition pos);

  /// \brief Move a Widget to the given Position
  /// The widget is being moved to the given position. If the widget is not in the WidgetPool or if there are no instances left nothing happens.
  /// The widget will always be visible at the new position.
  void MoveWidget(Widget * widget, PanelPosition pos, int index);

  /// \brief Move a Widget next to the given Widget
  /// The widget will be displayed next to the target Widget. The \ref PanelWidgetContainer decides what means "next".
  /// \param target Target Widget the other widget will be next this widget
  /// \param widget This widget will be displayed next the Target Widget
  void MoveNextTo(Widget* target, Widget* widget);

  /// \brief Remove a Widget out of a PanelBar
  /// The widget will be removed if it is in a PanelBar
  /// \param widget which should be removed
  void RemoveWidget(Widget * widget);

  /// \brief Toggle the visibility of all Panels
  void ToggleHide();
  /// \brief Toggle the visibility of the given Panel
  /// \param pos Panel which should be toggled
  void ToggleHide(PanelPosition pos);

private slots:
  void UpdateAddMenu();
  void UpdateCurrentMenu();
  void AddWidgetFromMenu(QAction * action);

private:
  PanelPosition GetParentPanel(Widget* widget);

  ThinSplitter* right_panel_splitter_;
  ThinSplitter* left_panel_splitter_;
  ThinSplitter* bottom_panel_splitter_;
  QMap<PanelPosition, PanelBar*> panels_;
  WidgetPool* pool_;
  QMenu* add_menu_;
  PanelPosition current_menu_;
};

}}

#endif

/// \example widget_example.py
///
/// Shows how to add a fancy python widget to the right side bar.
