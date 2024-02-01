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
#ifndef OST_GUI_PANEL_BAR_PANEL_BAR
#define OST_GUI_PANEL_BAR_PANEL_BAR


#include <ost/gui/module_config.hh>
#include <ost/gui/widget_pool.hh>
#include <ost/gui/widget.hh>
#include <ost/gui/thin_splitter.hh>

#include "button_box.hh"


#include <QStackedLayout>
#include <QWidget>
#include <QString>
#include <QAction>
#include <QActionGroup>
namespace ost { namespace gui {

class PanelWidgetContainer;

struct WidgetState
{
  WidgetState(Widget * w,bool v=true, QString n=""):
  widget(w), visible(v), name(n)
  { }
  Widget* widget;
  bool visible;
  QString name;
};

/// \brief Panel Bar which may contain other Widgets
/// A PanelBar might have different ViewModes which are \ref PanelWidgetContainer classes. Per default a PanelWidgetContainer has no ViewMode therefore first of all a view mode must be added. Otherwise nothing will be displayed.
/// The ViewMode is responsible for how the widget is displayed. This class keeps only track about the order of the added widgets and the ViewModes.
class DLLEXPORT_OST_GUI PanelBar : public Widget {
  Q_OBJECT
public:
  PanelBar(QWidget* parent);
  virtual ~PanelBar();

  virtual bool Save(const QString& prefix);
  virtual bool Restore(const QString& prefix);

  /// \brief Add a widget to this widget
  /// Puts a widget at the end of the current PanelWidgetContainer.
  /// \param widget Which should be added
  /// \param hidden if this is true widget will not be shown
  virtual void AddWidget(Widget* widget, bool hidden);

  // \brief Add a view mode to this widget
  // Adds a view mode to this class. The view mode will automaticaly added to the view mode menu see /ref PanelBar::AddGetViewModeMenu() for more information.
  // \param pwc View Mode to be added
  void AddViewMode(PanelWidgetContainer * pwc);
  /// \brief Get QMenu with view modes
  /// The returned menu can be used for exmaple in the menubar to allow the user choose the view mode.
  /// \return QMenu with all added view modes as actions.
  QMenu* GetViewModeMenu();

  /// \brief Get the index of a given widget
  /// The index of a widget is the position number. This method can be used to get the position of a widget within the PanelBar.
  /// \param widget get index of this widget
  /// \return index of the given widget -1 if Widget is not found.
  int GetIndex(Widget* widget);

  /// \brief Begin with dragging widgets
  /// This method is called by the \ref PanelManager class. Read more about this method at \ref PanelManager::StartDrag() .
  void DragEnter();
  /// \brief Stop dragging widgets
  /// This method is called by the \ref PanelManager class. Read more about this method at \ref PanelManager::EndDrag() .
  void EndDrag();

  /// \brief Get action which changes the visibility of this widget
  /// The returned action can be used for example in a MenuBar to allow hiding / showing this class.
  /// \return QAction that can change the visibility of this widget
  QAction* GetShowAction();

public slots:
  /// \brief This slot adds a given widget at the given position
  /// The Widget will be added at the given index. If the widget is already in this panel, it will be moved to the given index.
  /// \param widget Which should be inserted
  /// \param index the position where the widget will be inserted
  /// \param is_hidden marks if the widget will be inserted hidden or visible. true = widget is visible.
  void InsertWidget(Widget * widget, int index, bool is_hidden=false);
  /// \brief Change to the next view mode
  /// This method changes the view mode to the next view mode. If the view mode reaches the end of the list, the first one will be displayed.
  void ToggleViewMode();

  /// \brief Toggle the visibility of all Panels
  void ToggleHide();

  /// \brief Remove a Widget
  /// The widget will be removed if it is in this PanelBar
  /// \param widget which should be removed
  void RemoveWidget(Widget* widget);

private slots:
  void ChangeViewMode(QAction* action);
  void ToggleVisibility(Widget*);
  void ShowActionTrigger();
  void SetViewMode(int index);
private:
  void UpdateShow();

  QStackedLayout* layout_;
  QActionGroup* view_modes_;
  QMenu* view_mode_menu_;
  PanelWidgetContainer* current_view_mode_;
  QList<WidgetState> widget_states_;
  Widget* drop_box_;
  QAction* show_action_;
};

}}

#endif
