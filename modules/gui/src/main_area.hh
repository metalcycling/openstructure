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
#ifndef OST_GUI_MAIN_AREA_HH
#define OST_GUI_MAIN_AREA_HH


#include <ost/gui/module_config.hh>

#include <QMdiArea>
/*
  Author: Marco Biasini
 */

namespace ost { namespace gui {

/// \brief main window area
/// It is implemented as a MDI (multi document interface). This allows you to add custom widgets to it.
class DLLEXPORT_OST_GUI MainArea: public QMdiArea {
public:
  Q_OBJECT
public:
  MainArea(QWidget* parent=NULL);
  
  /// \brief switch between free window and tabbed window mode
  void EnableTabbedMode(bool tabbed_mode=true);
  
  /// \brief add volatile widget
  /// \param title string that is displayed in the gui
  /// \param widget is the widget to be added to the main area.
  void AddWidget(const QString& title, QWidget* widget);
  
  /// \brief add a widget whose geometry is preserved across application relaunches
  /// 
  /// For widgets that are volatile, use #AddWidget()
  /// 
  /// If tabbed mode is enabled, the widget geometry is ignored.
  /// 
  /// \param title string that is displayed in the gui
  /// \param name is the unique name (within the scope of the main area) for
  ///     the widget that is used to restore and save the widget geometry.
  /// \param widget is the widget to be added to the main area.
  /// \param width width of the widget inside the mdi
  /// \param height height of the widget inside the mdi
  /// \param x x position of the widget inside the mdi
  /// \param y y position of the widget inside the mdi
  /// \sa AddWidget
  void AddPersistentWidget(const QString& title, const QString& name, 
                           QWidget* widget, int width, int height, int x, int y);
  /// \brief add a widget with custom window states to the ui. The Geometry is preserved across application relaunches
  /// \sa AddPersistentWidget(const QString& title, const QString& name, QWidget* widget, int width, int height, int x, int y)
  /// \param title string that is displayed in the gui
  /// \param name is the unique name (within the scope of the main area) for
  ///     the widget that is used to restore and save the widget geometry.
  /// \param widget is the widget to be added to the main area.
  /// \param window_state custom window_state for the widget. See Qt Documentation to learn more about WindowStates.
  void AddPersistentWidget(const QString& title, const QString& name,
                             QWidget* widget, Qt::WindowStates window_state = Qt::WindowStates());
  
  /// \brief display the given widget inside the main area
  /// This method can be used to make a widget visible that has been added to the mdi area.
  /// This method should only be called if you are sure, that the widget has been added to the main area.
  /// Otherwise, there might be an unexpected behaviour!
  /// \param widget which you want to make visible
  void ShowSubWindow(QWidget* widget);
  
  /// \brief hides the given widget inside the main area
  /// This method can be used to hide a widget that has been added to this mdi area.
  /// This method should only be called if you are sure, that the widget has been added to the main area.
  /// Otherwise, there might be an unexpected behaviour!
  /// \param widget which you want to hide
  void HideSubWindow(QWidget* widget);

};

/// \example mdi_example.py
///
/// Shows how to add a custom python widget to the mdi-area of OpenStructure.

}}

#endif
