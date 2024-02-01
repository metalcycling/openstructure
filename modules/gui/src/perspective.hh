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
#ifndef OST_GUI_PERSPECTIVE_HH
#define OST_GUI_PERSPECTIVE_HH

/*
  Author: Marco Biasini, Stefan Scheuber
 */


#include <ost/gui/module_config.hh>

#include <QObject>

class QMenu;
class QMenuBar; 
class QStatusBar;
class QMainWindow;
class QTextEdit;
class QSlider;

namespace ost { namespace gui {


class SideBar;
class ThinSplitter;
class MainArea;
class PanelManager;
class MainArea;
class ButtonBar;

/// The perspective manages the layout of the widgets inside the main window.
/// It contains important classes which itself manages again a sub part of the whole layout.
class DLLEXPORT_OST_GUI Perspective : public QObject {
  Q_OBJECT  
public:
  Perspective(QMainWindow* parent);
  ~Perspective();

  /// \brief Returns Menubar
  /// \return The Menubar of the Application. Can be used to add some menupoints.
  QMenuBar* GetMenuBar();

  /// \brief get QMenu that corresponds to the given name
  /// \return The menu with the given name. If it does not exist, it will be created and returned.
  QMenu* GetMenu(const QString& name);
  
  /// \brief The PanelManager class organizes all the widgets in the side panels.
  /// \return The PanelManager instance which is used in this perspective.
  PanelManager* GetPanels();

  /// \brief Returns the MainArea
  /// \return The main area which is used in this perspective.
  MainArea* GetMainArea();

  /// \brief hide all bars
  /// Hides all side bars. Can be used if the MainArea should be expanded to full size.
  void HideAllBars();
  
  /// \brief read last perspective from QSettings
  /// \return Returns true if it was successful. Otherwise false
  bool Restore();

  /// \brief save settings (internal use)
  /// \return Returns true if it was successful. Otherwise false
  bool Save();

  void ChangeVerbositySlider(int value);

private slots:
  void AddSideBarWidget();
  void VerbosityChanged(int value);

public slots:
  /// \brief set status message
  /// \param m The message that will be displayed in the status bar.
  void StatusMessage(const String& m);


private:
  void SetupQuickAccessBar();
  QWidget*        central_;
  QMenuBar*        menu_bar_;
  MainArea*       main_area_;
  PanelManager*   panels_;
  QWidget*        quick_access_bar_;
  QStatusBar*     status_bar_;
  QSlider* verbosity_slider_;
};

}}


#endif
