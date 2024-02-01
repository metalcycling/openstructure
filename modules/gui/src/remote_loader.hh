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
#ifndef OST_GUI_REMOTE_LOADER_HH
#define OST_GUI_REMOTE_LOADER_HH


#include <ost/gfx/gfx_object.hh>
#include <ost/gui/module_config.hh>
#include <ost/gui/widget.hh>

#include <QToolButton>
#include <QLineEdit>
#include <QKeyEvent>
#include <QProgressBar>
#include <QActionGroup>
namespace ost { namespace gui {




/// \brief a remote file loader
class DLLEXPORT_OST_GUI RemoteLoader : public Widget {
  Q_OBJECT
public:
  RemoteLoader(QWidget* parent=NULL);
  
  virtual bool Save(const QString& prefix);
  virtual bool Restore(const QString& prefix);
  virtual ActionList GetActions();
protected:
  virtual void keyPressEvent(QKeyEvent* event);
private:
  void ShowProgressBar(bool visiblility);
  void BuildMenu(String active_loader=String());
private slots:
  void Clicked();
  void UrlClick();
  void DownloadFinished();
  void UpdateProgress(qint64 read, qint64 total);
  void RenameButton();

private:
  QLineEdit* line_edit_;
  QToolButton* button_;
  QProgressBar* progress_bar_;
  ActionList action_list_;
  QMenu* site_loader_menu_;
  QActionGroup* site_actions_;
  String selected_site_loader_;
  bool img_support_;
};

}}

#endif
