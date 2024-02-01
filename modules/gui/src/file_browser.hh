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
#ifndef OST_GUI_FILE_BROWSER_HH
#define OST_GUI_FILE_BROWSER_HH


#include <ost/gfx/gfx_object.hh>
#include <ost/gui/module_config.hh>
#include <ost/gui/widget.hh>

#include <QModelIndex>
#include <QDir>

class QDirModel;
class QListView;
class QComboBox;

namespace ost { namespace gui {




/// \brief a simplistic file browser
class FileBrowser : public Widget {
  Q_OBJECT
public:
  FileBrowser(QWidget* parent=NULL);
  FileBrowser(const QString& path, QWidget* parent=NULL);
  
  virtual ActionList GetActions();
  virtual bool Save(const QString& prefix);
  virtual bool Restore(const QString& prefix);

protected:
  virtual void keyPressEvent(QKeyEvent* event);
public slots:
  void Refresh();
private slots:
  void Init(const QString& path);
  void DoubleClicked(const QModelIndex& index);
  void ChangeToParentDirectory(int index);
  void Split();
  void ShowContextMenu(const QPoint& pos);
  void LoadCurrentObject();
  void LoadWithSystemEditor();
  void LoadWithSourceViewer();

private:
  void LoadObject(const QModelIndex& index);
  void UpdateMenu(const QString& path);
  void AddItem(const QDir& directory, const QString& mypath="");
  QComboBox* menu_;
  QDirModel* model_;
  QListView* view_;
  ActionList action_list_;
};

}}

#endif
