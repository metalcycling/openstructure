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
#ifndef OST_GUI_FILE_TYPE_DIALOG_HH
#define OST_GUI_FILE_TYPE_DIALOG_HH

/*
  Author: Stefan Scheuber, Marco Biasini
 */

#include <ost/gui/module_config.hh>

#ifndef Q_MOC_RUN

#include <ost/io/mol/entity_io_handler.hh>
#include <ost/io/seq/sequence_io_handler.hh>
#include <ost/io/mol/surface_io_handler.hh>
#include <ost/io/img/map_io_handler.hh>

#endif

#include <QDialog>
#include <QMetaType>

class QTableWidget;
class QLabel;
namespace ost { namespace gui {

/// \brief dialog to select a file type
class DLLEXPORT_OST_GUI FileTypeDialog : public QDialog {
  Q_OBJECT  
public:
  FileTypeDialog(const QString& file_name, QWidget* parent=NULL);
  io::EntityIOHandlerP GetEntityHandler();
  io::SequenceIOHandlerPtr GetSequenceHandler();
  io::SurfaceIOHandlerPtr GetSurfaceHandler();
  io::MapIOHandlerPtr GetMapHandler();
public slots:
  virtual void accept();
private:
  QTableWidget* list_;
  QLabel* label_;
  io::EntityIOHandlerP entity_handler_;
  io::SequenceIOHandlerPtr seq_handler_;
  io::SurfaceIOHandlerPtr surf_handler_;
  io::MapIOHandlerPtr map_handler_;

};

}}

Q_DECLARE_METATYPE(ost::io::EntityIOHandlerFactoryBaseP);
Q_DECLARE_METATYPE(ost::io::SequenceIOHandlerFactoryBasePtr);
Q_DECLARE_METATYPE(ost::io::SurfaceIOHandlerFactoryBasePtr);
Q_DECLARE_METATYPE(ost::io::MapIOHandlerFactoryBasePtr);

#endif
