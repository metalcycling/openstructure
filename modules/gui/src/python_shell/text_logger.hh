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
#ifndef TEXT_LOGGER_HH_
#define TEXT_LOGGER_HH_

/*
  Author: Andreas Schenk
 */

#include "python_interpreter.hh"

#include <QObject>
#include <QTextStream>

namespace ost{namespace gui{

class DLLEXPORT_OST_GUI TextLogger:public QObject,public QTextStream
{
  Q_OBJECT
public:
  TextLogger();
  TextLogger( QIODevice * device );
  TextLogger( FILE * fileHandle, QIODevice::OpenMode openMode = QIODevice::ReadWrite );
  TextLogger( QString * String, QIODevice::OpenMode openMode = QIODevice::ReadWrite );
  TextLogger( const QByteArray & array, QIODevice::OpenMode openMode = QIODevice::ReadOnly );
  virtual ~TextLogger();
  void SetCodeLogging(bool flag);
  void SetOutputLogging(bool flag);
  void SetErrorLogging(bool flag);
  bool GetCodeLogging();
  bool GetOutputLogging();
  bool GetErrorLogging();
  
public slots:
  void AppendCode(const QString& text); 
  void AppendOutput(unsigned int id,const QString& text); 
  
protected:
  bool log_code_;
  bool log_output_;
  bool log_error_;
};

}}//ns

#endif /*TEXT_LOGGER_HH_*/
