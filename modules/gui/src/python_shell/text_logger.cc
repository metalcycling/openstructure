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
#include "text_logger.hh"

/*
  Author: Andreas Schenk
 */

namespace ost{namespace gui{

TextLogger::TextLogger():
  QObject(),
  QTextStream(),
 log_code_(true),
 log_output_(true),
 log_error_(true)
{
}

TextLogger::TextLogger( QIODevice * device ):
  QObject(),
  QTextStream(device)
{
}

TextLogger::TextLogger( FILE * fileHandle, QIODevice::OpenMode openMode):
  QObject(),
  QTextStream(fileHandle,openMode)
{
}

TextLogger::TextLogger( QString * String, QIODevice::OpenMode openMode):
  QObject(),
  QTextStream(String,openMode)
{
}

TextLogger::TextLogger( const QByteArray & array, QIODevice::OpenMode openMode):
  QObject(),
  QTextStream(array,openMode)
{
}

TextLogger::~TextLogger()
{
}

void TextLogger::AppendCode(const QString& text)
{
  if(!log_code_){
    return;
  }
  operator<<(text);
  operator<<('\n');
  flush();
}
void TextLogger::AppendOutput(unsigned int id,const QString& text)
{
  operator<<(text);
  flush();
}
void TextLogger::SetCodeLogging(bool flag)
{
  log_code_=flag;
}

void TextLogger::SetOutputLogging(bool flag)
{
  log_output_=flag;  
}

void TextLogger::SetErrorLogging(bool flag)
{
  log_error_=flag;
}

bool TextLogger::GetCodeLogging()
{
  return log_code_;
}

bool TextLogger::GetOutputLogging()
{
  return log_output_;
}

bool TextLogger::GetErrorLogging()
{
  return log_error_;
}


}}//ns
