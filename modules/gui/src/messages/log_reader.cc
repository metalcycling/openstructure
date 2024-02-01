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
/*
 Author: Stefan Scheuber
 */
#include <ost/log.hh>

#include "message_widget.hh"
#include "log_reader.hh"

namespace ost { namespace gui {

LogReader::LogReader(QObject* parent) :
  QObject(parent) {
  ost::Logger& logger = ost::Logger::Instance();
  LogSinkPtr lsp = LogSinkPtr(this);
  logger.PushSink(lsp);
}

void LogReader::LogMessage(const String& message, int severity){
  if(this->parent()){
    if(MessageWidget* message_widget = qobject_cast<MessageWidget*>(this->parent())){
      QString q_message(message.c_str());
      if(q_message.endsWith("\n")){
        q_message.remove(q_message.size()-1,q_message.size()-1);
        if(!log_cache_.contains(severity)){
          message_widget->LogMessage(q_message,GetIconForSeverity(severity));
        }
        else{
          message_widget->LogMessage(log_cache_[severity]+q_message,GetIconForSeverity(severity));
          log_cache_.remove(severity);
        }
      }
      else{
        if(!log_cache_.contains(severity)){
          log_cache_[severity] = message.c_str();
        }
        else{
          log_cache_[severity].append(message.c_str());
        }
      }
    }
  }
  std::cerr << message;
}

QMessageBox::Icon LogReader::GetIconForSeverity(int severity){
  switch (severity) {
    case 1:
      return QMessageBox::Warning;
    case 0:
      return QMessageBox::Critical;
    default :
      return QMessageBox::Information;
  }
}

LogReader::~LogReader() {
}

}
} // ns
