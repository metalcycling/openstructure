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
#include <ost/log.hh>
#include <ost/log_sink.hh>
#include <iostream>
namespace ost {

MultiLogSink::MultiLogSink(){}

bool MultiLogSink::AddSink(LogSinkPtr& sink) {
  if ((std::find( sinks_.begin(), sinks_.end(), sink))==sinks_.end()) {
    sinks_.push_back(sink);
    return true;
  }
  return false;
}

bool MultiLogSink::RemoveSink(LogSinkPtr& sink){
  std::vector<LogSinkPtr>::iterator found=std::find(sinks_.begin(), 
                                                    sinks_.end(), sink);
  if (found!=sinks_.end() ){
    sinks_.erase(found);
    return true;
  }
  return false;
}

void MultiLogSink::LogMessage(const String& message, int severity){
  for (std::vector<LogSinkPtr>::const_iterator 
       i=sinks_.begin(), e=sinks_.end(); i!=e; ++i) {
    (*i)->LogMessage(message, severity);        
  }
}

} // ns

