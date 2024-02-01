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
#ifndef OST_LOG_SINK_HH
#define OST_LOG_SINK_HH

#include <ostream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <ost/module_config.hh>

namespace ost {

class DLLEXPORT LogSink {
public:
  LogSink(){};
  virtual ~LogSink() { }
  virtual void LogMessage(const String& message, int severity=0) {};
};

typedef boost::shared_ptr<LogSink> LogSinkPtr;

class DLLEXPORT StreamLogSink : public LogSink {
public:
  StreamLogSink(std::ostream& stream):stream_(stream){}
  virtual void LogMessage(const String& message, int severity){
    stream_ << message;
  }

private:
  std::ostream& stream_;
};

class DLLEXPORT StringLogSink : public LogSink {
public:
  StringLogSink():LogSink(),stream_(){}
  virtual void LogMessage(const String& message, int severity){
    stream_ << message;
  }
  String GetLog() const
  {
    return stream_.str();
  }

private:
  std::ostringstream stream_;
};

typedef boost::shared_ptr<StringLogSink> StringLogSinkPtr;

class DLLEXPORT FileLogSink : public LogSink {
public:
  FileLogSink(const String& file_name):stream_(file_name.c_str(), std::ios::out){}
  virtual void LogMessage(const String& message, int severity){
    if (stream_.is_open()){
      stream_ << message;
      stream_.flush();
    }
  }

  ~FileLogSink(){
    stream_.flush();
  }
private:
  std::ofstream stream_;
};

typedef boost::shared_ptr<FileLogSink> FileLogSinkPtr;


class DLLEXPORT_OST_BASE MultiLogSink : public LogSink {
public:
  MultiLogSink();
  bool AddSink(LogSinkPtr& observer);
  bool RemoveSink(LogSinkPtr& observer);
  void LogMessage(const String& message, int severity);
private:
  std::vector<LogSinkPtr> sinks_;
};

typedef boost::shared_ptr<MultiLogSink> MultiLogSinkPtr;

}
#endif
