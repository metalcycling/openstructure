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
#include <ost/log_sink.hh>
#include <ost/log.hh>

namespace ost {

Logger& Log = Logger::Instance();

Logger& Logger::Instance()
{
  static Logger inst;
  return inst;
}

Logger::Logger():
  level_(0),
  level_stack_(),
  sink_stack_()
{
  sink_stack_.push(LogSinkPtr(new StreamLogSink(std::cerr)));
}

Logger::Logger(const Logger&):
  level_(0),
  level_stack_(),
  sink_stack_()
{
  sink_stack_.push(LogSinkPtr(new StreamLogSink(std::cerr)));
}

Logger& Logger::operator=(const Logger&)
{
  return *this;
}

void Logger::PushVerbosityLevel(int l)
{
  level_stack_.push(l);
  level_=l;
}

void Logger::PopVerbosityLevel()
{
  if (!level_stack_.empty()) {
    level_stack_.pop();
  }
  if (level_stack_.empty()) {
    level_ = 0;
  } else {
    level_ = level_stack_.top();
  }
}

void Logger::PushFile(const String& fn)
{
  sink_stack_.push(LogSinkPtr(new FileLogSink(fn)));
}

void Logger::PopFile()
{
  WARN_DEPRECATED("Logger::PopFile is deprecated. Use Logger::PopSink instead");
  this->PopSink();
}

void Logger::PushSink(LogSinkPtr& sink)
{
  sink_stack_.push(sink);
}

void Logger::PopSink()
{
  if (sink_stack_.size() > 1) {
    sink_stack_.pop();
  } else {
    LOG_ERROR("Can't pop sink. There is only one sink left on the stack");
  }
}

} // ns

